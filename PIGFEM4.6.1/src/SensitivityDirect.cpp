/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated April 2017

##################################################################################
*/
#include "SensitivityDirect.h"
#include "Problem.h"
#include "Mesh.h"
#include "DofObject.h"
#include "SensitivityParameter.h"
#include "SensitivityFunction.h"
#include "NodalData.h"
#include "InternalVars.h"
#include "Solver.h"
#include "Utilities.h"
#include "Options.h"



// Solve the direct sensitivity problem for all of the functions and parameters
PetscErrorCode SensitivityDirect::solveProblem()
{
	Utilities::timer timer;
	PetscErrorCode ierr = 0;
	Mesh* mesh = _prob->get_mesh();
	PIGFEMPrint("COMPUTING SENSITIVITIES USING DIRECT METHOD\n");

	// Compute the function values
	for (id_type f=0; f<n_functions(); ++f)
		_function_vals[f] = _functions[f]->evaluate(_K, _U, _F_ext);
			

	// Solve the multiple right hand sides individully
	if (_individual) // No real testing there, just some standard
	{
		timer.start();
		assembleDerivatives();
		_sensitivity_assemble_time += timer.time_elapsed();

		// For every parameter, get the derivative of Uf and plug that into the direct equation
		for (id_type i=0; i<_parameters.size(); ++i)
		{
			// Create the Pseudo load
			timer.start();
			ierr = VecCopy(_delP_deld[i].f, _Pseudo);CHKERRQ(ierr);
			ierr = VecScale(_Pseudo, -1.0); CHKERRQ(ierr);
			_sensitivity_assemble_time += timer.time_elapsed();

			// Solve the direct system
			timer.start();
			ierr = KSPSolve(_ksp, _Pseudo, _dU_dd_vec.f);CHKERRQ(ierr);
			double solve_time = timer.time_elapsed();
			_sensitivity_solve_time += solve_time;

			// Plug into the direct equation for every function
			PetscScalar dot_ans;
			timer.start();
			for (id_type j=0; j<_functions.size(); ++j)
			{
				double sensitivity = 0.0;
				ierr = VecDot(_dU_dd_vec.f, _delf_delU[j].f, &dot_ans);CHKERRQ(ierr);
				sensitivity = dot_ans;
				sensitivity += _delf_deld[j][i];
				_sensitivities[j][i] = sensitivity;
			}
			_sensitivity_substitution_time += timer.time_elapsed();

			
			// Now store the U derivative field in a NodalData struct and use that to update the ISVs for all materials in the mesh
			if (_internal_vars->haveISVs())
			{
				_prob->get_solver()->store_solution(_dU_dd_vec, _dU_dd);
				_parameters[i]->updateSensitivityISVs(_prob, _dU_dd);
			}
		}
	}















	// Solve the system once using a direct solver and using MatMatSolve
	else
	{
		// Assemble the individual parameter derivatives
		timer.start();
		assembleDerivatives();

		// Assemble the individual parameter derivaives into the RHS matrix
		PetscInt low, high, n;
		n = _prob->get_dofs()->n_local_owned_free_dofs();
		ierr = VecGetOwnershipRange(_delP_deld[0].f, &low, &high);CHKERRQ(ierr);
		std::vector<PetscInt> rows(n);
		for (int i=0; i<n; ++i)
			rows[i] = i + low;
		for (id_type i=0; i<_parameters.size(); ++i)
		{
			ierr = VecScale(_delP_deld[i].f, -1.0); CHKERRQ(ierr);
			PetscScalar *array;
			ierr = VecGetArray(_delP_deld[i].f, &array);CHKERRQ(ierr);
			PetscInt cols = i;

			ierr = MatSetValues(_RHS, n, rows.data(), 1, &cols, array, INSERT_VALUES);CHKERRQ(ierr);
			VecRestoreArray(_delP_deld[i].f, &array);CHKERRQ(ierr);
		}
		ierr = MatAssemblyBegin(_RHS, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(_RHS, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		_sensitivity_assemble_time += timer.time_elapsed();

		// Create the matrix reordering index sets
		if (_step==0)
			ierr = MatGetOrdering((_K->ff), MATORDERINGND, &_rperm, &_cperm);CHKERRQ(ierr);
		// Factor the stiffness matrix (MATSOLVERSUPERLU_DIST MATSOLVERPETSC MATSOLVERMUMPS)
		// Solve the system using the RHS matrix
		Mat Fact;
		MatFactorInfo info;
		timer.start();
		ierr = MatGetFactor((_K->ff), MATSOLVERMUMPS, MAT_FACTOR_CHOLESKY, &Fact);CHKERRQ(ierr);
		info.fill = 5.0; // ????????????????????????????????????
		ierr = MatCholeskyFactorSymbolic(Fact, (_K->ff), _rperm, &info);CHKERRQ(ierr);
		ierr = MatCholeskyFactorNumeric(Fact,(_K->ff),&info);CHKERRQ(ierr);
		ierr = MatMatSolve(Fact, _RHS, _dUf_dd_block);CHKERRQ(ierr);
		ierr = MatAssemblyBegin(_dUf_dd_block, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(_dUf_dd_block, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatDestroy(&Fact);CHKERRQ(ierr);
		_sensitivity_solve_time += timer.time_elapsed();


		// Substitue into the direct sensitivity equations and distribute the solutions to all procs
		timer.start();
		for (id_type j=0; j<_functions.size(); ++j)
		{
			ierr = MatMultTranspose(_dUf_dd_block, _delf_delU[j].f, _Sens);CHKERRQ(ierr);
			PetscScalar* array;
			if (!mesh->serial())
			{
				ierr = VecScatterBegin(_scat_sens, _Sens, _seq_Sens, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
				ierr = VecScatterEnd(_scat_sens, _Sens, _seq_Sens, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
				ierr = VecGetArray(_seq_Sens, &array);CHKERRQ(ierr);
				for (id_type i=0; i<n_parameters(); ++i)
					_sensitivities[j][i] = array[i] + _delf_deld[j][i];
				ierr = VecRestoreArray(_seq_Sens, &array);
			}
			else
			{
				ierr = VecGetArray(_Sens, &array);CHKERRQ(ierr);
				for (id_type i=0; i<n_parameters(); ++i)
					_sensitivities[j][i] = array[i] + _delf_deld[j][i];
				ierr = VecRestoreArray(_Sens, &array);
			}
		}
		_sensitivity_substitution_time += timer.time_elapsed();
		

		// Now plug into the ISV update equations
		if (_internal_vars->haveISVs())
		{
			for (id_type i=0; i<n_parameters(); ++i)
			{
				ierr = MatGetColumnVector(_dUf_dd_block, _dU_dd_vec.f, i);CHKERRQ(ierr);
				_prob->get_solver()->store_solution(_dU_dd_vec, _dU_dd);
				_parameters[i]->updateSensitivityISVs(_prob, _dU_dd);
			}
		}
	}

	_step++;
	return ierr;
}






SensitivityDirect::SensitivityDirect()
	: _step(0), _individual(true)
{}


SensitivityDirect::~SensitivityDirect()
{
	// Destroy all things that existed no matter what
	if (_dU_dd != NULL)
		delete _dU_dd;

	if (_setup)
	{
		// _dU_dd_vec.destroy();
		// _dU_dd_vec.free();

		// Destroy the things that existed if we were running in individual mode
		if (_individual)
		{
			VecDestroy(&_Pseudo);
			KSPDestroy(&_ksp);
		}

		// Otherwise destroy the blocked variabled
		else
		{
			MatDestroy(&_RHS);
			MatDestroy(&_dUf_dd_block);
			ISDestroy(&_rperm);
			ISDestroy(&_cperm);
			VecDestroy(&_Sens);
			if (!_prob->get_mesh()->serial())
			{
				VecDestroy(&_seq_Sens);
				VecScatterDestroy(&_scat_sens);
			}
		}
	}
}
















PetscErrorCode SensitivityDirect::setup()
{
	PetscErrorCode ierr;
	Mesh* mesh = _prob->get_mesh();

	// Preallocate the things that will exist no matter what
	if (_setup)
		_dU_dd_vec.destroy();
	ierr = preallocateVector(&_dU_dd_vec.f, true);CHKERRQ(ierr);
	ierr = preallocateVector(&_dU_dd_vec.p, false);CHKERRQ(ierr);
	ierr = VecZeroEntries(_dU_dd_vec.p);

	if (n_parameters() <= 5)
	{
		if (_setup && _individual)
		{
			ierr = VecDestroy(&_Pseudo);CHKERRQ(ierr);
			ierr = KSPDestroy(&_ksp);CHKERRQ(ierr);
		}
		_individual = true;

		// Create the psudo load vector
		ierr = preallocateVector(&_Pseudo, true);CHKERRQ(ierr);

		// Create the ksp context
		ierr = _prob->get_solver()->initKSP(_ksp, true);
	}

	else
	{
		if (_setup && !_individual)
		{
			ierr = MatDestroy(&_RHS);CHKERRQ(ierr);
			ierr = MatDestroy(&_dUf_dd_block);CHKERRQ(ierr);
			ierr = VecDestroy(&_Sens);CHKERRQ(ierr);
			if (!mesh->serial())
			{
				ierr = VecDestroy(&_seq_Sens);CHKERRQ(ierr);
				ierr = VecScatterDestroy(&_scat_sens);CHKERRQ(ierr);
			}
		}
		_individual = false;

		// Create the dense matrices
		PetscInt n, N;
		n = _prob->get_dofs()->n_local_owned_free_dofs();
		N = _prob->get_dofs()->n_global_free_dofs();
		if (mesh->serial())
		{
			ierr = MatCreateSeqDense(mesh->get_comm(), N, n_parameters(), NULL, &_RHS);CHKERRQ(ierr);
			ierr = MatCreateSeqDense(mesh->get_comm(), N, n_parameters(), NULL, &_dUf_dd_block);CHKERRQ(ierr);
		}
		else
		{
			ierr = MatCreateDense(mesh->get_comm(), n, PETSC_DECIDE, N, n_parameters(), NULL, &_RHS);CHKERRQ(ierr);
			ierr = MatCreateDense(mesh->get_comm(), n, PETSC_DECIDE, N, n_parameters(), NULL, &_dUf_dd_block);CHKERRQ(ierr);
		}

		

		// Create the Vectors and scatter to get the vector of all parameter sensitivities on each process
		if (mesh->serial())
		{
			ierr = VecCreateSeq(mesh->get_comm(), n_parameters(), &_Sens);CHKERRQ(ierr);
		}
		else
		{
			ierr = VecCreateMPI(mesh->get_comm(), PETSC_DECIDE, n_parameters(), &_Sens);CHKERRQ(ierr);
			ierr = VecScatterCreateToAll(_Sens, &_scat_sens, &_seq_Sens);CHKERRQ(ierr);
		}
	}

	// The rest of initializations
	ierr = SensitivitySolver::setup();
	return ierr;
}

PetscErrorCode SensitivityDirect::init()
{
	_dU_dd = new NodalData( _prob->get_mesh() );
	_dU_dd->preallocate_storage( _prob->nndof() );

	return SensitivitySolver::init();
}