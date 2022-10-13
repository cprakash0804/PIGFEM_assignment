/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "SolverLinear.h"
#include "Problem.h"
#include "Assembler.h"
#include "Mesh.h"
#include "DofObject.h"
#include "NodalData.h"
#include "SubscaleModel.h"
#include "BodyLoad.h"
#include "Utilities.h"
#include "SensitivitySolver.h"
#include <iostream>



PETScLinearKSPSolver::PETScLinearKSPSolver()
{
}

PETScLinearKSPSolver::~PETScLinearKSPSolver()
{
}



// This function determines the local and global sizes of all vectors used in the global solution procedure
PetscErrorCode PETScLinearKSPSolver::preallocate_vectors()
{
	PetscErrorCode ierr;

	// If this FE system has already been set up before, I need to clear ou the old vectors
	if (_setup)
	{
		_U.destroy();
		_F_ext.destroy();
	}

	PetscInt n, N, n_const, N_const;
	DofObject* dofs = _prob->get_dofs();
	Mesh* mesh = _prob->get_mesh();
	n = dofs->n_local_owned_free_dofs();
	N = dofs->n_global_free_dofs();
	n_const = dofs->n_local_owned_const_dofs();
	N_const = dofs->n_global_const_dofs();

	// if(mesh->serial())
	// {
	// 	ierr = VecCreateSeq(mesh->get_comm(), N, _U.f);CHKERRQ(ierr);
	// 	ierr = VecCreateSeq(mesh->get_comm(), N, _F_ext.f);CHKERRQ(ierr);
	// 	ierr = VecCreateSeq(mesh->get_comm(), N_const, _U.p);CHKERRQ(ierr);
	// 	ierr = VecCreateSeq(mesh->get_comm(), N_const, _F_ext.p);CHKERRQ(ierr);
	// }
	// else
	// {
	// 	ierr = VecCreateMPI(mesh->get_comm(), n, N, _U.f);CHKERRQ(ierr);
	// 	ierr = VecCreateMPI(mesh->get_comm(), n, N, _F_ext.f);CHKERRQ(ierr);
	// 	ierr = VecCreateMPI(mesh->get_comm(), n_const, N_const, _U.p);CHKERRQ(ierr);
	// 	ierr = VecCreateMPI(mesh->get_comm(), n_const, N_const, _F_ext.p);CHKERRQ(ierr);
	// }
	if(mesh->serial())
	{
		ierr = VecCreateSeq(mesh->get_comm(), N, &_U.f);CHKERRQ(ierr);
		ierr = VecCreateSeq(mesh->get_comm(), N, &_F_ext.f);CHKERRQ(ierr);
		ierr = VecCreateSeq(mesh->get_comm(), N_const, &_U.p);CHKERRQ(ierr);
		ierr = VecCreateSeq(mesh->get_comm(), N_const, &_F_ext.p);CHKERRQ(ierr);
	}
	else
	{
		ierr = VecCreateMPI(mesh->get_comm(), n, N, &_U.f);CHKERRQ(ierr);
		ierr = VecCreateMPI(mesh->get_comm(), n, N, &_F_ext.f);CHKERRQ(ierr);
		ierr = VecCreateMPI(mesh->get_comm(), n_const, N_const, &_U.p);CHKERRQ(ierr);
		ierr = VecCreateMPI(mesh->get_comm(), n_const, N_const, &_F_ext.p);CHKERRQ(ierr);
	}

	return ierr;
}



// Solves the global system using PETSc.
// First assembles the matricies and vectors using assemble_linear_elasticity
// Then solves the partitioned system using the PETSc KSP context
PetscErrorCode PETScLinearKSPSolver::solve()
{
	//	Variables
	//
	//	K		-	Matrix that defines the linear system
	//	ksp		-	KSP context
	//	b, u	-	RHS, exact solution vectors

	KSP ksp;
	// AO ao;   				// Application Ordering object
	Vec temp1, temp2, Pf_store;
	PetscErrorCode ierr;
	Utilities::timer timer;

	// Print to screen
	PIGFEMPrint("\nSolving the Linear Problem...\n\n\n");

	// Preallocate all memory for matricies and vectors
	ierr = VecDuplicate(_F_ext.f, &Pf_store);CHKERRQ(ierr);
	ierr = VecDuplicate(_F_ext.f, &temp1);CHKERRQ(ierr);		  // Copy correct structure
	ierr = VecDuplicate(_F_ext.p, &temp2);CHKERRQ(ierr);
	ierr = VecCopy(_F_ext.f, Pf_store);CHKERRQ(ierr);      // Stores a copy of Ff for use later
	
	// Actually assemble the linear system
	Assembler* assem = _prob->get_assembler();
	timer.start();
	ierr = assem->assemble_new_load_step(_U.p, _F_ext, 1.0);CHKERRQ(ierr);
	if (getNonzeroStart())
		{ierr = VecAXPY(_U.p, 1.0, _Up_initial);CHKERRQ(ierr);}
	ierr = assem->assemble_linear(_K);CHKERRQ(ierr);
	_prob->_assemble_time = timer.time_elapsed();

	// Create the KSP context and set operators
	ierr = initKSP(ksp, true);CHKERRQ(ierr);

	// Start the solve process
	timer.start();
	ierr = MatMult(_K.fp, _U.p, temp1);CHKERRQ(ierr);    // Will store the product of Kfp*up in temp1
	ierr = VecAXPY(_F_ext.f, -1.0, temp1);CHKERRQ(ierr);  // Stored Pf - Kfp*Up in Pf. This will be te rhs of our solve
	
	// Actually Solve the system!!!!
	ierr = KSPSolve(ksp, _F_ext.f, _U.f);CHKERRQ(ierr);  // Solve the system
	
	// Compute the reaction forces
	ierr = MatMult(_K.pf, _U.f, temp2);CHKERRQ(ierr);
	ierr = MatMult(_K.pp, _U.p, _F_ext.p);CHKERRQ(ierr);
	ierr = VecAXPY(_F_ext.p, 1.0, temp2);CHKERRQ(ierr);
	_prob->_solve_time = timer.time_elapsed();
	
	// Stores the solution in the Problem object in an easier to access structure
	ierr = VecCopy(Pf_store, _F_ext.f);CHKERRQ(ierr);    // Restore the copy of Pf
	store_solution(_U, _prob->get_solution());
	store_solution(_F_ext, _prob->get_external_load());

	if (_prob->sensitivity())
	{
		timer.start();
		_prob->get_sensitivity_solver()->solveProblem();
		_prob->_sensitivity_time += timer.time_elapsed();
	}

	Output(0, 1.0);
	
	// Clean up local memory
	ierr = VecDestroy(&temp1);CHKERRQ(ierr);
	ierr = VecDestroy(&temp2);CHKERRQ(ierr);
	ierr = VecDestroy(&Pf_store);CHKERRQ(ierr);
	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

	_prob->get_solved() = true;
	return ierr;
}

