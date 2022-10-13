/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated July 2017

##################################################################################
*/
#include "SolverNonlinearKSP.h"
#include "Problem.h"
#include "Assembler.h"
#include "Options.h"
#include "Utilities.h"
#include <sstream>

SolverNonlinearKSP::SolverNonlinearKSP()
{
}
SolverNonlinearKSP::~SolverNonlinearKSP()
{
	KSPDestroy(&_ksp);
	VecDestroy(&_Res_f);
	VecDestroy(&_delta_Uf);
}




PetscErrorCode SolverNonlinearKSP::setup()
{
	PetscErrorCode ierr;

	// Clear any old structures
	if (_setup)
	{
		KSPDestroy(&_ksp);
		VecDestroy(&_delta_Uf);
	}

	// Create the basic matricies and vectors
	ierr = SolverNonlinear::setup();CHKERRQ(ierr);

	// Set up all specific structures here
	ierr = initKSP(_ksp, true);CHKERRQ(ierr);
	ierr = VecDuplicate(_U.f, &_delta_Uf);CHKERRQ(ierr);
	return ierr;
}






PetscErrorCode SolverNonlinearKSP::nonlinearStep(std::vector<std::vector<std::vector<double> > >& update_ISVs, bool& converged, double& Res_norm)
{
	// Some initializations
	PetscErrorCode ierr;
	double prev_Res_norm;
	id_type NR_iter = 0;
	id_type increase_count = 0; // Marks how many interations in a row the residual increases
	Assembler* assem = _prob->get_assembler();
	Utilities::timer timer;

	// Form the initial tangent stiffness matrix and internal force vector
	// For elastic problems this is a waste and will be exactly the same as the assembled system at the end
	//	of the previous load step. For inelastic systems, however, the updated ISVs are taken into account here
	// 	so I think it should actually be different. I suppose that is something I can test
	timer.start();
	ierr = assem->assemble_nonlinear(_K, _P_int, _delta_t, update_ISVs, _prob->get_solution(), true, true);CHKERRQ(ierr);
	// if (step_iter == 0)
	// {
	// 	ierr = MatSetOption(*_K.ff, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);CHKERRQ(ierr);
	// 	ierr = MatSetOption(*_K.fp, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);CHKERRQ(ierr);
	// 	ierr = MatSetOption(*_K.pf, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);CHKERRQ(ierr);
	// 	ierr = MatSetOption(*_K.pp, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);CHKERRQ(ierr);
	// }
	_prob->_assemble_time += timer.time_elapsed();

	// Form the initial residual
	ierr = VecCopy(_F_ext.f, _Res_f);CHKERRQ(ierr);
	ierr = VecAXPY(_Res_f, -1.0, _P_int.f);CHKERRQ(ierr);
	ierr = VecNorm(_Res_f, NORM_2, &Res_norm);CHKERRQ(ierr);
	double Res0_norm = Res_norm;

	// NEWTON-RAPHSON LOOP
	//---------------------------
	
	PIGFEMPrint("\tNR iteration: " << NR_iter << " Residual: " << Res_norm << std::endl);
	while(1)
	{
//		PIGFEMPrint("\tHEREREEE " << std::endl);
		// Solve for the change in displacement
		timer.start();
		ierr = KSPSolve(_ksp, _Res_f, _delta_Uf);CHKERRQ(ierr);	// Solve the system
		_prob->_solve_time += timer.time_elapsed();
		double iter_time = timer.time_elapsed();

		// Update the solution values
		ierr = VecAXPY(_U.f, 1.0, _delta_Uf);CHKERRQ(ierr);
		store_solution(_U, _prob->get_solution());

		// Assemble the new stiffness matrix and inernal force vectors
		timer.start();
		ierr = assem->assemble_nonlinear(_K, _P_int, _delta_t, update_ISVs, _prob->get_solution(), true, true);CHKERRQ(ierr);
		_prob->_assemble_time += timer.time_elapsed();

		// Compute the new residual vector
		//ierr = VecCopy(_Pf_int, Res_f);CHKERRQ(ierr);
		//ierr = VecAXPY(Res_f, -1.0, _Ff_ext);CHKERRQ(ierr);
		prev_Res_norm = Res_norm;
		ierr = VecCopy(_F_ext.f, _Res_f);CHKERRQ(ierr);
		ierr = VecAXPY(_Res_f, -1.0, _P_int.f);CHKERRQ(ierr);
		ierr = VecNorm(_Res_f, NORM_2, &Res_norm);CHKERRQ(ierr);

		// Increment the Newton-Raphson iteration counter
		NR_iter++;
		PIGFEMPrint("\tNR iteration: " << NR_iter  << " Iteration_Time: " << iter_time << " \tResidual: " << Res_norm << std::endl);

		// Check if convergence criteria is met (Converged)
		if (Res_norm < _rel_tol*Res0_norm
			|| Res_norm < _abs_tol)
		{

			converged = true;
			break;
		}

		// For some reason when a Nan is encountered KSP stores it somewhere and it shows up again the next time KSPSolve is ran...
		// Only way around this I've found is to destroy and reinitialize the KSP context
		if (!std::isfinite(Res_norm)) // NaN check
		{
			ierr = KSPDestroy(&_ksp);CHKERRQ(ierr);
			ierr = initKSP(_ksp, true);CHKERRQ(ierr);
			converged = false;
			break;
		}

		// Check the iteration limit counter (Not converged)
		if(NR_iter >= _max_iter)
		{
			converged = false;
			break;
		}

		// Final break condition. Check if the residual has increased for more than x iterations
		if (Res_norm > prev_Res_norm)
			increase_count++;
		else
			increase_count = 0; // reset
		if (increase_count >= 3)
		{
			converged = false;
			break;
		}
	}	// End NR loop --------------------------

	return ierr;
}
