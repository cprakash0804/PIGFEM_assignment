/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated April 2017

##################################################################################
*/
#include "SensitivityAdjoint.h"
#include "Problem.h"
#include "Mesh.h"
#include "Utilities.h"
#include "SensitivityParameter.h"
#include "SensitivityFunction.h"
#include "Solver.h"
#include "Options.h"




// Solve the direct sensitivity problem for all of the functions and parameters
PetscErrorCode SensitivityAdjoint::solveProblem()
{
	PetscErrorCode ierr = 0;
	PIGFEMPrint("COMPUTING SENSITIVITIES USING ADJOINT METHOD\n");

	// Assmble the current derivatives
	Utilities::timer timer;
	timer.start();
	assembleDerivatives();
	_sensitivity_assemble_time += timer.time_elapsed();

	// Compute the function values
	for (id_type f=0; f<_functions.size(); ++f)
		_function_vals[f] = _functions[f]->evaluate(_K, _U, _F_ext);

	// For every parameter, get the solution to the adjoint variable equation and plug it into the adjoint equation
	for (id_type i=0; i<_functions.size(); ++i)
	{
		// Solve the adjoint system (No point in building the pseudo load since its just one vector...)
		timer.start();
		ierr = KSPSolveTranspose(_ksp, _delf_delU[i].f, _Lambda);CHKERRQ(ierr);
		_sensitivity_solve_time += timer.time_elapsed();

		// Plug into the direct equation for every function
		PetscScalar dot_ans;
		timer.start();
		for (id_type j=0; j<_parameters.size(); ++j)
		{
			double sensitivity = 0.0;
			ierr = VecDot(_Lambda, _delP_deld[j].f, &dot_ans);CHKERRQ(ierr);
			sensitivity -= dot_ans;
			sensitivity += _delf_deld[i][j];
			_sensitivities[i][j] = sensitivity;
		}
		_sensitivity_substitution_time += timer.time_elapsed();
	}

	// Something about the interal state variables???????????/?

	return ierr;
}







SensitivityAdjoint::SensitivityAdjoint()
{}


SensitivityAdjoint::~SensitivityAdjoint()
{
	// Destroy all the PETSc objects
	if (_setup)
	{
		VecDestroy(&_Lambda);
		KSPDestroy(&_ksp);
	}
}
















PetscErrorCode SensitivityAdjoint::setup()
{
	PetscErrorCode ierr;

	// Preallocate the adjoint vector
	if (_setup)
	{
		ierr = VecDestroy(&_Lambda);CHKERRQ(ierr);
		ierr = KSPDestroy(&_ksp);CHKERRQ(ierr);
	}
	ierr = preallocateVector(&_Lambda, true);CHKERRQ(ierr);

	// Create the ksp context
	ierr = _prob->get_solver()->initKSP(_ksp, true);
	
	// The rest of initializations
	ierr = SensitivitySolver::setup();
	return ierr;
}