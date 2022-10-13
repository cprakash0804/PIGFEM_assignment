/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated July 2017

##################################################################################
*/
#include "SolverNonlinearSNES.h"
#include "Problem.h"
#include "Assembler.h"
#include "Mesh.h"
#include "Options.h"
#include "Utilities.h"
#include "NodalData.h"
#include <sstream>


typedef struct
{
	std::vector<std::vector<std::vector<double> > >* update_ISVs;
	SolverNonlinearSNES* solver;
	double delta_t;
} AppCtx;






SolverNonlinearSNES::SolverNonlinearSNES()
	: _misc_time(0.0), _assem_time(0.0)
{
}
SolverNonlinearSNES::~SolverNonlinearSNES()
{
	SNESDestroy(&_snes);
	VecDestroy(&_Res_f);
}




PetscErrorCode SolverNonlinearSNES::setup()
{
	PetscErrorCode ierr;

	// Clear any previous structures
	if (_setup)
	{
		SNESDestroy(&_snes);
		VecDestroy(&_Res_f);
	}

	// Create the basic matrixcies and vectors
	ierr = SolverNonlinear::setup();CHKERRQ(ierr);
	
	// Create the nonlinear solver
	ierr = SNESCreate(_prob->get_mesh()->get_comm(),&_snes);CHKERRQ(ierr);

	// Initialize the linear solver used within the nonlinear wrapper
	KSP ksp;
	ierr = SNESGetKSP(_snes, &ksp);
	ierr = initKSP(ksp, false);CHKERRQ(ierr);
	ierr = SNESSetType(_snes, SNESNEWTONLS);CHKERRQ(ierr);
	// ierr = SNESSetType(_snes, SNESNEWTONTR);CHKERRQ(ierr);

	/* Set nonlinear solver parameters
	 * Defaults:
	 *	atol: 1e-50
	 *	rtol: 1e-8
	 *	stol: 1e-8
	 *	maxits: 50
	 *	maxFs: 10000
	*/
	// ierr = SNESSetTolerances(_snes, _abs_tol, _rel_tol, PETSC_DEFAULT, _max_iter, PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = SNESSetTolerances(_snes, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 50, PETSC_DEFAULT);CHKERRQ(ierr);

	return ierr;
}







/*
 * Function to for the residual of the governing equations
 * Inputs:
 *	1. snes - SNES context (Don't know why I actually need this)
 *	2. x - current solution vector
 *	3. f - the residual vector that will be filled
 *	4. ctx - User context described below
 * Can get following directly:
 *	1. assembler object
 *	2. Current external force
 *	3. Prescribed displacement vector
 *	4. Current solution (I should actually use the solution provided here I think? ow to optimize assembler then?)
 * The user context contains:
 *	1. ISVs
 */
PetscErrorCode formResidual(SNES snes,Vec x,Vec f,void *ctx)
{

	Utilities::timer timer;
	timer.start();
	PetscErrorCode ierr = 0;
	AppCtx* user = (AppCtx*)ctx;
	Problem* prob = user->solver->_prob;
	Assembler* assem = prob->get_assembler();

	// Create a Nodal data struct with the current solution x
	NodalData curr_sol(prob->get_mesh());
	curr_sol.preallocate_storage(prob->nndof());
	user->solver->store_solution(x, user->solver->_U.p, &curr_sol);
	user->solver->_misc_time += timer.time_elapsed();
	timer.start();
	assem->assemble_nonlinear(user->solver->_K, user->solver->_P_int, user->delta_t,
							*(user->update_ISVs), &curr_sol, true, false);

	user->solver->_assem_time += timer.time_elapsed();

	timer.start();
	ierr = VecCopy((user->solver->_P_int.f), f);CHKERRQ(ierr);
	ierr = VecAXPY(f, -1.0, (user->solver->_F_ext.f));
	user->solver->_misc_time += timer.time_elapsed();
	return ierr;
}


/*
 * Function to for the jacobian of the governing equations
 */
PetscErrorCode formJacobian(SNES snes,Vec x,Mat Amat,Mat Pmat,void *ctx)
{
	Utilities::timer timer;
	timer.start();
	if (Amat != Pmat)
		err_message("When using SNES Solver, jacobian and preconditioner must be he same matrix");

	PetscErrorCode ierr = 0;
	AppCtx* user = (AppCtx*)ctx;
	Problem* prob = user->solver->_prob;
	Assembler* assem = prob->get_assembler();

	if (Amat != user->solver->_K.ff)
		err_message("In the SNES solver, the jacobian matrix doen't match the tangent matrix");

	// Create a Nodal data struct with the current solution x
	NodalData curr_sol(prob->get_mesh());
	curr_sol.preallocate_storage(prob->nndof());
	user->solver->store_solution(x, user->solver->_U.p, &curr_sol);
	user->solver->_misc_time += timer.time_elapsed();
	timer.start();
	assem->assemble_nonlinear(user->solver->_K, user->solver->_P_int, user->delta_t,
							*(user->update_ISVs), &curr_sol, false, true);
	user->solver->_assem_time += timer.time_elapsed();

	return ierr;
}




PetscErrorCode monitorFunction(SNES snes, PetscInt its, PetscReal norm, void *ctx)
{
	AppCtx* user = (AppCtx*)ctx;
	Utilities::timer timer;
	timer.start();
	PIGFEMPrint("\tNR iteration: " << its << " Residual: " << norm << std::endl);
	user->solver->_misc_time += timer.time_elapsed();
	return 0;
}






PetscErrorCode SolverNonlinearSNES::nonlinearStep(std::vector<std::vector<std::vector<double> > >& update_ISVs, bool& converged, double& Res_norm)
{
	// Some initializations
	PetscErrorCode ierr = 0;
	AppCtx user;
	user.update_ISVs = &update_ISVs;
	user.solver = this;
	user.delta_t = _delta_t;
	Utilities::timer timer;
	_misc_time = 0.0;
	_assem_time = 0.0;

	// Utilities::timer timer;

	// Set the residual function
	ierr = SNESSetFunction(_snes, _Res_f, formResidual, (void*)&user);CHKERRQ(ierr);
	// Set the jacobian function
	ierr = SNESSetJacobian(_snes, _K.ff, _K.ff, formJacobian, (void*)&user);CHKERRQ(ierr);

	// Set the progress monitor function
	ierr = SNESMonitorSet(_snes, monitorFunction, (void*)&user, NULL);CHKERRQ(ierr);

	// Solve the nonlinear system (Use a RHS of zero)
	timer.start();
	ierr = SNESSolve(_snes, NULL, _U.f);CHKERRQ(ierr);
	double total_time = timer.time_elapsed();
	_prob->_solve_time += (total_time - _misc_time - _assem_time);
	_prob->_assemble_time = _assem_time;
	store_solution(_U, _prob->get_solution());

	// Get the final converged? residual
	ierr = SNESGetFunctionNorm(_snes, &Res_norm);CHKERRQ(ierr);

	// Get whether or not the solver converged
	SNESConvergedReason reason;
	ierr = SNESGetConvergedReason(_snes, &reason);CHKERRQ(ierr);
	/* Convergence test
	 * 2 = absolute residual test
	 * 3 = relative residual test
	 * 4 = step size test
	 * 5 = max iter
	 * 7 = something to do with the trust region algorithm
	 * negative = diverged
	*/
	std::cout<<reason<<std::endl;
	if (reason==2 || reason==3 || reason==7 || reason==4 || reason==5)
		converged  = true;
	else
		converged = false;

	// If I am running sensitivity with this problem ten I need to update to the current jacobian
	if (_prob->sensitivity())
	{
		Assembler* assem = _prob->get_assembler();
		assem->assemble_nonlinear(_K, _P_int, _delta_t,
								  update_ISVs, _prob->get_solution(), false, true);
	}

	return ierr;
}



