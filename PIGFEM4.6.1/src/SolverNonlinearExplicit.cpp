/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated April 2017

##################################################################################
*/
#include "SolverNonlinearExplicit.h"
#include "Problem.h"
#include "Assembler.h"
#include "DofObject.h"
#include "NodalData.h"
#include "InternalVars.h"
#include "Writer.h"
#include "SubscaleModel.h"
#include "SensitivitySolver.h"
#include "elem.h"
#include "Options.h"
#include "Utilities.h"
#include <sstream>

SolverNonlinearExplicit::SolverNonlinearExplicit()
{
}

SolverNonlinearExplicit::~SolverNonlinearExplicit()
{
	KSPDestroy(&_ksp);
}




PetscErrorCode SolverNonlinearExplicit::setup()
{
	PetscErrorCode ierr;

	// Clear any old structures
	if (_setup)
	{
		KSPDestroy(&_ksp);
	}

	// Create the basic matricies and vectors
	ierr = SolverNonlinear::setup();CHKERRQ(ierr);

	// Set up all specific structures here
	ierr = initKSP(_ksp, true);CHKERRQ(ierr);
	return ierr;
}






PetscErrorCode SolverNonlinearExplicit::solve()
{
	// Get pointers to the main objects
	InternalVars* int_vars = _prob->get_internal_vars();

	// PETSc initializations
	PetscErrorCode ierr;
	ierr = VecZeroEntries(_U.f);CHKERRQ(ierr);
	ierr = VecZeroEntries(_U.p);CHKERRQ(ierr);
	ierr = VecZeroEntries(_F_ext.f);CHKERRQ(ierr);

	std::vector<std::vector<std::vector<double> > > update_ISVs;
	int_vars->preallocate_internal_vars_object(update_ISVs);

	id_type n_iterations = _T_final/_min_time_step;
	double current_t = 0.0;
	PIGFEMPrint("UTILIZING THE FULLY EXPLICIT SOLVER. ESTIMATED NUMBER OF TIME STEPS: " << n_iterations << std::endl);
	bool cont;
	explicit_solve(n_iterations, current_t, update_ISVs, _ksp, cont);

	return ierr;
}
