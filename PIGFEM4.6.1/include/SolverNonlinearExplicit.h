/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated July 2017

##################################################################################
*/
#ifndef _NLIN_EXPLICIT_SOLVER_H_
#define _NLIN_EXPLICIT_SOLVER_H_
#include "SolverNonlinear.h"
#include "petscksp.h"
#include "material.h"


class SolverNonlinearExplicit : public SolverNonlinear
{
	protected:
		KSP _ksp;

	public:
		SolverNonlinearExplicit();
		virtual ~SolverNonlinearExplicit();

		virtual PetscErrorCode setup();
	
		virtual PetscErrorCode solve();

	protected:

		virtual PetscErrorCode nonlinearStep() {err_message("No nonlinear steps in a fully explicit solve!");};
};


#endif