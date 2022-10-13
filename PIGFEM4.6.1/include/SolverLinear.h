/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _LIN_SOLVER_H_
#define _LIN_SOLVER_H_
#include "Solver.h"
#include "petscksp.h"
#include "Assembler.h"

class PETScLinearKSPSolver : public Solver
{
	public:
		PETScLinearKSPSolver();
		virtual ~PETScLinearKSPSolver();
	
		virtual PetscErrorCode solve();

	protected:
		virtual PetscErrorCode preallocate_vectors();
};


#endif