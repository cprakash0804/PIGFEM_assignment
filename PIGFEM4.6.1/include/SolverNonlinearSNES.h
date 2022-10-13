/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated July 2017

##################################################################################
*/
#ifndef _NLIN_SNES_SOLVER_H_
#define _NLIN_SNES_SOLVER_H_
#include "SolverNonlinear.h"
#include "petscsnes.h"
#include "material.h"


class SolverNonlinearSNES : public SolverNonlinear
{
	protected:
		SNES _snes;
		double _misc_time;
		double _assem_time;

		friend class SensitivityDirect;
		friend class SensitivityAdjoint;
		
	public:
		SolverNonlinearSNES();
		virtual ~SolverNonlinearSNES();

		virtual PetscErrorCode setup();

	protected:

		virtual PetscErrorCode nonlinearStep(std::vector<std::vector<std::vector<double> > >& update_ISVs, bool& converged, double& Res_norm);
		friend PetscErrorCode formResidual(SNES snes,Vec x,Vec f,void *ctx);
		friend PetscErrorCode formJacobian(SNES snes,Vec x,Mat Amat,Mat Pmat,void *ctx);
		friend PetscErrorCode monitorFunction(SNES snes,PetscInt its, PetscReal norm,void *mctx);
};


#endif
