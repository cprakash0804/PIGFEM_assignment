/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated July 2017

##################################################################################
*/
#ifndef _NLIN_KSP_SOLVER_H_
#define _NLIN_KSP_SOLVER_H_
#include "SolverNonlinear.h"
#include "petscksp.h"
#include "material.h"


class SolverNonlinearKSP : public SolverNonlinear
{
	protected:
		KSP _ksp;
		Vec _delta_Uf;

		friend class SensitivityDirect;
		friend class SensitivityAdjoint;

	public:
		SolverNonlinearKSP();
		virtual ~SolverNonlinearKSP();

		virtual PetscErrorCode setup();

	protected:

		virtual PetscErrorCode nonlinearStep(std::vector<std::vector<std::vector<double> > >& update_ISVs, bool& converged, double& Res_norm);
};


#endif
