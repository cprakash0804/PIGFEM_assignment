#ifndef _SENSITIVITY_ADJOINT_H_
#define _SENSITIVITY_ADJOINT_H_
#include "SensitivitySolver.h"

// A class that implements direct senstivity analysis
class SensitivityAdjoint : public SensitivitySolver
{
	public:
		
		SensitivityAdjoint();
		virtual ~SensitivityAdjoint();

		// Solve the adjoint sensitivity problem for all functions and parameters
		virtual PetscErrorCode solveProblem();

		// Initilize everything and make sure all structures are allocated
		virtual PetscErrorCode setup();

	private:

		// The PETSc objects that are preallocated no matter what
		Vec _Lambda;
		KSP _ksp;
};



#endif