#ifndef _SENSITIVITY_DIRECT_H_
#define _SENSITIVITY_DIRECT_H_
#include "SensitivitySolver.h"

class NodalData;

// A class that implements direct senstivity analysis
class SensitivityDirect : public SensitivitySolver
{
	public:
		
		SensitivityDirect();
		virtual ~SensitivityDirect();

		// Solve the direct sensitivity problem for all functions and parameters
		virtual PetscErrorCode solveProblem();

		// Initilize everything and make sure all structures are allocated
		virtual PetscErrorCode init();
		virtual PetscErrorCode setup();

	private:

		// The PETSc objects that are preallocated no matter what
		NodalData* _dU_dd;
		pvector _dU_dd_vec;

		// The PETSc objects that are preallocated if we are running in the "individual" mode
		KSP _ksp;
		Vec _Pseudo;

		// The PETSc objects that are preallocated if we are running in the "block" mode
		Mat _RHS, _dUf_dd_block;
		IS _rperm, _cperm;
		Vec _Sens, _seq_Sens;
		VecScatter _scat_sens;

		int _step;
		bool _individual;
};



#endif