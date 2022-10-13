#ifndef _RIGHT_LOAD_H_
#define _RIGHT_LOAD_H_
#include "SensitivityFunction.h"



/*
 * A function that give the load along the right side of a domain
 */
class SensitivityRightLoad : public SensitivityFunction
{
	private:

		// A Vector to pick out the appropriate dofs from the prescribed load vector
		Vec _L;
		bool _setup;
		double _area; // the area of the right hand side of the mesh (To convert force to stress)

		PetscErrorCode assembleLVector(Problem* prob);

	public:

		SensitivityRightLoad();
		virtual ~SensitivityRightLoad();

		virtual void setup();

		virtual SensitivityFunction* allocate() {return new SensitivityRightLoad;};

		virtual double evaluate(pmatrix* K, pvector* U, pvector* F_ext);

		virtual void assemble_dfdUf(Vec* dfdUf,
									pmatrix* K, pvector* U, pvector* F_ext);

		virtual void assemble_dfdUp(Vec* dfdUp,
									pmatrix* K, pvector* U, pvector* F_ext);

		virtual double get_parameter_partial(SensitivityParameter* param, SensitivitySolver* solver,
											 pmatrix* K, pvector* U, pvector* F_ext);
};


#endif