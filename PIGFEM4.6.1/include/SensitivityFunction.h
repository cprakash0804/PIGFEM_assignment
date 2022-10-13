#ifndef _SENSITIVITY_FUNCTION_H_
#define _SENSITIVITY_FUNCTION_H_
#include "petscmat.h"
#include "petscvec.h"
#include "common.h"

class Problem;
class SensitivityParameter;
class SensitivitySolver;


/*
 * A class that defines a general structure for all functions that
 * I would I like to tae the sensitivity of
 */
class SensitivityFunction
{
	protected:

		// The problem we're gonna be dealing with
		Problem* _prob;

	public:
		virtual ~SensitivityFunction();

		virtual void init() {};
		virtual void setup() {};

		SensitivityFunction* allocate_and_copy();
		virtual SensitivityFunction* allocate() = 0;
		virtual void copy(SensitivityFunction* new_func) {/* default to doing nothing */};

		virtual void attachProblem(Problem* prob) {_prob = prob;};

		virtual double evaluate(pmatrix* K, pvector* U, pvector* F_ext) = 0;

		virtual void assemble_dfdUf(Vec* dfdUf,
									pmatrix* K, pvector* U, pvector* F_ext) = 0;

		virtual void assemble_dfdUp(Vec* dfdUp,
									pmatrix* K, pvector* U, pvector* F_ext) = 0;

		virtual double get_parameter_partial(SensitivityParameter* param, SensitivitySolver* solver,
											 pmatrix* K, pvector* U, pvector* F_ext) = 0;
};



#endif