/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _LINEAR_PROBLEM_H_
#define _LINEAR_PROBLEM_H_
#include "Problem.h"


class ProblemLinear : public Problem
{
	public:

		/*
		 * Constructor, doesn't do anything
		 */
		ProblemLinear();

		/*
		 * Initializer. Makes sure all materials are linear
		 * Then calls basic Problem initializer
		 */
		virtual void init();

		virtual bool linear() {return true;};

	protected:

		/*
		 * Set the solver to the nonlinear solver
		 */
		virtual void generate_solver();
};


#endif
