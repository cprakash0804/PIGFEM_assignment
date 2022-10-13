/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _NONLINEAR_PROBLEM_H_
#define _NONLINEAR_PROBLEM_H_
#include "Problem.h"


class ProblemNonlinear : public Problem
{
	public:

		/*
		 * Constructor, doesn't do anything
		 */
		ProblemNonlinear();

		virtual bool linear() {return false;};

	protected:

		/*
		 * Set the solver to the nonlinear solver
		 */
		virtual void generate_solver();
};


#endif
