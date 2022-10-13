/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _LINEAR_THERMAL_H_
#define _LINEAR_THERMAL_H_

#include "ProblemLinear.h"

class ProblemLinearThermal : public ProblemLinear
{
	public:
		ProblemLinearThermal();

		// This function return the number of degrees of freedom present on each node
		virtual id_type nndof();

		// Returns the problem type enum for this problem
		virtual problem_type get_type();

		// Returns the physics of the problem that is being solved
		virtual classification get_classification() {return THERMAL;};

	protected:

		// This function will be overridden in derived classes to set the actual assembler and solver that will be used
		virtual void generate_assembler();

		virtual void fillBmat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& grad_x);
};


#endif
