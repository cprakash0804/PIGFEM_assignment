/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _NONLINEAR_STRUCTURAL_MULTISCALE_H_
#define _NONLINEAR_STRUCTURAL_MULTISCALE_H_

#include "ProblemNonlinear.h"

class ProblemNonlinearStructural_Multiscale : public ProblemNonlinear
{
	public:
		ProblemNonlinearStructural_Multiscale();

		// This function return the number of degrees of freedom present on each node
		virtual id_type nndof();

		// Returns the problem type enum for this problem
		virtual problem_type get_type();

		/*
		 * Function used to set a general single parameter
		 */
		virtual void set_parameter(std::string name, double val);
		virtual double get_parameter(std::string name);

		/*
		 * Function used to easily set a vector-valued parameter
		 */
		virtual void set_vec_parameter(std::string name, std::vector<double> val);
		virtual std::vector<double> get_vec_parameter(std::string name);


		// Returns the physics of the problem that is being solved
		virtual classification get_classification() {return STRUCTURAL;};
		
	protected:

		// This function will be overridden in derived classes to set the actual assembler and solver that will be used
		virtual void generate_assembler();

		virtual void fillBmat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& grad_x);

		bool _plane_strain;

		std::vector<double> _max_macro_strain;
		bool _strain_set;
};

#endif