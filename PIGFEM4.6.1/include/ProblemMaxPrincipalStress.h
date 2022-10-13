/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _NONLINEAR_STRUCTURAL_MPS_H_
#define _NONLINEAR_STRUCTURAL_MPS_H_

#include "ProblemNonlinear.h"
#include "mpi.h" // Used for the global reduction on the damage state



class ProblemMaxPrincipalStress : public ProblemNonlinear
{
	public:
		ProblemMaxPrincipalStress();

		// This function return the number of degrees of freedom present on each node
		virtual id_type nndof();

		// Returns the problem type enum for this problem
		virtual problem_type get_type();

		// Set any possible problem specific parameters
		virtual void set_parameter(std::string name, double val);
		virtual double get_parameter(std::string name);

		// Returns the physics of the problem that is being solved
		virtual classification get_classification() {return STRUCTURAL;};

		// A function to allow the problem the possiblity of repeating a time step based upon some converged solution criteria
		virtual bool repeat_time_step();

		// Function to update the time step and return whether or not the problem is controlling the time step
		virtual bool update_time_step(double& dt, double& rel_tol, double& abs_tol);

		// Initialization function
		virtual void init();
		
	protected:

		// This function will be overridden in derived classes to set the actual assembler and solver that will be used
		virtual void generate_assembler();

		virtual void fillBmat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& grad_x);

		bool _prev_step_had_damage;
		bool _plane_strain;
		double _sigma_max;
		double _solver_dt;
		double _solver_rel_tol;
		double _solver_abs_tol;
		unsigned char _convergence_case;
		id_type _max_order_decrease;	// How many orders of magnitude to decrease the stiffness of elements by
		id_type _order_decrease_step;	// The number of orders of magnitude to decrease the stiffness by each step
		std::vector<std::vector<id_type> > _curr_order_decrease;

		double _prob_rel_tol;
		double _prob_abs_tol;
};

#endif
