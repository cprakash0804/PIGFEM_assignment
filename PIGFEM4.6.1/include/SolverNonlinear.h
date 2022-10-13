/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated July 2017

##################################################################################
*/
#ifndef _NLIN_SOLVER_H_
#define _NLIN_SOLVER_H_
#include "Solver.h"
#include "petscksp.h"
#include "material.h"

class InternalVars;

class SolverNonlinear : public Solver
{
	protected:
		pvector _P_int;
		Vec _Res_f;
	public:
		SolverNonlinear();
		virtual ~SolverNonlinear();
		virtual PetscErrorCode solve();
	protected:
		virtual PetscErrorCode preallocate_vectors();

		void update_time_step(double& delta_t, double& curr_t, double final_t,
							  double min_dt, double max_dt, bool converged,
							  bool exceeds_rates, bool problem_repeat,
							  bool& problem_controlled_dt, double& rel_tol, double& abs_tol,
							  unsigned int curr_time_step_repeat);

		PetscErrorCode explicit_solve(id_type n_iterations, double& current_time, std::vector<std::vector<std::vector<double> > >& update_ISVs, KSP& ksp, bool& cont);

		virtual PetscErrorCode nonlinearStep(std::vector<std::vector<std::vector<double> > >& update_ISVs, bool& converged, double& Res_norm) = 0;



		virtual PetscErrorCode setup();



		void compute_int_var_avgs(InternalVars* vars,
								  std::map<std::string, std::vector<double> >& avg,
								  std::map<std::string, std::vector<double> >& maxes);

		bool check_rates(std::map<std::string, std::vector<double> >& avgs,
						 std::map<std::string, std::vector<double> >& avg_prev,
						 std::map<std::string, std::vector<double> >& maxes,
						 std::map<std::string, std::vector<double> >& rate_limits,
						 std::map<std::string, std::vector<double> >& check_limit);

		// Stores the global number of quadrature points associated with
		// each material in the mesh if the material has internal variables
		std::map<std::string, id_type> _nqp_global;
};



#endif
