/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#ifndef _SOLVER_H_
#define _SOLVER_H_
#include "petscksp.h"
#include "common.h"
#include <set>
#include <vector>

// Predeclaration of the problem class (not really sure why I have to do this here but oh wel
class Problem;
class NodalData;


struct solver_state
{
	Vec Uf;
	Vec Up;
	Vec Pf_ext;
	Vec Pp_ext;
	Vec Pf_int;
	Vec Pp_int;
	std::vector<std::vector<std::vector<double> > > internal_vars;
};

class Solver
{
	protected:
		// All of the PETSc matricies and vectors necessary to solve the system
		// They will be assembled by the assembler
		pmatrix _K;
		pvector _U;
		pvector _F_ext;
		Vec _Up_initial; // Stores the solution at the beginning of the problem. Gets added to the Dirichlet BC of the current run

		// Problem object and assembler object
		Problem* _prob;

		// Index sets for toe store solution function
		IS _global_f;
		IS _to_f;
		IS _global_p;
		IS _to_p;
		id_type _r_start; // The first row of the recieve vector that I own

		// Solver settings (Most of these are only used in the nonlinear solver)
		double _rel_tol, _abs_tol, _T_final, _max_time_step, _min_time_step, _delta_t;
		id_type _max_iter;
		bool _setup; // Stores whether or not this problem has been set up already
		bool _nonzero_start;


		friend class SensitivitySolver;
		friend class SensitivityAdjoint;
		friend class SensitivityDirect;


	public:

		Solver();
		virtual ~Solver();
		virtual PetscErrorCode setup();

		void attach_problem(Problem* problem) {_prob = problem;};

		virtual PetscErrorCode solve() = 0;
		PetscErrorCode store_solution(pvector& U, NodalData* data);
		PetscErrorCode store_solution(Vec& Uf, Vec& Up, NodalData* data);

		PetscErrorCode initialSolution(NodalData& starting_solution);

		/*
		 * Function used to set a general single parameter
		 */
		virtual void set_dparameter(std::string name, double val);
		virtual double get_dparameter(std::string name);
		virtual void set_iparameter(std::string name, id_type val);
		virtual id_type get_iparameter(std::string name);

		virtual double getDeltaT() {return 0.0;};
		bool getNonzeroStart() const {return _nonzero_start;};
//		std::vector<int> penetration_tags;

	protected:
		// Helper function for the preallocate matrices function
		void populate_local_neighbors(std::vector<std::set<id_type> >& dof_neighbors);

		//  Helper functions for dealing with the neighboring dofs for periodic bcs
		void combine_periodic_contributions(std::vector<std::set<id_type> >& dof_neighbors);

		// Preallocation functions
		virtual PetscErrorCode preallocate_vectors() = 0;
		PetscErrorCode preallocate_matrix(pmatrix& K);
		PetscErrorCode preallocate_index_sets();

		virtual PetscErrorCode initKSP(KSP& ksp, bool create);

		void Output(id_type step_iter, double current_t);
};



#endif
