#ifndef _SENSITVITY_SOLVER_H_
#define _SENSITVITY_SOLVER_H_
#include "petsc.h"
#include "petscmat.h"
#include "petscvec.h"
#include "petscksp.h"
#include "common.h"
#include <vector>



// Predeclarations
class Problem;
class SensitivityParameter;
class SensitivityFunction;
class InternalVars;



/*
 * A class that is used to solve a general sensitivity problem based on the
 * current state of the primal Problem stored in _prob
 * The sensitivity of various objective functions are found with respect to
 * several parameters
 */
class SensitivitySolver
{
	protected:

		// A pointer to the primal Problem
		Problem* _prob;

		// PETSc objects that describe the current state of the primal problem
		pmatrix* _K;
		pvector* _U;
		pvector* _F_ext;

		// The functions that we are finding the sensitivity of
		std::vector<SensitivityFunction*> _functions;

		// The parameters that we are taking the sensitivity with respect to
		std::vector<SensitivityParameter*> _parameters;

		// Vectors of the PETSc vectors associated with the current partial derivatives
		std::vector<Vec> _dUp_dd;			// Total derivative of the applied displacment wrt the parameters (Only filled for load parameters)
		std::vector<Vec> _dFf_dd;			// Total derivative of the applied external force wrt the parmeters (Only filled for load parameters)
		std::vector<pvector> _delP_deld;	// Partial derivative of the internal force wrt the parameters
		std::vector<pvector> _delf_delU;	// Partial derivative of the function wrt the solution vector
		std::vector<std::vector<double> > _delf_deld;	// Partial derivative of the objective functions wrt the parameters

		// Store the current sensitivities here
		std::vector<double> _function_vals;
		std::vector<std::vector<double> > _sensitivities;

		// Function to preallocate the appropriate vectors
		PetscErrorCode preallocateVector(Vec* vec, bool free_portion);

		// Timing variables
		double _sensitivity_solve_time, _sensitivity_substitution_time, _sensitivity_assemble_time;

		// Derivatives of the internal variables
		InternalVars* _internal_vars;

		// Check Whether or not this solver has been set up before
		bool _setup;

	public:

		// The Big 3
		SensitivitySolver();
		virtual ~SensitivitySolver();

		// Initilize everything and make sure all structures are allocated
		virtual PetscErrorCode init();
		virtual PetscErrorCode setup();

		// Attach a primal problem to the 
		void attachProblem(Problem* prob) {_prob = prob;};

		// Add a parameter to take sensitivity wrt
		void addParameter(SensitivityParameter* param);

		// Add a function to find the sensitivity of
		void addFunction(SensitivityFunction* func);

		// Function to call all of the appropriate assemblers to assemble the current derivatives
		PetscErrorCode assembleDerivatives();

		double get_function_val(id_type func);
		double get_sensitivity(id_type func, id_type param);

		id_type n_parameters() {return _parameters.size();};
		id_type n_functions() {return _functions.size();};

		double get_solve_time() {return _sensitivity_solve_time;};
		double get_assembly_time() {return _sensitivity_assemble_time;};
		double get_substitution_time() {return _sensitivity_substitution_time;};
		InternalVars* get_internal_vars() {return _internal_vars;};


		Vec* get_dPf(id_type param) {return &_delP_deld[param].f;};
		Vec* get_dPp(id_type param) {return &_delP_deld[param].p;};
		Vec* get_dUp(id_type param) {return &_dUp_dd[param];};
		Vec* get_dFf(id_type param) {return &_dFf_dd[param];};
		Vec* get_dfdUf(id_type func) {return &_delf_delU[func].f;};
		Vec* get_dfdUp(id_type func) {return &_delf_delU[func].p;};
		double get_dfdd(id_type func, id_type param) {return _delf_deld[func][param];};

		// Solve either the direct or adjoint sensitivity problem
		virtual PetscErrorCode solveProblem() = 0;

};





#endif