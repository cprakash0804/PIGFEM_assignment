/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#ifndef _PROBLEM_H_
#define _PROBLEM_H_
#include "common.h"
#include "material.h"
#include <vector>
#include <map>
#include <string>
#include <fstream>

class Assembler;
class BoundaryObject;
class DofObject;
class InternalVars;
class Mesh;
class NodalData;
class Solver;
class Elem;
class Writer;
class BodyLoad;
class SubscaleModel;
class SensitivitySolver;
class SensitivityMaterialParameter;
class SensitivityShapeParameter;
class SensitivityFunction;
class Options;





class Problem
{
	protected:

		// All of the major parts of the problem to be solved
		Mesh* _mesh;
		BoundaryObject* _boundary;
		DofObject* _dofs;
		Assembler* _assembler;
		Solver* _solver;
		SensitivitySolver* _sensitivity_solver;

		// INFORMATION ABOUT INTERNAL VARIABLES
		//----------------------------------------------------------------------------------------------------------------------------------------------
		/*
		 * object to store the internal variables
		 */
		InternalVars* _internal_vars;


		// INFORMATION ABOUT THE LOCAL SOLUTION
		//----------------------------------------------------------------------------------------------------------------------------------------------
		/*
		 * Objects that store the solution and external load at each node
		*/	
		NodalData* _solution;
		NodalData* _external_load;

		// Object used to output this problem to a file
		Writer* _writer;
		Writer* _coh_writer;
		Writer* _sensitivity_writer;
		Writer* _coh_failure_writer;
		Writer* _path_writer;
		Writer* _coh_openings_writer;

		// Boolean to store whether or not this problem has been solved yet
		bool _solved;

		// Stores what kind of method is used to enforce hanging node constraints
		enforcement_method _method;

		// Body loads applied to the body in question
		std::vector<BodyLoad*> _body_loads;

		// INFORMATION ABOUT THE SUBSCALE MODEL BEING SOLVED IF THIS IS A SUBSCALE SIMULATION
		//----------------------------------------------------------------------------------------------------------------------------------------------

		SubscaleModel* _subscale_model;

		// Stores whether or not the problem being run is a subscale simulation
		bool _subscale;

		bool _init;

		bool _setup;

		bool _sensitivity;

		// INFORMATION ABOUT THE SENSITIVITY PARAMETERS AND FUNCTIONS IF THEY EXIST IN THIS PROBLEM
		//----------------------------------------------------------------------------------------------------------------------------------------------

		// A vector of the sensitivity parameters
		std::vector<SensitivityParameter*> _parameters;

		// A vector of the functions this problem will take the sensitivity of
		std::vector<SensitivityFunction*> _sensitivity_functions;

		// INFORMATION ABOUT TIMING OF THE CODE
		//----------------------------------------------------------------------------------------------------------------------------------------------

		double _assemble_time, _solve_time, _misc_solve_time, _solver_write_time, _sensitivity_time;
		double _IGFEM_nodal_detection_time, _IGFEM_element_detection_time, _IGFEM_enrichments_time;

		// OPTIONS DATABASE. USED FOR STORING THINGS ABOUT THE SIMULATION BEING RUN
		//----------------------------------------------------------------------------------------------------------------------------------------------
		Options * _options;

		friend class Solver;
		friend class PETScLinearKSPSolver;
		friend class SolverNonlinear;
		friend class SolverNonlinearKSP;
		friend class SolverNonlinearSNES;
		friend class SolverNonlinearExplicit;


	public:

		Problem();
		virtual ~Problem();

		void attach_mesh(Mesh* mesh);

		void attach_subscale_model(SubscaleModel* model) {_subscale_model = model;};

		void add_body_load(BodyLoad* load);

		bool& get_solved() {return _solved;};

		bool& subscale() {return _subscale;};

		bool& sensitivity() {return _sensitivity;};

		virtual bool linear() = 0; 

		// MEMBER ACCESS OPERATORS
		// -----------------------------------------------------------------------------

		Mesh* get_mesh() const {return _mesh;};
		BoundaryObject* get_boundary() const {return _boundary;};
		DofObject* get_dofs() const {return _dofs;};
		Assembler* get_assembler() const {return _assembler;};
		Solver* get_solver() const {return _solver;};
		NodalData* get_solution() const {return _solution;};
		NodalData* get_external_load() const {return _external_load;};
		InternalVars* get_internal_vars() const {return _internal_vars;};
		Writer* get_writer() const {return _writer;};
		SubscaleModel* get_subscale_model() const;
		Writer* get_cohesive_writer() const;
		SensitivitySolver* get_sensitivity_solver() const {return _sensitivity_solver;};
		Options* getOptions() const {return _options;};
		Writer* getSensitivityWriter() const {return _sensitivity_writer;};
		Writer* getCohesiveFailureWriter() const {return _coh_failure_writer;};
		Writer* getPathWriter() const {return _path_writer;};
		Writer* getCohOpeningsWriter() const {return _coh_openings_writer;};


		// ITERATORS
		// -----------------------------------------------------------------------------
		typedef std::vector<BodyLoad*>::iterator body_load_iterator;
		body_load_iterator body_loads_begin() {return _body_loads.begin();};
		body_load_iterator body_loads_end() {return _body_loads.end();};





		// PURE VIRTUAL FUNCTIONS (NEED TO BE DEFINED BY CREATER OF NEW PROBLEMS)
		// -----------------------------------------------------------------------------


		// This function return the number of degrees of freedom present on each node
		virtual id_type nndof() = 0;

		// Returns the problem type enum for this problem
		virtual problem_type get_type() = 0;

		// A function to allow the problem the possiblity of repeating a time step based upon some converged solution criteria
		virtual bool repeat_time_step() {return false;};

		// Function to update the time step and return whether or not the problem is controlling the time step
		virtual bool update_time_step(double& dt, double& rel_tol, double& abs_tol) {return false;};

		virtual classification get_classification() = 0;


		// BIG FUNCTIONS THAT WILL BE CALLED BY USER
		// -----------------------------------------------------------------------------

		// Actually calls the solver to solve the problem at hand
		virtual void solve_problem();

		// Computes the L2 Norm of the error of the solution with respect to the analytical function provided as input
		double Compute_L2_Norm(std::vector<double> (*analytical_disp)(std::vector<double>));

		void writeStatistics();

		// Change how to enforce hanging node constraints (defaults to PI_MATRIX)
		enforcement_method& get_enforcement_method() {return _method;};

		// Computes all of the strain and stress values for all of the volumetric quadrature points in ana element
		void compute_stress_strain(std::vector<std::vector<double> >& strain, std::vector<std::vector<double> >& stress,
								   Elem* el);

		void stress_strain_kernel(std::vector<double>& strain, std::vector<double>& stress,
								  const std::vector<std::vector<double> >& shape_grad,
								  Material* mat, Material::input_params& input,
								  const std::vector<double>& elem_U_curr);

		/*
		 * Computes the actual solution value at the enrichment nodes
		 */
		void interpolate_nodal_enrichment(NodalData& solution, NodalData& external_load);

		/*
		 * Function used to set a general single parameter
		 */
		virtual void set_parameter(std::string name, double val) {};
		virtual double get_parameter(std::string name)  {err_message("Attemping to call get_parameter for a problem with no parameters");};


		/*
		 * Function used to easily set a vector-valued parameter
		 */
		virtual void set_vec_parameter(std::string name, std::vector<double> val) {};
		virtual std::vector<double> get_vec_parameter(std::string name) {err_message("Attemping to call get_vec_parameter for a problem with no parameters");};


		/*
		 *Function used to easily set a matrix-valued parameter
		 */
		virtual void set_mat_parameter(std::string name, DenseMatrix<double> val) {};
		virtual DenseMatrix<double> get_mat_parameter(std::string name) {err_message("Attemping to call get_mat_parameter for a problem with no parameters");};


		/*
		 * Functions to add parameters to take the sensitivity wrt
		 */
		void addMaterialSensitivityParameter(std::string mat_name, std::string param_name);
		void addShapeSensitivityParameter(id_type inclusion_id, std::string param_name);
		void addLoadSensitivityParameter(/* Not sure what parameters to place here yet */);

		id_type n_SensitivityFunctions() {return _sensitivity_functions.size();};
		id_type n_SensitivityParameters() {return _parameters.size();};

		void addSensitivityFunction(SensitivityFunction* func);
		

		// Initialization functions
		virtual void init();
		virtual void setup();

	protected:

		virtual void fillBmat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& grad_x) {};

		// This function will be overridden in derived classes to set the actual assembler that will be used
		virtual void generate_assembler() {};

		// This function will be overridden in derived classes to set the actual solver that will be used
		virtual void generate_solver() {};

		// Check that all of the materials in the mesh are the correct classification for the given problem type
		void check_materials();
};



#endif
