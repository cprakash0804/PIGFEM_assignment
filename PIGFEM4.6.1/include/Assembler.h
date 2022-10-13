/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _ASSEMBLER_H_
#define _ASSEMBLER_H_
#include "petscmat.h"
#include "petscvec.h"
#include "DenseMatrix.h"
#include "common.h"
#include "Utilities.h"
#include "material.h"
#include <vector>

// Forward declarations
class Problem;
struct solver_state;
class BodyLoad;
class Elem;

class Assembler
{
	protected:
		Problem* _prob;
		std::vector<std::vector<DenseMatrix<double> > > _Bmats;

	public:

		virtual ~Assembler();

		virtual void clear();

		void attach_problem(Problem* problem) {_prob = problem;};

		// General function for global assembly of a linear problem
		PetscErrorCode assemble_linear(pmatrix& K);

		// General function for global assembly of a nonlinear problem
		PetscErrorCode assemble_nonlinear(pmatrix& K, pvector& P_int,
										  double delta_t, std::vector<std::vector<std::vector<double> > >& init_internal_vars,
										  NodalData* solution, bool assembleFunc = true, bool assembleJac = true);

		// Function to assemble the presribed displacements and free loads at the beginning of each load step
		// (Or at the beginning of the analysis for linear)
		// This is made virtual so that calling it can update the assembler state at every load step is desired
		virtual PetscErrorCode assemble_new_load_step(Vec& Up, pvector& F_ext, double current_time_frac);


		/*
		 * Function used to set a general single parameter
		 */
		virtual void set_parameter(std::string name, double val) {};
		virtual double get_parameter(std::string name) {err_message("Attemping to call get_parameter for an assembler with no parameters");};


		/*
		 * Function used to easily set a vector-valued parameter
		 */
		virtual void set_vec_parameter(std::string name, std::vector<double> val) {};
		virtual std::vector<double> get_vec_parameter(std::string name) {err_message("Attemping to call get_vec_parameter for an assembler with no parameters");};


		/*
		 *Function used to easily set a matrix-valued parameter
		 */
		virtual void set_mat_parameter(std::string name, DenseMatrix<double> val) {};
		virtual DenseMatrix<double> get_mat_parameter(std::string name) {err_message("Attemping to call get_mat_parameter for an assembler with no parameters");};

		// Functions to do with storing the B-matrices for certain classes of problems
		virtual bool storedBmats() = 0;
		void storeBmats();

		friend class SensitivityParameter;
		friend class SensitivityShapeParameter;
		friend class SensitivityMaterialParameter;

	protected:

		// Assembles the elemental stiffness matrix for the element pointed to by el
		void assemble_elem_linear(Elem* el,
								  DenseMatrix<double>& K_el);

		// Assembles the elemental stiffness matrix and internal load vector for the element pointed to by el
		void assemble_elem_nonlinear(Elem* el,
									 DenseMatrix<double>& K_el, std::vector<double>& P_el_int,
									 double delta_t, std::vector<double>& elem_U_curr, std::vector<std::vector<double> >& init_elem_internal_vars,
									 bool assembleFunc, bool assembleJac);

		// Assembles the external force from body loads
		void assemble_elem_body_load(Elem* el, BodyLoad* load, std::vector<double>& P_el_ext);

		// Kernel function that actually does the computes the physics behind the specific problem
		// Computes the local contribution to the stiffness matrix and internal load vector
		// Also updates the internal variable storage in the input variable

		// B-matrix kernel variants
		virtual void KernelLinear(DenseMatrix<double>& K_el,
								  const std::vector<double>& shape, const DenseMatrix<double>& B,
								  Material* mat, Material::input_params& input) {};
		virtual void KernelNonlinear(DenseMatrix<double>& K_el, std::vector<double>& P_el_int,
									 const std::vector<double>& shape, const DenseMatrix<double>& B,
									 Material* mat, Material::input_params& input,
									 std::vector<double>& elem_U_curr, bool assembleFunc, bool assembleJac) {};
		virtual void KernelCohesive(DenseMatrix<double>& K_coh, std::vector<double>& P_int_coh,
									 const DenseMatrix<double>& N, const DenseMatrix<double>& rotation_matrix,
									 Material* mat, Material::input_params& input,
									 std::vector<double>& elem_U_curr, bool assembleFunc, bool assembleJac) {};
		
		// Shape function gradient kernel variants
		virtual void KernelLinear(DenseMatrix<double>& K_el,
								  const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
								  Material* mat, Material::input_params& input) {};
		virtual void KernelNonlinear(DenseMatrix<double>& K_el, std::vector<double>& P_el_int,
									 const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
									 Material* mat, Material::input_params& input,
									 std::vector<double>& elem_U_curr, bool assembleFunc, bool assembleJac) {};
		virtual void KernelCohesive(DenseMatrix<double>& K_coh, std::vector<double>& P_int_coh,
									 const std::vector<double>& shape, const DenseMatrix<double>& rotation_matrix,
									 Material* mat, Material::input_params& input,
									 std::vector<double>& elem_U_curr, bool assembleFunc, bool assembleJac) {};

		// Assembles the matrix used for the enforcement of hanging node constraints based on the penalty method
		// Then apply_penalty_constraint actually applies the constraint
		void assemble_penalty_matrix(DenseMatrix<double>& P, std::vector<double> constraint_vals);
		void apply_penalty_constraint(pmatrix& K, double max);

		// Sets the values for a matrix
		PetscErrorCode set_matrix(pmatrix& K, DenseMatrix<double>& mat, std::vector<id_type>& dofs);
		// Sets the values for a vector
		PetscErrorCode set_vec(pvector& P, std::vector<double>& vec, std::vector<id_type>& dofs);

		// Functions to do with storing the B-matrices for certain classes of problems
		virtual void fillBmat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& grad_x) {};
		DenseMatrix<double>& getBmat(id_type local_e, id_type qp);

};



#endif
