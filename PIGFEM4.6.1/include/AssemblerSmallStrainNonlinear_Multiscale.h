/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _SSNLMS_ASSEMBLER_H_
#define _SSNLMS_ASSEMBLER_H_
#include "Assembler.h"
#include "common.h"


class AssemblerSmallStrainNonlinear_Multiscale : public Assembler
{
	public:
		AssemblerSmallStrainNonlinear_Multiscale();
		virtual ~AssemblerSmallStrainNonlinear_Multiscale();

		virtual void clear();


		// Function to assemble the presribed displacements and free loads at the beginning of each load step
		// (Or at the beginning of the analysis for linear)
		// This is made virtual so that calling it can update the assembler state at every load step is desired
		virtual PetscErrorCode assemble_new_load_step(Vec& Up, pvector& F_ext, double current_time_frac);

		/*
		 * Function used to easily set a vector-valued parameter
		 */
		virtual void set_vec_parameter(std::string name, std::vector<double> val);
		virtual std::vector<double> get_vec_parameter(std::string name);

		virtual bool storedBmats();

	protected:

		// Kernel function that actually does the computes the physics behind the specific problem
		// Computes the local contribution to the stiffness matrix and internal load vector
		// Also updates the internal variable storage in the input variable
		virtual void KernelNonlinear(DenseMatrix<double>& K_el, std::vector<double>& P_el_int,
									 const std::vector<double>& shape, const DenseMatrix<double>& B,
									 Material* mat, Material::input_params& input,
									 std::vector<double>& elem_U_curr, bool assembleFunc, bool assembleJac);
		virtual void KernelNonlinear(DenseMatrix<double>& K_el, std::vector<double>& P_el_int,
									 const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
									 Material* mat, Material::input_params& input,
									 std::vector<double>& elem_U_curr, bool assembleFunc, bool assembleJac);
		virtual void KernelCohesive(DenseMatrix<double>& K_coh, std::vector<double>& P_int_coh,
									 const std::vector<double>& shape, const DenseMatrix<double>& rotation_matrix,
									 Material* mat, Material::input_params& input,
									 std::vector<double>& coh_U_curr, bool assembleFunc, bool assembleJac);

		virtual void fillBmat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& shape_grad);
	private:

		// The maximum macroscopic strain used in this problem (Should really be in problem class but the assembler needs access. Make them friends?)
		std::vector<double> _max_macro_strain;

		// The current time fraction of the simulation
		std::vector<double> _curr_macro_strain;

		// quick boolean check
		bool _strain_set;
};


#endif
