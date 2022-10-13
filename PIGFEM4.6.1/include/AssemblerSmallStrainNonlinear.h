/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _SSNL_ASSEMBLER_H_
#define _SSNL_ASSEMBLER_H_
#include "Assembler.h"
#include "common.h"


class AssemblerSmallStrainNonlinearStructural : public Assembler
{
	public:
		virtual ~AssemblerSmallStrainNonlinearStructural() {};

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
};


#endif
