/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "AssemblerSmallStrainNonlinear.h"
#include "material.h"
#include "Utilities.h"
#include <algorithm>




bool AssemblerSmallStrainNonlinearStructural::storedBmats()
{
	return true;
}


// Kernel function that actually does the computes the physics behind the specific problem
// Computes the local contribution to the stiffness matrix and internal load vector
// Also updates the internal variable storage in the input variable
//	B-matrix variant
void AssemblerSmallStrainNonlinearStructural::KernelNonlinear(DenseMatrix<double>& K_el, std::vector<double>& P_el_int,
													const std::vector<double>& shape, const DenseMatrix<double>& B,
													Material* mat, Material::input_params& input,
													std::vector<double>& elem_U_curr, bool assembleFunc, bool assembleJac)
{
	// Compute the current strain
	ProblemUtilities::SmallStrainFast(input.strain, B, elem_U_curr);

	// Compute the constitutive relations
	Material::output_params* output = mat->Constitutive(input);

	// Compute the stiffness matrix
	if (assembleJac)
		ProblemUtilities::SmallStrainKFast(K_el, B, output->Dmat);

	// Compute the internal force
	if (assembleFunc)
		ProblemUtilities::SmallStrainInternalForceFast(P_el_int, B, output->stress);
}


void AssemblerSmallStrainNonlinearStructural::KernelNonlinear(DenseMatrix<double>& K_el, std::vector<double>& P_el_int,
													const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
													Material* mat, Material::input_params& input,
													std::vector<double>& elem_U_curr, bool assembleFunc, bool assembleJac)
{
	// Compute the current strain
	ProblemUtilities::SmallStrainFast(input.strain, shape_grad, elem_U_curr);

	// Compute the constitutive relations
	Material::output_params* output = mat->Constitutive(input);

	// Compute the stiffness matrix
	if (assembleJac)
		ProblemUtilities::SmallStrainKFast(K_el, shape_grad, output->Dmat);
	
	// Compute the internal force
	if (assembleFunc)
		ProblemUtilities::SmallStrainInternalForceFast(P_el_int, shape_grad, output->stress);
}





void AssemblerSmallStrainNonlinearStructural::KernelCohesive(DenseMatrix<double>& K_coh, std::vector<double>& P_int_coh,
															 const std::vector<double>& shape, const DenseMatrix<double>& rotation_matrix,
															 Material* mat, Material::input_params& input,
															 std::vector<double>& coh_U_curr, bool assembleFunc, bool assembleJac)
{
	// Compute the current opening vector in the ttn coordinate system
	ProblemUtilities::cohesiveDeltaFast(input.delta, shape, coh_U_curr, rotation_matrix);

	// Compute the constitutive relations
	Material::output_params* output = mat->Constitutive(input);

	// Compute the stiffness matrix
	if (assembleJac)
		ProblemUtilities::cohesiveKFast(K_coh, shape, output->Dmat, rotation_matrix);

	// Compute the internal force
	if (assembleFunc)
		ProblemUtilities::cohesiveInternalForceFast(P_int_coh, shape, output->traction, rotation_matrix);	
}



// Select the appropriate B-matrix fill function
void AssemblerSmallStrainNonlinearStructural::fillBmat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& shape_grad)
{
	ProblemUtilities::Assemble_Small_Strain_B_Mat(B, shape_grad);
}
