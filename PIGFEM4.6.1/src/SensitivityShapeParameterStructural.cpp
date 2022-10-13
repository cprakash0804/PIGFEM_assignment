/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#include "SensitivityShapeParameterStructural.h"
#include "Utilities.h"
#include <algorithm>


// Kernel function that actually does the computes the physics behind the specific problem
// Computes the local contribution to the stiffness matrix and internal load vector
// Also updates the internal variable storage in the input variable
void SensitivityShapeParameterStructural::KernelVolumetric(std::vector<double>& dP_el,
											const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
											const std::vector<double>& dshape_dd, const std::vector<std::vector<double> >& dshape_grad_dd, double div_v,
											Material* mat, Material::input_params& input, Material::sensitivity_input_params& sens_input,
											const std::vector<double>& elem_U_curr, bool compute_intersected_terms)
{
	// Compute the current strain
	std::vector<double>& strain = sens_input.strain;
	ProblemUtilities::SmallStrainFast(strain, shape_grad, elem_U_curr);

	// This function call is necessary for either functionally graded materials or for materials with internal-state-variables
	Material::sensitivity_output_params* sens_output = mat->SensitivityConstitutive(sens_input);
	std::vector<double>& dstress = sens_output->dstress_dd;
	if (dstress.size() != strain.size()) // Called for a material that doesn't implement sensitivity
	{
		dstress.clear();
		dstress.resize(strain.size());
	}

	// Have to compute all of the shape function sensitivity terms
	if (compute_intersected_terms)
	{
		// Already have the strain, just have to copy it to my input object
		input.strain = strain;

		// Compute the constitutive relations for this material
		Material::output_params* output = mat->Constitutive(input);   // Computes the D-matrix for the current material and updates any internal variables

		std::vector<double> sig_div_v = Utilities::scale(output->stress, div_v);
		std::vector<double> dB_U;
		ProblemUtilities::SmallStrainFast(dB_U, dshape_grad_dd, elem_U_curr);
		std::vector<double> C_dB_U = (output->Dmat) * dB_U;
		std::vector<double> temp1 = Utilities::plus(dstress, C_dB_U);
		std::vector<double> temp2 = Utilities::plus(temp1, sig_div_v);
		std::vector<double> temp3;
		ProblemUtilities::SmallStrainInternalForceFast(temp3, shape_grad, temp2); // computes temp3 = B^T * temp2
		std::vector<double> temp4;
		ProblemUtilities::SmallStrainInternalForceFast(temp4, dshape_grad_dd, output->stress); // computes temp4 = B_star^T * sigma
		dP_el = Utilities::plus(temp3, temp4);
	}

	// This element only has ISV contributions
	else
		ProblemUtilities::SmallStrainInternalForceFast(dP_el, shape_grad, dstress);	
}

// B-matrix variant
void SensitivityShapeParameterStructural::KernelVolumetric(std::vector<double>& dP_el,
											const std::vector<double>& shape, const DenseMatrix<double>& B,
											const std::vector<double>& dshape_dd, const DenseMatrix<double>& dB_dd, double div_v,
											Material* mat, Material::input_params& input, Material::sensitivity_input_params& sens_input,
											const std::vector<double>& elem_U_curr, bool compute_intersected_terms)
{
	// Compute the current strain
	std::vector<double>& strain = sens_input.strain;
	ProblemUtilities::SmallStrainFast(strain, B, elem_U_curr);

	// This function call is necessary for either functionally graded materials or for materials with internal-state-variables
	Material::sensitivity_output_params* sens_output = mat->SensitivityConstitutive(sens_input);
	std::vector<double>& dstress = sens_output->dstress_dd;
	if (dstress.size() != strain.size()) // Called for a material that doesn't implement sensitivity
	{
		dstress.clear();
		dstress.resize(strain.size());
	}

	// Have to compute all of the shape function sensitivity terms
	if (compute_intersected_terms)
	{
		// Already have the strain, just have to copy it to my input object
		input.strain = strain;

		// Compute the constitutive relations for this material
		Material::output_params* output = mat->Constitutive(input);   // Computes the D-matrix for the current material and updates any internal variables

		std::vector<double> sig_div_v = Utilities::scale(output->stress, div_v);
		std::vector<double> dB_U;
		ProblemUtilities::SmallStrainFast(dB_U, dB_dd, elem_U_curr);
		std::vector<double> C_dB_U = (output->Dmat) * dB_U;
		std::vector<double> temp1 = Utilities::plus(dstress, C_dB_U);
		std::vector<double> temp2 = Utilities::plus(temp1, sig_div_v);
		std::vector<double> temp3;
		ProblemUtilities::SmallStrainInternalForceFast(temp3, B, temp2); // computes temp3 = B^T * temp2
		std::vector<double> temp4;
		ProblemUtilities::SmallStrainInternalForceFast(temp4, dB_dd, output->stress); // computes temp4 = B_star^T * sigma
		dP_el = Utilities::plus(temp3, temp4);
	}

	// This element only has ISV contributions
	else
		ProblemUtilities::SmallStrainInternalForceFast(dP_el, B, dstress);	
}


void SensitivityShapeParameterStructural::KernelVolumetricUpdate(const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
																 Material* mat, Material::sensitivity_input_params& input,
									 							 const std::vector<double>& elem_U_curr, std::vector<double>& delem_U_curr)
{
	// Compute the current strain
	ProblemUtilities::SmallStrainFast(input.strain, shape_grad, elem_U_curr);

	// Compute the current strain derivative
	ProblemUtilities::SmallStrainFast(input.dstrain_dd, shape_grad, delem_U_curr);

	
	mat->updateSensitivityISVs(input);
}

// B-matrix variant
void SensitivityShapeParameterStructural::KernelVolumetricUpdate(const std::vector<double>& shape, const DenseMatrix<double>& B,
																 Material* mat, Material::sensitivity_input_params& input,
									 							 const std::vector<double>& elem_U_curr, std::vector<double>& delem_U_curr)
{
	// Compute the current strain
	ProblemUtilities::SmallStrainFast(input.strain, B, elem_U_curr);

	// Compute the current strain derivative
	ProblemUtilities::SmallStrainFast(input.dstrain_dd, B, delem_U_curr);

	
	mat->updateSensitivityISVs(input);
}















void SensitivityShapeParameterStructural::KernelCohesive(std::vector<double>& dP_el_coh,
											const std::vector<double>& shape, const DenseMatrix<double>& rotation_matrix,
											const std::vector<double>& dshape_dd, const DenseMatrix<double>& drotation_matrix_dd, double div_v,
											Material* mat, Material::input_params& input, Material::sensitivity_input_params& sens_input,
											const std::vector<double>& coh_U_curr, bool compute_intersected_terms)
{
	// Compute the current opening vector in the ttn coordinate system
	std::vector<double>& delta = sens_input.delta;
	ProblemUtilities::cohesiveDeltaFast(delta, shape, coh_U_curr, rotation_matrix);

	// Utilize the material constitutive function to get Dcr for any type of cohesive material
	Material::sensitivity_output_params* sens_output = mat->SensitivityConstitutive(sens_input);   // Computes the A-matrix for the current material
	std::vector<double>& dtraction = sens_output->dtraction_dd;
	if (dtraction.size() != delta.size()) // Called for a material that doesn't implement sensitivity
	{
		dtraction.clear();
		dtraction.resize(delta.size());
	}

	// Already have the delta, just have to copy it to my input object
	input.delta = delta;
	
	if (compute_intersected_terms)
	{
		// Compute the constitutive relations for this material
		Material::output_params* output = mat->Constitutive(input);   // Computes the D-matrix for the current material and updates any internal variables

		std::vector<double> t_div_v = Utilities::scale(output->traction, div_v);
		std::vector<double> dR_NU;
		ProblemUtilities::cohesiveDeltaFast(dR_NU, shape, coh_U_curr, drotation_matrix_dd);
		std::vector<double> A_dR_NU = (output->Dmat) * dR_NU;
		std::vector<double> temp1 = Utilities::plus(dtraction, A_dR_NU);
		std::vector<double> temp2 = Utilities::plus(temp1, t_div_v);
		std::vector<double> temp3;
		ProblemUtilities::cohesiveInternalForceFast(temp3, shape, temp2, rotation_matrix); // computes Nc^T * R^T * temp2
		std::vector<double> temp4;
		ProblemUtilities::cohesiveInternalForceFast(temp4, shape, output->traction, drotation_matrix_dd); // computes Nc^T * R_star^T * traction
		dP_el_coh = Utilities::plus(temp3, temp4);
	}

	// This element only has ISV contributions
	else
		ProblemUtilities::cohesiveInternalForceFast(dP_el_coh, shape, dtraction, rotation_matrix); // computes Nc^T * R^T * temp2	
}




void SensitivityShapeParameterStructural::KernelCohesiveUpdate(const std::vector<double>& shape, const DenseMatrix<double>& rotation_matrix,
										 					   Material* mat, Material::sensitivity_input_params& input,
															   const std::vector<double>& coh_U_curr, std::vector<double>& dcoh_U_curr)
{
	// Compute the current opening vector in the ntt coordinate system
	ProblemUtilities::cohesiveDeltaFast(input.delta, shape, coh_U_curr, rotation_matrix);

	// Compute the derivative of the opening
	ProblemUtilities::cohesiveDeltaFast(input.ddelta_dd, shape, dcoh_U_curr, rotation_matrix);
	
	// Update the internal variable sensitivities
	mat->updateSensitivityISVs(input);
}
