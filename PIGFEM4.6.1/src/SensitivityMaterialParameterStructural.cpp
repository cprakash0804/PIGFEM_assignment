/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#include "SensitivityMaterialParameterStructural.h"
#include "Utilities.h"
#include <algorithm>


// Kernel function that actually does the computes the physics behind the specific problem
// Computes the local contribution to the stiffness matrix and internal load vector
// Also updates the internal variable storage in the input variable
void SensitivityMaterialParameterStructural::KernelVolumetric(std::vector<double>& dP_el,
											const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
											Material* mat, Material::sensitivity_input_params& input,
											const std::vector<double>& elem_U_curr)
{
	// Compute the current strain
	std::vector<double>& strain = input.strain;
	ProblemUtilities::SmallStrainFast(strain, shape_grad, elem_U_curr);

	// Compute the constitutive relations for this quadrature point
	Material::sensitivity_output_params* output = mat->SensitivityConstitutive(input);
	std::vector<double>& dstress = output->dstress_dd;
	if (dstress.size() != strain.size()) // Called for a material that doesn't implement sensitivity
	{
		dstress.clear();
		dstress.resize(strain.size());
	}

	// Compute the internal force
	ProblemUtilities::SmallStrainInternalForceFast(dP_el, shape_grad, dstress);
}

// B-matrix variant
void SensitivityMaterialParameterStructural::KernelVolumetric(std::vector<double>& dP_el,
											const std::vector<double>& shape, const DenseMatrix<double>& B,
											Material* mat, Material::sensitivity_input_params& input,
											const std::vector<double>& elem_U_curr)
{
	// Compute the current strain
	std::vector<double>& strain = input.strain;
	ProblemUtilities::SmallStrainFast(strain, B, elem_U_curr);

	// Compute the constitutive relations for this quadrature point
	Material::sensitivity_output_params* output = mat->SensitivityConstitutive(input);
	std::vector<double>& dstress = output->dstress_dd;
	if (dstress.size() != strain.size()) // Called for a material that doesn't implement sensitivity
	{
		dstress.clear();
		dstress.resize(strain.size());
	}

	// Compute the internal force
	ProblemUtilities::SmallStrainInternalForceFast(dP_el, B, dstress);
}


void SensitivityMaterialParameterStructural::KernelVolumetricUpdate(const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
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
void SensitivityMaterialParameterStructural::KernelVolumetricUpdate(const std::vector<double>& shape, const DenseMatrix<double>& B,
																	Material* mat, Material::sensitivity_input_params& input,
									 								const std::vector<double>& elem_U_curr, std::vector<double>& delem_U_curr)
{
	// Compute the current strain
	ProblemUtilities::SmallStrainFast(input.strain, B, elem_U_curr);

	// Compute the current strain derivative
	ProblemUtilities::SmallStrainFast(input.dstrain_dd, B, delem_U_curr);

	
	mat->updateSensitivityISVs(input);
}














void SensitivityMaterialParameterStructural::KernelCohesive(std::vector<double>& dP_el_coh,
											const std::vector<double>& shape, const DenseMatrix<double>& rotation_matrix,
											Material* mat, Material::sensitivity_input_params& input,
											const std::vector<double>& coh_U_curr)
{
	// Compute the current opening vector in the ntt coordinate system
	std::vector<double>& delta = input.delta;
	ProblemUtilities::cohesiveDeltaFast(delta, shape, coh_U_curr, rotation_matrix);

	// Utilize the material constitutive function to get Dcr for any type of cohesive material
	Material::sensitivity_output_params* output = mat->SensitivityConstitutive(input);   // Computes the A-matrix for the current material
	std::vector<double>& dtraction = output->dtraction_dd;
	if (dtraction.size() != delta.size()) // Called for a material that doesn't implement sensitivity
	{
		dtraction.clear();
		dtraction.resize(delta.size());
	}

	// Compute the partial derivative of the internal force
	ProblemUtilities::cohesiveInternalForceFast(dP_el_coh, shape, dtraction, rotation_matrix);
}




void SensitivityMaterialParameterStructural::KernelCohesiveUpdate(const std::vector<double>& shape, const DenseMatrix<double>& rotation_matrix,
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
