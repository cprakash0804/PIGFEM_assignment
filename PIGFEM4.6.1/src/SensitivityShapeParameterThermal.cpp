/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#include "SensitivityShapeParameterThermal.h"
#include "Utilities.h"
#include <algorithm>


// Kernel function that actually does the computes the physics behind the specific problem
// Computes the local contribution to the stiffness matrix and internal load vector
// Also updates the internal variable storage in the input variable
void SensitivityShapeParameterThermal::KernelVolumetric(std::vector<double>& dP_el,
											const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
											const std::vector<double>& dshape_dd, const std::vector<std::vector<double> >& dshape_grad_dd, double div_v,
											Material* mat, Material::input_params& input, Material::sensitivity_input_params& sens_input,
											const std::vector<double>& elem_U_curr, bool compute_intersected_terms)
{
	err_message("Thermal shape sensitivity analysis not implemented yet!");
}

// B-matrix variant
void SensitivityShapeParameterThermal::KernelVolumetric(std::vector<double>& dP_el,
											const std::vector<double>& shape, const DenseMatrix<double>& B,
											const std::vector<double>& dshape_dd, const DenseMatrix<double>& dB_dd, double div_v,
											Material* mat, Material::input_params& input, Material::sensitivity_input_params& sens_input,
											const std::vector<double>& elem_U_curr, bool compute_intersected_terms)
{
	err_message("Thermal shape sensitivity analysis not implemented yet!");
}


void SensitivityShapeParameterThermal::KernelVolumetricUpdate(const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
															  Material* mat, Material::sensitivity_input_params& input,
									 						  const std::vector<double>& elem_U_curr, std::vector<double>& delem_U_curr)
{
	err_message("Thermal shape sensitivity analysis not implemented yet!");
}

// B-matrix variant
void SensitivityShapeParameterThermal::KernelVolumetricUpdate(const std::vector<double>& shape, const DenseMatrix<double>& B,
															  Material* mat, Material::sensitivity_input_params& input,
									 						  const std::vector<double>& elem_U_curr, std::vector<double>& delem_U_curr)
{
	err_message("Thermal shape sensitivity analysis not implemented yet!");
}














void SensitivityShapeParameterThermal::KernelCohesive(std::vector<double>& dP_el_coh,
											const std::vector<double>& shape, const DenseMatrix<double>& rotation_matrix,
											const std::vector<double>& dshape_dd, const DenseMatrix<double>& drotation_matrix_dd, double div_v,
											Material* mat, Material::input_params& input, Material::sensitivity_input_params& sens_input,
											const std::vector<double>& coh_U_curr, bool compute_intersected_terms)
{
	err_message("Thermal shape sensitivity analysis not implemented yet!");
}




void SensitivityShapeParameterThermal::KernelCohesiveUpdate(const std::vector<double>& shape, const DenseMatrix<double>& rotation_matrix,
										 					Material* mat, Material::sensitivity_input_params& input,
															const std::vector<double>& coh_U_curr, std::vector<double>& dcoh_U_curr)
{
	err_message("Thermal shape sensitivity analysis not implemented yet!");
}
