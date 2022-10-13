/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#ifndef _SENSITIVITY_MATERIAL_PARAMETER_STRUCTURAL_H_
#define _SENSITIVITY_MATERIAL_PARAMETER_STRUCTURAL_H_
#include "SensitivityMaterialParameter.h"


class SensitivityMaterialParameterStructural : public SensitivityMaterialParameter
{
	public:

		// The functions that need to be implemented to compute the sensitivities
		virtual void KernelVolumetric(std::vector<double>& dP_el,
									  const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
									  Material* mat, Material::sensitivity_input_params& input,
									  const std::vector<double>& elem_U_curr);
		// B-matrix variant
		virtual void KernelVolumetric(std::vector<double>& dP_el,
									  const std::vector<double>& shape, const DenseMatrix<double>& B,
									  Material* mat, Material::sensitivity_input_params& input,
									  const std::vector<double>& elem_U_curr);
		virtual void KernelCohesive(std::vector<double>& dP_el_coh,
									const std::vector<double>& shape, const DenseMatrix<double>& rotation_matrix,
									Material* mat, Material::sensitivity_input_params& input,
									const std::vector<double>& coh_U_curr);

		// The functions that need to be implemented to update internal variable senitivities
		virtual void KernelVolumetricUpdate(const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
											Material* mat, Material::sensitivity_input_params& input,
									 		const std::vector<double>& elem_U_curr, std::vector<double>& delem_U_curr);
		// B-matrix variant
		virtual void KernelVolumetricUpdate(const std::vector<double>& shape, const DenseMatrix<double>& B,
											Material* mat, Material::sensitivity_input_params& input,
									 		const std::vector<double>& elem_U_curr, std::vector<double>& delem_U_curr);
		virtual void KernelCohesiveUpdate(const std::vector<double>& shape, const DenseMatrix<double>& rotation_matrix,
										  Material* mat, Material::sensitivity_input_params& input,
										  const std::vector<double>& coh_U_curr, std::vector<double>& dcoh_U_curr);
};


#endif