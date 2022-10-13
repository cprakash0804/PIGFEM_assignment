/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#ifndef _SENSITIVITY_MATERIAL_PARAMETER_H_
#define _SENSITIVITY_MATERIAL_PARAMETER_H_
#include "SensitivityParameter.h"
#include <string>
#include <vector>
#include "material.h"


class Elem;

class SensitivityMaterialParameter : public SensitivityParameter
{
	protected:

		// The name fo the material this parameter came from
		std::string _mat_name;

		// The name of the paramter itself
		std::string _param_name;

		// The functions that need to be implemented to compute the sensitivities
		virtual void KernelVolumetric(std::vector<double>& dP_el,
									  const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
									  Material* mat, Material::sensitivity_input_params& input,
									  const std::vector<double>& elem_U_curr) = 0;
		// B-matrix variant
		virtual void KernelVolumetric(std::vector<double>& dP_el,
									  const std::vector<double>& shape, const DenseMatrix<double>& shape_grad,
									  Material* mat, Material::sensitivity_input_params& input,
									  const std::vector<double>& elem_U_curr) = 0;
		virtual void KernelCohesive(std::vector<double>& dP_el_coh,
									const std::vector<double>& shape, const DenseMatrix<double>& rotation_matrix,
									Material* mat, Material::sensitivity_input_params& input,
									const std::vector<double>& coh_U_curr) = 0;

		virtual void assemble_elem_dPdd(Elem* el, std::vector<double>& dP_el, std::vector<double>& elem_sol, Problem* prob);

	public:

		void set_mat_name(std::string mat_name) {_mat_name = mat_name;};
		void set_param_name(std::string param_name) {_param_name = param_name;};
		std::string get_mat_name() {return _mat_name;};
		std::string get_param_name() {return _param_name;};
		virtual void init();


		virtual sensitivity_parameter_type get_type() {return MATERIAL;};
};




#endif