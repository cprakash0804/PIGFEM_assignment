#ifndef _MACRO_STRAIN_LINEAR_H_
#define _MACRO_STRAIN_LINEAR_H_
#include "BodyLoad.h"


class BodyLoadMacroStrainLinear : public BodyLoad
{
	private:
		/*
		 * The macro strain that is to be applied everywhere
		 */
		std::vector<double> _macro_strain;


	public:

		/*
		 * Main function that each Boy force will have to implement
		 * Takes the (x, y, z) coordinates and returns the 
		 */
		virtual void BF_Kernel(std::vector<double>& P_el, const std::vector<double>& points,
							   const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
							   const double& J, const double& W,
							   Material* mat, Material::input_params& input);


		/*
		 * Function used to easily set a vector-valued parameter
		 */
		virtual void set_vec_parameter(std::string name, std::vector<double> val);
		virtual std::vector<double> get_vec_parameter(std::string name);


		/*
		 * Differentiates between the types of body loads that are being applied
		 */
		virtual body_load_type get_type() {return MULTISCALE_LINEAR;};
};



#endif
