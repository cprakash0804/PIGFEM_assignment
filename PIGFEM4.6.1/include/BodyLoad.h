#ifndef _BODY_FORCE_H_
#define _BODY_FORCE_H_
#include <string>
#include <vector>
#include "common.h"
#include "material.h"
#include "DenseMatrix.h"



class BodyLoad
{
	public:

		/*
		 * Main function that each Boy force will have to implement
		 * Takes the (x, y, z) coordinates and returns the 
		 */
		virtual void BF_Kernel(std::vector<double>& P_el, const std::vector<double>& points,
							   const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
							   const double& J, const double& W,
							   Material* mat, Material::input_params& input) = 0;


		/*
		 * Function used to set a general single parameter
		 */
		virtual void set_parameter(std::string name, double val) {};
		virtual double get_parameter(std::string name) {err_message("Attemping to call get_parameter for a body load with no parameters");};


		/*
		 * Function used to easily set a vector-valued parameter
		 */
		virtual void set_vec_parameter(std::string name, std::vector<double> val) {};
		virtual std::vector<double> get_vec_parameter(std::string name) {err_message("Attemping to call get_vec_parameter for a body load with no parameters");};


		/*
		 *Function used to easily set a matrix-valued parameter
		 */
		virtual void set_mat_parameter(std::string name, DenseMatrix<double> val) {};
		virtual DenseMatrix<double> get_mat_parameter(std::string name) {err_message("Attemping to call get_mat_parameter for a body load with no parameters");};


		/*
		 * Differentiates between the types of body loads that are being applied
		 */
		virtual body_load_type get_type() = 0;

};



#endif
