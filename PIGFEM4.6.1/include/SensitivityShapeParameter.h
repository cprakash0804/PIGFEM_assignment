/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#ifndef _SENSITIVITY_SHAPE_PARAMETER_H_
#define _SENSITIVITY_SHAPE_PARAMETER_H_
#include "SensitivityParameter.h"
#include "DenseMatrix.h"
#include <string>

class NodalData;

class SensitivityShapeParameter : public SensitivityParameter
{
	protected:

		// The inclusion id for this shape parameter
		id_type _inclusion_id;

		// The name of the paramter itself
		std::string _param_name;

		// Shape function derivatives and velocity
		std::vector<std::vector<std::vector<double> > > _dshape;
		std::vector<std::vector<std::vector<std::vector<double> > > > _dshape_grad;
		std::vector<std::vector<double> > _div_v;
		std::vector<std::vector<DenseMatrix<double> > > _dRot;
		std::vector<std::vector<DenseMatrix<double> > > _dBmats;
		NodalData* _velocity;

		// The functions that need to be implemented to compute the sensitivities
		virtual void KernelVolumetric(std::vector<double>& dP_el,
									  const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
									  const std::vector<double>& dshape_dd, const std::vector<std::vector<double> >& dshape_grad_dd, double div_v,
									  Material* mat, Material::input_params& input, Material::sensitivity_input_params& sens_input,
									  const std::vector<double>& elem_U_curr, bool compute_intersected_terms) {};
		// B-matrix variant
		virtual void KernelVolumetric(std::vector<double>& dP_el,
									  const std::vector<double>& shape, const DenseMatrix<double>& shape_grad,
									  const std::vector<double>& dshape_dd, const DenseMatrix<double>& dshape_grad_dd, double div_v,
									  Material* mat, Material::input_params& input, Material::sensitivity_input_params& sens_input,
									  const std::vector<double>& elem_U_curr, bool compute_intersected_terms) {};
		virtual void KernelCohesive(std::vector<double>& dP_el_coh,
									const std::vector<double>& shape, const DenseMatrix<double>& rotation_matrix,
									const std::vector<double>& dshape_dd, const DenseMatrix<double>& drotation_matrix_dd, double div_v,
									Material* mat, Material::input_params& input, Material::sensitivity_input_params& sens_input,
									const std::vector<double>& coh_U_curr, bool compute_intersected_terms) {};

		// Function to handle all of the precomputaion of the nodal velocities and shape function derivatives
		void compute_velocity(Problem* prob);

		// Handle the elemental assembly
		virtual void assemble_elem_dPdd(Elem* el, std::vector<double>& dP_el, std::vector<double>& elem_sol, Problem* prob);

	public:

		SensitivityShapeParameter();
		virtual ~SensitivityShapeParameter();

		virtual void init();
		

		
		double& get_div_v(id_type local_e, id_type qp);
		std::vector<double>& get_shape(id_type local_e, id_type qp);
		std::vector<std::vector<double> >& get_shape_grad(id_type local_e, id_type qp);
		DenseMatrix<double>& get_rot_matrix(id_type local_e, id_type qp);
		DenseMatrix<double>& getBmat(id_type local_e, id_type qp);
		NodalData* get_velocity() const {return _velocity;};

		void set_inc_id(id_type inc_id) {_inclusion_id = inc_id;};
		void set_param_name(std::string param_name) {_param_name = param_name;};
		id_type get_inc_id() {return _inclusion_id;};
		std::string get_param_name() {return _param_name;};

		virtual sensitivity_parameter_type get_type() {return SHAPE;};
};




#endif