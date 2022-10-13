/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include "common.h"
#include "DenseMatrix.h"
#include <string.h>

#define MAX_MATERIAL_SIZE 1024 // Max size of a material in bytes (can be changed)



// Predeclaratinos
class SensitivityMaterialParameter;


class SensitivityParameter;

class Material
{
	public:
		// Input and output structs used in the call to Constitutive
		// Can add params to each for new types of problems
		struct input_params
		{
			id_type dim;
			std::vector<double> strain;
			std::vector<double> strain_rate;
			std::vector<double>* internal_vars;
			std::vector<double> delta; // Opening used for cohesive materials
			std::vector<double> location;
			double delta_t;
			double temp_change; // Used for calculation of residual strain
			bool plane_strain; // boolean as to whether or not this problem is in plane strain

			input_params()
				: dim(1), delta_t(0.0), temp_change(0.0), plane_strain(false)
			{}
		};

		struct sensitivity_input_params
		{
			id_type dim;
			std::vector<double> strain;
			std::vector<double> dstrain_dd;
			std::vector<double> strain_rate;
			std::vector<double>* internal_vars;
			std::vector<double>* internal_vars_sensitivity;
			std::vector<double> delta; // Opening used for cohesive materials
			std::vector<double> ddelta_dd;
			double delta_t;
			double temp_change; // Used for calculation of residual strain
			bool plane_strain; // boolean as to whether or not this problem is in plane strain
			std::string sensitivity_mat_name;
			std::string sensitivity_param_name;
			id_type parameter_id;

			sensitivity_input_params()
				: dim(1), delta_t(0.0), temp_change(0.0), plane_strain(false)
			{}
		};

		struct output_params
		{
			DenseMatrix<double> Dmat;
			std::vector<double> stress;
			std::vector<double> traction; // For cohesive materials
		};

		struct sensitivity_output_params
		{
			std::vector<double> dstress_dd;
			std::vector<double> dtraction_dd; // For cohesive materials
		};

	protected:
		id_type _id;   // A global material id
		std::string _name;
		bool _name_set;

		// curr_dim is the dimension associated with the current output matrix
		id_type _curr_dim;

		// The output data structure that will be returned
		Material::output_params _output;
		Material::sensitivity_output_params _sensitivity_output;

		// Function to compute the linear D matrix for the given dimension
		virtual void computeLinearDmat(Material::input_params& input) {};

	public:

		Material();
		virtual ~Material() {};


		virtual Material::output_params* Constitutive(Material::input_params& input) = 0;
		virtual material_type get_type() = 0;
		virtual classification get_classification() = 0;
		virtual bool linear() = 0;
		virtual void set_parameter(std::string name, double val) = 0;
		virtual double get_parameter(std::string name) = 0;

		// Internal variable functions (Defaults to having no internal variables)
		virtual id_type n_internal_vars() {return 0;};
		virtual std::vector<double> init_internal_vars() {return std::vector<double>();};
		virtual std::vector<double> max_internal_vars() {return std::vector<double>();};
		virtual std::vector<std::string> internal_vars_name() {return std::vector<std::string>();};
		virtual std::vector<bool> internal_vars_print() {return std::vector<bool>();};

		// Functions specific to cohesive materials
		virtual bool is_cohesive() {return false;};
		virtual id_type get_cohesive_damage_zone(std::vector<double>& delta) {return 9999999;};
		virtual double getNormalizedCohesiveFailure(std::vector<double>& delta_ntt) {return 0.0;};
		
		void set_id(id_type id) {_id = id;};
		id_type get_id() { return _id;};
		void set_name(std::string name);
		std::string get_name() {return _name;};
		bool name_set() {return _name_set;};
		
		virtual void pack(char* buf) {err_message("There are issues with packing materials currently. Please contact the developer about this issue if you need to pack a material");};
		static Material* unpack(char* buf);
		Material* allocate_and_copy();
		virtual Material* allocate() = 0;              // Allocates a new derived type
		virtual void copy(Material* other_mat) = 0;    // Copies new derived type





		// Functions related to material sensitivity
		virtual SensitivityMaterialParameter* getSensitivityParameter(std::string name) {err_message("Attempting to get the sensitivity parameter of a material that does not support sensitivity!");};
		virtual Material::sensitivity_output_params* SensitivityConstitutive(Material::sensitivity_input_params& input) {return &_sensitivity_output;};
		virtual void updateSensitivityISVs(Material::sensitivity_input_params& input) {/* Default to doing nothing for materials with no internal variables */};

};


#endif
