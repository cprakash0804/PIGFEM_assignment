/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#ifndef _NU_XN_COH_MATERIAL_H_
#define _NU_XN_COH_MATERIAL_H_

#include "material.h"


/*
 * Cohesive Law developed by Xu and Needleman (1993)
 * Separate shear and normal traction components
 */
class XNCohesiveNUMaterial : public Material
{
	protected:
		// Linear elastic parameters
		double _sigma_cn;
		double _delta_cn;
		double _delta_ct;
		double _q;
		bool _sn_set;
		bool _dn_set;
		bool _dt_set;
		bool _q_set;
		bool _ready; // Stores whether or not this material is ready for calls to constitutive
		
		void set_scn(double scn);
		double get_scn();
		void set_dcn(double dcn);
		double get_dcn();
		void set_dct(double dct);
		double get_dct();
		void set_q(double q);
		double get_q();

		void getDnDt(std::vector<double>& delta_ttn, double& Dn, double& Dt);

		virtual void Compute_linear_D_Mat(Material::input_params& input) {};


	public:
		XNCohesiveNUMaterial();
		virtual ~XNCohesiveNUMaterial() {};

		// Main constitutive subroutine
		virtual Material::output_params* Constitutive(Material::input_params& input);

		// Helper functions		
		virtual material_type get_type() {return XN_COHESIVE_NO_UNLOADING;};
		virtual classification get_classification() {return STRUCTURAL;};
		virtual bool linear() {return false;};
		virtual void set_parameter(std::string name, double val);
		virtual double get_parameter(std::string name);

		// Functions specific to cohesive materials
		virtual bool is_cohesive() {return true;};
		virtual id_type get_cohesive_damage_zone(std::vector<double>& delta);
		virtual double getNormalizedCohesiveFailure(std::vector<double>& delta_ntt);
		
		// Annoying functions needed for copying a material
		virtual Material* allocate();			// virtual function
		virtual void copy(Material* other_mat);	// virtual function

		// Functions related to material sensitivity
		virtual SensitivityMaterialParameter* getSensitivityParameter(std::string name);
		virtual Material::sensitivity_output_params* SensitivityConstitutive(Material::sensitivity_input_params& input);
};

#endif