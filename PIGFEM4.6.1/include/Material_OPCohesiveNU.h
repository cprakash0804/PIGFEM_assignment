/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#ifndef _NU_OP_COH_MATERIAL_H_
#define _NU_OP_COH_MATERIAL_H_

#include "material.h"

/*
 * Cohesive law described by Ortiz and Pandolfi (1999)
 * Tractions described by an "effective" opening which
 *	is a combination of normal and tangential components
 */
class OPCohesiveNUMaterial : public Material
{
	protected:
		// Linear elastic parameters
		double _sigma_c;
		double _delta_c;
		double _beta2;
		bool _s_set;
		bool _d_set;
		bool _b_set;
		bool ready; // Stores whether or not this material is ready for calls to constitutive
		double _visc;
		bool _visc_set;

		void set_sc(double sc);
		double get_sc();
		void set_dc(double dc);
		double get_dc();
		void set_b(double b);
		double get_b();
		void set_visc(double visc);
		double get_visc();

		double getDeltaEff(std::vector<double>& delta_ttn);

	public:
		OPCohesiveNUMaterial();
		virtual ~OPCohesiveNUMaterial() {};

		// Main constitutive subroutine
		virtual Material::output_params* Constitutive(Material::input_params& input);

		// Helper functions
		virtual material_type get_type() {return OP_COHESIVE_NO_UNLOADING;};
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
