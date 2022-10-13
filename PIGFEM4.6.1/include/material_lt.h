/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _LTMATERIAL_H_
#define _LTMATERIAL_H_

#include "material.h"

// This material is simlar to the LEI material, but instead of a constituitive matrix based on the C_ijkl tensor,
// the constiutive matrix is just an identity multiplied by the the conductivity
class LinearThermalMaterial : public Material
{
	protected:
		double _kappa;
		bool kappa_set;
		
		void set_kappa(double E);
		double get_kappa();

		virtual void computeLinearDmat(Material::input_params& input);
		
	public:
		LinearThermalMaterial();
		virtual ~LinearThermalMaterial() {};
		
		// Main constitutive subroutine
		virtual Material::output_params* Constitutive(Material::input_params& input);

		// Helper functions
		virtual material_type get_type() {return LINEAR_THERMAL;};
		virtual classification get_classification() {return THERMAL;};
		virtual bool linear() {return true;};
		virtual void set_parameter(std::string name, double val);
		virtual double get_parameter(std::string name);
		
		// Annoying functions needed for communicating a material
		virtual Material* allocate();			// virtual function
		virtual void copy(Material* other_mat);	// virtual function
};

#endif
