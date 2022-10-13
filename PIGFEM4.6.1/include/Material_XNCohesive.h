/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#ifndef _XN_COH_MATERIAL_H_
#define _XN_COH_MATERIAL_H_

#include "Material_XNCohesiveNU.h"


/*
 * Cohesive Law developed by Xu and Needleman (1993)
 * Separate shear and normal traction components
 * This class implments linear unloading in both normal
 *	and tangential directions on top of the no-unloading class
 */
class XNCohesiveMaterial : public XNCohesiveNUMaterial
{
	public:
		virtual ~XNCohesiveMaterial() {};

		// Main constitutive subroutine
		virtual Material::output_params* Constitutive(Material::input_params& input);

		// Helper functions
		virtual material_type get_type() {return XN_COHESIVE;};

		// Internal variable functions
		virtual id_type n_internal_vars() {return 2;};
		virtual std::vector<double> init_internal_vars() {return std::vector<double>(2, 0.0);};
		virtual std::vector<double> max_internal_vars() {return std::vector<double>(2, 1.0);};
		virtual std::vector<std::string> internal_vars_name();
		virtual std::vector<bool> internal_vars_print();
		
		// Annoying functions needed for copying a material
		virtual Material* allocate();			// virtual function

		// Functions related to material sensitivity
		virtual Material::sensitivity_output_params* SensitivityConstitutive(Material::sensitivity_input_params& input);
		virtual void updateSensitivityISVs(Material::sensitivity_input_params& input);
};

#endif