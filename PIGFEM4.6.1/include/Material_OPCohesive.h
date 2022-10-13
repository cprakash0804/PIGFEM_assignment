/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#ifndef _OP_COH_MATERIAL_H_
#define _OP_COH_MATERIAL_H_

#include "Material_OPCohesiveNU.h"
#include "Solver.h"

/*
 * Cohesive law described by Ortiz and Pandolfi (1999)
 * Tractions described by an "effective" opening which
 *	is a combination of normal and tangential components
 * This class implements linear unloading on top of the
 *	no-unloading class
 */
class OPCohesiveMaterial : public OPCohesiveNUMaterial
{
	public:
		virtual ~OPCohesiveMaterial() {};

		// Main constitutive subroutine
		virtual Material::output_params* Constitutive(Material::input_params& input);

		// Helper functions
		virtual material_type get_type() {return OP_COHESIVE;};

		// Internal variable functions
		virtual id_type n_internal_vars() {return 1;};
		virtual std::vector<double> init_internal_vars() {return std::vector<double>(1,0.0);};
		virtual std::vector<double> max_internal_vars() {return std::vector<double>(1,1.0);};
		virtual std::vector<std::string> internal_vars_name();
		virtual std::vector<bool> internal_vars_print();
		
		// Annoying functions needed for copying a material
		virtual Material* allocate();			// virtual function

		// Functions related to material sensitivity
		virtual Material::sensitivity_output_params* SensitivityConstitutive(Material::sensitivity_input_params& input);
		virtual void updateSensitivityISVs(Material::sensitivity_input_params& input);
};

#endif
