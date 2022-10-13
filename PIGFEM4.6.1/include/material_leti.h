/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _LETIMATERIAL_H_
#define _LETIMATERIAL_H_

#include "material.h"

class LinearElasticTransverselyIsotropicMaterial : public Material
{
	protected:
		double _E1;
		double _E2;
		double _G12;
		double _G23;
		double _nu12;
		double _alpha;	// First rotation angle (about current z-axis)
		double _beta;	// Secnd rotation angle (about current y-axis)
		int n_set;
		bool E1_set;
		bool E2_set;
		bool G12_set;
		bool G23_set;
		bool nu12_set;
		
		void set_E1(double E1);
		double get_E1();
		void set_E2(double E2);
		double get_E2();
		void set_G12(double G12);
		double get_G12();
		void set_G23(double G23);
		double get_G23();
		void set_nu12(double nu12);
		double get_nu12();
		void set_alpha(double alpha);
		double get_alpha();
		void set_beta(double beta);
		double get_beta();

		// Thermal expansion parameters
		double _thermal_exp1;
		bool _thermal_exp1_set;
		double _thermal_exp2;
		bool _thermal_exp2_set;
		void set_thermal_exp(double t, int dir);
		double get_thermal_exp(int dir);

		// Function to compute the D matrix for the given dimension
		virtual void computeLinearDmat(Material::input_params& input);
		
	public:
		LinearElasticTransverselyIsotropicMaterial();
		virtual ~LinearElasticTransverselyIsotropicMaterial() {};

		// Main constitutive subroutine
		virtual Material::output_params* Constitutive(Material::input_params& input);

		// Helper functions
		virtual material_type get_type() {return LINEAR_ELASTIC_TRANSVERSELY_ISOTROPIC;};
		virtual classification get_classification() {return STRUCTURAL;};
		virtual bool linear() {return true;};
		virtual void set_parameter(std::string name, double val);
		virtual double get_parameter(std::string name);
		
		// Annoying functions needed for communicating a material
		virtual Material* allocate();			// virtual function
		virtual void copy(Material* other_mat);	// virtual function

};

#endif
