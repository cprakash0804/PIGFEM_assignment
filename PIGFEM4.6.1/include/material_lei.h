/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _LEIMATERIAL_H_
#define _LEIMATERIAL_H_

#include "material.h"

class LinearElasticIsotropicMaterial : public Material
{
	protected:
		// Elastic parameters
		double _E;
		double _nu;
		double _lambda;
		double _mu;
		double _K;
		int n_set;
		bool E_set;
		bool nu_set;
		bool lambda_set;
		bool mu_set;
		bool K_set;
		
		void set_E(double E);
		double get_E();
		void set_nu(double nu);
		double get_nu();
		void set_lambda(double lambda);
		double get_lambda();
		void set_mu(double mu);
		double get_mu();
		void set_K(double K);
		double get_K();
		void set_thermal_shrink(double t);
		double get_thermal_shrink();

		// Thermal expansion parameters
		double _thermal_exp;
		bool _thermal_exp_set;
		void set_thermal_exp(double t);
		double get_thermal_exp();

		// Function to compute the D matrix for the given dimension
		virtual void computeLinearDmat(Material::input_params& input);

		// A vector fo the sensitivities of the D matrix
		void computeDmatSensitivities(id_type dim, bool plane_strain);

		// A function to efficiently compute the stress state and avoid zero multiplies
		void computeStress(std::vector<double>& stress, const DenseMatrix<double>& Dmat, const std::vector<double>& strain);

		std::vector<DenseMatrix<double> > _Dmat_sensitivities;
		bool _sensitivities_computed;
		
	public:
		LinearElasticIsotropicMaterial();
		virtual ~LinearElasticIsotropicMaterial() {};

		// Main constitutive subroutine
		virtual Material::output_params* Constitutive(Material::input_params& input);

		// Helper functions
		virtual material_type get_type() {return LINEAR_ELASTIC_ISOTROPIC;};
		virtual classification get_classification() {return STRUCTURAL;};
		virtual bool linear() {return true;};
		virtual void set_parameter(std::string name, double val);
		virtual double get_parameter(std::string name);
		
		// Annoying functions needed for communicating a material
		virtual Material* allocate();			// virtual function
		virtual void copy(Material* other_mat);	// virtual function

		// Functions related to material sensitivity
		virtual SensitivityMaterialParameter* getSensitivityParameter(std::string name);
		virtual Material::sensitivity_output_params* SensitivityConstitutive(Material::sensitivity_input_params& input); 
};

#endif
