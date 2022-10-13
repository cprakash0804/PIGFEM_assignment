/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _LEIPCDMATERIAL_H_
#define _LEIPCDMATERIAL_H_

#include "material.h"
#include "Utilities.h"

#define max_max_p_stress 1e308 // Pretty much the largest number I can store

class LinearElasticIsotropicProblemControlledDamageMaterial : public Material
{
	protected:
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

		// Thermal expansion parameters
		double _thermal_exp;
		bool _thermal_exp_set;
		void set_thermal_exp(double t);
		double get_thermal_exp();

		// The linear constitutive matrix
		DenseMatrix<double> _linear_Dmat;

		// Function to compute the D matrix for the given dimension
		virtual void computeLinearDmat(Material::input_params& input);

		void apply_damage(DenseMatrix<double>& D, double damage);
		void compute_stress(std::vector<double>& stress, const DenseMatrix<double>& D, const std::vector<double>& strain);
		
	public:
		LinearElasticIsotropicProblemControlledDamageMaterial();
		virtual ~LinearElasticIsotropicProblemControlledDamageMaterial() {};

		// Main constitutive subroutine
		virtual Material::output_params* Constitutive(Material::input_params& input);

		// Helper functions
		virtual material_type get_type() {return LINEAR_ELASTIC_ISOTROPIC_PROBLEM_CONTROLLED_DAMAGE;};
		virtual classification get_classification() {return STRUCTURAL;};
		virtual bool linear() {return true;};
		virtual void set_parameter(std::string name, double val);
		virtual double get_parameter(std::string name);

		// Internal variable functions
		virtual id_type n_internal_vars() {return 1;};
		virtual std::vector<double> init_internal_vars() {return std::vector<double>(1,0.0);};
		virtual std::vector<double> max_internal_vars() {return std::vector<double>(1,999999999999999999999.9);}; // Actual maximum is 1.0 but I don't want this to get caught by the rate limiter
		virtual std::vector<std::string> internal_vars_name();
		virtual std::vector<bool> internal_vars_print();
		
		// Annoying functions needed for communicating a material
		virtual Material* allocate();			// virtual function
		virtual void copy(Material* other_mat);	// virtual function
};

#endif
