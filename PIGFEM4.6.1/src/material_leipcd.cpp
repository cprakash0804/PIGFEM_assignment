/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "material_leipcd.h"
#include "Utilities.h"
#include <cmath>
#include <algorithm>
#include <iostream>



Material* LinearElasticIsotropicProblemControlledDamageMaterial::allocate()
{
	return new LinearElasticIsotropicProblemControlledDamageMaterial;
}
// This function assumes the current material is set up correctly. Which it should't be possible for it not to be
void LinearElasticIsotropicProblemControlledDamageMaterial::copy(Material* other_mat)
{
	other_mat->set_name(_name);
	other_mat->set_id(_id);
	if(E_set)
		other_mat->set_parameter("E", _E);
	if(nu_set)
		other_mat->set_parameter("nu", _nu);
	if(lambda_set)
		other_mat->set_parameter("lambda", _lambda);
	if(mu_set)
		other_mat->set_parameter("mu", _mu);
	if(K_set)
		other_mat->set_parameter("K", _K);
	if (_thermal_exp_set)
		other_mat->set_parameter("THERMAL_EXPANSION", _thermal_exp);
}


LinearElasticIsotropicProblemControlledDamageMaterial::LinearElasticIsotropicProblemControlledDamageMaterial()
	: n_set(0), E_set(false), nu_set(false), lambda_set(false), mu_set(false), K_set(false),
	  _thermal_exp(0.0), _thermal_exp_set(false)
{
}

std::vector<std::string> LinearElasticIsotropicProblemControlledDamageMaterial::internal_vars_name()
{
	std::vector<std::string> vec = {"damage"};
	return vec;
}

std::vector<bool> LinearElasticIsotropicProblemControlledDamageMaterial::internal_vars_print()
{
	std::vector<bool> vec = {true};
	return vec;
}




void LinearElasticIsotropicProblemControlledDamageMaterial::computeLinearDmat(Material::input_params& input)
{
	if(n_set<2)
		err_message("Please set two independant parameters for a linear elastic isotropic material.");

	id_type dim = input.dim;
	bool plane_strain = input.plane_strain;

	// Compute the constitutive matrix
	std::vector<id_type> voigt;
	if(dim==1)
	{
		_linear_Dmat.resize(1,1);
		_linear_Dmat(0,0) = _E;
	}
	else if(dim==2)
	{
		_linear_Dmat.clear();
		_linear_Dmat.resize(3, 3);
		if(plane_strain)
		{
			double coef = _E / ((1.0+_nu) * (1.0-2.0*_nu));
			_linear_Dmat(0,0) = coef * (1.0-_nu);
			_linear_Dmat(0,1) = coef * (_nu);
			_linear_Dmat(1,0) = coef * (_nu);
			_linear_Dmat(1,1) = coef * (1.0-_nu);
			_linear_Dmat(2,2) = coef * (1.0-2.0*_nu)/2.0;
		}
		else
		{
			double coef = _E / (1.0-_nu*_nu);
			_linear_Dmat(0,0) = coef;
			_linear_Dmat(0,1) = coef * (_nu);
			_linear_Dmat(1,0) = coef * (_nu);
			_linear_Dmat(1,1) = coef;
			_linear_Dmat(2,2) = coef * (1.0-_nu)/2.0;	// Division by 2 is to account for the engineering strain we use
		}
	}
	else if(dim==3)
	{
		_linear_Dmat.clear();
		_linear_Dmat.resize(6, 6);
		double coef = _E / ((1.0+_nu) * (1.0-2.0*_nu));
		_linear_Dmat(0,0) = coef * (1.0-_nu);
		_linear_Dmat(0,1) = coef * (_nu);
		_linear_Dmat(0,2) = coef * (_nu);
		_linear_Dmat(1,0) = coef * (_nu);
		_linear_Dmat(1,1) = coef * (1.0-_nu);
		_linear_Dmat(1,2) = coef * (_nu);
		_linear_Dmat(2,0) = coef * (_nu);
		_linear_Dmat(2,1) = coef * (_nu);
		_linear_Dmat(2,2) = coef * (1.0-_nu);
		_linear_Dmat(3,3) = coef * (1.0-2.0*_nu)/2;
		_linear_Dmat(4,4) = coef * (1.0-2.0*_nu)/2;
		_linear_Dmat(5,5) = coef * (1.0-2.0*_nu)/2;
	}
	else
		err_message("Dimensions for the D matrix must be less than or equal to 3.");
}

void LinearElasticIsotropicProblemControlledDamageMaterial::apply_damage(DenseMatrix<double>& D, double damage)
{
	unsigned char nrows = D.n_rows();

	if (damage != 0.0)
	{
		switch (nrows)
		{
			case 3: // Probably the most common case (2D)
				for (unsigned char i=0; i<2; ++i)
					for (unsigned char j=0; j<2; ++j)
						D(i,j) = (1.0 - damage) * _linear_Dmat(i, j);
				D(2,2) = (1.0 - damage) * _linear_Dmat(2,2);
				break;
			case 6:	// (3D)
				for (unsigned char i=0; i<3; ++i)
					for (unsigned char j=0; j<2; ++j)
						D(i,j) = (1.0 - damage) * _linear_Dmat(i, j);
				D(3,3) = (1.0 - damage) * _linear_Dmat(3,3);
				D(4,4) = (1.0 - damage) * _linear_Dmat(4,4);
				D(5,5) = (1.0 - damage) * _linear_Dmat(5,5);
				break;
			case 1: // (1D)
				D(0,0) = (1.0 - damage) * _linear_Dmat(0,0);
				break;
			default:
				err_message("Invalid number of rows in the apply damage function.");
		}
	}
	else
	{
		switch (nrows)
		{
			case 3: // Probably the most common case (2D)
				for (unsigned char i=0; i<2; ++i)
					for (unsigned char j=0; j<2; ++j)
						D(i,j) = _linear_Dmat(i, j);
				D(2,2) = _linear_Dmat(2,2);
				break;
			case 6:	// (3D)
				for (unsigned char i=0; i<3; ++i)
					for (unsigned char j=0; j<2; ++j)
						D(i,j) = _linear_Dmat(i, j);
				D(3,3) = _linear_Dmat(3,3);
				D(4,4) = _linear_Dmat(4,4);
				D(5,5) = _linear_Dmat(5,5);
				break;
			case 1: // (1D)
				D(0,0) = _linear_Dmat(0,0);
				break;
			default:
				err_message("Invalid number of rows in the apply damage function.");
		}
	}
}
void LinearElasticIsotropicProblemControlledDamageMaterial::compute_stress(std::vector<double>& stress, const DenseMatrix<double>& D, const std::vector<double>& strain)
{
	unsigned char nrows = D.n_rows();

	switch (nrows)
	{
		case 3: // Probably the most common case (2D)
			stress[0] = D(0,0)*strain[0] + D(0,1)*strain[1];
			stress[1] = D(1,0)*strain[0] + D(1,1)*strain[1];
			stress[2] = D(2,2)*strain[2];
			break;


		case 6:	// (3D)
			stress[0] = D(0,0)*strain[0] + D(0,1)*strain[1] + D(0,2)*strain[2];
			stress[1] = D(1,0)*strain[0] + D(1,1)*strain[1] + D(1,2)*strain[2];
			stress[2] = D(2,0)*strain[0] + D(2,1)*strain[1] + D(2,2)*strain[2];
			stress[3] = D(3,3)*strain[3];
			stress[4] = D(4,4)*strain[4];
			stress[5] = D(5,5)*strain[5];
			break;
	

		case 1: // (1D)
			stress[0] = D(0,0) * strain[0];
			break;
		

		default:
			err_message("Invalid number of rows in the apply damage function.");
	}
}



Material::output_params* LinearElasticIsotropicProblemControlledDamageMaterial::Constitutive(Material::input_params& input)
{
	// Get references to the output matrix and stress
	DenseMatrix<double>& Dmat_bar = _output.Dmat;
	std::vector<double>& stress = _output.stress;

	// Update the D matrix for 
	if(input.dim != _curr_dim)
	{
		computeLinearDmat(input);
		_curr_dim = input.dim;
		Dmat_bar = _linear_Dmat;
		stress.resize(3);
	}

	// Get a reference to the damage and strain
	if ((*input.internal_vars).size() != 1)
		err_message("Invalid number of internal variables");
	double& curr_damage = (*input.internal_vars)[0];
	std::vector<double>& strain = input.strain;
	if (_thermal_exp_set)
	{
		std::vector<double> thermal_strain(strain.size(), 0.0);
		double strain_val = -1.0 * input.temp_change * _thermal_exp; // Get the volumetric strain component from the material thermal expansion ratio
		switch (input.dim)
		{
			case 1:
				thermal_strain[0] = strain_val;
				break;
			case 2:
				std::fill(thermal_strain.begin(), thermal_strain.begin()+2, strain_val);
				break;
			case 3:
				std::fill(thermal_strain.begin(), thermal_strain.begin()+3, strain_val);
				break;
			default:
				err_message("Wrong dimension somehow...");
		}
		Utilities::_VecAXPY(1.0, thermal_strain, strain); // strain = 1.0*thermal_strain + strain
	}

	apply_damage(Dmat_bar, curr_damage);
	compute_stress(stress, Dmat_bar, strain);
	//Dmat_bar = (1.0 - curr_damage) * _linear_Dmat;
	//stress = Dmat_bar*strain;

	return &_output;
}



void LinearElasticIsotropicProblemControlledDamageMaterial::set_parameter(std::string name, double val)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="E" || name=="YOUNG'S MODULUS" || name=="YOUNGS MODULUS")
		set_E(val);
	else if(name=="NU" || name=="POISSON'S RATIO" || name=="POISSONS RATIO")
		set_nu(val);
	else if(name=="LAMBDA" || name=="FIRST LAME CONSTANT")
		set_lambda(val);
	else if(name=="MU" || name=="G" || name=="SECOND LAME CONSTANT" || name=="SHEAR MODULUS")
		set_mu(val);
	else if(name=="K" || name=="BULK MODULUS")
		set_K(val);
	else if (name=="THERMAL_EXPANSION" || name=="THERMAL EXPANSION" || name=="THERMAL_SHINKAGE" || name=="THERMAL SHRINKAGE")
		set_thermal_exp(val);
	else
		err_message(name << " is not a valid parameter name");
}
double LinearElasticIsotropicProblemControlledDamageMaterial::get_parameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="E" || name=="YOUNG'S MODULUS" || name=="YOUNGS MODULUS")
		return get_E();
	else if(name=="NU" || name=="POISSON'S RATIO" || name=="POISSONS RATIO")
		return get_nu();
	else if(name=="LAMBDA" || name=="FIRST LAME CONSTANT")
		return get_lambda();
	else if(name=="MU" || name=="G" || name=="SECOND LAME CONSTANT" || name=="SHEAR MODULUS")
		return get_mu();
	else if(name=="K" || name=="BULK MODULUS")
		return get_K();
	else if (name=="THERMAL_EXPANSION" || name=="THERMAL EXPANSION" || name=="THERMAL_SHINKAGE" || name=="THERMAL SHRINKAGE")
		return get_thermal_exp();
	else
		err_message(name << " is not a valid parameter name");
}



void LinearElasticIsotropicProblemControlledDamageMaterial::set_E(double E)
{
	if(n_set>=2 && !E_set)
		err_message("Can only set 2 independant parameters in a linear elastic isotropic material.");
	
	_E = E;
	if(!E_set) // setting value for the first time
	{
		n_set = n_set + 1;
		E_set = true;
	}
	if(n_set==2)
	{
		if(nu_set) // Parameters E and nu were set
		{
			_lambda = (_E*_nu)/((1.0+_nu)*(1.0-2.0*_nu));
			_mu = _E/(2.0*(1.0+_nu));
			_K = _E/(3.0*(1-2.0*_nu));
		}
		else if(lambda_set) // Parameters E and lambda were set
		{
			double R = sqrt(_E*_E + 9.0*_lambda*_lambda + 2.0*_E*_lambda);
			_nu = 2.0*_lambda/(_E+_lambda+R);
			_mu = (_E-3.0*_lambda+R)/4.0;
			_K = (_E+3.0*_lambda+R)/6.0;
		}
		else if(mu_set)
		{
			_lambda = _mu*(_E-2.0*_mu)/(3.0*_mu-_E);
			_nu = _E/(2.0*_mu) - 1.0;
			_K = _E*_mu/(3.0*(3.0*_mu-_E));
		}
		else if(K_set)
		{
			_lambda = 3.0*_K*(3.0*_K-_E)/(9.0*_K-_E);
			_nu = (3.0*_K-_E)/(6.0*_K);
			_mu = 3.0*_K*_E/(9.0*_K-_E);
		}
	}
}
double LinearElasticIsotropicProblemControlledDamageMaterial::get_E()
{
	if(n_set>=2 || E_set)
		return _E;
	else
		err_message("Young's modulus could not be determined.");
}




void LinearElasticIsotropicProblemControlledDamageMaterial::set_nu(double nu)
{
	if(n_set>=2 && !nu_set)
		err_message("Can only set 2 independant parameters in a linear elastic isotropic material.");
	
	_nu = nu;
	if(!nu_set) // setting value for the first time
	{
		n_set = n_set + 1;
		nu_set = true;
	}
	if(n_set==2)
	{
		if(E_set) // Parameters E and nu were set
		{
			_lambda = (_E*_nu)/((1.0+_nu)*(1.0-2.0*_nu));
			_mu = _E/(2.0*(1.0+_nu));
			_K = _E/(3.0*(1-2.0*_nu));
		}
		else if(lambda_set) // Parameters nu and lambda were set
		{
			_E = _lambda*(1.0+_nu)*(1.0-2.0*_nu)/_nu;
			_mu = _lambda*(1.0-2.0*_nu)/(2.0*_nu);
			_K = _lambda*(1.0+_nu)/(3.0*_nu);
		}
		else if(mu_set)
		{
			_lambda = 2.0*_mu*_nu/(1.0-2.0*_nu);
			_E = 2.0*_mu*(1.0+_nu);
			_K = 2.0*_mu*(1.0+_nu)/(3.0*(1.0-2.0*_nu));
		}
		else if(K_set)
		{
			_lambda = 3.0*_K*_nu/(1+_nu);
			_E = 3.0*_K*(1.0-2.0*_nu);
			_mu = 3.0*_K*(1.0-2.0*_nu)/(2.0*(1.0+_nu));
		}
	}
}
double LinearElasticIsotropicProblemControlledDamageMaterial::get_nu()
{
	if(n_set>=2 || nu_set)
		return _nu;
	else
		err_message("Poisson's ratio could not be determined.");
}




void LinearElasticIsotropicProblemControlledDamageMaterial::set_lambda(double lambda)
{
	if(n_set>=2 && !lambda_set)
		err_message("Can only set 2 independant parameters in a linear elastic isotropic material.");
	
	_lambda = lambda;
	if(!lambda_set) // setting value for the first time
	{
		n_set = n_set + 1;
		lambda_set = true;
	}
	if(n_set==2)
	{
		if(E_set) // Parameters E and lambda were set
		{
			double R = sqrt(_E*_E + 9.0*_lambda*_lambda + 2.0*_E*_lambda);
			_nu = 2.0*_lambda/(_E+_lambda+R);
			_mu = (_E-3.0*_lambda+R)/4.0;
			_K = (_E+3.0*_lambda+R)/6.0;
		}
		else if(nu_set) // Parameters nu and lambda were set
		{
			_E = _lambda*(1.0+_nu)*(1.0-2.0*_nu)/_nu;
			_mu = _lambda*(1.0-2.0*_nu)/(2.0*_nu);
			_K = _lambda*(1.0+_nu)/(3.0*_nu);
		}
		else if(mu_set)
		{
			_nu = _lambda/(2.0*(_lambda+_mu));
			_E = _mu*(3.0*_lambda+2.0*_mu)/(_lambda+_mu);
			_K = _lambda+2.0*_mu/3.0;
		}
		else if(K_set)
		{
			_nu = _lambda/(3.0*_K-_lambda);
			_E = 9.0*_K*(_K-_lambda)/(3.0*_K-_lambda);
			_mu = 3.0*(_K-_lambda)/2.0;
		}
	}
}
double LinearElasticIsotropicProblemControlledDamageMaterial::get_lambda()
{
	if(n_set>=2 || lambda_set)
		return _lambda;
	else
		err_message("First Lame constant could not be determined.");
}




void LinearElasticIsotropicProblemControlledDamageMaterial::set_mu(double mu)
{
	if(n_set>=2 && !mu_set)
		err_message("Can only set 2 independant parameters in a linear elastic isotropic material.");
	
	_mu = mu;
	if(!mu_set) // setting value for the first time
	{
		n_set = n_set + 1;
		mu_set = true;
	}
	if(n_set==2)
	{
		if(E_set) // Parameters E and mu were set
		{
			_lambda = _mu*(_E-2.0*_mu)/(3.0*_mu-_E);
			_nu = _E/(2.0*_mu) - 1.0;
			_K = _E*_mu/(3.0*(3.0*_mu-_E));
		}
		else if(nu_set) // Parameters nu and lambda were set
		{
			_lambda = 2.0*_mu*_nu/(1.0-2.0*_nu);
			_E = 2.0*_mu*(1.0+_nu);
			_K = 2.0*_mu*(1.0+_nu)/(3.0*(1.0-2.0*_nu));
		}
		else if(lambda_set)
		{
			_nu = _lambda/(2.0*(_lambda+_mu));
			_E = _mu*(3.0*_lambda+2.0*_mu)/(_lambda+_mu);
			_K = _lambda+2.0*_mu/3.0;
		}
		else if(K_set)
		{
			_nu = (3.0*_K-2.0*_mu)/(2.0*(3.0*_K+_mu));
			_E = 9.0*_K*_mu/(3.0*_K+_mu);
			_lambda = _K-2.0*_mu/3.0;
		}
	}
}
double LinearElasticIsotropicProblemControlledDamageMaterial::get_mu()
{
	if(n_set>=2 || mu_set)
		return _mu;
	else
		err_message("Shear modulus could not be determined.");
}




void LinearElasticIsotropicProblemControlledDamageMaterial::set_K(double K)
{
	if(n_set>=2 && !K_set)
		err_message("Can only set 2 independant parameters in a linear elastic isotropic material.");
	
	_K = K;
	if(!K_set) // setting value for the first time
	{
		n_set = n_set + 1;
		K_set = true;
	}
	if(n_set==2)
	{
		if(E_set) // Parameters E and K were set
		{
			_lambda = 3.0*_K*(3.0*_K-_E)/(9.0*_K-_E);
			_nu = (3.0*_K-_E)/(6.0*_K);
			_mu = 3.0*_K*_E/(9.0*_K-_E);
		}
		else if(nu_set) // Parameters nu and K were set
		{
			_lambda = 3.0*_K*_nu/(1+_nu);
			_E = 3.0*_K*(1.0-2.0*_nu);
			_mu = 3.0*_K*(1.0-2.0*_nu)/(2.0*(1.0+_nu));
		}
		else if(lambda_set)
		{
			_nu = _lambda/(3.0*_K-_lambda);
			_E = 9.0*_K*(_K-_lambda)/(3.0*_K-_lambda);
			_mu = 3.0*(_K-_lambda)/2.0;
		}
		else if(mu_set)
		{
			_nu = (3.0*_K-2.0*_mu)/(2.0*(3.0*_K+_mu));
			_E = 9.0*_K*_mu/(3.0*_K+_mu);
			_lambda = _K-2.0*_mu/3.0;
		}
	}
}
double LinearElasticIsotropicProblemControlledDamageMaterial::get_K()
{
	if(n_set>=2 || K_set)
		return _K;
	else
		err_message("Bulk modulus could not be determined.");
}


void LinearElasticIsotropicProblemControlledDamageMaterial::set_thermal_exp(double val)
{
	_thermal_exp = val;
	_thermal_exp_set = true;
}
double LinearElasticIsotropicProblemControlledDamageMaterial::get_thermal_exp()
{
	if (_thermal_exp_set)
		return _thermal_exp;
	else
		err_message("Thermal expansion ratio could not be determined.");
}
