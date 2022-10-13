/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "material_lei.h"
#include "SensitivityMaterialParameterStructural.h"
#include "Utilities.h"
#include <cmath>
#include <algorithm>
#include <iostream>



Material* LinearElasticIsotropicMaterial::allocate()
{
	return new LinearElasticIsotropicMaterial;
}
// This function assumes the current material is set up correctly. Which it should't be possible for it not to be
void LinearElasticIsotropicMaterial::copy(Material* other_mat)
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


LinearElasticIsotropicMaterial::LinearElasticIsotropicMaterial()
	: n_set(0), E_set(false), nu_set(false), lambda_set(false), mu_set(false), K_set(false),
	  _thermal_exp(0.0), _thermal_exp_set(false), _sensitivities_computed(false)
{
}




void LinearElasticIsotropicMaterial::computeLinearDmat(Material::input_params& input)
{
	if(n_set<2)
		err_message("Please set two independant parameters for a linear elastic isotropic material.");

	id_type dim = input.dim;
	bool plane_strain = input.plane_strain;

	// Compute the constitutive matrix
	if(dim==1)
	{
		_output.Dmat.resize(1,1);
		_output.Dmat(0,0) = _E;
	}
	else if(dim==2)
	{
		_output.Dmat.clear();
		_output.Dmat.resize(3, 3);
		if(plane_strain)    // If the use set the plane strain boolean to true
		{
			double coef = _E / ((1.0+_nu) * (1.0-2.0*_nu));
			_output.Dmat(0,0) = coef * (1.0-_nu);		_output.Dmat(0,1) = coef * (_nu);
			_output.Dmat(1,0) = coef * (_nu);			_output.Dmat(1,1) = coef * (1.0-_nu);
			_output.Dmat(2,2) = coef * (1.0-2.0*_nu)/2.0;
		}
		else	// Defaults to true
		{
			double coef = _E / (1.0-_nu*_nu);
			_output.Dmat(0,0) = coef;			_output.Dmat(0,1) = coef * (_nu);
			_output.Dmat(1,0) = coef * (_nu);	_output.Dmat(1,1) = coef;
			_output.Dmat(2,2) = coef * (1.0-_nu)/2.0;	// Division by 2 is to account for the engineering strain we use
		}
	}
	else if(dim==3)
	{
		_output.Dmat.clear();
		_output.Dmat.resize(6, 6);
		double coef = _E / ((1.0+_nu) * (1.0-2.0*_nu));
		_output.Dmat(0,0) = coef * (1.0-_nu);			_output.Dmat(0,1) = coef * (_nu);			_output.Dmat(0,2) = coef * (_nu);
		_output.Dmat(1,0) = coef * (_nu);				_output.Dmat(1,1) = coef * (1.0-_nu);		_output.Dmat(1,2) = coef * (_nu);
		_output.Dmat(2,0) = coef * (_nu);				_output.Dmat(2,1) = coef * (_nu);			_output.Dmat(2,2) = coef * (1.0-_nu);
		_output.Dmat(3,3) = coef * (1.0-2.0*_nu)/2;		_output.Dmat(4,4) = coef * (1.0-2.0*_nu)/2;	_output.Dmat(5,5) = coef * (1.0-2.0*_nu)/2;
	}
	else
		err_message("Dimensions for the D matrix must be less than or equal to 3.");
}







// A vector fo the sensitivities of the D matrix
void LinearElasticIsotropicMaterial::computeDmatSensitivities(id_type dim, bool plane_strain)
{
	_Dmat_sensitivities.resize(5); // 5 parameters

	// Do E, G, lambda, and K all with the same algorithm, just nu is different
	std::vector<double> coefs(4);
	if (dim == 1)
	{
		for (id_type p=0; p<5; ++p)
		{
			_Dmat_sensitivities[p].clear();
			_Dmat_sensitivities[p].resize(1,1);
		}	
		coefs[0] = 1.0;								// E
		coefs[1] = (1.0+_nu) * (1.0-2.0*_nu) / _nu; // Lambda
		coefs[2] = 2.0 * (1.0+_nu);					// G
		coefs[3] = 3.0 * (1.0-2.0*_nu);				// K
	}

	else if (dim == 2)
	{
		for (id_type p=0; p<5; ++p)
		{
			_Dmat_sensitivities[p].clear();
			_Dmat_sensitivities[p].resize(3,3);
		}	
		if (plane_strain)
		{
			coefs[0] = 1.0 / ((1.0+_nu) * (1.0-2.0*_nu));	// E
			coefs[1] = 1.0 / _nu; 							// Lambda
			coefs[2] = 2.0 / (1.0-2.0*_nu);					// G
			coefs[3] = 3.0 / (1.0+_nu);						// K
		}
		else
		{
			coefs[0] = 1.0 / (1.0-_nu*_nu);					// E
			coefs[1] = (1.0-2.0*_nu) / (_nu*(1.0-_nu)); 	// Lambda
			coefs[2] = 2.0 / (1.0-_nu);						// G
			coefs[3] = 3.0 * (1.0-2.0*_nu) / (1.0-_nu*_nu);	// K
		}
	}

	else if (dim == 3)
	{
		for (id_type p=0; p<5; ++p)
		{
			_Dmat_sensitivities[p].clear();
			_Dmat_sensitivities[p].resize(6,6);
		}	
		coefs[0] = 1.0 / ((1.0+_nu) * (1.0-2.0*_nu));	// E
		coefs[1] = 1.0 / _nu; 							// Lambda
		coefs[2] = 2.0 / (1.0-2.0*_nu);					// G
		coefs[3] =3.0 / (1.0+_nu);						// K
	}

	for (id_type p=0; p<4; ++p)
	{
		if (dim == 1)
			_Dmat_sensitivities[p](0,0) = coefs[p];

		else if (dim == 2)
		{
			if (plane_strain)    // If the use set the plane strain boolean to true
			{
				_Dmat_sensitivities[p](0,0) = coefs[p] * (1.0-_nu);			_Dmat_sensitivities[p](0,1) = coefs[p] * (_nu);
				_Dmat_sensitivities[p](1,0) = coefs[p] * (_nu);				_Dmat_sensitivities[p](1,1) = coefs[p] * (1.0-_nu);
				_Dmat_sensitivities[p](2,2) = coefs[p] * (1.0-2.0*_nu)/2.0;
			}
			else	// Defaults to true
			{
				_Dmat_sensitivities[p](0,0) = coefs[p];					_Dmat_sensitivities[p](0,1) = coefs[p] * (_nu);
				_Dmat_sensitivities[p](1,0) = coefs[p] * (_nu);			_Dmat_sensitivities[p](1,1) = coefs[p];
				_Dmat_sensitivities[p](2,2) = coefs[p] * (1.0-_nu)/2.0;	// Division by 2 is to account for the engineering strain we use
			}
		}

		else if (dim == 3)
		{
			_Dmat_sensitivities[p](0,0) = coefs[p] * (1.0-_nu);			_Dmat_sensitivities[p](0,1) = coefs[p] * (_nu);				_Dmat_sensitivities[p](0,2) = coefs[p] * (_nu);
			_Dmat_sensitivities[p](1,0) = coefs[p] * (_nu);				_Dmat_sensitivities[p](1,1) = coefs[p] * (1.0-_nu);			_Dmat_sensitivities[p](1,2) = coefs[p] * (_nu);
			_Dmat_sensitivities[p](2,0) = coefs[p] * (_nu);				_Dmat_sensitivities[p](2,1) = coefs[p] * (_nu);				_Dmat_sensitivities[p](2,2) = coefs[p] * (1.0-_nu);
			_Dmat_sensitivities[p](3,3) = coefs[p] * (1.0-2.0*_nu)/2;	_Dmat_sensitivities[p](4,4) = coefs[p] * (1.0-2.0*_nu)/2;	_Dmat_sensitivities[p](5,5) = coefs[p] * (1.0-2.0*_nu)/2;
		}
	}

	// Now do poisson ratio sensitivity since its weird
	if (dim == 1)
		_Dmat_sensitivities[4](0,0) = 0.0;

	else if (dim == 2)
	{
		if (plane_strain)    // If the use set the plane strain boolean to true
		{
			double coef = _E / ((1.0+_nu) * (1.0-2.0*_nu));
			double dcoef = _E * (1.0+4.0*_nu) / pow(1.0-_nu-2.0*_nu*_nu, 2);
			_Dmat_sensitivities[4](0,0) = dcoef*(1.0-_nu) - coef;			_Dmat_sensitivities[4](0,1) = dcoef*_nu + coef;
			_Dmat_sensitivities[4](1,0) = dcoef*_nu + coef;				_Dmat_sensitivities[4](1,1) = dcoef*(1.0-_nu) - coef;
			_Dmat_sensitivities[4](2,2) = dcoef*(1.0-2.0*_nu)/2.0 - coef;
		}
		else	// Defaults to true
		{
			double coef = _E / (1.0-_nu*_nu);
			double dcoef = 2.0*_E*_nu / pow(1.0-_nu*_nu, 2);
			_Dmat_sensitivities[4](0,0) = dcoef;				_Dmat_sensitivities[4](0,1) = dcoef*_nu + coef;
			_Dmat_sensitivities[4](1,0) = dcoef*_nu + coef;	_Dmat_sensitivities[4](1,1) = dcoef;
			_Dmat_sensitivities[4](2,2) = dcoef*(1.0-_nu)/2.0 - 0.5*coef;	// Division by 2 is to account for the engineering strain we use
		}
	}

	else if (dim == 3)
	{
		double coef = _E / ((1.0+_nu) * (1.0-2.0*_nu));
		double dcoef = _E * (1.0+4.0*_nu) / pow(1.0-_nu-2.0*_nu*_nu, 2);
		_Dmat_sensitivities[4](0,0) = dcoef*(1.0-_nu) - coef;			_Dmat_sensitivities[4](0,1) = dcoef*_nu + coef;					_Dmat_sensitivities[4](0,2) = dcoef*_nu + coef;
		_Dmat_sensitivities[4](1,0) = dcoef*_nu + coef;					_Dmat_sensitivities[4](1,1) = dcoef*(1.0-_nu) - coef;			_Dmat_sensitivities[4](1,2) = dcoef*_nu + coef;
		_Dmat_sensitivities[4](2,0) = dcoef*_nu + coef;					_Dmat_sensitivities[4](2,1) = dcoef*_nu + coef;					_Dmat_sensitivities[4](2,2) = dcoef*(1.0-_nu) - coef;
		_Dmat_sensitivities[4](3,3) = dcoef*(1.0-2.0*_nu)/2.0 - coef;	_Dmat_sensitivities[4](4,4) = dcoef*(1.0-2.0*_nu)/2.0 - coef;	_Dmat_sensitivities[4](5,5) = dcoef*(1.0-2.0*_nu)/2.0 - coef;
	}

	_sensitivities_computed = true;
}




// A function to efficiently compute the stress state and avoid zero multiplies
void LinearElasticIsotropicMaterial::computeStress(std::vector<double>& stress, const DenseMatrix<double>& Dmat, const std::vector<double>& strain)
{
	if (Dmat.n_rows() == 1) // 1-dimensional
	{
		stress.resize(1);
		stress[0] = Dmat(0,0) * strain[0];
	}

	else if (Dmat.n_rows() == 3) // 2-dimensional
	{
		stress.resize(3);
		stress[0] = Dmat(0,0)*strain[0] + Dmat(0,1)*strain[1];
		stress[1] = Dmat(1,0)*strain[0] + Dmat(1,1)*strain[1];
		stress[2] = Dmat(2,2)*strain[2];
	}

	else if (Dmat.n_rows() == 6) // 3-dimensional
	{
		stress.resize(6);
		stress[0] = Dmat(0,0)*strain[0] + Dmat(0,1)*strain[1] + Dmat(0,2)*strain[2];
		stress[1] = Dmat(1,0)*strain[0] + Dmat(1,1)*strain[1] + Dmat(1,2)*strain[2];
		stress[2] = Dmat(2,0)*strain[0] + Dmat(2,1)*strain[1] + Dmat(2,2)*strain[2];
		stress[3] = Dmat(3,3)*strain[3];
		stress[4] = Dmat(4,4)*strain[4];
		stress[5] = Dmat(5,5)*strain[5];
	}

	else
		err_message("Invalid dimension of the D matrix for a linear elastic isotropic material!");
}



Material::output_params* LinearElasticIsotropicMaterial::Constitutive(Material::input_params& input)
{
	if(input.dim != _curr_dim)
	{
		computeLinearDmat(input);
		_curr_dim = input.dim;
	}
	if (_thermal_exp_set)
	{
		std::vector<double>& strain = input.strain;
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
		computeStress(_output.stress, _output.Dmat, strain);
	}
	else
		computeStress(_output.stress, _output.Dmat, input.strain);
	
	
	return &_output;
}





Material::sensitivity_output_params* LinearElasticIsotropicMaterial::SensitivityConstitutive(Material::sensitivity_input_params& input)
{
	// Get the internal variables and references to the output
	std::vector<double>& dstress_dd = _sensitivity_output.dstress_dd;

// 1a
	// Get the dimension of the problem and the input opening vector (int ntt coordinates hopefully)
	// ntt stands for normal-tangential-tangential
	id_type dim = input.dim;
	std::vector<double>& strain = input.strain;
	bool plane_strain = input.plane_strain;

	

	// Traction vector partial derivatives
	dstress_dd.clear();
	dstress_dd.resize(strain.size());

	// Partial derivative is only non-zero if the material I am taing the derivative wrt is the same as the current material I am in
	if (input.sensitivity_mat_name == _name) 
	{
		// Actually compute the partials of the Dmatrix
		if (!_sensitivities_computed)
			computeDmatSensitivities(dim, plane_strain);

		// Young's modulus sensitivity
		if (input.sensitivity_param_name == "E")
			computeStress(dstress_dd, _Dmat_sensitivities[0], strain);

		// Lambda sensitivity
		else if (input.sensitivity_param_name == "LAMBDA")
			computeStress(dstress_dd, _Dmat_sensitivities[1], strain);

		// Shear Modulus sensitivity
		else if (input.sensitivity_param_name == "MU")
			computeStress(dstress_dd, _Dmat_sensitivities[2], strain);

		// Bulk Modulus sensitivity
		else if (input.sensitivity_param_name == "K")
			computeStress(dstress_dd, _Dmat_sensitivities[3], strain);

		// Poisson's ratio sensitivity
		else if (input.sensitivity_param_name == "NU")
			computeStress(dstress_dd, _Dmat_sensitivities[4], strain);

		else
			err_message(input.sensitivity_param_name << " is an invalid parameter name for a linear elastic isotropic material!");
	}

	// Return a pointer to the output structure
	return &_sensitivity_output;
}












SensitivityMaterialParameter* LinearElasticIsotropicMaterial::getSensitivityParameter(std::string param_name)
{
	SensitivityMaterialParameter* ret = new SensitivityMaterialParameterStructural;
	ret->set_mat_name(_name);
	std::transform(param_name.begin(), param_name.end(), param_name.begin(), ::toupper); // Capitilize the name
	if (param_name=="E" || param_name=="YOUNG'S MODULUS" || param_name=="YOUNGS MODULUS")
		ret->set_param_name("E");
	else if (param_name=="NU" || param_name=="POISSON'S RATIO" || param_name=="POISSONS RATIO")
		ret->set_param_name("NU");
	else if (param_name=="LAMBDA" || param_name=="FIRST LAME CONSTANT")
		ret->set_param_name("LAMBDA");
	else if (param_name=="MU" || param_name=="G" || param_name=="SECOND LAME CONSTANT" || param_name=="SHEAR MODULUS")
		ret->set_param_name("MU");
	else if (param_name=="K" || param_name=="BULK MODULUS")
		ret->set_param_name("K");
	else if (param_name=="THERMAL_EXPANSION" || param_name=="THERMAL EXPANSION" || param_name=="THERMAL_SHINKAGE" || param_name=="THERMAL SHRINKAGE")
		ret->set_param_name("THERMAL_EXPANSION");
	else
	{
		delete ret;
		err_message(param_name << " is not a valid parameter name");
	}

	return ret;
}







void LinearElasticIsotropicMaterial::set_parameter(std::string name, double val)
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
double LinearElasticIsotropicMaterial::get_parameter(std::string name)
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



void LinearElasticIsotropicMaterial::set_E(double E)
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
double LinearElasticIsotropicMaterial::get_E()
{
	if(n_set>=2 || E_set)
		return _E;
	else
		err_message("Young's modulus could not be determined.");
}




void LinearElasticIsotropicMaterial::set_nu(double nu)
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
double LinearElasticIsotropicMaterial::get_nu()
{
	if(n_set>=2 || nu_set)
		return _nu;
	else
		err_message("Poisson's ratio could not be determined.");
}




void LinearElasticIsotropicMaterial::set_lambda(double lambda)
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
double LinearElasticIsotropicMaterial::get_lambda()
{
	if(n_set>=2 || lambda_set)
		return _lambda;
	else
		err_message("First Lame constant could not be determined.");
}




void LinearElasticIsotropicMaterial::set_mu(double mu)
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
double LinearElasticIsotropicMaterial::get_mu()
{
	if(n_set>=2 || mu_set)
		return _mu;
	else
		err_message("Shear modulus could not be determined.");
}




void LinearElasticIsotropicMaterial::set_K(double K)
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
double LinearElasticIsotropicMaterial::get_K()
{
	if(n_set>=2 || K_set)
		return _K;
	else
		err_message("Bulk modulus could not be determined.");
}


void LinearElasticIsotropicMaterial::set_thermal_exp(double val)
{
	_thermal_exp = val;
	_thermal_exp_set = true;
}
double LinearElasticIsotropicMaterial::get_thermal_exp()
{
	if (_thermal_exp_set)
		return _thermal_exp;
	else
		err_message("Thermal expansion ratio could not be determined.");
}