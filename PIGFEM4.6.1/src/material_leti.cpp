/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "material_leti.h"
#include "Utilities.h"
#include <cmath>
#include <algorithm>
#include <iostream>


Material* LinearElasticTransverselyIsotropicMaterial::allocate()
{
	return new LinearElasticTransverselyIsotropicMaterial;
}
// This function assumes the current material is set up correctly. Which it should't be possible for it not to be
void LinearElasticTransverselyIsotropicMaterial::copy(Material* other_mat)
{
	other_mat->set_name(_name);
	other_mat->set_id(_id);
	if(E1_set)
		other_mat->set_parameter("E1", _E1);
	if(E2_set)
		other_mat->set_parameter("E2", _E2);
	if(G12_set)
		other_mat->set_parameter("G12", _G12);
	if(G23_set)
		other_mat->set_parameter("G23", _G23);
	if(nu12_set)
		other_mat->set_parameter("nu12", _nu12);
	if (_thermal_exp1_set)
		other_mat->set_parameter("THERMAL_EXPANSION1", _thermal_exp1);
	if (_thermal_exp2_set)
		other_mat->set_parameter("THERMAL_EXPANSION2", _thermal_exp2);
}


LinearElasticTransverselyIsotropicMaterial::LinearElasticTransverselyIsotropicMaterial()
	: _alpha(0.0), _beta(0.0), n_set(0), E1_set(false), E2_set(false), G12_set(false), G23_set(false), nu12_set(false),
	  _thermal_exp1(0.0), _thermal_exp1_set(false), _thermal_exp2(0.0), _thermal_exp2_set(false)
{
}




void LinearElasticTransverselyIsotropicMaterial::computeLinearDmat(Material::input_params& input)
{
	if(n_set<5)
		err_message("Please set all five parameters for a linear elastic transversely isotropic material.");

	if (_alpha != 0.0 || _beta != 0.0)
		err_message("Rotation of a transversely isotropic material is currently not supported");

	id_type dim = input.dim;
	bool plane_strain = input.plane_strain;

	// Storage of the dependant parameters
	double _E3 = _E2;
	double _G13 = _G12;
	double _nu13 = _nu12;
	double _nu23 = _E2 / (2.0*_G23) - 1.0; // Nt sure if its E2 or E3 but it doesn't matter for this case

	// Computation of some stiffness constants
	double delta = 1.0 - 2.0*(_E3/_E1)*_nu12*_nu13*_nu23 - (_E3/_E1)*_nu13*_nu13 - 
					(_E3/_E2)*_nu23*_nu23 - (_E2/_E1)*_nu12*_nu12;
	double C11 = _E1 * (1.0-(_E3/_E2)*_nu23*_nu23) / delta;
	double C22 = _E2 * (1.0-(_E3/_E1)*_nu13*_nu13) / delta;
	double C33 = _E3 * (1.0-(_E2/_E1)*_nu12*_nu12) / delta;
	double C12 = (_E2*_nu12 + _E3*_nu13*_nu23) / delta;
	double C13 = _E3 * (_nu12*_nu23+_nu13) / delta;
	double C23 = (_E3/_E1) * (_E1*_nu23+_E2*_nu12*_nu13) / delta;
	double C44 = _G23;
	double C55 = _G13;
	double C66 = _G12;

	// Compute the constitutive matrix
	if(dim==1)
	{
		_output.Dmat.resize(1,1);
		_output.Dmat(0,0) = _E1;
	}
	else if(dim==2)
	{
		_output.Dmat.clear();
		_output.Dmat.resize(3, 3);
		if(plane_strain)    // If the use set the plane strain boolean to true
		{
			_output.Dmat(0,0) = C11;
			_output.Dmat(0,1) = C12;
			_output.Dmat(1,0) = C12;
			_output.Dmat(1,1) = C22;
			_output.Dmat(2,2) = C66;
		}
		else	// Defaults to true. NOTE: This will be much more complicated once we allow rotations
		{
			_output.Dmat(0,0) = C11 - C13/C33;
			_output.Dmat(0,1) = C12 - C23/C33;
			_output.Dmat(1,0) = C12 - C23/C33;
			_output.Dmat(1,1) = C22 - C23/C33;
			_output.Dmat(2,2) = C66;
		}
	}
	else if(dim==3)
	{
		_output.Dmat.clear();
		_output.Dmat.resize(6, 6);
		_output.Dmat(0,0) = C11;
		_output.Dmat(0,1) = C12;
		_output.Dmat(0,2) = C13;
		_output.Dmat(1,0) = C12;
		_output.Dmat(1,1) = C22;
		_output.Dmat(1,2) = C23;
		_output.Dmat(2,0) = C13;
		_output.Dmat(2,1) = C23;
		_output.Dmat(2,2) = C33;
		_output.Dmat(3,3) = C44;
		_output.Dmat(4,4) = C55;
		_output.Dmat(5,5) = C66;
	}
	else
		err_message("Dimensions for the D matrix must be less than or equal to 3.");
}



Material::output_params* LinearElasticTransverselyIsotropicMaterial::Constitutive(Material::input_params& input)
{
	if(input.dim != _curr_dim)
	{
		computeLinearDmat(input);
		_curr_dim = input.dim;
	}
	if (_thermal_exp1_set && _thermal_exp2_set)
	{
		std::vector<double>& strain = input.strain;
		std::vector<double> thermal_strain(strain.size(), 0.0);
		switch (input.dim)
		{
			case 1:
				thermal_strain[0] = -1.0 * input.temp_change * _thermal_exp1;
				break;
			case 2:
				thermal_strain[0] = -1.0 * input.temp_change * _thermal_exp1;
				thermal_strain[1] = -1.0 * input.temp_change * _thermal_exp2;
				break;
			case 3:
				thermal_strain[0] = -1.0 * input.temp_change * _thermal_exp1;
				thermal_strain[1] = -1.0 * input.temp_change * _thermal_exp2;
				thermal_strain[2] = -1.0 * input.temp_change * _thermal_exp2;
				break;
			default:
				err_message("Wrong dimension somehow...");
		}
		Utilities::_VecAXPY(1.0, thermal_strain, strain); // strain = 1.0*thermal_strain + strain
		_output.stress = _output.Dmat * strain;
	}
	else
		_output.stress = _output.Dmat * input.strain;

	return &_output;
}



void LinearElasticTransverselyIsotropicMaterial::set_parameter(std::string name, double val)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="E1" || name=="LONGITUDINAL MODULUS")
		set_E1(val);
	else if(name=="E2" || name=="TRANSVERSE MODULUS")
		set_E2(val);
	else if(name=="G12")
		set_G12(val);
	else if(name=="G23")
		set_G23(val);
	else if(name=="NU12")
		set_nu12(val);
	else if(name=="ALPHA")
		set_alpha(val);
	else if (name=="THERMAL_EXPANSION1" || name=="THERMAL EXPANSION1" || name=="THERMAL_SHINKAGE1" || name=="THERMAL SHRINKAGE1")
		set_thermal_exp(val, 1);
	else if (name=="THERMAL_EXPANSION2" || name=="THERMAL EXPANSION2" || name=="THERMAL_SHINKAGE2" || name=="THERMAL SHRINKAGE2")
		set_thermal_exp(val, 2);
	else
		err_message(name << " is not a valid parameter name");
}
double LinearElasticTransverselyIsotropicMaterial::get_parameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="E1" || name=="LONGITUDINAL MODULUS")
		return get_E1();
	else if(name=="E2" || name=="TRANSVERSE MODULUS")
		return get_E2();
	else if(name=="G12")
		return get_G12();
	else if(name=="G23")
		return get_G23();
	else if(name=="NU12")
		return get_nu12();
	else if(name=="ALPHA")
		return get_alpha();
	else if(name=="BETA")
		return get_beta();
	else if (name=="THERMAL_EXPANSION1" || name=="THERMAL EXPANSION1" || name=="THERMAL_SHINKAGE1" || name=="THERMAL SHRINKAGE1")
		return get_thermal_exp(1);
	else if (name=="THERMAL_EXPANSION2" || name=="THERMAL EXPANSION2" || name=="THERMAL_SHINKAGE2" || name=="THERMAL SHRINKAGE2")
		return get_thermal_exp(2);
	else
		err_message(name << " is not a valid parameter name");
}



void LinearElasticTransverselyIsotropicMaterial::set_E1(double E1)
{
	if(!E1_set)
		n_set++;
	
	_E1 = E1;
	E1_set = true;
}
double LinearElasticTransverselyIsotropicMaterial::get_E1()
{
	if(E1_set)
		return _E1;
	else
		err_message("Longintudinal young's modulus could not be determined.");
}


void LinearElasticTransverselyIsotropicMaterial::set_E2(double E2)
{
	if(!E2_set)
		n_set++;
	
	_E2 = E2;
	E2_set = true;
}
double LinearElasticTransverselyIsotropicMaterial::get_E2()
{
	if(E2_set)
		return _E2;
	else
		err_message("Transverse young's modulus could not be determined.");
}


void LinearElasticTransverselyIsotropicMaterial::set_G12(double G12)
{
	if(!G12_set)
		n_set++;
	
	_G12 = G12;
	G12_set = true;
}
double LinearElasticTransverselyIsotropicMaterial::get_G12()
{
	if(G12_set)
		return _G12;
	else
		err_message("G12 could not be determined.");
}


void LinearElasticTransverselyIsotropicMaterial::set_G23(double G23)
{
	if(!G23_set)
		n_set++;
	
	_G23 = G23;
	G23_set = true;
}
double LinearElasticTransverselyIsotropicMaterial::get_G23()
{
	if(G23_set)
		return _G23;
	else
		err_message("G23 could not be determined.");
}


void LinearElasticTransverselyIsotropicMaterial::set_nu12(double nu12)
{
	if(!nu12_set)
		n_set++;
	
	_nu12 = nu12;
	nu12_set = true;
}
double LinearElasticTransverselyIsotropicMaterial::get_nu12()
{
	if(nu12_set)
		return _nu12;
	else
		err_message("nu12 could not be determined.");
}


void LinearElasticTransverselyIsotropicMaterial::set_alpha(double alpha)
{
	_alpha = alpha;
}
double LinearElasticTransverselyIsotropicMaterial::get_alpha()
{
	return _alpha;
}


void LinearElasticTransverselyIsotropicMaterial::set_beta(double beta)
{
	_beta = beta;
}
double LinearElasticTransverselyIsotropicMaterial::get_beta()
{
	return _beta;
}


void LinearElasticTransverselyIsotropicMaterial::set_thermal_exp(double val, int dir)
{
	if (dir==1)
	{
		_thermal_exp1 = val;
		_thermal_exp1_set = true;
	}
	else if (dir==2)
	{
		_thermal_exp2 = val;
		_thermal_exp2_set = true;
	}
	else
		err_message("Invalid direction for thermal expansion ratio");
}
double LinearElasticTransverselyIsotropicMaterial::get_thermal_exp(int dir)
{
	if (dir==1)
	{
		if (_thermal_exp1_set)
			return _thermal_exp1;
		else
			err_message("Thermal expansion ratio could not be determined.");
	}
	else if (dir==2)
	{
		if (_thermal_exp2_set)
			return _thermal_exp2;
		else
			err_message("Thermal expansion ratio could not be determined.");
	}
	else
		err_message("Invalid direction for thermal expansion ratio");
}
