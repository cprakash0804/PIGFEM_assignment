/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "material_lt.h"
#include <cmath>
#include <algorithm>
#include <iostream>


Material* LinearThermalMaterial::allocate()
{
	return new LinearThermalMaterial;
}
// This function assumes the current material is set up correctly. Which it should't be possible for it not to be
void LinearThermalMaterial::copy(Material* other_mat)
{
	other_mat->set_name(_name);
	other_mat->set_id(_id);
	if(kappa_set)
		other_mat->set_parameter("kappa", _kappa);
}

LinearThermalMaterial::LinearThermalMaterial()
	: _kappa(0.0), kappa_set(false)
{
}


// FIXME: Need to be able to do plane stress solutions as well
void LinearThermalMaterial::computeLinearDmat(Material::input_params& input)
{
	if(!kappa_set)
		err_message("Please set a value for conductivity.");

	id_type dim = input.dim;

	_output.Dmat.clear();
	_output.Dmat.resize(dim,dim);
	for(id_type i=0; i<dim; ++i)
		_output.Dmat(i,i) = _kappa;
}



Material::output_params* LinearThermalMaterial::Constitutive(Material::input_params& input)
{
	if(input.dim != _curr_dim)
	{
		computeLinearDmat(input);
		_curr_dim = input.dim;
	}

	return &_output;
}



void LinearThermalMaterial::set_parameter(std::string name, double val)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="K" || name=="KAPPA" || name=="CONDUCTIVITY")
		set_kappa(val);
	else
		err_message(name << " is not a valid parameter name");
}
double LinearThermalMaterial::get_parameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="K" || name=="KAPPA" || name=="CONDUCTIVITY")
		return get_kappa();
	else
		err_message(name << " is not a valid parameter name");
}



void LinearThermalMaterial::set_kappa(double kappa)
{
	_kappa = kappa;
	kappa_set = true;
}
double LinearThermalMaterial::get_kappa()
{
	if(kappa_set)
		return _kappa;
	else
		err_message("Please set a value for kappa before attempting ot retrieve it.");
}
