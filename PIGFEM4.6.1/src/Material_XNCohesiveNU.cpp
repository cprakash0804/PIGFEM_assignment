/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#include "Material_XNCohesiveNU.h"
#include "SensitivityMaterialParameterStructural.h"
#include "Utilities.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>


Material* XNCohesiveNUMaterial::allocate()
{
	return new XNCohesiveNUMaterial;
}
// This function assumes the current material is set up correctly. Which it should't be possible for it not to be
void XNCohesiveNUMaterial::copy(Material* other_mat)
{
	other_mat->set_name(_name);
	other_mat->set_id(_id);
	if (_sn_set)
		other_mat->set_parameter("SIGMA_N", _sigma_cn);
	if (_dn_set)
		other_mat->set_parameter("DELTA_N", _delta_cn);
	if (_dt_set)
		other_mat->set_parameter("DELTA_T", _delta_ct);
	if (_q_set)
		other_mat->set_parameter("Q", sqrt(_q));
}


XNCohesiveNUMaterial::XNCohesiveNUMaterial()
	: _sigma_cn(0), _delta_cn(0), _delta_ct(0), _q(0),
	  _sn_set(false), _dn_set(false), _dt_set(false), _q_set(false), _ready(false)
{
}


void XNCohesiveNUMaterial::getDnDt(std::vector<double>& delta_ttn, double& Dn, double& Dt)
{
	// Replacing delta_ttn with a small number to prevent blowup at the beginning of loading
	double norm = 0;
	for(id_type i=0; i<delta_ttn.size(); ++i)
		norm += delta_ttn[i]*delta_ttn[i];
	if(norm < 1e-40) // Since this is norm squared
	{
		for (id_type i=0; i< delta_ttn.size(); ++i)
			delta_ttn[i] = 1e-20;
	}

	// Compute normal and transverse componenets of the opening
	Dn = delta_ttn[delta_ttn.size()-1];
	Dt = 0;
	if (delta_ttn.size() == 2)
		Dt = delta_ttn[0];
	else if (delta_ttn.size() == 3)
		Dt = sqrt(pow(delta_ttn[0], 2) + pow(delta_ttn[1], 2));
}


Material::output_params* XNCohesiveNUMaterial::Constitutive(Material::input_params& input)
{
	if(!_ready)
		err_message("Not all material parameters have been set for the Ortiz-Pandolfi Cohesive material.");

	// Get the internal variables and references to the output
	DenseMatrix<double>& Dmat = _output.Dmat;
	std::vector<double>& traction_ttn = _output.traction;

	// Get the dimension of the problem and the input opening vector (int ttn coordinates hopefully)
	// ttn stands for normal-tangential-tangential
	id_type dim = input.dim;
	std::vector<double>& delta_ttn = input.delta;

	// Get the normal and tangential compoennet of the opening
	double Dn, Dt;
	getDnDt(delta_ttn, Dn, Dt);

	// Uncoupled Form
	// -------------------------------------------------
	// Some temporary variables
	double shear_rat = Dt/_delta_ct;
	double shear_rat_exp = exp( -pow(shear_rat, 2) );
	double norm_rat = Dn/_delta_cn;
	double norm_rat_exp = exp( -1.0*norm_rat );
	double phi_n = _sigma_cn * _delta_cn * exp(1);

	// Scalar traction components
	double Tn = (phi_n/_delta_cn) * norm_rat_exp * norm_rat;
	double Tt = (2.0*phi_n*Dt/pow(_delta_ct,2)) * _q * shear_rat_exp;

	// Traction vector
	traction_ttn.resize(dim);
	traction_ttn[dim-1] = Tn;
	if (dim==2)
		traction_ttn[0] = Tt;
	else if (dim==3)
	{
		traction_ttn[0] = Tt * (delta_ttn[0]/Dt);
		traction_ttn[1] = Tt * (delta_ttn[1]/Dt);
	}

	// Tangent matrix
	if (dim==1)
	{
		double dTndDn = (phi_n/_delta_cn) * norm_rat_exp * ((1-norm_rat)/_delta_cn);
		Dmat.resize(1,1);
		Dmat(0,0) = dTndDn;
	}
	else if (dim==2)
	{
		double dTndDn = (phi_n/_delta_cn) * norm_rat_exp * ((1-norm_rat)/_delta_cn);
		double dTndDt = 0.0;
		double dTtdDn = 0.0;
		double dTtdDt = (2.0*phi_n/pow(_delta_ct,2)) * _q * (1-2.0*pow(shear_rat,2)) * shear_rat_exp;
		Dmat.resize(2,2);
		Dmat(0,0) = dTtdDt;
		Dmat(0,1) = dTtdDn;
		Dmat(1,0) = dTndDt;
		Dmat(1,1) = dTndDn;
	}
	else if (dim==3)
	{
		double dTndDn = (phi_n/_delta_cn) * norm_rat_exp * ((1-norm_rat)/_delta_cn);
		double dTndDt = 0.0;
		double dTtdDn = 0.0;
		double dTtdDt = (2.0*phi_n/pow(_delta_ct,2)) * _q * (1-2.0*pow(shear_rat,2)) * shear_rat_exp;
		Dmat.resize(3,3);
		Dmat(0,0) = dTtdDt * pow(delta_ttn[0]/Dt,2) + (Tt/Dt) * (1-pow(delta_ttn[0]/Dt,2));
		Dmat(0,1) = (dTtdDt - Tt/Dt) * (delta_ttn[0] * delta_ttn[1] / (Dt*Dt));
		Dmat(0,2) = dTtdDn * (delta_ttn[0]/Dt);
		Dmat(1,0) = (dTtdDt - Tt/Dt) * (delta_ttn[0] * delta_ttn[1] / (Dt*Dt));
		Dmat(1,1) = dTtdDt * pow(delta_ttn[1]/Dt,2) + (Tt/Dt) * (1-pow(delta_ttn[1]/Dt,2));
		Dmat(1,2) = dTtdDn * (delta_ttn[1]/Dt);
		Dmat(2,0) = dTndDt * (delta_ttn[0]/Dt);
		Dmat(2,1) = dTndDt * (delta_ttn[1]/Dt);
		Dmat(2,2) = dTndDn;
	}
	

	return &_output;
}




Material::sensitivity_output_params* XNCohesiveNUMaterial::SensitivityConstitutive(Material::sensitivity_input_params& input)
{
	if(!_ready)
		err_message("Not all material parameters have been set for the Continuum Damage Material.");

	// Get the internal variables and references to the output
	std::vector<double>& dtraction_ttn_dd = _sensitivity_output.dtraction_dd;

	// Get the dimension of the problem and the input opening vector
	// ttn stands for normal-tangential-tangential
	id_type dim = input.dim;
	std::vector<double>& delta_ttn = input.delta;

	// Traction vector partial derivatives
	dtraction_ttn_dd.clear();
	dtraction_ttn_dd.resize(delta_ttn.size());

	// Partial derivative is only non-zero if the material I am taing the derivative wrt is the same as the current material I am in
	if (input.sensitivity_mat_name == _name) 
	{
		// Get the normal and tangential compoennet of the opening
		double Dn, Dt;
		getDnDt(delta_ttn, Dn, Dt);

		double dTn(0.0), dTt(0.0);

		// Sigma_cn sensitivity
		if (input.sensitivity_param_name == "SIGMA_CN"){
			dTn = (Dn/_delta_cn) * exp(1.0-Dn/_delta_cn );
			dTt = (2.0*_delta_cn*_q*Dt/pow(_delta_ct,2)) * exp(1.0 - pow(Dt/_delta_ct, 2));
		}

		// Delta_cn sensitivity
		else if (input.sensitivity_param_name == "DELTA_CN"){
			dTn = (_sigma_cn*Dn/pow(_delta_cn,2)) * (Dn/_delta_cn-1.0) * exp(1.0-Dn/_delta_cn);
			dTt = (2.0*_sigma_cn*Dt*_q/pow(_delta_ct,2)) * exp(1.0-pow(Dt/_delta_ct,2));
		}

		// Delta_ct sensitivity
		else if (input.sensitivity_param_name == "DELTA_CT"){
			dTn = 0.0;
			dTt = (4.0*_sigma_cn*_delta_cn*_q*Dt/pow(_delta_ct,3)) * (pow(Dt/_delta_ct, 2)-1.0) * exp(1.0 - pow(Dt/_delta_ct, 2));
		}

		// Q sensitivity (WTF is going on here? It only matches FD w/o the 2.0)
		else if (input.sensitivity_param_name == "Q"){
			dTn = 0.0;
			// dTt = (2.0*_sigma_cn*_delta_cn*Dt/pow(_delta_ct,2)) * exp(1.0 - pow(Dt/_delta_ct, 2));
			dTt = (_sigma_cn*_delta_cn*Dt/pow(_delta_ct,2)) * exp(1.0 - pow(Dt/_delta_ct, 2));
		}

		// Combined delta_c sensitivity (if _delta_cn and _delta_ct are actually equal - this is checked in the sensitivity parameter constructor)
		else if (input.sensitivity_param_name == "DELTA_C"){
			dTn = (_sigma_cn*Dn/pow(_delta_cn,2)) * (Dn/_delta_cn-1.0) * exp(1.0-Dn/_delta_cn);
			dTt = (2.0*_sigma_cn*_q*Dt/pow(_delta_ct,2)) * (2.0*pow(Dt/_delta_ct, 2)-1.0) * exp(1.0 - pow(Dt/_delta_ct, 2));
		}

		else
			err_message("Unknown sensitivity parameter for the Xu-Needleman cohesive material!");

		// Update the sensitivity of the traction vector
		dtraction_ttn_dd[dtraction_ttn_dd.size()-1] = dTn;
		if (dim==2)
			dtraction_ttn_dd[0] = dTt;
		else if (dim==3)
		{
			dtraction_ttn_dd[0] = dTt * (delta_ttn[0]/Dt);
			dtraction_ttn_dd[1] = dTt * (delta_ttn[1]/Dt);
		}
	}

	return &_sensitivity_output;
}




SensitivityMaterialParameter* XNCohesiveNUMaterial::getSensitivityParameter(std::string param_name)
{
	SensitivityMaterialParameter* ret = new SensitivityMaterialParameterStructural;
	ret->set_mat_name(_name);
	std::transform(param_name.begin(), param_name.end(), param_name.begin(), ::toupper); // Capitilize the name
	if(param_name=="SIGMA" || param_name=="SIGMA_C" || param_name=="SIGMA C" || param_name=="SIGMA_N" || param_name=="SIGMA_CN" || param_name=="SIGMA CN")
		ret->set_param_name("SIGMA_CN");
	else if(param_name=="DELTA_N" || param_name=="DELTA_CN" || param_name=="DELTA CN" || param_name=="DELTA_CRITICAL_N" || param_name=="DELTA CRITICAL N")
		ret->set_param_name("DELTA_CN");
	else if(param_name=="Q")
		ret->set_param_name("Q");
	else if(param_name=="DELTA_T" || param_name=="DELTA_CT" || param_name=="DELTA CT" || param_name=="DELTA_CRITICAL_T" || param_name=="DELTA CRITICAL T")
		ret->set_param_name("DELTA_CT");
	else if(param_name=="DELTA" || param_name=="DELTA_C" || param_name=="DELTA C" || param_name=="DELTA_CRITICAL" || param_name=="DELTA CRITICAL")
	{
		if (std::fabs(_delta_cn/_delta_ct - 1.0) < 1e-14)
			ret->set_param_name("DELTA_C");
		else
			err_message("To take combined delta_c sensitivity, delta_cn and delta_c must be equal!");
	}
	else
	{
		delete ret;
		err_message(param_name << " is not a valid parameter name");
	}

	return ret;
}





double XNCohesiveNUMaterial::getNormalizedCohesiveFailure(std::vector<double>& delta_ttn)
{
	double Dn, Dt;
	getDnDt(delta_ttn, Dn, Dt);
	Dn = Dn / _delta_cn;
	Dt = Dt / _delta_ct;
	if (Dn > Dt)
		return Dn;
	else
		return Dt;
}
id_type XNCohesiveNUMaterial::get_cohesive_damage_zone(std::vector<double>& delta)
{
	double norm = getNormalizedCohesiveFailure( delta );

	if (norm <= 0)
		return 0;
	else if (norm <= 1)
		return 1;
	else if (norm <= 5.0) // Kinda arbitrary selection
		return 2;
	else
		return 3;
}




void XNCohesiveNUMaterial::set_parameter(std::string name, double val)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="SIGMA_N" || name=="SIGMA_C" || name=="SIGMA_CN" || name=="SIGMA CN" || name=="SIGMA_N_CRITICAL" || name=="SIGMA_N CRITICAL")
		set_scn(val);
	else if(name=="DELTA_N" || name=="DELTA_CN" || name=="DELTA CN" || name=="DELTA_N_CRITICAL" || name=="DELTA_N CRITICAL")
		set_dcn(val);
	else if(name=="DELTA_T" || name=="DELTA_CT" || name=="DELTA CT" || name=="DELTA_T_CRITICAL" || name=="DELTA_T CRITICAL")
		set_dct(val);
	else if(name=="DELTA" || name=="DELTA_C" || name=="DELTA C" || name=="DELTA_CRITICAL" || name=="DELTA CRITICAL"){ // Option for equal critical openings
		set_dcn(val);
		set_dct(val);
	}
	else if(name=="Q")
		set_q(val);
	else
		err_message(name << " is not a valid parameter name");
}
double XNCohesiveNUMaterial::get_parameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="SIGMA_N" || name=="SIGMA_CN" || name=="SIGMA CN" || name=="SIGMA_N_CRITICAL" || name=="SIGMA_N CRITICAL")
		return get_scn();
	else if(name=="DELTA_N" || name=="DELTA_CN" || name=="DELTA CN" || name=="DELTA_N_CRITICAL" || name=="DELTA_N CRITICAL")
		return get_dcn();
	else if(name=="DELTA_T" || name=="DELTA_CT" || name=="DELTA CT" || name=="DELTA_T_CRITICAL" || name=="DELTA_T CRITICAL")
		return get_dct();
	else if(name=="DELTA" || name=="DELTA_C" || name=="DELTA C" || name=="DELTA_CRITICAL" || name=="DELTA CRITICAL") // Default to the normal opening
		return get_dcn();
	else if(name=="Q")
		return get_q();
	else
		err_message(name << " is not a valid parameter name");
}




void XNCohesiveNUMaterial::set_scn(double scn)
{
	_sigma_cn = scn;
	_sn_set = true;
	if (_sn_set && _dn_set && _dn_set && _q_set)
		_ready = true;
}
double XNCohesiveNUMaterial::get_scn()
{
	if(_sn_set)
		return _sigma_cn;
	else
		err_message("Sigma critical normal could not be determined.");
}




void XNCohesiveNUMaterial::set_dcn(double dcn)
{
	_delta_cn = dcn;
	_dn_set = true;
	if (_sn_set && _dn_set && _dn_set && _q_set)
		_ready = true;
}
double XNCohesiveNUMaterial::get_dcn()
{
	if(_dn_set)
		return _delta_cn;
	else
		err_message("Delta critical normal could not be determined.");
}




void XNCohesiveNUMaterial::set_dct(double dct)
{
	_delta_ct = dct;
	_dt_set = true;
	if (_sn_set && _dn_set && _dn_set && _q_set)
		_ready = true;
}
double XNCohesiveNUMaterial::get_dct()
{
	if(_dt_set)
		return _delta_ct;
	else
		err_message("Delta critical tangential could not be determined.");
}




void XNCohesiveNUMaterial::set_q(double q)
{
	_q = q;
	_q_set = true;
	if (_sn_set && _dn_set && _dn_set && _q_set)
		_ready = true;
}
double XNCohesiveNUMaterial::get_q()
{
	if(_q_set || _ready)
		return _q;
	else
		err_message("q could not be determined.");
}