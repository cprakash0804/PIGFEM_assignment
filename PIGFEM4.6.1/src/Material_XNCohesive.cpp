/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#include "Material_XNCohesive.h"
#include "Utilities.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>



Material* XNCohesiveMaterial::allocate()
{
	return new XNCohesiveMaterial;
}


std::vector<std::string> XNCohesiveMaterial::internal_vars_name()
{
	std::vector<std::string> vec = {"norm_max_opening", "shear_max_opening"};
	return vec;
}

std::vector<bool> XNCohesiveMaterial::internal_vars_print()
{
	std::vector<bool> vec = {false, false};
	return vec;
}



Material::output_params* XNCohesiveMaterial::Constitutive(Material::input_params& input)
{
	if(!_ready)
		err_message("Not all material parameters have been set for the Ortiz-Pandolfi Cohesive material.");

	// Get the internal variables and references to the output
	double& Dn_max = (*input.internal_vars)[0];
	double& Dt_max = (*input.internal_vars)[1];
	DenseMatrix<double>& Dmat = _output.Dmat;
	std::vector<double>& traction_ttn = _output.traction;\

	// Get the dimension of the problem and the input opening vector (int ttn coordinates hopefully)
	// ttn stands for tangential-tangential-normal
	id_type dim = input.dim;
	std::vector<double>& delta_ttn = input.delta;

	// Get the normal and tangential compoennet of the opening
	double Dn, Dt;
	getDnDt(delta_ttn, Dn, Dt);

	// Uncoupled Form
	// -------------------------------------------------
	// Some temporary variables
	double phi_n = _sigma_cn * _delta_cn * exp(1);

	// Scalar traction components
	double Tn, dTndDn, dTndDt, Tt, dTtdDn, dTtdDt;
	if (Dn >= Dn_max || Dn <= 0) // Following the exponential curve
	{
		double norm_rat = Dn/_delta_cn;
		double norm_rat_exp = exp( -1.0*norm_rat );
		Tn = (phi_n/_delta_cn) * norm_rat_exp * norm_rat;
		dTndDn = (phi_n/_delta_cn) * norm_rat_exp * ((1-norm_rat)/_delta_cn);
		dTndDt = 0.0;
		Dn_max = Dn; // Update the ISV
	}
	else // Linear Unloading regime
	{
		double norm_rat_max = Dn_max/_delta_cn;
		double norm_rat_exp_max = exp( -1.0*norm_rat_max );
		double Tn_max = (phi_n/_delta_cn) * norm_rat_exp_max * norm_rat_max;
		Tn = (Dn/Dn_max) * Tn_max;
		dTndDn = Tn_max / Dn_max;
		dTndDt = 0.0;
	}
	if (std::fabs(Dt) >= Dt_max) // Following the exponential curve
	{
		double shear_rat = Dt/_delta_ct;
		double shear_rat_exp = exp( -pow(shear_rat, 2) );
		Tt = (2.0*phi_n*_q/_delta_ct) * shear_rat_exp * shear_rat;
		dTtdDn = 0.0;
		dTtdDt = (2.0*phi_n/pow(_delta_ct,2)) * _q * (1-2.0*pow(shear_rat,2)) * shear_rat_exp;
		Dt_max = std::fabs(Dt); // Update the ISV
	}
	else // Linear Unloading regime
	{
		double shear_rat_max = Dt_max/_delta_ct;
		double shear_rat_exp_max = exp( -pow(shear_rat_max, 2) );
		double Tt_max = (2.0*phi_n*_q/_delta_ct) * shear_rat_exp_max * shear_rat_max;
		Tt = (Dt/Dt_max) * Tt_max;
		dTtdDn = 0.0;
		dTtdDt = Tt_max / Dt_max;
	}


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
		Dmat.resize(1,1);
		Dmat(0,0) = dTndDn;
	}
	else if (dim==2)
	{
		Dmat.resize(2,2);
		Dmat(0,0) = dTtdDt;
		Dmat(0,1) = dTtdDn;
		Dmat(1,0) = dTndDt;
		Dmat(1,1) = dTndDn;
	}
	else if (dim==3)
	{
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







Material::sensitivity_output_params* XNCohesiveMaterial::SensitivityConstitutive(Material::sensitivity_input_params& input)
{
	if(!_ready)
		err_message("Not all material parameters have been set for the Continuum Damage Material.");

	// Get the internal variables and references to the output
	double& Dn_max = (*input.internal_vars)[0];
	double& Dt_max = (*input.internal_vars)[1];
	id_type param_id = input.parameter_id;
	double& dDn_max_dd = (*input.internal_vars_sensitivity)[param_id * n_internal_vars()];
	double& dDt_max_dd = (*input.internal_vars_sensitivity)[param_id * n_internal_vars() + 1];
	std::vector<double>& dtraction_ttn_dd = _sensitivity_output.dtraction_dd;

	// Get the dimension of the problem and the input opening vector
	// ttn stands for normal-tangential-tangential
	id_type dim = input.dim;
	std::vector<double>& delta_ttn = input.delta;

	// Get the normal and tangential compoennet of the opening
	double Dn, Dt;
	getDnDt(delta_ttn, Dn, Dt);

	// Traction vector partial derivatives
	dtraction_ttn_dd.clear();
	dtraction_ttn_dd.resize(delta_ttn.size());
	double phi_n = _sigma_cn * _delta_cn * exp(1);

	// Partial derivative is only non-zero if the material I am taing the derivative wrt is the same as the current material I am in
	if (input.sensitivity_mat_name == _name) 
	{
		double dTn(0.0), dTt(0.0);

		// Sigma_cn sensitivity
		if (input.sensitivity_param_name == "SIGMA_CN"){
			if (Dn >= Dn_max || Dn <= 0)
				dTn = (Dn/_delta_cn) * exp(1.0-Dn/_delta_cn );
			else
				dTn = (Dn/_delta_cn) * exp(1.0-Dn_max/_delta_cn );
			if (std::fabs(Dt) >= Dt_max)
				dTt = (2.0*_delta_cn*_q*Dt/pow(_delta_ct,2)) * exp(1.0 - pow(Dt/_delta_ct, 2));
			else
				dTt = (2.0*_delta_cn*_q*Dt/pow(_delta_ct,2)) * exp(1.0 - pow(Dt_max/_delta_ct, 2));
		}

		// Delta_cn sensitivity
		else if (input.sensitivity_param_name == "DELTA_CN"){
			if (Dn >= Dn_max || Dn <= 0)
				dTn = (_sigma_cn*Dn/pow(_delta_cn,2)) * (Dn/_delta_cn-1.0) * exp(1.0-Dn/_delta_cn);
			else
				dTn = (_sigma_cn*Dn/pow(_delta_cn,2)) * (Dn_max/_delta_cn-1.0) * exp(1.0-Dn_max/_delta_cn);
			if (std::fabs(Dt) >= Dt_max)
				dTt = (2.0*_sigma_cn*Dt*_q/pow(_delta_ct,2)) * exp(1.0-pow(Dt/_delta_ct,2));
			else
				dTt = (2.0*_sigma_cn*Dt*_q/pow(_delta_ct,2)) * exp(1.0-pow(Dt_max/_delta_ct, 2));
		}

		// Delta_ct sensitivity
		else if (input.sensitivity_param_name == "DELTA_CT"){
			if (Dn >= Dn_max || Dn <= 0)
				dTn = 0.0;
			else
				dTn = 0.0;
			if (std::fabs(Dt) >= Dt_max)
				dTt = (4.0*_sigma_cn*_delta_cn*_q*Dt/pow(_delta_ct,3)) * (pow(Dt/_delta_ct, 2)-1.0) * exp(1.0 - pow(Dt/_delta_ct, 2));
			else
				dTt = (4.0*_sigma_cn*_delta_cn*_q*Dt/pow(_delta_ct,3)) * (pow(Dt_max/_delta_ct, 2)-1.0) * exp(1.0 - pow(Dt_max/_delta_ct, 2));
		}

		// Q sensitivity (WTF is going on here? It only matches FD w/o the 2.0)
		else if (input.sensitivity_param_name == "Q"){
			if (Dn >= Dn_max || Dn <= 0)
				dTn = 0.0;
			else
				dTn = 0.0;
			// if (std::fabs(Dt) >= Dt_max)
			// 	dTt = (2.0*_sigma_cn*_delta_cn*Dt/pow(_delta_ct,2)) * exp(1.0-pow(Dt/_delta_ct, 2));
			// else
			// 	dTt = (2.0*_sigma_cn*_delta_cn*Dt/pow(_delta_ct,2)) * exp(1.0-pow(Dt_max/_delta_ct, 2));
			if (std::fabs(Dt) >= Dt_max)
				dTt = (_sigma_cn*_delta_cn*Dt/pow(_delta_ct,2)) * exp(1.0-pow(Dt/_delta_ct, 2));
			else
				dTt = (_sigma_cn*_delta_cn*Dt/pow(_delta_ct,2)) * exp(1.0-pow(Dt_max/_delta_ct, 2));
		}

		// Combined delta_c sensitivity (if _delta_cn and _delta_ct are actually equal - this is checked in the sensitivity parameter constructor)
		else if (input.sensitivity_param_name == "DELTA_C"){
			if (Dn >= Dn_max || Dn <= 0)
				dTn = (_sigma_cn*Dn/pow(_delta_cn,2)) * (Dn/_delta_cn-1.0) * exp(1.0-Dn/_delta_cn);
			else
				dTn = (_sigma_cn*Dn/pow(_delta_cn,2)) * (Dn_max/_delta_cn-1.0) * exp(1.0-Dn_max/_delta_cn);
			if (std::fabs(Dt) >= Dt_max)
				dTt = (2.0*_sigma_cn*_q*Dt/pow(_delta_ct,2)) * (2.0*pow(Dt/_delta_ct, 2)-1.0) * exp(1.0 - pow(Dt/_delta_ct, 2));
			else
				dTt = (2.0*_sigma_cn*_q*Dt/pow(_delta_ct,2)) * (2.0*pow(Dt_max/_delta_ct, 2)-1.0) * exp(1.0 - pow(Dt_max/_delta_ct, 2));
		}

		// Tn = (_sigma_cn * Dn/_delta_cn) * exp(1.0-Dn_max/_delta_cn);
		// Tt = (2.0*_sigma_cn*_delta_cn*_q*Dt/pow(_delta_ct,2)) * exp(1.0-pow(Dt_max/_delta_ct, 2));

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

	// Total derivative wrt the internal state variable term is present regardless of whether or not
	// I'm taking sensitivities wrt this material's parameters
	if (0.0 < Dn && Dn < Dn_max) // unloading in the normal direction
	{
		double norm_rat_max = Dn_max/_delta_cn;
		double norm_rat_exp_max = exp( -1.0*norm_rat_max );
		double Tn_max = (phi_n/_delta_cn) * norm_rat_exp_max * norm_rat_max;
		double Tn = (Dn/Dn_max) * Tn_max;
		double dTn_dDn_max = Tn * (-1.0/_delta_cn);
		dtraction_ttn_dd[dtraction_ttn_dd.size()-1] += dTn_dDn_max * dDn_max_dd;
	}
	if (std::fabs(Dt) < Dt_max) // Unloading in the tangential direction
	{
		double shear_rat_max = Dt_max/_delta_ct;
		double shear_rat_exp_max = exp( -pow(shear_rat_max, 2) );
		// double Tt_max = (2.0*phi_n*_q/_delta_ct) * shear_rat_exp_max * shear_rat_max;
		// double Tt = (Dt/Dt_max) * Tt_max;
		// double dTt_dDt_max = Tt * (-2.0*Dt_max / pow(_delta_ct,2));
		double dTt_dDt_max = (2.0*phi_n*_q*Dt/pow(_delta_ct,2)) * shear_rat_exp_max * (-2.0*Dt_max/pow(_delta_ct,2));
		double dTt_dd = dTt_dDt_max * dDt_max_dd;

		if (dim==2)
			dtraction_ttn_dd[0] += dTt_dd;
		else if (dim==3)
		{
			dtraction_ttn_dd[0] += dTt_dd * (delta_ttn[0]/Dt);
			dtraction_ttn_dd[1] += dTt_dd * (delta_ttn[1]/Dt);
		}
	}

	return &_sensitivity_output;
}



void XNCohesiveMaterial::updateSensitivityISVs(Material::sensitivity_input_params& input)
{
	if(!_ready)
		err_message("Not all material parameters have been set for the Continuum Damage Material.");

	// Get the internal variables and references to the output
	double& Dn_max = (*input.internal_vars)[0];
	double& Dt_max = (*input.internal_vars)[1];
	id_type param_id = input.parameter_id;
	double& dDn_max_dd = (*input.internal_vars_sensitivity)[param_id * n_internal_vars()];
	double& dDt_max_dd = (*input.internal_vars_sensitivity)[param_id * n_internal_vars() + 1];

// 1a
	// Get the dimension of the problem and the input opening vector (int ttn coordinates hopefully)
	// ttn stands for tangential-tangential-normal
	id_type dim = input.dim;
	std::vector<double>& delta_ttn = input.delta;
	std::vector<double>& ddelta_ttn_dd = input.ddelta_dd;

	// Get the normal and tangential compoennet of the opening
	double Dn, Dt;
	getDnDt(delta_ttn, Dn, Dt);

	// If we are continuing in the loading regime in either direction we must update the ISV derivatives
	if (Dn >= Dn_max)
		dDn_max_dd = ddelta_ttn_dd[ddelta_ttn_dd.size() - 1]; // Very simple update

	if (std::fabs(Dt) >= Dt_max)
	{
		if (dim==2)
			dDt_max_dd = ddelta_ttn_dd[0] * std::fabs(Dt)/Dt;	// ddelta_t_dd * sign(Dt)
		else if (dim==3)
			dDt_max_dd = (1.0/Dt) * (delta_ttn[0]*ddelta_ttn_dd[0] + delta_ttn[1]*ddelta_ttn_dd[1]);
	}
}
