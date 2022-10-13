/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "material_cdm.h"
#include "Utilities.h"
#include "SensitivityMaterialParameterStructural.h"
#include "Utilities.h"
#include <cmath>
#include <algorithm>
#include <iostream>



Material* ContinuumDamageModelMaterial::allocate()
{
	return new ContinuumDamageModelMaterial;
}
// This function assumes the current material is set up correctly. Which it should't be possible for it not to be
void ContinuumDamageModelMaterial::copy(Material* other_mat)
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
	if(P1_set)
		other_mat->set_parameter("P1", _P1);
	if(P1_set)
		other_mat->set_parameter("P2", _P2);
	if(P1_set)
		other_mat->set_parameter("Yin", _Yin);
	if(P1_set)
		other_mat->set_parameter("mu_visc", _mu_visc);
	if (_thermal_exp_set)
		other_mat->set_parameter("THERMAL_EXPANSION", _thermal_exp);
	if(sigmay_t_set)
		other_mat->set_parameter("sigmay_t", _sigmay_t);
	if(sigmay_c_set)
		other_mat->set_parameter("sigmay_c", _sigmay_c);
	if(epsilon_t_set)
		other_mat->set_parameter("epsilon_t", _epsilon_t);
	if(epsilon_c_set)
		other_mat->set_parameter("epsilon_c", _epsilon_c);	
	if(H_ro_set)
		other_mat->set_parameter("H_ro", _H_ro);	
	if(n_ro_set)
		other_mat->set_parameter("n_ro", _n_ro);	
}


ContinuumDamageModelMaterial::ContinuumDamageModelMaterial()
	: n_set(0), E_set(false), nu_set(false), lambda_set(false), mu_set(false), K_set(false),
	  P1_set(false), P2_set(false), Yin_set(false), mu_visc_set(false), ready(false),
	  _thermal_exp(0.0), _thermal_exp_set(false), sigmay_c_set(false), sigmay_t_set(false),
	  epsilon_t_set(false), epsilon_c_set(false), H_ro_set(false), n_ro_set(false)
{
}



std::vector<std::string> ContinuumDamageModelMaterial::internal_vars_name()
{
	std::vector<std::string> vec = {"damage", "damage_threshold"};
	return vec;
}

std::vector<bool> ContinuumDamageModelMaterial::internal_vars_print()
{
	std::vector<bool> vec = {true, false};
	return vec;
}



void ContinuumDamageModelMaterial::computeLinearDmat(Material::input_params& input)
{
	if(n_set<2)
		err_message("Please set two independant parameters for a linear elastic isotropic material.");

	id_type dim = input.dim;
	bool plane_strain = input.plane_strain;

	// Compute the constitutive matrix
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


void ContinuumDamageModelMaterial::apply_damage(DenseMatrix<double>& D, double damage)
{
	if (damage != 0.0)
		D = _linear_Dmat * (1.0-damage);
	else
		D = _linear_Dmat;
	/*
	unsigned char nrows = D.n_rows();

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
	}*/
}
void ContinuumDamageModelMaterial::compute_stress(std::vector<double>& stress, const DenseMatrix<double>& D, const std::vector<double>& strain)
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




Material::output_params* ContinuumDamageModelMaterial::Constitutive(Material::input_params& input)
{
		
		if(!ready)
			err_message("Not all material parameters have been set for the Continuum Damage Material.");

		// Get the internal variables and references to the output
		bool penetration_iter = Utilities::penetration_store(0);
		//input.internal_vars->resize(8);
		double& curr_damage = (*input.internal_vars)[0];
		double& curr_damage_thresh = (*input.internal_vars)[1];
		double dt = input.delta_t;
		double dmu = dt*_mu_visc;
		DenseMatrix<double>& Dmat_bar = _output.Dmat;
		std::vector<double>& stress = _output.stress;
		double fake_curr_damage;
		double fake_curr_damage_thresh;
		
		
		// Saving stresses (I wish there was a better way!), only works for 2D at the moment.
		double& prev_stress_1 = (*input.internal_vars)[2];
		double& prev_stress_2 = (*input.internal_vars)[3];
		double& prev_stress_3 = (*input.internal_vars)[4];
		
		std::vector<double> prev_stress;
		prev_stress.push_back(prev_stress_1);
		prev_stress.push_back(prev_stress_2);
		prev_stress.push_back(prev_stress_3);
		//-----------------------------------------------------------------------------------------
		// Saving strains (I wish there was a better way!), only works for 2D at the moment.
		double& prev_strain_1 = (*input.internal_vars)[5];
		double& prev_strain_2 = (*input.internal_vars)[6];
		double& prev_strain_3 = (*input.internal_vars)[7];
		
		std::vector<double> prev_strain;
		prev_strain.push_back(prev_strain_1);
		prev_strain.push_back(prev_strain_2);
		prev_strain.push_back(prev_strain_3); 
		//------------------------------------------------------------------------------------------
		double& initial_hardening_von = (*input.internal_vars)[8];
		double& initial_damage_threshold = (*input.internal_vars)[9];
		
		// Compute the linear D matrix for this material
		Dmat_bar = _linear_Dmat;

		if(input.dim != _curr_dim)
			{
				computeLinearDmat(input);
				_curr_dim = input.dim;
				Dmat_bar = _linear_Dmat;
			}

		// Compute the strain
		std::vector<double> strain = input.strain;
	//	if (_thermal_exp_set)
		//	strain = Utilities::plus(strain, ComputeThermalStrain(input));

		// Compute the Strain energy (0.5*\epsilon'*D*\epsilon) (Known as Y_bar)
		std::vector<double> dYbar_depsilon = _linear_Dmat*prev_strain;
		double strain_energy = Utilities::dot(dYbar_depsilon, prev_strain);
		strain_energy *= 0.5;
		double tau_bar = sqrt(2*strain_energy);
		
		double J_2_strain = (pow(prev_strain[0]-prev_strain[1],2) + pow(prev_strain[0],2) + pow(prev_strain[1],2))/6 + pow(prev_strain[2],2);
		double I_1_strain = prev_strain[0]+prev_strain[1];
		double damage_criterion = 6*J_2_strain+2*I_1_strain*(_epsilon_c-_epsilon_t)-2*_epsilon_c*_epsilon_t;
		
		bool new_damage = false;
		//double tmp = 0.0;
		double damage_func = 0.0;

		// If there is new damage present in the material, do all of the updating
		if(damage_criterion>0)
		{
			
			if (penetration_iter==0)
				if (initial_damage_threshold == 0)
					initial_damage_threshold=tau_bar*0.999;
			
			double tau_bar_0 = initial_damage_threshold;
			//double tau_bar_0 = 0;
			
			damage_func = 1.0 - tau_bar_0*(1-_P1)/tau_bar-_P1*exp(_P2*(tau_bar_0-tau_bar));
			
			// Update the damage parameters
			if (penetration_iter==0)
			{
				
				double prev_damage = curr_damage;
				//curr_damage += (dmu/(1.0+dmu)) * (tau_bar-curr_damage_thresh)*dG_dYbar;
				curr_damage += (dmu/(1.0+dmu)) * (damage_func-curr_damage_thresh);
				if (curr_damage<prev_damage)
					curr_damage = prev_damage;
				
			//	if (curr_damage>.99)
			//		curr_damage = .99;
				
				double prev_damage_thresh = curr_damage_thresh;
				curr_damage_thresh = (curr_damage_thresh+dmu*damage_func) / (1.0+dmu);
				if (curr_damage_thresh<prev_damage_thresh)
					curr_damage_thresh=prev_damage_thresh;
			}
		}
		
		stress = prev_stress;
		//curr_damage =0;		
		double J_2 = (pow(stress[0]-stress[1],2) + pow(stress[0],2) + pow(stress[1],2))/6 + pow(stress[2],2);
		double I_1 = stress[0] + stress[1];
		double von_Mises = sqrt(0.5*(pow(stress[0]-stress[1],2) + pow(stress[0],2) + pow(stress[1],2))+3*pow(stress[2],2));
		double yield_criterion = 6*J_2+2*I_1*(_sigmay_c-_sigmay_t)-2*_sigmay_c*_sigmay_t;

		//Gradient of the flow rule
		std::vector<double> dfdsigma;
		dfdsigma.push_back(1/(2*von_Mises)*(2*stress[0]-stress[1]));
		dfdsigma.push_back(1/(2*von_Mises)*(2*stress[1]-stress[0]));
		dfdsigma.push_back(1/(2*von_Mises)*(6*stress[2]));
		
		//define strain increment vector
		std::vector<double> dstrain;
		dstrain.push_back(strain[0]-prev_strain[0]);
		dstrain.push_back(strain[1]-prev_strain[1]);	
		dstrain.push_back(strain[2]-prev_strain[2]);
		
		//define stress increment vector
		std::vector<double> dstress;
		if (yield_criterion>0)
		{
			if (penetration_iter==0)
				if (initial_hardening_von==0)
					initial_hardening_von = von_Mises;
			
			_H_ro *= pow((initial_hardening_von/von_Mises),_n_ro);
			//_H_ro=0;
			
			DenseMatrix<double> temp_tensor = Dmat_bar;
			for(id_type i=0; i<dfdsigma.size(); ++i)
				for(id_type j=0; j<dfdsigma.size(); ++j)
					temp_tensor(i,j) = dfdsigma[i]*dfdsigma[j];


			std::vector<double> temp_vector = Dmat_bar*dfdsigma;	
			
			double temp_scalar=0;			
			for(id_type i=0; i<dfdsigma.size(); ++i)
					temp_scalar += dfdsigma[i]*temp_vector[i];	
			
			DenseMatrix<double> Dmat_bar_increment = (Dmat_bar*temp_tensor)*Dmat_bar;
			
			
			//std::cout <<von_Mises<<std::endl;
			for(id_type i=0; i<dfdsigma.size(); ++i)
				for(id_type j=0; j<dfdsigma.size(); ++j)			
					Dmat_bar(i,j) = Dmat_bar(i,j) + Dmat_bar_increment(i,j)*(-1/(_H_ro+temp_scalar));
				
				
			//undamaged stress 
			dstress = Dmat_bar*dstrain;
			for(id_type i=0; i<stress.size(); ++i)
				stress[i] = prev_stress[i] + dstress[i];	
				//stress[i] = stress[i] + dstress[i];				
			//update the state variable
			if (penetration_iter==0)
			{
				prev_stress_1 = stress[0];
				prev_stress_2 = stress[1];
				prev_stress_3 = stress[2];
			}
			//apply the damage
			for(id_type i=0; i<stress.size(); ++i)
				stress[i] = (1.0-curr_damage)*stress[i];
			Dmat_bar = (1.0-curr_damage)*Dmat_bar;
			//apply_damage(Dmat_bar,curr_damage);
		}else{
			//undamaged stress
			
			dstress = Dmat_bar*dstrain;
			//stress = Dmat_bar*strain;
			for(id_type i=0; i<stress.size(); ++i)
				stress[i] = prev_stress[i] + dstress[i];
				//stress[i] = stress[i] + dstress[i];
			//update the state variable
			if (penetration_iter==0)
			{
				prev_stress_1 = stress[0];
				prev_stress_2 = stress[1];
				prev_stress_3 = stress[2];
			}
			// Compute he current D matrix
			for(id_type i=0; i<stress.size(); ++i)
				stress[i] = (1.0-curr_damage)*stress[i];
			Dmat_bar = (1.0-curr_damage)*Dmat_bar;
			//apply_damage(Dmat_bar,curr_damage);

		}
		//updade the state variables storing strains
		if (penetration_iter==0)
		{
			prev_strain_1 = strain [0];
			prev_strain_2 = strain [1];
			prev_strain_3 = strain [2]; 
		}
		

		// Return a pointer to the output structure
		return &_output;
}






Material::sensitivity_output_params* ContinuumDamageModelMaterial::SensitivityConstitutive(Material::sensitivity_input_params& input)
{
	err_message("Sensitivity analysis is not implemented for a continuum damage material yet!");

	// Return a pointer to the output structure
	return &_sensitivity_output;
}












SensitivityMaterialParameter* ContinuumDamageModelMaterial::getSensitivityParameter(std::string param_name)
{
	SensitivityMaterialParameter* ret = new SensitivityMaterialParameterStructural;
	ret->set_mat_name(_name);
	std::transform(param_name.begin(), param_name.end(), param_name.begin(), ::toupper); // Capitilize the name
	
	err_message("Sensitivity analysis is not implemented for a continuum damae material yet!");

	return ret;
}


















void ContinuumDamageModelMaterial::updateSensitivityISVs(Material::sensitivity_input_params& input)
{
	// if (!ready)
	// 	err_message("Not all material parameters have been set for the Continuum Damage Material.");

	// // Get the internal variables and references to the output
	// double& curr_damage = (*input.internal_vars)[0];
	// double& curr_damage_thresh = (*input.internal_vars)[1];
	// double& ddamage_dd = (*input.internal_vars_sensitivity)[param_id * n_internal_vars()];
	// double& dthresh_dd = (*input.internal_vars_sensitivity)[param_id * n_internal_vars() + 1];
	// double dt = input.delta_t;
	// double dmu = dt*_mu_visc;

	// // Get the dimension of the problem
	// id_type dim = input.dim;
	// std::vector<double>& strain = input.strain;

	// // Compute the linear D matrix for this material
	// if(dim != _curr_dim)
	// {
	// 	computeLinearDmat(input);
	// 	_curr_dim = dim;
	// 	Dmat_bar = _linear_Dmat;
	// 	stress.resize(3);
	// }

	// // Compute the Strain energy (0.5*\epsilon'*D*\epsilon) (Known as Y_bar)
	// std::vector<double> dYbar_depsilon = strain*_linear_Dmat;
	// double strain_energy = Utilities::dot(dYbar_depsilon, strain);
	// strain_energy *= 0.5;

	// // Check the criteria to see if there is new damage present
	// // First Check: If the strain energy is greater than the initial threshold (Yin)
	// // Second Check: If the damage function exceeds the curent damage threshold (internal variable)
	// bool new_damage = false;
	// double tmp = 0.0;
	// double damage_func = 0.0;
	// if(strain_energy >= _Yin)
	// {
	// 	tmp = (strain_energy - _Yin)/(_P1*_Yin);    // This is purely to make the code more readable. This value is used several times. Could probably some up with a better name for it I guess..
	// 	damage_func = 1.0 - exp(-1.0*pow(tmp, _P2));
	// 	if(damage_func > curr_damage_thresh)
	// 		new_damage = true;
	// 	else
	// 		new_damage = false;
	// }
	// else
	// 	new_damage = false;
}







void ContinuumDamageModelMaterial::set_parameter(std::string name, double val)
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
	else if(name=="P1")
		set_P1(val);
	else if(name=="P2")
		set_P2(val);
	else if(name=="YIN" || name=="DAMAGETHRESHOLD" || name=="DAMAGE THRESHOLD")
		set_Yin(val);
	else if(name=="MU_VISC" || name=="MU_VISCOUS" || name=="MU VISC" || name=="MU VISCOUS")
		set_mu_visc(val);
	else if (name=="THERMAL_EXPANSION" || name=="THERMAL EXPANSION" || name=="THERMAL_SHINKAGE" || name=="THERMAL SHRINKAGE")
		set_thermal_exp(val);
	else if (name=="SIGMAY_C")
		set_sigmay_c(val);
	else if (name=="SIGMAY_T")
		set_sigmay_t(val);
	else if (name=="EPSILON_T")
		set_epsilon_t(val);	
	else if (name=="EPSILON_C")
		set_epsilon_c(val);	
	else if (name=="H_RO")
		set_H_ro(val);	
	else if (name=="N_RO")
		set_n_ro(val);	
	else
		err_message(name << " is not a valid parameter name");
}
double ContinuumDamageModelMaterial::get_parameter(std::string name)
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
	else if(name=="P1")
		return get_P1();
	else if(name=="P2")
		return get_P2();
	else if(name=="YIN" || name=="DAMAGETHRESHOLD" || name=="DAMAGE THRESHOLD")
		return get_Yin();
	else if(name=="MU_VISC" || name=="MU_VISCOUS" || name=="MU VISC" || name=="MU VISCOUS")
		return get_mu_visc();
	else if (name=="THERMAL_EXPANSION" || name=="THERMAL EXPANSION" || name=="THERMAL_SHINKAGE" || name=="THERMAL SHRINKAGE")
		return get_thermal_exp();
	else if (name=="SIGMAY_C")
		return get_sigmay_c();
	else if (name=="SIGMAY_T")
		return get_sigmay_t();
	else if (name=="EPSILON_T")
		return get_epsilon_t();	
	else if (name=="EPSILON_C")
		return get_epsilon_c();	
	else if (name=="H_RO")
		return get_H_ro();	
	else if (name=="N_RO")
		return get_n_ro();
	else
		err_message(name << " is not a valid parameter name");
}


void ContinuumDamageModelMaterial::set_E(double E)
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
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set && sigmay_t_set && sigmay_c_set && epsilon_t_set && epsilon_c_set && H_ro_set && n_ro_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_E()
{
	if(n_set>=2 || E_set)
		return _E;
	else
		err_message("Young's modulus could not be determined.");
}




void ContinuumDamageModelMaterial::set_nu(double nu)
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
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set && sigmay_t_set && sigmay_c_set && epsilon_t_set && epsilon_c_set && H_ro_set && n_ro_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_nu()
{
	if(n_set>=2 || nu_set)
		return _nu;
	else
		err_message("Poisson's ratio could not be determined.");
}




void ContinuumDamageModelMaterial::set_lambda(double lambda)
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
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set && sigmay_t_set && sigmay_c_set && epsilon_t_set && epsilon_c_set && H_ro_set && n_ro_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_lambda()
{
	if(n_set>=2 || lambda_set)
		return _lambda;
	else
		err_message("First Lame constant could not be determined.");
}




void ContinuumDamageModelMaterial::set_mu(double mu)
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
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set && sigmay_t_set && sigmay_c_set && epsilon_t_set && epsilon_c_set && H_ro_set && n_ro_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_mu()
{
	if(n_set>=2 || mu_set)
		return _mu;
	else
		err_message("Shear modulus could not be determined.");
}




void ContinuumDamageModelMaterial::set_K(double K)
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
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set && sigmay_t_set && sigmay_c_set && epsilon_t_set && epsilon_c_set && H_ro_set && n_ro_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_K()
{
	if(n_set>=2 || K_set)
		return _K;
	else
		err_message("Bulk modulus could not be determined.");
}




void ContinuumDamageModelMaterial::set_P1(double P1)
{
	_P1 = P1;
	P1_set = true;
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set && sigmay_t_set && sigmay_c_set && epsilon_t_set && epsilon_c_set && H_ro_set && n_ro_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_P1()
{
	if(P1_set)
		return _P1;
	else
		err_message("P1 could not be determined.");
}




void ContinuumDamageModelMaterial::set_P2(double P2)
{
	_P2 = P2;
	P2_set = true;
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set && sigmay_t_set && sigmay_c_set && epsilon_t_set && epsilon_c_set && H_ro_set && n_ro_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_P2()
{
	if(P2_set)
		return _P2;
	else
		err_message("P2 could not be determined.");
}




void ContinuumDamageModelMaterial::set_Yin(double Yin)
{
	_Yin = Yin;
	Yin_set = true;
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set && sigmay_t_set && sigmay_c_set && epsilon_t_set && epsilon_c_set && H_ro_set && n_ro_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_Yin()
{
	if(Yin_set)
		return _Yin;
	else
		err_message("Yin could not be determined.");
}




void ContinuumDamageModelMaterial::set_mu_visc(double mu_visc)
{
	_mu_visc = mu_visc;
	mu_visc_set = true;
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set && sigmay_t_set && sigmay_c_set && epsilon_t_set && epsilon_c_set && H_ro_set && n_ro_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_mu_visc()
{
	if(mu_visc_set)
		return _mu_visc;
	else
		err_message("mu_visc could not be determined.");
}



void ContinuumDamageModelMaterial::set_sigmay_t(double sigmay_t)
{
	_sigmay_t = sigmay_t;
	sigmay_t_set = true;
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set && sigmay_t_set && sigmay_c_set && epsilon_t_set && epsilon_c_set && H_ro_set && n_ro_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_sigmay_t()
{
	if(sigmay_t_set)
		return _sigmay_t;
	else
		err_message("sigmay_t could not be determined.");
}


void ContinuumDamageModelMaterial::set_sigmay_c(double sigmay_c)
{
	_sigmay_c = sigmay_c;
	sigmay_c_set = true;
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set && sigmay_t_set && sigmay_c_set && epsilon_t_set && epsilon_c_set && H_ro_set && n_ro_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_sigmay_c()
{
	if(sigmay_t_set)
		return _sigmay_c;
	else
		err_message("sigmay_t could not be determined.");
}


void ContinuumDamageModelMaterial::set_epsilon_c(double epsilon_c)
{
	_epsilon_c = epsilon_c;
	epsilon_c_set = true;
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set && sigmay_t_set && sigmay_c_set && epsilon_t_set && epsilon_c_set && H_ro_set && n_ro_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_epsilon_c()
{
	if(epsilon_t_set)
		return _epsilon_c;
	else
		err_message("epsilon_t could not be determined.");
}

void ContinuumDamageModelMaterial::set_epsilon_t(double epsilon_t)
{
	_epsilon_t = epsilon_t;
	epsilon_t_set = true;
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set && sigmay_t_set && sigmay_c_set && epsilon_t_set && epsilon_c_set && H_ro_set && n_ro_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_epsilon_t()
{
	if(epsilon_t_set)
		return _epsilon_t;
	else
		err_message("epsilon_t could not be determined.");
}

void ContinuumDamageModelMaterial::set_H_ro(double H_ro)
{
	_H_ro = H_ro;
	H_ro_set = true;
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set && sigmay_t_set && sigmay_c_set && epsilon_t_set && epsilon_c_set && H_ro_set && n_ro_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_H_ro()
{
	if(H_ro_set)
		return _H_ro;
	else
		err_message("H_ro could not be determined.");
}

void ContinuumDamageModelMaterial::set_n_ro(double n_ro)
{
	_n_ro = n_ro;
	n_ro_set = true;
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set && sigmay_t_set && sigmay_c_set && epsilon_t_set && epsilon_c_set && H_ro_set && n_ro_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_n_ro()
{
	if(n_ro_set)
		return _n_ro;
	else
		err_message("n_ro could not be determined.");
}






void ContinuumDamageModelMaterial::set_thermal_exp(double val)
{
	_thermal_exp = val;
	_thermal_exp_set = true;
}
double ContinuumDamageModelMaterial::get_thermal_exp()
{
	if (_thermal_exp_set)
		return _thermal_exp;
	else
		err_message("Thermal expansion ratio could not be determined.");
}
