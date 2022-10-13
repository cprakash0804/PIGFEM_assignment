/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#include "Material_OPCohesive.h"
#include "Utilities.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>


Material* OPCohesiveMaterial::allocate()
{
	return new OPCohesiveMaterial;
}


std::vector<std::string> OPCohesiveMaterial::internal_vars_name()
{
	std::vector<std::string> vec = {"maximum_opening"};
	return vec;
}

std::vector<bool> OPCohesiveMaterial::internal_vars_print()
{
	std::vector<bool> vec = {false};
	return vec;
}



Material::output_params* OPCohesiveMaterial::Constitutive(Material::input_params& input)
{
	if(!ready)
		err_message("Not all material parameters have been set for the Ortiz-Pandolfi Cohesive material.");

	// Get the internal variables and references to the output
	if((*input.internal_vars).size() != 1)
		err_message("Invalid number of internal variables for a Ortiz-Pandolfi Cohesive material.");
	double& delta_max = (*input.internal_vars)[0];
	DenseMatrix<double>& Dmat = _output.Dmat;
	std::vector<double>& traction_ttn = _output.traction;
	double dt = input.delta_t;

	// Get the dimension of the problem and the input opening vector (int ttn coordinates hopefully)
	// ttn stands for tangential-tangential-normal
	id_type dim = input.dim;
	std::vector<double>& delta_ttn = input.delta;

double& penetrationn = (*input.internal_vars)[1];
double& penet_tan_dir = (*input.internal_vars)[2];
bool penetration_iter = Utilities::penetration_store(0);
bool penetration;

/*if (delta_max==0)
{
	penetration = 0;
}else{
	if (penetration_iter==1)
	{
		if(delta_ttn[delta_ttn.size()-1]>=0)
		{
			penetration = 0;
		}else{
			penetration = 1;
		}
	}else{
		if (delta_max>0)
		{
			penetration = 0;
		}else{
			penetration = 1;
		}
	}

}*/



if (delta_max==0)
{
	penetration = 0;
	penetrationn = penetration;
	penet_tan_dir = 1;
}else{
	if (penetration_iter==1)
	{
		penet_tan_dir=1;
		if(delta_ttn[delta_ttn.size()-1]>=0)
		{
			penetration = 0;
		}else{
			penetration = 1;
		}
		penetrationn = penetration;
	}

}



//
//if (delta_max<0.000000000001 && delta_max>=0)
//{
//	if (delta_ttn[delta_ttn.size()-1]>=0)
//	{
//		penetration = 0;
//	}else{
//		penetration = 1;
//	}
//}else{
//	if(delta_max>0)
//	{
//		penetration = 0;
//	}else{
//		penetration=1;
//	}
//}
penetrationn=0;
double Dtraction_Ddelta = 0;

//Bilienar
double delta_c_bi = 1e-4;
double delta_0_bi = 5e-4;

//Trapezoidal
double delta_0_1_tr = 3e-4;
double delta_0_2_tr = 5e-4;


if (penetrationn==0)
{
	// Get the effective opening
	double delta_eff = getDeltaEff(delta_ttn);

	// Effective traction based on loading (Eq 28) or unloading (Eq 30)
	double trac_eff = 0;
	//double Dtraction_Ddelta = 0;

	if(delta_eff>=delta_max || delta_eff<0.0) // Loading, continue along curve and update ISV
	{
	//if (penet_tan_dir==1)
	//{
		
	// Exponential Behavior	----------------------------------------------------------------------------------------------------------
		trac_eff =  _sigma_c * (delta_eff/_delta_c) * exp(1.0-delta_eff/_delta_c)+_visc/_sigma_c*(delta_eff-delta_max)/dt/_delta_c;
		Dtraction_Ddelta =  (_sigma_c/_delta_c) * (1.0-delta_eff/_delta_c) * exp(1.0-delta_eff/_delta_c)+_visc/_sigma_c/dt/_delta_c;
		

		
	/*	if (trac_eff<1 && delta_eff>_delta_c)
		{
			trac_eff=1;
			Dtraction_Ddelta = 0.0;
			
		} */
		
		if (penetration_iter==0)
			delta_max = delta_eff; 
		
		/*if (penetration_iter = 1)
			if (delta_eff>_delta_c && delta_eff<2*_delta_c)
				delta_eff=1;*/
			
		
		
		
	// Bilinear Behavior-----------------------------------------------------------------------
	/*	if (delta_eff<=delta_c_bi)
		{
			trac_eff = _sigma_c/delta_c_bi*delta_eff;
			Dtraction_Ddelta =_sigma_c/delta_c_bi;
		}
		else if (delta_eff>delta_c_bi && delta_eff<=delta_0_bi)
		{
			trac_eff = -_sigma_c/(delta_0_bi-delta_c_bi)*(delta_eff-delta_c_bi)+_sigma_c;
			Dtraction_Ddelta =-_sigma_c/(delta_0_bi-delta_c_bi);
		}
		else
		{
			trac_eff = 0;
			Dtraction_Ddelta =0;
		}
		
		if (trac_eff<1 && delta_eff>delta_c_bi)
		{
			trac_eff=1;
			Dtraction_Ddelta = 0.0;
			
		} 
		
		if (penetration_iter==0)
			delta_max = delta_eff; 	*/	
			

	/*	{
			delta_max = delta_eff; // Update the internal state variable!
			if(delta_eff>=delta_max)
			{
				penet_tan_dir = 1;
			}else{
				penet_tan_dir = 0;
			}
		}	*/
	// Trapezoidal Behavior-----------------------------------------------------------------------	
	/*	if (delta_eff<=_delta_c)
		{
			trac_eff = _sigma_c/_delta_c*delta_eff;
			Dtraction_Ddelta =_sigma_c/_delta_c;
		}
		else if(delta_eff>_delta_c && delta_eff<=delta_0_1_tr)
		{
			trac_eff = _sigma_c;
			Dtraction_Ddelta =0;
		}
		else if (delta_eff>delta_0_1_tr && delta_eff<=delta_0_2_tr)
		{
			trac_eff = -_sigma_c/(delta_0_2_tr-delta_0_1_tr)*(delta_eff-delta_0_1_tr)+_sigma_c;
			Dtraction_Ddelta = -_sigma_c/(delta_0_2_tr-delta_0_1_tr);
		}
		else
		{
			trac_eff = 0;
			Dtraction_Ddelta =0;
		}
		
		if (trac_eff<1 && delta_eff>delta_c_bi)
		{
			trac_eff=1;
			Dtraction_Ddelta = 0.0;
			
		} 
		
		if (penetration_iter==0)
			delta_max = delta_eff; */
	
	}
	else	// Unloading. Unload to origin and don't update ISV 
	{
		//
		//Exponential Behavior-----------------------------------------------------------
		double max_traction = _sigma_c * (delta_max/_delta_c) * exp(1.0-delta_max/_delta_c);//+_visc/_sigma_c*(delta_eff-delta_max)/dt/_delta_c;  // Traction at max opening
		//if (max_traction<1 && delta_max>_delta_c)
		//	max_traction=1;
		
		trac_eff = max_traction * (delta_eff/delta_max);
		Dtraction_Ddelta = max_traction/delta_max; 
		
		//Bilinear Behavior------------------------------------------------------------
	/*	double max_traction;
		if (delta_max<=delta_c_bi)
		{
			max_traction = _sigma_c/delta_c_bi*delta_max;
		}
		else if (delta_max>delta_c_bi && delta_max<=delta_0_bi)
		{
			max_traction = -_sigma_c/(delta_0_bi-delta_c_bi)*(delta_max-delta_c_bi)+_sigma_c;
		}
		else
		{
			max_traction = 0;
		} 
		
		trac_eff = max_traction * (delta_eff/delta_max);
		Dtraction_Ddelta = max_traction/delta_max; */
		
		//Trapezoidal Behavior------------------------------------------------------------
	/*	double max_traction;
		if (delta_max<=_delta_c)
		{
			max_traction = _sigma_c/_delta_c*delta_max;
		}
		else if(delta_max>_delta_c && delta_max<=delta_0_1_tr)
		{
			max_traction = _sigma_c;
		}
		else if (delta_max>delta_0_1_tr && delta_max<=delta_0_2_tr)
		{
			max_traction = -_sigma_c/(delta_0_2_tr-delta_0_1_tr)*(delta_max-delta_0_1_tr)+_sigma_c;
		}
		else
		{
			max_traction = 0;
		}
		
		trac_eff = max_traction * (delta_eff/delta_max);
		Dtraction_Ddelta = max_traction/delta_max; */
		
	} 

	// Actual traction vector
	std::vector<double> bracket_term(delta_ttn.size(), 0.0); // Term in brackets of Eq (43)
	traction_ttn.resize(delta_ttn.size());
	for(id_type i=0; i<delta_ttn.size(); ++i)
		bracket_term[i] = _beta2 * delta_ttn[i];
	bracket_term[dim-1] += (1.0-_beta2) * delta_ttn[dim-1]; // Add the second term of the bracketed expression (in the ttn coord system the dot product is just equal to the normal component of delta_ttn)
	for(id_type i=0; i<delta_ttn.size(); ++i)
		traction_ttn[i] = (trac_eff/delta_eff) * bracket_term[i];

	// Derivatives of delta_eff wrt delta_ttn
	std::vector<double> Ddeltaeff_Ddeltattn(dim);
	Ddeltaeff_Ddeltattn[dim-1] = delta_ttn[dim-1]/delta_eff;
	for(id_type i=0; i<(dim-1); i++)
		Ddeltaeff_Ddeltattn[i] = _beta2 * delta_ttn[i]/delta_eff;

	// Computation of the D matrix (Derivative of Eq 43 wrt delta_ttn)
	Dmat.resize(dim, dim);
	for(id_type i=0; i<dim; ++i)
		for(id_type k=0; k<dim; ++k)
			Dmat(i, k) = bracket_term[i] * (Dtraction_Ddelta*delta_eff - trac_eff)*(Ddeltaeff_Ddeltattn[k]/pow(delta_eff,2)) + 
						(trac_eff/delta_eff) * (_beta2*Utilities::kron(i,k)+(1.0-_beta2)*Utilities::kron(i,dim-1)*Utilities::kron(k,dim-1)); // NOTE: this is assuming that the normal vector is [0,0,1] (3D)
}
else
{

	double traction_t;
	double EE = 100000000;
	double traction_n = EE*delta_ttn[1];
	double unloading_check = sqrt(pow(delta_ttn[0],2))-sqrt(pow(delta_max,2));
	if (penet_tan_dir==1)
	{
	//if (unloading_check>=0)
	//{
		traction_t = _sigma_c * (delta_ttn[0]/_delta_c) * exp(1.0-sqrt(pow(delta_ttn[0],2))/_delta_c);
		Dtraction_Ddelta = (_sigma_c/_delta_c) * (1.0-sqrt(pow(delta_ttn[0],2))/_delta_c) * exp(1.0-sqrt(pow(delta_ttn[0],2))/_delta_c);

		if (penetration_iter==0)
			delta_max = delta_ttn[0];
		{
			
			if(unloading_check>=0)
			{
				penet_tan_dir = 1;
			}else{
				penet_tan_dir = 0;
			}

		} 
	}else{ // Unloading Case
		double max_traction = _sigma_c * (delta_max/_delta_c) * exp(1.0-delta_max/_delta_c);
		traction_t = max_traction * (delta_ttn[0]/delta_max);
		Dtraction_Ddelta = max_traction/delta_max;

	} 




	if (Dtraction_Ddelta==0)
	{
		Dtraction_Ddelta = 0.000000001;
	}

	Dmat.resize(dim, dim);
	Dmat(0,0)=Dtraction_Ddelta;
	Dmat(0,1)=0;
	Dmat(1,0)=0;
	Dmat(1,1)=EE;


	/* double EE = 10000000;
	double traction_n = EE*delta_ttn[1];
	double traction_t = _sigma_c * (delta_ttn[0]/_delta_c) * exp(1.0-delta_ttn[0]/_delta_c);
	delta_max = (-1);//*delta_ttn[0];
	double temp = traction_t/delta_ttn[0];
	if (delta_ttn[0]==0)
		temp = 0.000001;

	Dmat.resize(dim, dim);
	Dmat(0,0)=temp;
	Dmat(0,1)=0;
	Dmat(1,0)=0;
	Dmat(1,1)=EE; */

	traction_ttn.resize(delta_ttn.size());
	traction_ttn[0] = traction_t;
	traction_ttn[1] = traction_n;
}

	// Return a pointer to the output structure
	return &_output;
}





Material::sensitivity_output_params* OPCohesiveMaterial::SensitivityConstitutive(Material::sensitivity_input_params& input)
{
	if(!ready)
		err_message("Not all material parameters have been set for the Continuum Damage Material.");

	// Get the internal variables and references to the output
	double& delta_max = (*input.internal_vars)[0];
	id_type param_id = input.parameter_id;
	double& ddelta_max_dd = (*input.internal_vars_sensitivity)[param_id * n_internal_vars()];
	std::vector<double>& dtraction_ttn_dd = _sensitivity_output.dtraction_dd;

	// Get the effective opening
	id_type dim = input.dim;
	std::vector<double>& delta_ttn = input.delta;
	double delta_eff = getDeltaEff(delta_ttn);

	// Traction vector partial derivatives
	dtraction_ttn_dd.clear();
	dtraction_ttn_dd.resize(delta_ttn.size());

	// Used a lot so I'll just compute it here
	std::vector<double> bracket_term(delta_ttn.size(), delta_ttn[dim-1]); // Term in brackets of Eq (43)
	for(id_type i=0; i<(dim-1); ++i)
		bracket_term[i] = _beta2 * delta_ttn[i];

	// just store boolean here for ease of use
	bool loading = (delta_eff >= delta_max);

	// Partial derivative is only non-zero if the material I am taing the derivative wrt is the same as the current material I am in
	if (input.sensitivity_mat_name == _name) 
	{
		// Sigma_c sensitivity
		if (input.sensitivity_param_name == "SIGMA_C")
		{
			// 3b PARTIAL derivative
			double dtrac_eff_dsigmac;
			if (loading)
				dtrac_eff_dsigmac = delta_eff/_delta_c*exp(1-delta_eff/_delta_c);
			else
				dtrac_eff_dsigmac = delta_eff/_delta_c*exp(1-delta_max/_delta_c);

			for(id_type i=0; i<delta_ttn.size(); ++i)
				dtraction_ttn_dd[i] = (dtrac_eff_dsigmac/delta_eff) * bracket_term[i];
		}

		// Delta_c sensitivity
		else if (input.sensitivity_param_name == "DELTA_C")
		{
			// 3c PARTIAL derivative
			double dtrac_eff_ddeltac;
			if (loading)
				dtrac_eff_ddeltac = (_sigma_c*delta_eff/(_delta_c*_delta_c)*exp(1-delta_eff/_delta_c)) * (delta_eff/_delta_c - 1);
			else
				dtrac_eff_ddeltac = (_sigma_c*delta_eff/(_delta_c*_delta_c)*exp(1-delta_max/_delta_c)) * (delta_max/_delta_c - 1);

			for(id_type i=0; i<delta_ttn.size(); ++i)
				dtraction_ttn_dd[i] = (dtrac_eff_ddeltac/delta_eff) * bracket_term[i];
		}

		// Beta sensitivity
		else if (input.sensitivity_param_name == "BETA")
		{
			err_message("Ortiz-Pandolphi Beta sensitivity has not been derived yet!");
		}
	}

	// Total derivative wrt the internal state variable term is present regardless of whether or not I'm in the same material
	if (!loading)
	{
		double dtrac_eff_ddeltamax = (_sigma_c * delta_eff / _delta_c) * exp(1-delta_max/_delta_c) * (-1.0/_delta_c);
		for(id_type i=0; i<delta_ttn.size(); ++i)
			dtraction_ttn_dd[i] += ddelta_max_dd * (dtrac_eff_ddeltamax/delta_eff) * bracket_term[i];
	}

	return &_sensitivity_output;
}



void OPCohesiveMaterial::updateSensitivityISVs(Material::sensitivity_input_params& input)
{
	if(!ready)
		err_message("Not all material parameters have been set for the Continuum Damage Material.");

	// Get the internal variables and references to the output
	double& delta_max = (*input.internal_vars)[0];
	id_type param_id = input.parameter_id;
	double& ddelta_max_dd = (*input.internal_vars_sensitivity)[param_id * n_internal_vars()];

	// Get the effective opening
	id_type dim = input.dim;
	std::vector<double>& delta_ttn = input.delta;
	double delta_eff = getDeltaEff(delta_ttn);

	// Loading, continue along curve and update ISV
	if(delta_eff >= delta_max)
	{
		std::vector<double>& ddelta_ttn_dd = input.ddelta_dd;

		// 2b
		double ddelta_eff_dd = 2.0*delta_ttn[dim-1]*ddelta_ttn_dd[dim-1];
		for (id_type i=0; i<(dim-1); ++i)
			ddelta_eff_dd += 2.0*_beta2*delta_ttn[i]*ddelta_ttn_dd[i];
		ddelta_eff_dd *= 1.0/(2.0*delta_eff);

		// Update the sensitivity internal state variables!
		ddelta_max_dd = ddelta_eff_dd;
	}
}
