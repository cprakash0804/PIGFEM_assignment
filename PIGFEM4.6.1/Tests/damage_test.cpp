#include <iostream>
#include <vector>
#include <fstream>
#include "material.h"
#include "material_cdm.h"
using namespace std;

#define E 3.5e9
#define nu 0.3
#define P1 0.1
#define P2 1
#define mu_visc 20
#define Yin 75e6
#define e_max 0.5
#define t_max 1



int main(void)
{
	Material* mat = new ContinuumDamageModelMaterial;
	std::vector<double> int_vars = {0,0};
	std::vector<double> strain(1,0);
	std::vector<double> strain_hist;
	std::vector<double> stress_hist;

	mat->set_name("Epoxy");
	mat->set_parameter("E", E);
	mat->set_parameter("nu", nu);
	mat->set_parameter("P1", P1);
	mat->set_parameter("P2", P2);
	mat->set_parameter("mu_visc", mu_visc);
	mat->set_parameter("Yin", Yin);

	int n_pts = 1001;
	double delta_t = double(t_max)/double(n_pts-1);


	Material::input_params input;
	input.dim = 1;
	input.delta_t = delta_t;
	input.internal_vars = &(int_vars);

	for(int i=1; i<=n_pts; ++i)
	{
		strain[0] = e_max*double(i-1)/double(n_pts-1);
		input.strain = strain;

		Material::output_params* output = mat->Constitutive(input);
		strain_hist.push_back(strain[0]);
		stress_hist.push_back(output->stress[0]);
	}

	std::ofstream myfile;
	myfile.open("Output/damage_output.txt", std::ofstream::out);

	for(unsigned int i=0; i<strain_hist.size(); ++i)
	{
		myfile << strain_hist[i] << ", " << stress_hist[i] << endl;
	}


	myfile.close();

	delete mat;


}
