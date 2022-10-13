#include <iostream>
#include <fstream>
#include <cstddef>
#include <sys/time.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <set>
#include <map>
#include <string>
#include <cmath>
#include "Mesh.h"
#include "Problem.h"
#include "ProblemMaxPrincipalStress.h"
#include "ProblemNonlinearStructural.h"
#include "ProblemNonlinearStructural_Multiscale.h"
#include "ProblemLinearElasticity.h"
#include "ProblemLinearThermal.h"
#include "BoundaryObject.h"
#include "DofObject.h"
#include "Inclusion.h"
#include "InclusionPlane.h"
#include "InclusionEllipse.h"
#include "InclusionEllipsoid.h"
#include "material.h"
#include "material_cdm.h"
#include "material_lei.h"
#include "material_leti.h"
#include "material_leimps.h"
#include "material_leipcd.h"
#include "material_lt.h"
#include "Material_OPCohesive.h"
#include "BodyLoad.h"
#include "BodyLoadMacroStrainLinear.h"
#include "SubscaleModelLinear.h"
#include "SubscaleModelNonlinear.h"
#include "Writer.h"
#include "common.h"
#include "mpi.h"

using namespace std;

// Path to the mpiexec executable
// ./petsc-3.6.2/arch-linux2-c-opt/bin/mpiexec


#define pi 3.141592653589799323

// Define the x-axial loading strain
#define strain 0.05

// Parameters to control size and number of elements in the domain
#define n_elem_per_inclusion 20	// Define the number of elements to put across the diameter of the smallest fiber
#define yexpansion_plies 0.2	// Increase the size of the domain (add matrix) by this much on each side of the domain (0.05=5%)
#define yexpansion 0.02
#define xexpansion 0.02
#define bounding_ply_limit 0.17 // the transversely isotropic bounding plies cover top and bottom 17% of the domain

// Set some flags here
#define cohesive_flag true
#define max_principal_stress_flag false
#define continuum_damage_flag true
#define bounding_plies_flag false

// Define matrix properties (Taken from Scott's ppt)
#define E1 2.38e9
#define nu1 0.43
#define sigma_max 85e6 // From http://www.epoxyworktops.com/epoxy-resin/mech-properties.html

// Define carbon fiber transverse properties (Taken from Scott's ppt)
#define E2 19.5e9
#define nu2 0.45

// Define cohesive properties
//#define sigmac 61e6
//#define deltac 1.8394e-6 // Taken from balancing the G_c in Scott's ppt with the integral of Ortiz-Pandolfi curve for the given sigmac
#define sigmac 61e5
#define deltac 1.8394e-7

// Define transversely isotrpoic properties (Taken from Scott's ppt)
#define Eti1 36.2e9
#define Eti2 7.11e9
#define G12  2.32e9
#define G23  2.17e9
#define nu12 0.335

// Define damage properties
#define P1 0.1
#define P2 1
#define mu_visc 20
#define Yin 75e3




template <typename T>
void print_vec(const std::vector<T> vec)
{
	if(vec.size() >= 1)
	{
		std::cout << "{";
		for(unsigned int i=0; i<(vec.size()-1); ++i)
			std::cout << vec[i] << ", ";
		std::cout << vec[vec.size()-1] << "}";
	}
}


void read_inclusions(Mesh* mesh, std::string filename, std::vector<double>& maxes, std::vector<double>& mins, double& min_rad)
{
	LinearElasticIsotropicMaterial fiber;
	fiber.set_parameter("E", E2);
	fiber.set_parameter("nu", nu2);
	fiber.set_name("Glass");

	ifstream myfile;
	myfile.open( filename );
	if (!myfile.good())
		err_message("Invalid input file name!");

	std::string dummy;
	char c;

	// Read in the circles (x, y, r)
	double x, y, r;
	std::vector<double> center;
	Ellipse_Inclusion ellipse;
	ellipse.set_material(&fiber);
	while (myfile.peek() != EOF)
	{
		// Read the data from the file
		myfile >> x >> c >> y >> c >> r;
		std::getline(myfile, dummy); // Read in the endline
		x = x / 1000; // To convert to mm
		y = y / 1000; // To convert to mm
		r = r / 1000; // To convert to mm
		maxes[0] = std::max(maxes[0], x+r);
		maxes[1] = std::max(maxes[1], y+r);
		mins[0] = std::min(mins[0], x-r);
		mins[1] = std::min(mins[1], y-r);
		min_rad = std::min(min_rad, r);

		// Set the ellipse parameters
		center = {x, y};
		ellipse.set_vec_parameter("center", center);
		ellipse.set_parameter("a", r);
		ellipse.set_parameter("b", r);

		// Add the inclusion to the mesh
		mesh->add_inclusion(&ellipse);
	}

	myfile.close();
}



// The main function of this test is to test the mesh adaptivity functions
int main (int argc, char* argv[])
{
	// Make sure that a file name was actualy passed in first thing
	if (argc <= 1)
		err_message("Please input a filename to read circular inclusions from.");

	SubscaleModel* model;
	if (cohesive_flag || max_principal_stress_flag || continuum_damage_flag)
		model = new SubscaleModelNonlinear(&argc, &argv);
	else
		model = new SubscaleModelLinear(&argc, &argv);

	Mesh* mesh = model->get_mesh();
	Problem* prob = model->get_problem();
	int rank = mesh->get_rank();
	int rank_see = 0;

	// Read in the inclusions (Add fiber material in this function)
	if(rank==rank_see)
		cout << "Reading inclusions...\n";
	std::vector<double> maxes(2, -99999999999999999999999999.0); // Maximum domain limits (x and y)
	std::vector<double> mins(2, 99999999999999999999999999.0);   // Minimum domain limits (x and y)
	double min_rad = 99999999999999999999999999.0; 				 // Minimum fiber radius (So I know how many elements across to use)
	std::string filename(argv[1]);
	filename = "Input/" + filename;
	read_inclusions(mesh, filename, maxes, mins, min_rad);


	// Generate a mesh
	if(rank==rank_see)
		cout << "Generating the mesh...";
	double x_lims = maxes[0] - mins[0];
	double y_lims = maxes[1] - mins[1];
	mins[0] = mins[0] - x_lims*xexpansion;
	mins[1] = mins[1] - y_lims*yexpansion;
	maxes[0] = maxes[0] + x_lims*xexpansion;
	maxes[1] = maxes[1] + y_lims*yexpansion;
	id_type x_elem = std::ceil((maxes[0]-mins[0]) / (min_rad*2)) * n_elem_per_inclusion;
	id_type y_elem = std::ceil((maxes[1]-mins[1]) / (min_rad*2)) * n_elem_per_inclusion;
	mesh->generate_mesh(QUAD4, mins[0], maxes[0], x_elem, mins[1], maxes[1], y_elem); // Actually generate the mesh

	// Add the matrix to every element (linear elastic and fails at a max priniciple stress)
	if (max_principal_stress_flag)
	{
		LinearElasticIsotropicProblemControlledDamageMaterial matrix;
		matrix.set_parameter("E", E1);
		matrix.set_parameter("nu", nu1);
		matrix.set_name("Epoxy");
		mesh->add_material(&matrix);
		mesh->set_material("Epoxy");
	}
	else if (continuum_damage_flag)
	{
		ContinuumDamageModelMaterial matrix;
		matrix.set_parameter("E", E1);
		matrix.set_parameter("nu", nu1);
		matrix.set_parameter("P1", P1);
		matrix.set_parameter("P2", P2);
		matrix.set_parameter("mu_visc", mu_visc);
		matrix.set_parameter("Yin", Yin);
		matrix.set_name("Epoxy");
		mesh->add_material(&matrix);
		mesh->set_material("Epoxy");
	}
	else // Linear elastic
	{
		LinearElasticIsotropicMaterial matrix;
		matrix.set_parameter("E", E1);
		matrix.set_parameter("nu", nu1);
		matrix.set_name("Epoxy");
		mesh->add_material(&matrix);
		mesh->set_material("Epoxy");
	}


	if (cohesive_flag)
	{
		// Create the cohesive material (Will automatically be applied to any material boundary
		//	that doesn't have another cohesive material specified for the touching material pair)
		OPCohesiveMaterial coh_mat;
		coh_mat.set_parameter("sigma", sigmac);
		coh_mat.set_parameter("delta", deltac);
		coh_mat.set_parameter("beta", 1.0);
		coh_mat.set_name("Epoxy-Glass");
		mesh->add_material(&coh_mat);
	}

	// Add the inclusinos for this test problem
	
	LinearElasticIsotropicMaterial fiber;
	fiber.set_parameter("E", E2);
	fiber.set_parameter("nu", nu2);
	fiber.set_name("Glass");

	std::vector<double> center;
	Ellipse_Inclusion ellipse;
	ellipse.set_material(&fiber);
	center = {0.32, 0.19};
	ellipse.set_vec_parameter("center", center);
	ellipse.set_parameter("a", 0.11);
	ellipse.set_parameter("b", 0.11);
	mesh->add_inclusion(&ellipse);
	center = {0.60, 0.57};
	ellipse.set_vec_parameter("center", center);
	ellipse.set_parameter("a", 0.24);
	ellipse.set_parameter("b", 0.24);
	mesh->add_inclusion(&ellipse);

	// Define the macroscopic strain
	std::vector<double> strain_load = {strain, 0.0, 0.0};
	model->set_macro_strain(strain_load);

	// Solve the model
	if(rank==rank_see)
		cout << "\nSolving the model...";
	model->solve();
	
	// Output results to file
	if(rank==rank_see)
		cout << "\nOutputting results to file...";
	char buf[50];
	sprintf(buf, "Output/IGFEM_out_%i", mesh->n_inclusions());
	prob->get_writer()->write( buf );
	
	
	// Output more detailed results
	if(rank==rank_see)
	{
		cout << "\n\nPARTITION " << rank_see << "\n--------------------------------------------------\n";
		cout << "Number of Elements: " << mesh->n_global_elem() << endl;
		cout << "Number of Nodes: " << mesh->n_global_nodes()+mesh->n_global_enrich_nodes() << endl;
		cout << "Number of dofs: " << prob->get_dofs()->n_global_dofs() << endl;
	}


	delete model;
}
