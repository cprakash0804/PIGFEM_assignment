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
#include "Material_OPCohesiveNU.h"
#include "Writer.h"
#include "common.h"
#include "mpi.h"

using namespace std;

// Path to the mpiexec executable
// ./petsc-3.6.2/arch-linux2-c-opt/bin/mpiexec


#define pi 3.141592653589799323

// Define the x-axial loading strain
#define strain 0.03

// Parameters to control size and number of elements in the domain
#define n_elem_per_inclusion 5	// Define the number of elements to put across the diameter of the smallest fiber
#define yexpansion_plies 0.2	// Increase the size of the domain (add matrix) by this much on each side of the domain (0.05=5%)
#define yexpansion 0.05
#define xexpansion 0.05
#define bounding_ply_limit 0.17 // the transversely isotropic bounding plies cover top and bottom 17% of the domain

// Set some flags here
#define cohesive_flag true
#define output_coh true
#define max_principal_stress_flag false
#define continuum_damage_flag false
#define bounding_plies_flag false

// Define matrix properties (Taken from Scott's ppt)
#define E1 2380
#define nu1 0.43
#define sigma_max 85 // From http://www.epoxyworktops.com/epoxy-resin/mech-properties.html

// Define carbon fiber transverse properties (Taken from Scott's ppt)
#define E2 19.5e3
#define nu2 0.45

// Define cohesive properties
#define sigmac 61
//#define deltac 1.016666666666667e-05 // Scott's, not taking into account fracture energy...just matching deltac from bi-linear
#define deltac 1e-4 // Taken from balancing the G_c in Scott's ppt with the integral of Ortiz-Pandolfi curve for the given sigmac

// Define transversely isotropic properties (Taken from Scott's ppt)
#define Eti1 36.2e3
#define Eti2 7.11e3
#define G12  2.32e3
#define G23  2.17e3
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
		//x = x / 1000; // To convert to m
		//y = y / 1000; // To convert to m
		//r = r / 1000; // To convert to m
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

	Mesh mesh(&argc, &argv);
	int rank = mesh.get_rank();
	int rank_see = 0;

	// Read in the inclusions (Add fiber material in this function)
	if(rank==rank_see)
		cout << "Reading inclusions...\n";
	std::vector<double> maxes(2, -99999999999999999999999999.0); // Maximum domain limits (x and y)
	std::vector<double> mins(2, 99999999999999999999999999.0);   // Minimum domain limits (x and y)
	double min_rad = 99999999999999999999999999.0; 				 // Minimum fiber radius (So I know how many elements across to use)
	std::string filename(argv[1]);
	filename = "Input/" + filename;
	read_inclusions(&mesh, filename, maxes, mins, min_rad);

	// Generate a mesh
	if(rank==rank_see)
		cout << "Generating the mesh...";
	double x_lims = maxes[0] - mins[0];
	double y_lims = maxes[1] - mins[1];
	mins[0] = mins[0] - x_lims*xexpansion;
	maxes[0] = maxes[0] + x_lims*xexpansion;
	id_type x_elem = std::ceil((maxes[0]-mins[0]) / (min_rad*2)) * n_elem_per_inclusion;
	id_type y_elem;
	if (bounding_plies_flag)
	{
		mins[1] = mins[1] - y_lims*yexpansion_plies;
		maxes[1] = maxes[1] + y_lims*yexpansion_plies;
		y_elem = std::ceil((maxes[1]-mins[1]) / (min_rad*2)) * n_elem_per_inclusion;
	}
	else
	{
		mins[1] = mins[1] - y_lims*yexpansion;
		maxes[1] = maxes[1] + y_lims*yexpansion;
		y_elem = std::ceil((maxes[1]-mins[1]) / (min_rad*2)) * n_elem_per_inclusion;
	}
	mesh.generate_mesh(TRI3, mins[0], maxes[0], x_elem, mins[1], maxes[1], y_elem); // Actually generate the mesh


	// Add the matrix to every element (linear elastic and fails at a max priniciple stress)
	if (max_principal_stress_flag)
	{
		LinearElasticIsotropicProblemControlledDamageMaterial matrix;
		matrix.set_parameter("E", E1);
		matrix.set_parameter("nu", nu1);
		matrix.set_name("Epoxy");
		mesh.add_material(&matrix);
		mesh.set_material("Epoxy");
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
		mesh.add_material(&matrix);
		mesh.set_material("Epoxy");
	}
	else // Linear elastic
	{
		LinearElasticIsotropicMaterial matrix;
		matrix.set_parameter("E", E1);
		matrix.set_parameter("nu", nu1);
		matrix.set_name("Epoxy");
		mesh.add_material(&matrix);
		mesh.set_material("Epoxy");
	}

	// If we are doing bounding plies, add them here
	if (bounding_plies_flag)
	{
		// Create the element sets that are defined as the top and bottom ?% of the mesh
		std::set<id_type> top_ply, bottom_ply;
		double top_limit = maxes[1] - y_lims*bounding_ply_limit;
		double bottom_limit = mins[1] + y_lims*bounding_ply_limit;
		// Loop over every active element and see if all of its nodes lie within the top or bottom x% of the domain
		for (Mesh::element_iterator it=mesh.active_elements_begin(), end=mesh.active_elements_end(); it!=end; ++it)
		{
			bool in_top_ply = true;
			bool in_bottom_ply = true;
			for (id_type n=0; n<(*it)->n_nodes(); ++n)
			{
				std::vector<double> coords = (*it)->get_node(n)->get_coords();
				double y = coords[1];
				if (y < top_limit)
					in_top_ply = false;
				if(y > bottom_limit)
					in_bottom_ply = false;
			}
			if (in_top_ply)
				top_ply.insert((*it)->get_id());
			if (in_bottom_ply)
				bottom_ply.insert((*it)->get_id());
		}
		mesh.add_elemset("top_ply", top_ply);
		mesh.add_elemset("bottom_ply", bottom_ply);

		// Create the transversely isotropic material to the top and the bottom (Overwrites the matrix set earlier)
		LinearElasticTransverselyIsotropicMaterial bounds;
		bounds.set_parameter("E1", Eti1);
		bounds.set_parameter("E2", Eti2);
		bounds.set_parameter("G12", G12);
		bounds.set_parameter("G23", G23);
		bounds.set_parameter("nu12", nu12);
		bounds.set_name("0_degree_plies");
		mesh.add_material(&bounds);
		mesh.set_material_from_elemset("0_degree_plies", "top_ply");
		mesh.set_material_from_elemset("0_degree_plies", "bottom_ply");
	}

	if (cohesive_flag)
	{
		// Create the cohesive material (Will automatically be applied to any material boundary
		//	that doesn't have another cohesive material specified for the touching material pair)
		OPCohesiveNUMaterial coh_mat; // no unloading.
		coh_mat.set_parameter("sigma", sigmac);
		coh_mat.set_parameter("delta", deltac);
		coh_mat.set_parameter("beta", 1.0);
		coh_mat.set_name("Epoxy-Glass");
		mesh.add_material(&coh_mat);
	}

	// Assign BCs
	double app_disp = (maxes[0]-mins[0]) * strain;
	BoundaryObject* boundary = mesh.get_boundary();
	boundary->set_dirichlet_bcs_from_nodeset("right", 0, app_disp);
	//boundary->set_dirichlet_bcs_from_nodeset("right", 1, x_disp);
	boundary->set_dirichlet_bcs_from_nodeset("left", 0, 0.0);
	boundary->set_dirichlet_bcs_from_nodeset("bottom", 1, 0.0);
	for(unsigned int d=0; d<mesh.dim(); ++d)
		boundary->set_dirichlet_bc(0, d, 0);





	// Define the problem
	Problem * problem;
	if (max_principal_stress_flag)
	{
		problem = new ProblemMaxPrincipalStress;
		problem->set_parameter("sigma_max", sigma_max);
	}
	else if (continuum_damage_flag || cohesive_flag)
		problem = new ProblemNonlinearStructural;
	else // linear elastic
		problem = new ProblemLinearElasticity;

	// Solve the problem
	if(rank==rank_see)
		cout << "\nSolving the problem...";
	problem->output_cohesive() = output_coh;
	problem->attach_mesh(&mesh);
	problem->solve_problem();
	
	// Output results to file
	if(rank==rank_see)
		cout << "\nOutputting results to file...";
	char buf[50];
	sprintf(buf, "Output/IGFEM_out_%i", problem->get_mesh()->n_inclusions());
	problem->get_writer()->write( buf );
	
	
	// Output more detailed results
	if(rank==rank_see)
	{
		cout << "\n\nPARTITION " << rank_see << "\n--------------------------------------------------\n";
		cout << "Number of Elements: " << problem->get_mesh()->n_global_elem() << endl;
		cout << "Number of Nodes: " << problem->get_mesh()->n_global_nodes()+problem->get_mesh()->n_global_enrich_nodes() << endl;
		cout << "Number of dofs: " << problem->get_dofs()->n_global_dofs() << endl;
	}


	problem->Finalize();
	delete problem;
}
