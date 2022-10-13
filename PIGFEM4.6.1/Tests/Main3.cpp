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
#include "Writer.h"
#include "common.h"
#include "mpi.h"

using namespace std;


#define pi 3.141592653589799323

// Define the x-axial loading strain
#define strain 0.02

// Define matrix properties (Taken from Scott's ppt)
#define E1 2.38e9
#define nu1 0.43
#define sigma_max 85e6 // From http://www.epoxyworktops.com/epoxy-resin/mech-properties.html

// Define carbon fiber transverse properties (Taken from Scott's ppt)
#define E2 19.5e9
#define nu2 0.45

// Define cohesive properties
#define sigmac 61e6
#define deltac 1.8394e-6 // Taken from balancing the G_c in Scott's ppt with the integral of Ortiz-Pandolfi curve for the given sigmac

// Define transversely isotrpoic properties (Taken from Scott's ppt)
#define Eti1 36.2e9
#define Eti2 7.11e9
#define G12  2.32e9
#define G23  2.17e9
#define nu12 0.335

// Parameters to control size and number of elements in the domain
#define n_elem_per_inclusion 5	// Define the number of elements to put across the diameter of the smallest fiber
#define yexpansion 0.2	// Increase the size of the domain (add matrix) by this much on each side of the domain (0.05=5%)
#define xexpansion 0.05
#define bounding_ply_limit 0.17 // the transversely isotropic bounding plies cover top and bottom 3% of the domain
//#define xexpansion 0.35
//#define yexpansion 0.35
//#define bounding_ply_limit 0.2

// Define damage properties (not currently used here)
#define P1 0.1
#define P2 1
#define mu_visc 20
#define Yin 75e3
// The size of the square domain in the positive and negative x and y directions (centers at (0,0))
// ~/projects/petsc-3.6.2/arch-linux2-c-opt/bin/mpiexec

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
	mins[1] = mins[1] - y_lims*yexpansion;
	maxes[0] = maxes[0] + x_lims*xexpansion;
	maxes[1] = maxes[1] + y_lims*yexpansion;
	id_type x_elem = std::ceil((maxes[0]-mins[0]) / (min_rad*2)) * n_elem_per_inclusion;
	id_type y_elem = std::ceil((maxes[1]-mins[1]) / (min_rad*2)) * n_elem_per_inclusion;
	mesh.generate_mesh(QUAD4, mins[0], maxes[0], x_elem, mins[1], maxes[1], y_elem); // Actually generate the mesh

	// Add the matrix to every element (linear elastic and fails at a max priniciple stress)
	LinearElasticIsotropicMaterial matrix;
	matrix.set_parameter("E", E1);
	matrix.set_parameter("nu", nu1);
	matrix.set_name("Epoxy");
	mesh.add_material(&matrix);
	mesh.set_material("Epoxy");

	// Assign BCs
	double app_disp = (maxes[0]-mins[0]) * strain;
	BoundaryObject* boundary = mesh.get_boundary();
	boundary->set_dirichlet_bcs_from_nodeset("right", 0, app_disp);
	//boundary->set_dirichlet_bcs_from_nodeset("right", 1, x_disp);
	boundary->set_dirichlet_bcs_from_nodeset("left", 0, 0.0);
	for(unsigned int d=0; d<mesh.dim(); ++d)
		boundary->set_dirichlet_bc(0, d, 0);

	// Create the cohesive material (Will automatically be applied to any material boundary
	//	that doesn't have another cohesive material specified for the touching material pair)
/*	OPCohesiveMaterial coh_mat;
	coh_mat.set_parameter("sigma", sigmac);
	coh_mat.set_parameter("delta", deltac);
	coh_mat.set_parameter("beta", 1.0);
	coh_mat.set_name("Epoxy-Glass");
	mesh.add_material(&coh_mat); */

	// Define the problem
	ProblemLinearElasticity problem;
	//ProblemNonlinearStructural problem;
	problem.attach_mesh(&mesh);
	
	// Solve the problem!
	if(rank==rank_see)
		cout << "\nSolving the problem...";
	problem.solve_problem();

	// Output results to file
	if(rank==rank_see)
		cout << "\nOutputting results to file...";
	char buf[50];
	sprintf(buf, "Output/IGFEM_out_%i", mesh.n_inclusions());
	problem.get_writer()->write( buf );
	
	// Output more detailed results
	if(rank==rank_see)
	{
		cout << "\n\nPARTITION " << rank_see << "\n--------------------------------------------------\n";
		cout << "Number of Elements: " << problem.get_mesh()->n_global_active_elem() << endl;
		cout << "Number of Nodes: " << problem.get_mesh()->n_global_nodes()+problem.get_mesh()->n_global_enrich_nodes() << endl;
		cout << "Number of dofs: " << problem.get_dofs()->n_global_dofs() << endl;
	}


	problem.Finalize();
}
