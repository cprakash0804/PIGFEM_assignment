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
#include "ProblemNonlinearStructural.h"
#include "ProblemLinearElasticity.h"
#include "ProblemLinearThermal.h"
#include "BoundaryObject.h"
#include "DofObject.h"
#include "Inclusion.h"
#include "InclusionPlane.h"
#include "InclusionEllipse.h"
#include "InclusionEllipsoid.h"
#include "elem.h"
#include "point1.h"
#include "edge2.h"
#include "tri3.h"
#include "quad4.h"
#include "tet4.h"
#include "hex8.h"
#include "material.h"
#include "material_cdm.h"
#include "material_lei.h"
#include "material_lt.h"
#include "Material_OPCohesive.h"
#include "mpi.h"

using namespace std;


#define pi 3.141592653589799323

#define a 1.06
#define d_size 4
//#define z_size 5e-2
#define z_size 4
#define nelem 9
#define strain 0.01

#define E1 2.5e9
#define E2 70e9
#define nu1 0
#define nu2 0

#define sigmac 10e6
#define deltac 1.5e-7

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
			cout << vec[i] << ", ";
		cout << vec[vec.size()-1] << "}";
	}
}


void read_inclusions(Mesh* mesh, std::string filename, std::vector<double>& maxes, std::vector<double>& mins, double& min_rad)
{
	LinearElasticIsotropicMaterial mat2;
	mat2.set_parameter("E", E2);
	mat2.set_parameter("nu", nu2);
	mat2.set_name("Glass");

	ifstream myfile;
	myfile.open( filename );
	if (!myfile.good())
		std::cerr << "Invalid input file name!." << std::endl;

	std::string dummy;
	char c;

	// Read in the circles (x, y, r)
	double x, y, r;
	std::vector<double> center;
	Ellipse_Inclusion ellipse;
	ellipse.set_material(&mat2);
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
	Mesh mesh(&argc, &argv);
	int rank = mesh.get_rank();
	int rank_see = 0;

	// Read in the inclusions
	if(rank==rank_see)
		cout << "Reading inclusions...\n";
	std::vector<double> maxes(2, -99999999999999999999999999.0);
	std::vector<double> mins(2, 99999999999999999999999999.0);
	double min_rad = 99999999999999999999999999.0; // Really large number
	read_inclusions(&mesh, "Input/39_fibers.txt", maxes, mins, min_rad);

	// Generate a mesh
	if(rank==rank_see)
		cout << "Generating the mesh...";
	id_type n_elem_per_inclusion = 10;
	double expansion = 0.05;
	double x_lims = maxes[0] - mins[0];
	double y_lims = maxes[1] - mins[1];
	mins[0] = mins[0] - x_lims*expansion;
	mins[1] = mins[1] - y_lims*expansion;
	maxes[0] = maxes[0] + x_lims*expansion;
	maxes[1] = maxes[1] + y_lims*expansion;
	id_type x_elem = std::ceil((maxes[0]-mins[0]) / (min_rad*2)) * n_elem_per_inclusion;
	id_type y_elem = std::ceil((maxes[1]-mins[1]) / (min_rad*2)) * n_elem_per_inclusion;
	mesh.generate_mesh(QUAD4, mins[0], maxes[0], x_elem, mins[1], maxes[1], y_elem);

	// Add a material to the mesh
	//ContinuumDamageModelMaterial material;
	LinearElasticIsotropicMaterial material;
	material.set_parameter("E", E1);
	material.set_parameter("nu", nu1);
	material.set_name("Epoxy");
	mesh.add_material(&material);
	mesh.set_material("Epoxy");

	// Cohesive material

	OPCohesiveMaterial coh_mat;
	coh_mat.set_parameter("sigma", sigmac);
	coh_mat.set_parameter("delta", deltac);
	coh_mat.set_parameter("beta", 1.0);
	coh_mat.set_name("Epoxy-Glass");
	mesh.add_material(&coh_mat);


	// Assign BCs
	double app_disp = (maxes[0]-mins[0]) * strain;
	BoundaryObject* boundary = mesh.get_boundary();
	boundary->set_dirichlet_bcs_from_nodeset("right", 0, app_disp);
	//boundary->set_dirichlet_bcs_from_nodeset("right", 1, x_disp);
	boundary->set_dirichlet_bcs_from_nodeset("left", 0, 0.0);
	for(unsigned int d=0; d<mesh.dim(); ++d)
		boundary->set_dirichlet_bc(0, d, 0);





	// Define the problem
	ProblemNonlinearStructural problem;
	//ProblemLinearElasticity problem;
	problem.attach_mesh(&mesh);
	// Solve the problem!
	if(rank==rank_see)
		cout << "\nSolving the problem...";
	problem.solve_problem();

	// Compute the L2 Norm of the error
	if(rank==rank_see)
		cout << "\nComputing the L2 Norm of the error...";
	double L2 = 0.0;

	// Output results to file
	if(rank==rank_see)
		cout << "\nOutputting results to file...";
	char buf[50];
	sprintf(buf, "Output/IGFEM5_out_%i", mesh.n_inclusions());
	problem.writeToVTK( buf );
	
	// Output more detailed results
	if(rank==rank_see)
	{
		cout << "\n\nPARTITION " << rank_see << "\n--------------------------------------------------\n";
		cout << "Number of Elements: " << problem.get_mesh()->n_global_elem() << endl;
		cout << "Number of Nodes: " << problem.get_mesh()->n_global_nodes() << endl;
		cout << "Number of dofs: " << problem.get_dofs()->n_global_dofs() << endl;
		cout << "L2 Norm of the error: " << L2 << endl << endl;
	}


	problem.Finalize();
}
