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
#include "mpi.h"

using namespace std;

#define pi 3.141592653589799323
#define Force 50e6
#define a 0.76
#define d_size 4
//#define z_size 5e-2
#define z_size 4
#define nelem 4
#define x_disp 1e0

#define mu1 1.13e9
#define mu2 31.15e9
#define nu1 0.33
#define nu2 0.22

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




int main (int argc, char* argv[])
{
	Mesh mesh(&argc, &argv);
	int rank = mesh.get_rank();
	int rank_see = 0;

	// Generate a mesh
	if(rank==rank_see)
		cout << "Generating the mesh...";
	//mesh.generate_mesh(QUAD4, -d_size, d_size, nelem, -d_size, d_size, nelem);
	mesh.generate_mesh(TET4, -d_size, d_size, nelem, -d_size, d_size, nelem, -d_size, d_size, nelem);

	// Add a material to the mesh
	ContinuumDamageModelMaterial material;
	material.set_parameter("mu", mu1);
	material.set_parameter("nu", nu1);
	material.set_parameter("P1", P1);
	material.set_parameter("P2", P2);
	material.set_parameter("Yin", Yin);
	material.set_parameter("mu_visc", mu_visc);
	material.set_name("Epoxy");
	mesh.add_material(&material);
	//mesh.set_material("Epoxy"); // This is only for non-IGFEM meshes


	// IGFEM Materials
	mesh.set_matrix_material("Epoxy");     // This is the first thing that is really inherent to IGFEM
	LinearElasticIsotropicMaterial material2;
	material2.set_parameter("mu", mu2);
	material2.set_parameter("nu", nu2);
	material2.set_name("Glass");
	mesh.add_material(&material2); // Test adding inclusion with and without this line
	if(mesh.dim()==2)
	{
		Ellipse_Inclusion ellipse;
		std::vector<double> center = {0, 0};
		ellipse.set_vec_parameter("center", center);
		ellipse.set_parameter("alpha", 0.0);
		ellipse.set_parameter("a", a);
		ellipse.set_parameter("b", a);
		ellipse.set_material(&material2);
		mesh.add_inclusion(&ellipse);	// This will add the material to the mesh is line 596? is commented out
	}
	else if(mesh.dim()==3)
	{
		Ellipsoid_Inclusion ellipsoid;
		std::vector<double> center = {0, 0, 0};
		ellipsoid.set_vec_parameter("center", center);
		ellipsoid.set_parameter("alpha", 0.0);
		ellipsoid.set_parameter("beta", 0.0);
		ellipsoid.set_parameter("gamma", 0.0);
		ellipsoid.set_parameter("a", a);
		ellipsoid.set_parameter("b", a);
		ellipsoid.set_parameter("c", a);
		ellipsoid.set_material(&material2);
		mesh.add_inclusion(&ellipsoid);	// This will add the material to the mesh is line 596? is commented out
	}
	else
		err_message("Unknown mesh dimension");


	// Assign BCs
	BoundaryObject* boundary = mesh.get_boundary();
	boundary->set_dirichlet_bcs_from_nodeset("right", 0, x_disp);
	boundary->set_dirichlet_bcs_from_nodeset("left", 0, 0.0);
	for(unsigned int d=0; d<mesh.dim(); ++d)
		boundary->set_dirichlet_bc(0, d, 0);
	if(mesh.dim()==3)			// Extra stupid BC necessary
		boundary->set_dirichlet_bc(nelem*(nelem+1), 2, 0); // back left node
	
	// Define the problem
	ProblemNonlinearStructural problem;
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
	sprintf(buf, "Output/IGFEM3_out_%i.txt", nelem);
	problem.print_solution( std::string(buf) );
	
	// Output more detailed results
	if(rank==rank_see)
	{
		cout << "\n\nPARTITION " << rank_see << "\n--------------------------------------------------\n";
		cout << "Number of Elements: " << problem.get_mesh()->n_elem() << endl;
		cout << "Number of Nodes: " << problem.get_mesh()->n_nodes() << endl;
		cout << "Number of dofs: " << problem.get_dofs()->n_global_dofs() << endl;
		cout << "L2 Norm of the error: " << L2 << endl << endl;

		// Output the solution
/*
		cout << "\nSolution:\n\n";
		for(unsigned int n=0; n<(mesh.n_local_nodes()+mesh.n_local_enrich_nodes()); ++n)
		{
			unsigned int g_node = mesh.get_node_local(n)->get_id();
			for(unsigned int d=0; d<mesh.dim(); ++d)
			{
				if(n<mesh.n_local_nodes())
					std::cout << "Node(";
				else
					std::cout << "Enrichment Node(";

				std::cout << g_node << ") dof(" << d << "): " << problem.get_solution_local(n, d) << std::endl;
			}
		}
*/
	}


	problem.Finalize();
}
