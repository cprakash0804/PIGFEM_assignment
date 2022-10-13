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
#define Force 50e6
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





// The main function of this test is to test the mesh adaptivity functions
int main (int argc, char* argv[])
{
	Mesh mesh(&argc, &argv);
	int rank = mesh.get_rank();
	int rank_see = 0;

	// Generate a mesh
	mesh.generate_mesh(TRI3, 0.0, 1.0, 10, 0.0, 1.0, 10);

	// Add a material to the mesh
	//ContinuumDamageModelMaterial material;
	LinearElasticIsotropicMaterial material;
	material.set_parameter("E", E1);
	material.set_parameter("nu", nu1);
	material.set_name("Epoxy");
	mesh.add_material(&material);
	mesh.set_material("Epoxy");

	// Cohesive material
/*
	OPCohesiveMaterial coh_mat;
	coh_mat.set_parameter("sigma", sigmac);
	coh_mat.set_parameter("delta", deltac);
	coh_mat.set_parameter("beta", 1.0);
	coh_mat.set_name("Epoxy-Glass");
	mesh.add_material(&coh_mat);
*/

	// Assign BCs
	double app_disp = 1.0 * strain;
	BoundaryObject* boundary = mesh.get_boundary();
	boundary->set_dirichlet_bcs_from_nodeset("right", 0, app_disp);
	//boundary->set_dirichlet_bcs_from_nodeset("right", 1, x_disp);
	boundary->set_dirichlet_bcs_from_nodeset("left", 0, 0.0);
	for(unsigned int d=0; d<mesh.dim(); ++d)
		boundary->set_dirichlet_bc(0, d, 0);





	// Define the problem
	//ProblemNonlinearStructural problem;
	ProblemLinearElasticity problem;
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
