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
#include "Material_OPCohesiveNU.h"
#include "Writer.h"
#include "common.h"
#include "mpi.h"

using namespace std;


#define pi 3.141592653589799323
#define x_disp 0.04

// elastic properties
#define E1 2500
#define nu1 0.0

// The size of the square domain in the positive and negative x and y directions (centers at (0,0))
//    ___________________   _
//    |        |        |    |
//    |        |        |    |
//    |        |        |    H
//    |        |        |    |
//    |________|________|   _|
//     
//    |---L1---|---L2---|

#define L1 2.0
#define L2 2.0
#define H 2.0

// cohesive properties
#define sigma1 5.0
#define delta1 (0.04/3.5)


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

	//  Axial Tension Case
	mesh.generate_mesh(QUAD4, 0.0, L1+L2, 1, 0.0, H, 1);

	// Add a matrix material to the mesh
	LinearElasticIsotropicMaterial material;
	material.set_parameter("E", E1);
	material.set_parameter("nu", nu1);
	material.set_name("Epoxy");
	mesh.add_material(&material);
	mesh.set_material("Epoxy"); // This is only for non-IGFEM meshes

	// Create Inclusions
	Plane_Inclusion plane;
	std::vector<double> vec = {L1, 0.0};
	plane.set_vec_parameter("Point", vec);
	vec = {1.0, 0.0};
	plane.set_vec_parameter("Normal", vec);
	material.set_name("Epoxy1");
	plane.set_material(&material);
	mesh.add_inclusion(&plane);

	// Cohesive material
	OPCohesiveNUMaterial coh_mat;
	coh_mat.set_parameter("sigma", sigma1);
	coh_mat.set_parameter("delta", delta1);
	coh_mat.set_parameter("beta", 1.0);
	coh_mat.set_name("Epoxy-Epoxy1");
	mesh.add_material(&coh_mat);
	mesh.add_cohesive_material_pair("Epoxy-Epoxy1", "Epoxy", "Epoxy1");

	// Assign BCs
	BoundaryObject* boundary = mesh.get_boundary();
	boundary->set_dirichlet_bcs_from_nodeset("right", 0, x_disp);
	boundary->set_dirichlet_bcs_from_nodeset("bottom", 1, 0.0);
	boundary->set_dirichlet_bcs_from_nodeset("left", 0, 0.0);
	for(unsigned int d=0; d<mesh.dim(); ++d)
		boundary->set_dirichlet_bc(0, d, 0);


	// Define the problem
	Problem * problem;
	problem = new ProblemNonlinearStructural;

	// Solve the problem
	if(rank==rank_see)
		cout << "\nSolving the problem...";
	problem->attach_mesh(&mesh);
	problem->solve_problem();


	// Output results to file
	if(rank==rank_see)
		cout << "\nOutputting results to file...";
	// problem.print_solution( "Output/IGFEM9_out.txt" );
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
