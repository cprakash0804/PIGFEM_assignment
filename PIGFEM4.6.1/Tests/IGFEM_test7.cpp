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
#include "mpi.h"

using namespace std;


#define pi 3.141592653589799323
#define x_disp 0.15

#define E1 1e9
#define nu1 0

#define L1 0.2
#define L2 0.1
#define L3 0.25
#define H 0.1
#define sigma1 50e6
#define sigma2 60e6
#define delta1 3e-2
#define delta2 3e-2
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

	//  Axial Tension Case
	mesh.generate_mesh(QUAD4, 0.0, L1+L2+L3, 1, 0.0, H, 1);

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
	// Second plane
	vec = {L1+L2, 0.0};
	plane.set_vec_parameter("Point", vec);
	vec = {-1.0, 0.0};
	plane.set_vec_parameter("Normal", vec);
	material.set_name("Epoxy2");
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
	// Second cohesive material
	coh_mat.set_parameter("sigma", sigma2);
	coh_mat.set_parameter("delta", delta2);
	coh_mat.set_parameter("beta", 1.0);
	coh_mat.set_name("Epoxy-Epoxy2");
	mesh.add_material(&coh_mat);
	mesh.add_cohesive_material_pair("Epoxy-Epoxy2", "Epoxy", "Epoxy2");

	// Assign BCs
	BoundaryObject* boundary = mesh.get_boundary();
	boundary->set_dirichlet_bcs_from_nodeset("right", 0, x_disp);
	//boundary->set_dirichlet_bcs_from_nodeset("right", 1, x_disp);
	boundary->set_dirichlet_bcs_from_nodeset("left", 0, 0.0);
	for(unsigned int d=0; d<mesh.dim(); ++d)
		boundary->set_dirichlet_bc(0, d, 0);


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

	// Output results to file
	if(rank==rank_see)
		cout << "\nOutputting results to file...";
	problem.print_solution( "Output/IGFEM7_out.txt" );
	
	// Output more detailed results
	if(rank==rank_see)
	{
		cout << "\n\nPARTITION " << rank_see << "\n--------------------------------------------------\n";
		cout << "Number of Elements: " << problem.get_mesh()->n_elem() << endl;
		cout << "Number of Nodes: " << problem.get_mesh()->n_nodes() << endl;
		cout << "Number of dofs: " << problem.get_dofs()->n_global_dofs() << endl;
		cout << "L2 Norm of the error: " << 0.0 << endl << endl;
	}


	problem.Finalize();
}
