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
#include "InclusionCircle.h"
#include "InclusionEllipse.h"
#include "InclusionSphere.h"
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
#include "Material_XNCohesive.h"
#include "Material_XNCohesiveNU.h"
#include "SensitivityRightLoad.h"
#include "Solver.h"
#include "Writer.h"
#include "mpi.h"

using namespace std;


#define pi 3.141592653589799323
#define Force 50e6
#define a 2.23
#define d_size 4
//#define z_size 5e-2
#define z_size 4
#define nelem 1
#define x_disp 0.1

#define H 0.1
#define L1 0.2
#define L2 0.1
#define L3 0.25
#define Uapp 0.15

#define E1 1e9 * 1.001
#define E2 1e9
#define nu1 0.3
#define nu2 0

#define P1 0.1
#define P2 1
#define mu_visc 20
#define Yin 75e3

// Solver and damage settings
#define max_damage_order 4
#define damage_order_step 2
#define plane_strain 1
#define rel_tol 1e-8
#define abs_tol 1e-8
#define prob_rel_tol 1e-12
#define prob_abs_tol 1e-12
#define T_final 1
#define max_time_step 1e-2
#define min_time_step 1e-4
#define max_iter 20
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

	//  One Plane Axial Tension Case
	mesh.generate_mesh(TRI3, 0, 1, 1, 0, 1, 1);
	// std::string filename(argv[1]);
	// filename = "Input/" + filename;
	// mesh.read_mesh(filename);

	// Add a material to the mesh
	LinearElasticIsotropicMaterial material;
	material.set_parameter("E", E1*1.5);
	material.set_parameter("nu", nu1);
	material.set_name("Epoxy");
	mesh.add_material(&material);
	mesh.set_material("Epoxy"); // This is only for non-IGFEM meshes

	// IGFEM Materials
	Plane_Inclusion plane;
	std::vector<double> vec = {0.8201, 0.7};
	plane.set_vec_parameter("point", vec);
	vec = {-1.0,-1.0,0.0};
	plane.set_vec_parameter("normal", vec);
	material.set_name("Epoxy1");
	plane.set_material(&material); // I don't actualy care about the material difference, just the cohesive behavoir
	mesh.add_inclusion(&plane);

	vec = {0.43, 0.0};
	plane.set_vec_parameter("point", vec);
	vec = {1.0,0.1,0.0};
	plane.set_vec_parameter("normal", vec);
	material.set_name("Epoxy2");
	plane.set_material(&material); // I don't actualy care about the material difference, just the cohesive behavoir
	mesh.add_inclusion(&plane);
	// Ellipse_Inclusion ellipse;
	// // std::vector<double> vec = {-0.001, 0.5};
	// // std::vector<double> vec = {0, 0.5};
	// std::vector<double> vec = {0.001, 0.5};
	// ellipse.set_vec_parameter("center", vec);
	// ellipse.set_parameter("a", 0.7);
	// ellipse.set_parameter("b", 0.55);
	// material.set_name("Epoxy1");
	// ellipse.set_material(&material); // I don't actualy care about the material difference, just the cohesive behavoir
	// mesh.add_inclusion(&ellipse);

	// Cohesive material
	// XNCohesiveMaterial coh_mat;
	// coh_mat.set_parameter("sigma_cn", 4e6);
	// coh_mat.set_parameter("delta_cn", x_disp/5);
	// coh_mat.set_parameter("delta_ct", x_disp/5);
	// coh_mat.set_parameter("r", 0.0);
	// coh_mat.set_parameter("q", 1.0);
	OPCohesiveMaterial coh_mat;
	coh_mat.set_parameter("sigma", 4e6);
	coh_mat.set_parameter("delta", x_disp/5);
	coh_mat.set_parameter("beta", 1.0);
	coh_mat.set_name("Epoxy-Epoxy1");
	mesh.add_material(&coh_mat);

	// Assign BCs
	BoundaryObject* boundary = mesh.get_boundary();
	boundary->set_dirichlet_bcs_from_nodeset("right", 0, x_disp);
	boundary->set_dirichlet_bcs_from_nodeset("left", 0, 0.0);
	boundary->set_dirichlet_bcs_from_nodeset("left", 1, 0.0);
	// boundary->lock_down_direction(0);









	// //  Two PlaneAxial Tension Case
	// mesh.generate_mesh(QUAD4, 0, L1+L2+L3, 1, 0, H, 1);
	// // std::string filename(argv[1]);
	// // filename = "Input/" + filename;
	// // mesh.read_mesh(filename);

	// // Add a material to the mesh
	// LinearElasticIsotropicMaterial material;
	// material.set_parameter("E", E1);
	// material.set_parameter("nu", nu1);
	// material.set_name("Epoxy");
	// mesh.add_material(&material);
	// mesh.set_material("Epoxy"); // This is only for non-IGFEM meshes

	// // IGFEM Materials
	// Plane_Inclusion plane;
	// std::vector<double> vec = {L1, 0.0};
	// plane.set_vec_parameter("Point", vec);
	// vec = {1.0, 0.0};
	// plane.set_vec_parameter("Normal", vec);
	// material.set_name("Epoxy1");
	// plane.set_material(&material); // I don't actualy care about the material difference, just the cohesive behavoir
	// mesh.add_inclusion(&plane);
	// // Second plane
	// Plane_Inclusion plane2;
	// vec = {L1+L2, 0.0};
	// plane2.set_vec_parameter("Point", vec);
	// vec = {-1.0, 0.0};
	// plane2.set_vec_parameter("Normal", vec);
	// material.set_name("Epoxy2");
	// plane2.set_material(&material); // I don't actualy care about the material difference, just the cohesive behavoir
	// mesh.add_inclusion(&plane2);

	// // Cohesive material
	// OPCohesiveNUMaterial coh_mat;
	// coh_mat.set_parameter("sigma", 50e6);
	// coh_mat.set_parameter("delta", 0.03);
	// coh_mat.set_parameter("beta", 1.0);
	// coh_mat.set_name("Epoxy-Epoxy1");
	// mesh.add_material(&coh_mat);
	// mesh.add_cohesive_material_pair("Epoxy-Epoxy1", "Epoxy", "Epoxy1");
	// // Second cohesive material
	// coh_mat.set_parameter("sigma", 60e6);
	// coh_mat.set_parameter("delta", 0.03);
	// coh_mat.set_parameter("beta", 1.0);
	// coh_mat.set_name("Epoxy-Epoxy2");
	// mesh.add_material(&coh_mat);
	// mesh.add_cohesive_material_pair("Epoxy-Epoxy2", "Epoxy", "Epoxy2");

	// // Assign BCs
	// BoundaryObject* boundary = mesh.get_boundary();
	// boundary->set_dirichlet_bcs_from_nodeset("right", 0, Uapp);
	// boundary->set_dirichlet_bcs_from_nodeset("left", 0, 0.0);
	// boundary->set_dirichlet_bcs_from_nodeset("left", 1, 0.0);









	// Define the problem
	//ProblemLinearElasticity problem;
	ProblemNonlinearStructural problem;
	problem.attach_mesh(&mesh);
	problem.set_parameter("plane_strain", plane_strain);

	// // Add Sensitivity info
	// SensitivityRightLoad function;
	// problem.addSensitivityFunction(&function);
	// // problem.addMaterialSensitivityParameter("Epoxy-Epoxy", "sigma_c");
	// // problem.addMaterialSensitivityParameter("Epoxy-Epoxy", "delta_c");
	// // problem.addMaterialSensitivityParameter("Epoxy", "E");
	// // problem.addMaterialSensitivityParameter("Epoxy", "nu");
	// // problem.addMaterialSensitivityParameter("Epoxy", "mu");
	// // problem.addMaterialSensitivityParameter("Epoxy", "lambda");
	// // problem.addMaterialSensitivityParameter("Epoxy", "K");
	// problem.addShapeSensitivityParameter(0, "X");
	problem.init();

	// Set solver parameters
	Solver* solver = problem.get_solver();
	solver->set_dparameter("rel_tol", rel_tol);
	solver->set_dparameter("abs_tol", abs_tol);
	solver->set_dparameter("max_time_step", max_time_step);
	solver->set_dparameter("min_time_step", min_time_step);
	solver->set_dparameter("T_final", T_final);
	solver->set_iparameter("max_iter", max_iter);

	// Finally actually solve the problem
	if(rank==rank_see)
		cout << "\nSolving the problem...";

	// Output some problem information to the screen
	if(rank==rank_see)
	{
		cout << "\n\nPROBLEM DETAILS\n---------------------------------------------------\n";
		cout << "Problem Type: " << problem_type_names[problem.get_type()] << endl;
		cout << "Number of Mesh Partitions: " << mesh.n_ranks() << endl;
		cout << "Number of Inclusions: " << mesh.n_inclusions() << endl;
		cout << "Number of Elements: " << problem.get_mesh()->n_global_elem() << endl;
		cout << "Number of Nodes: " << problem.get_mesh()->n_global_nodes()+problem.get_mesh()->n_global_enrich_nodes()  << " (" << problem.get_mesh()->n_global_enrich_nodes() << " enriched)" << endl;
		cout << "Number of dofs: " << problem.get_dofs()->n_global_dofs()  << " (" << problem.get_dofs()->n_global_free_dofs() << " free)" << endl;
	}

	problem.solve_problem();


	problem.Finalize();
}