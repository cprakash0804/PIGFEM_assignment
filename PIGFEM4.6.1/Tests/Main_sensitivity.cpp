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
#include "elem.h"
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
#include "InclusionPolyhedron.h"
#include "InclusionCircle.h"
#include "InclusionSphere.h"
#include "material.h"
#include "material_cdm.h"
#include "material_lei.h"
#include "material_leti.h"
#include "material_leimps.h"
#include "material_leipcd.h"
#include "material_lt.h"
#include "Material_OPCohesive.h"
#include "Material_OPCohesiveNU.h"
#include "Material_XNCohesiveNU.h"
#include "Material_XNCohesive.h"
#include "Writer.h"
#include "Solver.h"
#include "STL_Reader.h"
#include "SubscaleModel.h"
#include "SubscaleModelLinear.h"
#include "SubscaleModelNonlinear.h"
#include "SensitivityRightLoad.h"
#include "common.h"
#include "mpi.h"

using namespace std;

// Path to the mpiexec executable
// ./petsc-3.6.2/arch-linux2-c-opt/bin/mpiexec


#define pi 3.141592653589799323

// Define the x-axial loading strain
#define strain 0.025

// Parameters to control size and number of elements in the domain
#define n_elem_per_inclusion 10	// Define the number of elements to put across the diameter of the smallest fiber
#define xexpansion 0.02
#define yexpansion 0.02
#define zexpansion 0.02
#define yexpansion_plies 0.21	// Increase the size of the domain (add matrix) by this much on each side of the domain (0.05=5%)
#define bounding_ply_limit 0.20 // the transversely isotropic bounding plies cover top and bottom 17% of the domain

// Set some flags here
#define cohesive_flag true
#define output_coh false
#define max_principal_stress_flag false
#define continuum_damage_flag false
#define bounding_plies_flag false
#define subscale_flag false

// Define matrix properties (Taken from Scott's ppt)
#define E1 2380 // MPa
#define nu1 0.43
#define sigma_max 70 // From http://www.epoxyworktops.com/epoxy-resin/mech-properties.html (MPa)

// Define carbon fiber transverse properties (Taken from Scott's ppt)
#define E2 19.5e3 // MPa
#define nu2 0.45

// Define cohesive properties
#define sigmac 10 // MPa
//#define deltac 1.8394e-3 // (mm) Taken from balancing the G_c in Scott's ppt with the integral of Ortiz-Pandolfi curve for the given sigmac
//#define deltac (61e-3)/6000 // mm
#define deltac 100e-6

// Define transversely isotrpoic properties (Taken from Scott's ppt)
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

// Solver and damage settings
#define max_damage_order 4
#define damage_order_step 2
#define plane_strain 1
#define rel_tol 1e-8
#define abs_tol 1e-8
#define prob_rel_tol 1e-6
#define prob_abs_tol 1e-6
#define T_final 1
#define max_time_step 1e-2
#define min_time_step 1e-5
#define max_iter 20






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

void add_circle(Mesh* mesh, Inclusion* circle, double x, double y, double r)
{
	std::vector<double> center = {x, y};
	circle->set_vec_parameter("center", center);
	circle->set_parameter("r", r);
	mesh->add_inclusion(circle);
}

void add_sphere(Mesh* mesh, Inclusion* sphere, double x, double y, double z, double r)
{
	std::vector<double> center = {x, y, z};
	sphere->set_vec_parameter("center", center);
	sphere->set_parameter("r", r);
	mesh->add_inclusion(sphere);
}





void read_circles(Mesh* mesh, std::string filename, std::vector<double>& maxes, std::vector<double>& mins, double& min_rad, id_type& dimension)
{
	// Currently only support reading in 2D fiber geometries
	dimension = 2;
	maxes = {-99999999999999999999999999., -99999999999999999999999999.};
	mins = {99999999999999999999999999., 99999999999999999999999999.};
	min_rad = 99999999999999999999999999.;

	LinearElasticIsotropicMaterial fiber;
	fiber.set_parameter("E", E2);
	fiber.set_parameter("nu", nu2);
	fiber.set_name("Carbon");

	ifstream myfile;
	myfile.open( filename );
	if (!myfile.good())
		err_message("Invalid input file name!");

	std::string dummy;
	char c;

	// Read in the circles (x, y, r)
	double x, y, z, r;
	Circle_Inclusion circle;
	Sphere_Inclusion sphere;
	circle.set_material(&fiber);
	sphere.set_material(&fiber);
	while (myfile.peek() != EOF)
	{
		// Read the data from the file
		myfile >> x >> c >> y >> c >> r;
		myfile.get(c);
		if (!(c == '\n' || c == '\r') && myfile.peek() != EOF) // This is reading in spheres instead of circles (so c should be a comma currently) (if its the end of a circle file c will still be a comma so need to check EOF)
		{
			z = r;
			myfile >> r;
			add_sphere(mesh, &sphere, x, y, z, r);
			dimension = 3;
			if (maxes.size() == 2)
			{
				maxes.push_back(-99999999999999999999999999.);
				mins.push_back(99999999999999999999999999.);
			}
			maxes[2] = std::max(maxes[2], z+r);
			mins[2] = std::min(mins[2], z-r);

			std::getline(myfile, dummy); // Read in the endline
		}
		else
			add_circle(mesh, &circle, x, y, r); // Don't need to read in the endline already got it

		// Update the x and y maximal and minimal coordinates
		maxes[0] = std::max(maxes[0], x+r);
		maxes[1] = std::max(maxes[1], y+r);
		mins[0] = std::min(mins[0], x-r);
		mins[1] = std::min(mins[1], y-r);
		min_rad = std::min(min_rad, r);
	}

	myfile.close();
}







void read_stl(Mesh* mesh, std::string filename, std::vector<double>& maxes, std::vector<double>& mins, double& min_rad, id_type& dimension)
{
	// Actually read the file
	STL_Reader reader;
	reader.read(filename);

	// Create the Polyhedral inclusion
	LinearElasticIsotropicMaterial poly;
	poly.set_parameter("E", E2);
	poly.set_parameter("nu", nu2);
	poly.set_name("Carbon");
	Polyhedron_Inclusion polyhedron;
	polyhedron.set_material(&poly);
	polyhedron.set_facets(reader.get_facets());
	polyhedron.set_normals(reader.get_normals());
	mesh->add_inclusion(&polyhedron);

	std::vector<double> bounds = polyhedron.get_domain_lims();
	maxes = {bounds[3], bounds[4], bounds[5]};
	mins = {bounds[0], bounds[1], bounds[2]};
	std::vector<double> sizes = {bounds[3]-bounds[0], bounds[4]-bounds[1], bounds[5]-bounds[2]};
	min_rad = *std::min_element(sizes.begin(), sizes.end());
	min_rad /= 2;
	dimension = 3;
}





















// The main function of this test is to test the mesh adaptivity functions
int main (int argc, char* argv[])
{
	// Make sure that a file name was actualy passed in first thing
	// if (argc <= 1)
	// 	err_message("Please input a filename to read inclusions from.");

	if (!subscale_flag)
	{
		Mesh mesh(&argc, &argv);
		int rank = mesh.get_rank();
		int rank_see = 0;

		// Read in the inclusions (Add fiber material in this function)
		if(rank==rank_see)
			cout << "Reading inclusions...\n";
		std::vector<double> maxes;	// Maximum domain limits (x and y)
		std::vector<double> mins;	// Minimum domain limits (x and y)
		double min_rad;				// Minimum fiber radius (So I know how many elements across to use)
		id_type dimension;
		std::string filename(argv[1]);
		filename = "Input/" + filename;

		// Check what kind of file I'm reading in
		if (filename.find(".csv", filename.length()-4) != std::string::npos ||
			filename.find(".txt", filename.length()-4) != std::string::npos)
		{
			read_circles(&mesh, filename, maxes, mins, min_rad, dimension);
		}
		else if (filename.find(".stl", filename.length()-4) != std::string::npos)
		{
			read_stl(&mesh, filename, maxes, mins, min_rad, dimension);
		}
		else
			err_message("Unkonwn file type!");
		

		// Generate a mesh
		if(rank==rank_see)
			cout << "Generating the mesh...";

		double x_lims = maxes[0] - mins[0];
		double y_lims = maxes[1] - mins[1];
		if (dimension == 2)
		{
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
		}
		else if (dimension ==  3)
		{
			double z_lims = maxes[2] - mins[2];
			mins[0] = mins[0] - x_lims*xexpansion;
			maxes[0] = maxes[0] + x_lims*xexpansion;
			mins[2] = mins[2] - z_lims*zexpansion;
			maxes[2] = maxes[2] + z_lims*zexpansion;
			id_type x_elem = std::ceil((maxes[0]-mins[0]) / (min_rad*2)) * n_elem_per_inclusion;
			id_type z_elem = std::ceil((maxes[2]-mins[2]) / (min_rad*2)) * n_elem_per_inclusion;
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
			mesh.generate_mesh(TET4, mins[0], maxes[0], x_elem, mins[1], maxes[1], y_elem, mins[2], maxes[2], z_elem); // Actually generate the mesh
		}
		else
			err_message("Invalid mesh dimension somehow.");


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
			OPCohesiveMaterial coh_mat;
			coh_mat.set_parameter("sigma", sigmac);
			coh_mat.set_parameter("delta", deltac);
			coh_mat.set_parameter("beta", 1.0);

			// XNCohesiveMaterial coh_mat;
			// coh_mat.set_parameter("sigma_n", sigmac);
			// // coh_mat.set_parameter("delta_n", deltac);
			// // coh_mat.set_parameter("delta_t", deltac);
			// coh_mat.set_parameter("delta_c", deltac);
			// coh_mat.set_parameter("q", 1.0);

			coh_mat.set_name("Epoxy-Carbon");
			mesh.add_material(&coh_mat);
		}

		// Assign BCs
		double app_disp = (maxes[0]-mins[0]) * strain;
		BoundaryObject* boundary = mesh.get_boundary();
		boundary->set_dirichlet_bcs_from_nodeset("right", 0, app_disp);
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
			problem->set_parameter("max_damage_order", max_damage_order);
			problem->set_parameter("damage_step", damage_order_step);
			problem->set_parameter("rel_tol", prob_rel_tol);
			problem->set_parameter("abs_tol", prob_abs_tol);
		}
		else if (continuum_damage_flag || cohesive_flag)
			problem = new ProblemNonlinearStructural;
		else // linear elastic
			problem = new ProblemLinearElasticity;
		problem->set_parameter("plane_strain", plane_strain);

		// Set up the problem
		if(rank==rank_see)
			cout << "\nSetting up the problem...";
		problem->output_cohesive() = output_coh;
		problem->attach_mesh(&mesh);

		// Add Sensitivity info
		SensitivityRightLoad function;
		problem->addSensitivityFunction(&function);
		problem->addMaterialSensitivityParameter("Epoxy-Carbon", "sigma_c");
		problem->addMaterialSensitivityParameter("Epoxy-Carbon", "delta_c");
		// problem->addShapeSensitivityParameter(0, "X");
		problem->init();

		// Set solver parameters
		Solver* solver = problem->get_solver();
		solver->set_dparameter("rel_tol", rel_tol);
		solver->set_dparameter("abs_tol", abs_tol);
		solver->set_dparameter("max_time_step", max_time_step);
		solver->set_dparameter("min_time_step", min_time_step);
		solver->set_dparameter("T_final", T_final);
		solver->set_iparameter("max_iter", max_iter);

		// Output some problem information to the screen
		if(rank==rank_see)
		{
			cout << "\n\nPROBLEM DETAILS\n---------------------------------------------------\n";
			cout << "Problem Type: " << problem_type_names[problem->get_type()] << endl;
			cout << "Number of Mesh Partitions: " << mesh.n_ranks() << endl;
			cout << "Number of Inclusions: " << mesh.n_inclusions() << endl;
			cout << "Number of Elements: " << problem->get_mesh()->n_global_elem() << endl;
			cout << "Number of Nodes: " << problem->get_mesh()->n_global_nodes()+problem->get_mesh()->n_global_enrich_nodes()  << " (" << problem->get_mesh()->n_global_enrich_nodes() << " enriched)" << endl;
			cout << "Number of dofs: " << problem->get_dofs()->n_global_dofs()  << " (" << problem->get_dofs()->n_global_free_dofs() << " free)" << endl;
		}

		// Finally actually solve the problem
		problem->solve_problem();

		// Cleanup
		problem->Finalize();
		delete problem;
	}













	// Running a subscale model
	else
	{
		if (max_principal_stress_flag)
			err_message("Matrix damage is not allowed in a subscale model!");

		SubscaleModel* model;
		if (cohesive_flag || continuum_damage_flag)
			model = new SubscaleModelNonlinear(&argc, &argv);
		else
			model = new SubscaleModelLinear(&argc, &argv);

		Mesh* mesh = model->get_mesh();
		int rank = mesh->get_rank();
		int rank_see = 0;

		// Read in the inclusions (Add fiber material in this function)
		if(rank==rank_see)
			cout << "Reading inclusions...\n";
		std::vector<double> maxes;	// Maximum domain limits (x and y)
		std::vector<double> mins;	// Minimum domain limits (x and y)
		double min_rad;				// Minimum fiber radius (So I know how many elements across to use)
		id_type dimension;
		std::string filename(argv[1]);
		filename = "Input/" + filename;

		// Check what kind of file I'm reading in
		if (filename.find(".csv", filename.length()-4) != std::string::npos ||
			filename.find(".txt", filename.length()-4) != std::string::npos)
		{
			read_circles(mesh, filename, maxes, mins, min_rad, dimension);
		}
		else if (filename.find(".stl", filename.length()-4) != std::string::npos)
		{
			read_stl(mesh, filename, maxes, mins, min_rad, dimension);
		}
		else
			err_message("Unkonwn file type!");
		

		// Generate a mesh
		if(rank==rank_see)
			cout << "Generating the mesh...";

		double x_lims = maxes[0] - mins[0];
		double y_lims = maxes[1] - mins[1];
		if (dimension == 2)
		{
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
			mesh->generate_mesh(TRI3, mins[0], maxes[0], x_elem, mins[1], maxes[1], y_elem); // Actually generate the mesh
		}
		else if (dimension ==  3)
		{
			double z_lims = maxes[2] - mins[2];
			mins[0] = mins[0] - x_lims*xexpansion;
			maxes[0] = maxes[0] + x_lims*xexpansion;
			mins[2] = mins[2] - z_lims*zexpansion;
			maxes[2] = maxes[2] + z_lims*zexpansion;
			id_type x_elem = std::ceil((maxes[0]-mins[0]) / (min_rad*2)) * n_elem_per_inclusion;
			id_type z_elem = std::ceil((maxes[2]-mins[2]) / (min_rad*2)) * n_elem_per_inclusion;
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
			mesh->generate_mesh(TET4, mins[0], maxes[0], x_elem, mins[1], maxes[1], y_elem, mins[2], maxes[2], z_elem); // Actually generate the mesh
		}
		else
			err_message("Invalid mesh dimension somehow.");


		// Add the matrix to every element (linear elastic and fails at a max priniciple stress)
		if (continuum_damage_flag)
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

		// If we are doing bounding plies, add them here
		if (bounding_plies_flag)
		{
			// Create the element sets that are defined as the top and bottom ?% of the mesh
			std::set<id_type> top_ply, bottom_ply;
			double top_limit = maxes[1] - y_lims*bounding_ply_limit;
			double bottom_limit = mins[1] + y_lims*bounding_ply_limit;
			// Loop over every active element and see if all of its nodes lie within the top or bottom x% of the domain
			for (Mesh::element_iterator it=mesh->active_elements_begin(), end=mesh->active_elements_end(); it!=end; ++it)
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
			mesh->add_elemset("top_ply", top_ply);
			mesh->add_elemset("bottom_ply", bottom_ply);

			// Create the transversely isotropic material to the top and the bottom (Overwrites the matrix set earlier)
			LinearElasticTransverselyIsotropicMaterial bounds;
			bounds.set_parameter("E1", Eti1);
			bounds.set_parameter("E2", Eti2);
			bounds.set_parameter("G12", G12);
			bounds.set_parameter("G23", G23);
			bounds.set_parameter("nu12", nu12);
			bounds.set_name("0_degree_plies");
			mesh->add_material(&bounds);
			mesh->set_material_from_elemset("0_degree_plies", "top_ply");
			mesh->set_material_from_elemset("0_degree_plies", "bottom_ply");
		}

		if (cohesive_flag)
		{
			// Create the cohesive material (Will automatically be applied to any material boundary
			//	that doesn't have another cohesive material specified for the touching material pair)
			// OPCohesiveMaterial coh_mat;
			// coh_mat.set_parameter("sigma", sigmac);
			// coh_mat.set_parameter("delta", deltac);
			// coh_mat.set_parameter("beta", 1.0);

			XNCohesiveMaterial coh_mat;
			coh_mat.set_parameter("sigma_n", sigmac);
			coh_mat.set_parameter("delta_n", deltac);
			coh_mat.set_parameter("delta_t", deltac);
			coh_mat.set_parameter("q", 1.0);
			// coh_mat.set_parameter("tau_c", sigma_c);

			coh_mat.set_name("Epoxy-Glass");
			mesh->add_material(&coh_mat);
		}

		// Assign BCs
		// No need to assign BCs here, that's done within the subscale model (Do I do anything about the zero node?)





		// Define the problem
		Problem* problem = model->get_problem();
		problem->set_parameter("plane_strain", plane_strain);

		// Set up the problem
		if(rank==rank_see)
			cout << "\nSetting up the problem...";
		problem->output_cohesive() = output_coh;
		model->init();

		// Set solver parameters
		Solver* solver = problem->get_solver();
		solver->set_dparameter("rel_tol", rel_tol);
		solver->set_dparameter("abs_tol", abs_tol);
		solver->set_dparameter("max_time_step", max_time_step);
		solver->set_dparameter("min_time_step", min_time_step);
		solver->set_dparameter("T_final", T_final);
		solver->set_iparameter("max_iter", max_iter);

		// Define the macroscopic strain
		if (dimension == 2)
		{
			std::vector<double> strain_load = {strain, 0.0, 0.0};
			model->set_macro_strain(strain_load);
		}
		else
		{
			std::vector<double> strain_load = {strain, 0.0, 0.0, 0.0, 0.0, 0.0};
			model->set_macro_strain(strain_load);
		}

		// Output some problem information to the screen
		if(rank==rank_see)
		{
			cout << "\n\nPROBLEM DETAILS\n---------------------------------------------------\n";
			cout << "Problem Type: " << problem_type_names[problem->get_type()] << endl;
			cout << "Number of Mesh Partitions: " << mesh->n_ranks() << endl;
			cout << "Number of Inclusions: " << mesh->n_inclusions() << endl;
			cout << "Number of Elements: " << problem->get_mesh()->n_global_elem() << endl;
			cout << "Number of Nodes: " << problem->get_mesh()->n_global_nodes()+problem->get_mesh()->n_global_enrich_nodes()  << " (" << problem->get_mesh()->n_global_enrich_nodes() << " enriched)" << endl;
			cout << "Number of dofs: " << problem->get_dofs()->n_global_dofs()  << " (" << problem->get_dofs()->n_global_free_dofs() << " free)" << endl;
		}


		// Finally actually solve the problem
		model->solve();

		// Cleanup
		delete model;
	}
}
