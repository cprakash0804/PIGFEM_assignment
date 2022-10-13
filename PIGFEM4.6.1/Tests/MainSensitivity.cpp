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
#include "Options.h"
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
#include "Utilities.h"
#include "mpi.h"
#include <algorithm>
#include <cctype>

using namespace std;

// Path to the mpiexec executable
// ./petsc-3.6.2/arch-linux2-c-opt/bin/mpiexec


#define pi 3.141592653589799323

// Define the x-axial loading strain
#define strain 0.01

// Parameters to control size and number of elements in the domain
#define n_elem_per_inclusion 10	// Define the number of elements to put across the diameter of the smallest fiber
#define xexpansion 0.02
#define yexpansion 0.02
#define zexpansion 0.02
#define yexpansion_plies 0.21	// Increase the size of the domain (add matrix) by this much on each side of the domain (0.05=5%)
#define bounding_ply_limit 0.20 // the transversely isotropic bounding plies cover top and bottom 17% of the domain

// Set some flags here
#define cohesive_flag true
#define max_principal_stress_flag false
#define continuum_damage_flag false
#define bounding_plies_flag false
#define subscale_flag false

// Define matrix properties (Taken from Scott's ppt)
#define E1 2380 // MPa
#define nu1 0.43
#define sigma_max 70 // From Sottos

// Define carbon fiber transverse properties (Taken from Scott's ppt)
#define E2 19.5e3 // MPa
#define nu2 0.45

// Define cohesive properties
#define sigmac 10 // MPa
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
#define rel_tol 1e-12
#define abs_tol 1e-12
#define prob_rel_tol 1e-6
#define prob_abs_tol 1e-6
#define T_final 1
#define max_time_step 2e-2
#define min_time_step 1e-5
#define max_iter 20







void read_csv(Mesh* mesh, std::string filename, std::vector<double>& maxes, std::vector<double>& mins, double& min_rad, int& dimension,
			  std::vector<double>& sigma_c, std::vector<double>& delta_c, bool& individual);

void read_stl(Mesh* mesh, std::string filename, std::vector<double>& maxes, std::vector<double>& mins, double& min_rad, int& dimension,
			  std::vector<double>& sigma_c, std::vector<double>& delta_c, bool& individual);



// The main function of this test is to test the mesh adaptivity functions
int main (int argc, char* argv[])
{
	// Define the mesh
	Mesh mesh(&argc, &argv);
	int rank = mesh.get_rank();
	int rank_see = 0;

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
	problem->attach_mesh(&mesh);


	// Parse the command line input for all command line arguments
	problem->getOptions()->parseInput(argc, argv);
	std::string in_file = problem->getOptions()->getOption("-IN");
	if (in_file == "")
	{
		delete problem;
		err_message("Must provide an input file with the -in option");
	}


	// Read in the inclusions (Add fiber material in this function)
	bool write_to_screen = problem->getOptions()->getBoolOption("-WRITE_TO_SCREEN");
	if(rank==rank_see && write_to_screen)
		cout << "Reading inclusions...\n";
	std::vector<double> maxes;	// Maximum domain limits (x, y, and z)
	std::vector<double> mins;	// Minimum domain limits (x, y, and z)
	double min_rad;				// Minimum fiber radius (So I know how many elements across to use)
	int dimension;
	bool individual;
	bool subscale = problem->getOptions()->getBoolOption("-SUBSCALE");
	std::vector<double> sigma_c, delta_c; // only used if we're reading in indivdual fiber properties

	// Check what kind of file I'm reading in
	if (in_file.find(".csv", in_file.length()-4) != std::string::npos ||
		in_file.find(".txt", in_file.length()-4) != std::string::npos)
		read_csv(&mesh, in_file, maxes, mins, min_rad, dimension, sigma_c, delta_c, individual);
	else if (in_file.find(".stl", in_file.length()-4) != std::string::npos)
		read_stl(&mesh, in_file, maxes, mins, min_rad, dimension, sigma_c, delta_c, individual);
	else
		err_message("Unkonwn file type!");
	


	// Generate a mesh
	{
		if(rank==rank_see && write_to_screen)
			cout << "Generating the mesh...";

		if (problem->getOptions()->getOption("-MESH") == "") // No mesh file provided. Generate a mesh
		{
			double x_lims = maxes[0] - mins[0];
			mins[0] = mins[0] - x_lims*xexpansion;
			maxes[0] = maxes[0] + x_lims*xexpansion;
			id_type x_elem = std::ceil((maxes[0]-mins[0]) / (min_rad*2)) * n_elem_per_inclusion;
			double y_lims = maxes[1] - mins[1];
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
			maxes = {0.00925, 0.008689};
			mins = {-0.00375, -0.003739}; // Fix the boundaries for the shape sensitivity
			if (dimension == 2)
				mesh.generate_mesh(TRI3, mins[0], maxes[0], x_elem, mins[1], maxes[1], y_elem); // Actually generate the mesh
			
			else if (dimension ==  3)
			{
				double z_lims = maxes[2] - mins[2];
				mins[2] = mins[2] - z_lims*zexpansion;
				maxes[2] = maxes[2] + z_lims*zexpansion;
				id_type z_elem = std::ceil((maxes[2]-mins[2]) / (min_rad*2)) * n_elem_per_inclusion;
				mesh.generate_mesh(TET4, mins[0], maxes[0], x_elem, mins[1], maxes[1], y_elem, mins[2], maxes[2], z_elem); // Actually generate the mesh
			}
			else
				err_message("Invalid mesh dimension somehow.");

			// If we are doing bounding plies, add them here (Can only do this if we are generating a mesh. If done with ana abaqus meshI guess we could do it with a plane inclusion)
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
			}
		}
		else
		{
			std::string abaqus_filename(argv[2]);
			abaqus_filename = "Input/" + abaqus_filename;
			mesh.read_mesh(abaqus_filename);
		}
	}



	// Assign materials
	{
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
		if (bounding_plies_flag)
		{
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
			// OPCohesiveMaterial coh_mat;
			// coh_mat.set_parameter("sigma", sigmac);
			// coh_mat.set_parameter("delta", deltac);
			// coh_mat.set_parameter("beta", 1.0);

			XNCohesiveMaterial coh_mat;
			coh_mat.set_parameter("SIGMA_CN", sigmac);
			coh_mat.set_parameter("DELTA_C", deltac);
			// coh_mat.set_parameter("delta_n", deltac);
			// coh_mat.set_parameter("delta_t", deltac);
			coh_mat.set_parameter("q", 1.0);	// Ratio of normal to tangential strengths

			if (!individual)
			{
				coh_mat.set_name("Epoxy-Carbon");
				mesh.add_material(&coh_mat);
			}
			else
			{
				std::stringstream ss;
				for (id_type m=0; m<sigma_c.size(); ++m)
				{
					// Set the individual cohesive name
					ss.str("");
					ss << m;
					std::string carbon_name = "Carbon" + ss.str();
					std::string coh_mat_name = "Epoxy-" + carbon_name;
					coh_mat.set_name(coh_mat_name);

					// Set the individual properties
					coh_mat.set_parameter("SIGMA_CN", sigma_c[m]);
					coh_mat.set_parameter("DELTA_C", delta_c[m]);
					// if (coh_mat.get_type()==OP_COHESIVE || coh_mat.get_type()==OP_COHESIVE_NO_UNLOADING)
					// 	coh_mat.set_parameter("DELTA_C", delta_c[m]);
					// else if (coh_mat.get_type()==XN_COHESIVE || coh_mat.get_type()==XN_COHESIVE_NO_UNLOADING)
					// {
					// 	coh_mat.set_parameter("DELTA_CN", delta_c[m]);
					// 	coh_mat.set_parameter("DELTA_CT", delta_c[m]);
					// }
					// else
					// 	err_message("Unknown cohesive material!");

					// Add info to the mesh
					mesh.add_material(&coh_mat);
					mesh.add_cohesive_material_pair(coh_mat_name, "Epoxy", carbon_name);
				}
			}
		}
	}



	// Assign BCs
	{
		if (subscale)
		{
			// std::vector<int> voigt = {1,3,6};
			// std::vector<double> strain_load(voigt[dimension-1]);
			// strain_load[0] = strain;
			// 	model->set_macro_strain(strain_load);
		}
		else
		{
			double app_disp = (maxes[0]-mins[0]) * strain;
			BoundaryObject* boundary = mesh.get_boundary();
			boundary->set_dirichlet_bcs_from_nodeset("right", 0, app_disp);
			//boundary->set_dirichlet_bcs_from_nodeset("right", 1, x_disp);
			boundary->set_dirichlet_bcs_from_nodeset("left", 0, 0.0);
			boundary->set_dirichlet_bcs_from_nodeset("bottom", 1, 0.0);
			for(unsigned int d=0; d<mesh.dim(); ++d)
				boundary->set_dirichlet_bc(0, d, 0);
		}
	}



	// Add Sensitivity info
	{
		if (problem->getOptions()->getBoolOption("-SENSITIVITY"))
		{
			SensitivityRightLoad function;
			problem->addSensitivityFunction(&function);
			// if (individual)
			// {
			// 	std::stringstream ss;
			// 	for (id_type m=0; m<sigma_c.size(); ++m)
			// 	{
			// 		// Set the individual cohesive name
			// 		ss.str("");
			// 		ss << m;
			// 		std::string coh_mat_name = "Epoxy-Carbon" + ss.str();
			// 		problem->addMaterialSensitivityParameter(coh_mat_name, "sigma");
			// 		// problem->addMaterialSensitivityParameter(coh_mat_name, "delta");
			// 	}
			// }
			// else
			// {
			// 	problem->addMaterialSensitivityParameter("Epoxy-Carbon", "sigma");
			// 	// problem->addMaterialSensitivityParameter("Epoxy-Carbon", "delta");
			// }

			// Shape sensitivity parameters
			for (id_type i=0; i<mesh.n_inclusions(); ++i)
			{
				problem->addShapeSensitivityParameter(i, "X");
				problem->addShapeSensitivityParameter(i, "Y");
				problem->addShapeSensitivityParameter(i, "R");
			}
		}
	}

	// Set up the problem
	if(rank==rank_see && write_to_screen)
		cout << "\nSetting up the problem...";
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
	if(rank==rank_see && write_to_screen)
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
	if (subscale)
	{
		// model->solve();
		// delete model;
	}
	else
	{
		problem->solve_problem();
		delete problem;
	}




















































	// // Running a subscale model
	// else
	// {
	// 	if (max_principal_stress_flag)
	// 		err_message("Matrix damage is not allowed in a subscale model!");

	// 	SubscaleModel* model;
	// 	if (cohesive_flag || continuum_damage_flag)
	// 		model = new SubscaleModelNonlinear(&argc, &argv);
	// 	else
	// 		model = new SubscaleModelLinear(&argc, &argv);

	// 	Mesh* mesh = model->get_mesh();
	// 	int rank = mesh->get_rank();
	// 	int rank_see = 0;

	// 	// Read in the inclusions (Add fiber material in this function)
	// 	if(rank==rank_see)
	// 		cout << "Reading inclusions...\n";
	// 	std::vector<double> maxes;	// Maximum domain limits (x and y)
	// 	std::vector<double> mins;	// Minimum domain limits (x and y)
	// 	double min_rad;				// Minimum fiber radius (So I know how many elements across to use)
	// 	int dimension;
	// 	std::string in_file(argv[1]);
	// 	in_file = "Input/" + in_file;

	// 	// Check what kind of file I'm reading in
	// 	std::vector<double> sigma_c, delta_c; // Only sed if we're using individualized properties
	// 	bool individual;
	// 	if (in_file.find(".csv", in_file.length()-4) != std::string::npos ||
	// 		in_file.find(".txt", in_file.length()-4) != std::string::npos)
	// 		read_csv(mesh, in_file, maxes, mins, min_rad, dimension, sigma_c, delta_c, individual);
	// 	else if (in_file.find(".stl", in_file.length()-4) != std::string::npos)
	// 		read_stl(mesh, in_file, maxes, mins, min_rad, dimension, sigma_c, delta_c, individual);
	// 	else
	// 		err_message("Unkonwn file type!");
		

	// 	if (argc == 2)
	// 	{
	// 		// Generate a mesh
	// 		if(rank==rank_see)
	// 			cout << "Generating the mesh...";

	// 		double x_lims = maxes[0] - mins[0];
	// 		mins[0] = mins[0] - x_lims*xexpansion;
	// 		maxes[0] = maxes[0] + x_lims*xexpansion;
	// 		id_type x_elem = std::ceil((maxes[0]-mins[0]) / (min_rad*2)) * n_elem_per_inclusion;
	// 		double y_lims = maxes[1] - mins[1];
	// 		id_type y_elem;
	// 		if (bounding_plies_flag)
	// 		{
	// 			mins[1] = mins[1] - y_lims*yexpansion_plies;
	// 			maxes[1] = maxes[1] + y_lims*yexpansion_plies;
	// 			y_elem = std::ceil((maxes[1]-mins[1]) / (min_rad*2)) * n_elem_per_inclusion;
	// 		}
	// 		else
	// 		{
	// 			mins[1] = mins[1] - y_lims*yexpansion;
	// 			maxes[1] = maxes[1] + y_lims*yexpansion;
	// 			y_elem = std::ceil((maxes[1]-mins[1]) / (min_rad*2)) * n_elem_per_inclusion;
	// 		}
	// 		if (dimension == 2)
	// 			mesh->generate_mesh(TRI3, mins[0], maxes[0], x_elem, mins[1], maxes[1], y_elem); // Actually generate the mesh
			
	// 		else if (dimension ==  3)
	// 		{
	// 			double z_lims = maxes[2] - mins[2];
	// 			mins[2] = mins[2] - z_lims*zexpansion;
	// 			maxes[2] = maxes[2] + z_lims*zexpansion;
	// 			id_type z_elem = std::ceil((maxes[2]-mins[2]) / (min_rad*2)) * n_elem_per_inclusion;
	// 			mesh->generate_mesh(TET4, mins[0], maxes[0], x_elem, mins[1], maxes[1], y_elem, mins[2], maxes[2], z_elem); // Actually generate the mesh
	// 		}
	// 		else
	// 			err_message("Invalid mesh dimension somehow.");

	// 		// If we are doing bounding plies, add them here (Can only do this if we are generating a mesh. If done with ana abaqus meshI guess we could do it with a plane inclusion)
	// 		if (bounding_plies_flag)
	// 		{
	// 			// Create the element sets that are defined as the top and bottom ?% of the mesh
	// 			std::set<id_type> top_ply, bottom_ply;
	// 			double top_limit = maxes[1] - y_lims*bounding_ply_limit;
	// 			double bottom_limit = mins[1] + y_lims*bounding_ply_limit;
	// 			// Loop over every active element and see if all of its nodes lie within the top or bottom x% of the domain
	// 			for (Mesh::element_iterator it=mesh->active_elements_begin(), end=mesh->active_elements_end(); it!=end; ++it)
	// 			{
	// 				bool in_top_ply = true;
	// 				bool in_bottom_ply = true;
	// 				for (id_type n=0; n<(*it)->n_nodes(); ++n)
	// 				{
	// 					std::vector<double> coords = (*it)->get_node(n)->get_coords();
	// 					double y = coords[1];
	// 					if (y < top_limit)
	// 						in_top_ply = false;
	// 					if(y > bottom_limit)
	// 						in_bottom_ply = false;
	// 				}
	// 				if (in_top_ply)
	// 					top_ply.insert((*it)->get_id());
	// 				if (in_bottom_ply)
	// 					bottom_ply.insert((*it)->get_id());
	// 			}
	// 			mesh->add_elemset("top_ply", top_ply);
	// 			mesh->add_elemset("bottom_ply", bottom_ply);

	// 			// Create the transversely isotropic material to the top and the bottom (Overwrites the matrix set earlier)
	// 			LinearElasticTransverselyIsotropicMaterial bounds;
	// 			bounds.set_parameter("E1", Eti1);
	// 			bounds.set_parameter("E2", Eti2);
	// 			bounds.set_parameter("G12", G12);
	// 			bounds.set_parameter("G23", G23);
	// 			bounds.set_parameter("nu12", nu12);
	// 			bounds.set_name("0_degree_plies");
	// 			mesh->add_material(&bounds);
	// 			mesh->set_material_from_elemset("0_degree_plies", "top_ply");
	// 			mesh->set_material_from_elemset("0_degree_plies", "bottom_ply");
	// 		}
	// 	}
	// 	else if (argc == 3)
	// 		err_message("Subscale models do not support Abaqus mesh files!");
	// 	else
	// 		err_message("Unknown number of command line parameters!");


	// 	// Add the matrix to every element (linear elastic and fails at a max priniciple stress)
	// 	if (continuum_damage_flag)
	// 	{
	// 		ContinuumDamageModelMaterial matrix;
	// 		matrix.set_parameter("E", E1);
	// 		matrix.set_parameter("nu", nu1);
	// 		matrix.set_parameter("P1", P1);
	// 		matrix.set_parameter("P2", P2);
	// 		matrix.set_parameter("mu_visc", mu_visc);
	// 		matrix.set_parameter("Yin", Yin);
	// 		matrix.set_name("Epoxy");
	// 		mesh->add_material(&matrix);
	// 		mesh->set_material("Epoxy");
	// 	}
	// 	else // Linear elastic
	// 	{
	// 		LinearElasticIsotropicMaterial matrix;
	// 		matrix.set_parameter("E", E1);
	// 		matrix.set_parameter("nu", nu1);
	// 		matrix.set_name("Epoxy");
	// 		mesh->add_material(&matrix);
	// 		mesh->set_material("Epoxy");
	// 	}

	// 	if (cohesive_flag)
	// 	{
	// 		// Create the cohesive material (Will automatically be applied to any material boundary
	// 		//	that doesn't have another cohesive material specified for the touching material pair)
	// 		OPCohesiveMaterial coh_mat;
	// 		coh_mat.set_parameter("sigma", sigmac);
	// 		coh_mat.set_parameter("delta", deltac);
	// 		coh_mat.set_parameter("beta", 1.0);
	// 		coh_mat.set_name("Epoxy-Glass");
	// 		mesh->add_material(&coh_mat);


	// 		if (!individual)
	// 		{
	// 			coh_mat.set_name("Epoxy-Carbon");
	// 			mesh->add_material(&coh_mat);
	// 		}
	// 		else
	// 		{
	// 			std::stringstream ss;
	// 			for (id_type m=0; m<sigma_c.size(); ++m)
	// 			{
	// 				// Set the individual cohesive name
	// 				ss << m;
	// 				std::string carbon_name = "Carbon" + ss.str();
	// 				std::string coh_mat_name = "Epoxy-" + carbon_name;
	// 				coh_mat.set_name(coh_mat_name);

	// 				// Set the individual properties
	// 				coh_mat.set_parameter("SIGMA_C", sigma_c[m]);
	// 				if (coh_mat.get_type()==OP_COHESIVE || coh_mat.get_type()==OP_COHESIVE_NO_UNLOADING)
	// 					coh_mat.set_parameter("DELTA_C", delta_c[m]);
	// 				else if (coh_mat.get_type()==XN_COHESIVE || coh_mat.get_type()==XN_COHESIVE_NO_UNLOADING)
	// 				{
	// 					coh_mat.set_parameter("DELTA_CN", delta_c[m]);
	// 					coh_mat.set_parameter("DELTA_CT", delta_c[m]);
	// 				}
	// 				else
	// 					err_message("Unknown cohesive material!");

	// 				// Add info to the mesh
	// 				mesh->add_material(&coh_mat);
	// 				mesh->add_cohesive_material_pair(coh_mat_name, "Epoxy", carbon_name);
	// 			}
	// 		}
	// 	}

	// 	// Assign BCs
	// 	// No need to assign BCs here, that's done within the subscale model (Do I do anything about the zero node?)





	// 	// Define the problem
	// 	Problem* problem = model->get_problem();
	// 	problem->set_parameter("plane_strain", plane_strain);

	// 	// Set up the problem
	// 	if(rank==rank_see)
	// 		cout << "\nSetting up the problem...";
	// 	model->init();

	// 	// Set solver parameters
	// 	Solver* solver = problem->get_solver();
	// 	solver->set_dparameter("rel_tol", rel_tol);
	// 	solver->set_dparameter("abs_tol", abs_tol);
	// 	solver->set_dparameter("max_time_step", max_time_step);
	// 	solver->set_dparameter("min_time_step", min_time_step);
	// 	solver->set_dparameter("T_final", T_final);
	// 	solver->set_iparameter("max_iter", max_iter);

	// 	// Define the macroscopic strain
	// 	if (dimension == 2)
	// 	{
	// 		std::vector<double> strain_load = {strain, 0.0, 0.0};
	// 		model->set_macro_strain(strain_load);
	// 	}
	// 	else
	// 	{
	// 		std::vector<double> strain_load = {strain, 0.0, 0.0, 0.0, 0.0, 0.0};
	// 		model->set_macro_strain(strain_load);
	// 	}

	// 	// Output some problem information to the screen
	// 	if(rank==rank_see)
	// 	{
	// 		cout << "\n\nPROBLEM DETAILS\n---------------------------------------------------\n";
	// 		cout << "Problem Type: " << problem_type_names[problem->get_type()] << endl;
	// 		cout << "Number of Mesh Partitions: " << mesh->n_ranks() << endl;
	// 		cout << "Number of Inclusions: " << mesh->n_inclusions() << endl;
	// 		cout << "Number of Elements: " << problem->get_mesh()->n_global_elem() << endl;
	// 		cout << "Number of Nodes: " << problem->get_mesh()->n_global_nodes()+problem->get_mesh()->n_global_enrich_nodes()  << " (" << problem->get_mesh()->n_global_enrich_nodes() << " enriched)" << endl;
	// 		cout << "Number of dofs: " << problem->get_dofs()->n_global_dofs()  << " (" << problem->get_dofs()->n_global_free_dofs() << " free)" << endl;
	// 	}


	// 	// Finally actually solve the problem
	// 	model->solve();

	// 	// Cleanup
	// 	delete model;
	// }
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


void read_csv(Mesh* mesh, std::string filename, std::vector<double>& maxes, std::vector<double>& mins, double& min_rad, int& dimension,
			  std::vector<double>& sigma_c, std::vector<double>& delta_c, bool& individual)
{
	// Open the file
	ifstream myfile;
	myfile.open( filename );
	if (!myfile.good())
		err_message("Invalid input file name!");

	// Set some initial values
	individual = false; // Will be reset if it is individual properties
	dimension = -1;
	maxes = {-99999999999999999999999999., -99999999999999999999999999., -99999999999999999999999999.};
	mins = {99999999999999999999999999., 99999999999999999999999999., 99999999999999999999999999.};
	min_rad = 99999999999999999999999999.;

	// Prepare a material for these inclusions
	LinearElasticIsotropicMaterial fiber;
	fiber.set_parameter("E", E2);
	fiber.set_parameter("nu", nu2);
	std::string mat_name = "Carbon";
	fiber.set_name(mat_name);

	// Create the inclusion objects
	Circle_Inclusion circle;
	Sphere_Inclusion sphere;
	circle.set_material(&fiber);
	sphere.set_material(&fiber);

	// Read lines of comma-separated values until the end of the file
	id_type count = 0;
	std::stringstream ss;
	for (std::string line; std::getline(myfile, line); )
    {
        // Remove whitespace
        line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());

        // Parse into comma-separated values
        std::vector<std::string> vals = Utilities::splitString(line, ",");

        // If the individual material properties are being specified, set the individual material name here
        if (vals.size()==5 || vals.size()==6)
        {
        	if (!individual && count!=0)	 // We've read in thing's already and they havent had sigma and delta specified...
        		err_message("Input file with multiple input formats!");
        	individual = true;
			ss.str("");
			ss << count;
			std::string mat_name = "Carbon" + ss.str();
			fiber.set_name(mat_name); // Set the fiber material name
		}


        // Depending on how many entries are on this line, we will add circles or spheres with or without uniques cohesive parameters
        double x, y, z, r, sig, del;
        z = 0.0;
        sig = 0.0;
        del = 0.0;
        if (vals.size()==3 || vals.size()==5) // circle
        {
        	if (dimension == 3)
				err_message("Input file with multiple input formats!");
			dimension = 2;
			x = std::stod(vals[0]);
			y = std::stod(vals[1]);
			r = std::stod(vals[2]);
			add_circle(mesh, &circle, x, y, r); // Fiber material name has already been set if this is an individual properties sim

			// If this is an individual properties simulation, 
			if (vals.size()==3)
			{
				if (individual)
        			err_message("Input file with multiple input formats!");
        	}
        	else
        	{
        		sig = std::stod(vals[3]);
				del = std::stod(vals[4]);
        	}
        }
        else if (vals.size()==4 || vals.size()==6) // sphere
        {
        	if (dimension == 2)
				err_message("Input file with multiple input formats!");
			dimension = 3;
			x = std::stod(vals[0]);
			y = std::stod(vals[1]);
			z = std::stod(vals[2]);
			r = std::stod(vals[3]);
			add_sphere(mesh, &sphere, x, y, z, r); // Fiber material name has already been set if this is an individual properties sim

			// If this is an individual properties simulation, 
			if (vals.size()==4)
			{
				if (individual)
        			err_message("Input file with multiple input formats!");
        	}
        	else
        	{
        		sig = std::stod(vals[4]);
				del = std::stod(vals[5]);
        	}
        }
        else
        	err_message("Unknown structure of input csv file!");

        // Set the maximum dimension values
        min_rad = std::min(min_rad, r);
        maxes[0] = std::max(maxes[0], x+r);
        mins[0] = std::min(mins[0], x-r);
        maxes[1] = std::max(maxes[1], y+r);
        mins[1] = std::min(mins[1], y-r);
        if (dimension == 3)
        {
        	maxes[2] = std::max(maxes[2], z+r);
        	mins[2] = std::min(mins[2], z-r);
        }

        // If this is an individual properties simulation, add the sigma and delta to a vector
        if (individual)
        {
        	sigma_c.push_back(sig);
        	delta_c.push_back(del);
        }

        count++;
    }

    // truncate the maxes and mins if this is a 2D sim
    if (dimension == 2)
    {
    	maxes.resize(2);
    	mins.resize(2);
    }

    // Finally, close the file
    myfile.close();
}


void read_stl(Mesh* mesh, std::string filename, std::vector<double>& maxes, std::vector<double>& mins, double& min_rad, int& dimension,
			  std::vector<double>& sigma_c, std::vector<double>& delta_c, bool& individual)
{
	individual = false;

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









