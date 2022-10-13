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
#include "ProblemNonlinearStructural_Shrinkage.h"
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
#include <iomanip>

using namespace std;

// Path to the mpiexec executable
// ./petsc-3.6.2/arch-linux2-c-opt/bin/mpiexec



#define pi 3.141592653589799323

// Define the x-axial loading strain
#define strain 0.03

// Define the initial temperature change
#define temp_change -150

// Parameters to control size and number of elements in the domain
#define n_elem_per_inclusion 10	// Define the number of elements to put across the diameter of the smallest fiber
#define xexpansion 0.005
#define yexpansion 0.02
#define zexpansion 0.02
#define yexpansion_plies 0.15	// Increase the size of the domain (add matrix) by this much on each side of the domain (0.05=5%)
#define bounding_ply_limit 0.135 // the transversely isotropic bounding plies cover top and bottom 17% of the domain

// Set some flags here
#define cohesive_flag true
#define max_principal_stress_flag false
#define continuum_damage_flag false
#define bounding_plies_flag false
#define subscale_flag false
#define thermal_shrinkage_flag false

// Define matrix properties (Taken from Scott's ppt)
#define E1 3760 // MPa
#define nu1 0.39
#define alpha1 0.00009897 // From Chris's data
#define sigma_max 70 // From Sottos

//#define E1 2380 // MPa
//#define nu1 0.43
//#define alpha1 0.00009897 // From Chris's data
//#define sigma_max 70 // From Sottos


//#define E1 2500 // MPa
//#define nu1 0.3
//#define alpha1 0.00009897 // From Chris's data
//#define sigma_max 70 // From Sottos

// Define carbon fiber transverse properties (Taken from Scott's ppt)
#define E2 19.5e3 // MPa
#define nu2 0.45
#define alpha2 7.5e-6 // From Pradere and Sauder (2008)

//#define E2 70e3 // MPa
//#define nu2 0.23
//#define alpha2 7.5e-6 // From Pradere and Sauder (2008)

// Define cohesive properties
#define sigmac 10 // MPa
#define deltac 1.5e-4
#define viscosity 0.5

// Define transversely isotrpoic properties (Taken from Scott's ppt)
// glass fiber cte: long: 5e-6
#define Eti1 36.2e3
#define Eti2 7.11e3
#define G12  2.32e3
#define G23  2.17e3
#define nu12 0.335
#define alphat1 7.169e-6 //calculated from equations at http://www.substech.com/dokuwiki/doku.php?id=estimations_of_composite_materials_properties
#define alphat2 6.10e-5 //????????????????????????

// Define damage properties
#define P1 0.95
#define P2 2
#define mu_visc 10 
#define Yin 0
#define epsilon_t 0.035
#define epsilon_c 0.35
#define sigmay_t 62
#define sigmay_c 79
#define H_ro 2e4
#define n_ro 12

// Solver and damage settings
#define max_damage_order 2
#define damage_order_step 1
#define plane_strain 1
#define rel_tol 1e-8
#define abs_tol 1e-8
#define prob_rel_tol 1e-6
#define prob_abs_tol 1e-6
#define T_final 10
#define max_time_step 2e-1
#define min_time_step 1e-6
#define max_iter 40







void read_csv(Mesh* mesh, std::string filename, std::vector<double>& maxes, std::vector<double>& mins, double& min_rad, int& dimension,
			  std::vector<double>& sigma_c, std::vector<double>& delta_c, bool& individual);

void read_stl(Mesh* mesh, std::string filename, std::vector<double>& maxes, std::vector<double>& mins, double& min_rad, int& dimension,
			  std::vector<double>& sigma_c, std::vector<double>& delta_c, bool& individual);
void runProblem(Problem* problem);


// The main function of this test is to test the mesh adaptivity functions
int main (int argc, char* argv[])
{	
	// Define the mesh
	Mesh mesh(&argc, &argv);

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
	if (thermal_shrinkage_flag)
	{
		problem = new ProblemNonlinearStructural_Shrinkage;
		problem->set_parameter("DELTA_T", temp_change);
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
	PIGFEMPrint("Reading inclusions...\n");
	std::vector<double> maxes(2,1.0);	// Maximum domain limits (x, y, and z)
	std::vector<double> mins(2,0.0);	// Minimum domain limits (x, y, and z)
	double min_rad = 1.0;				// Minimum fiber radius (So I know how many elements across to use)
	int dimension=2;
	bool individual = false;
	bool subscale = problem->getOptions()->getBoolOption("-SUBSCALE");
	std::vector<double> sigma_c, delta_c; // only used if we're reading in indivdual fiber properties

	//in_file=input_address;
	// Check what kind of file I'm reading in
	if (in_file.find(".csv", in_file.length()-4) != std::string::npos ||
		in_file.find(".txt", in_file.length()-4) != std::string::npos)
		read_csv(&mesh, in_file, maxes, mins, min_rad, dimension, sigma_c, delta_c, individual);
	else if (in_file.find(".stl", in_file.length()-4) != std::string::npos)
		read_stl(&mesh, in_file, maxes, mins, min_rad, dimension, sigma_c, delta_c, individual);
	else if (in_file == "none")
	{
		PIGFEMPrint("Simulation with no incusions\n");
	}
	else
		err_message("Unkonwn file type!");



	// Generate a mesh
	id_type x_elem(0), y_elem(0), z_elem(0);
	{
		PIGFEMPrint("Generating the mesh...");

		if (problem->getOptions()->getOption("-MESH") == "") // No mesh file provided. Generate a mesh
		{
			double x_lims = maxes[0] - mins[0];
			mins[0] = mins[0] - x_lims*xexpansion;
			maxes[0] = maxes[0] + x_lims*xexpansion;
			x_elem = std::ceil((maxes[0]-mins[0]) / (min_rad*2)) * n_elem_per_inclusion;
			double y_lims = maxes[1] - mins[1];
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
			if (dimension == 2)
				mesh.generate_mesh(TRI3, mins[0], maxes[0], x_elem, mins[1], maxes[1], y_elem); // Actually generate the mesh

			else if (dimension ==  3)
			{
				double z_lims = maxes[2] - mins[2];
				mins[2] = mins[2] - z_lims*zexpansion;
				maxes[2] = maxes[2] + z_lims*zexpansion;
				z_elem = std::ceil((maxes[2]-mins[2]) / (min_rad*2)) * n_elem_per_inclusion;
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
			std::string abaqus_filename = problem->getOptions()->getOption("-MESH");;
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
			if (thermal_shrinkage_flag)
				matrix.set_parameter("THERMAL_EXPANSION", alpha1);
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
			matrix.set_parameter("epsilon_t", epsilon_t);
			matrix.set_parameter("epsilon_c", epsilon_c);
			matrix.set_parameter("sigmay_c", sigmay_c);
			matrix.set_parameter("sigmay_t", sigmay_t);
			matrix.set_parameter("H_ro", H_ro);
			matrix.set_parameter("n_ro", n_ro);
			if (thermal_shrinkage_flag)
				matrix.set_parameter("THERMAL_EXPANSION", alpha1);
			matrix.set_name("Epoxy");
			mesh.add_material(&matrix);
			mesh.set_material("Epoxy");
		}
		else // Linear elastic
		{
			LinearElasticIsotropicMaterial matrix;
			matrix.set_parameter("E", E1);
			matrix.set_parameter("nu", nu1);
			if (thermal_shrinkage_flag)
				matrix.set_parameter("THERMAL_EXPANSION", alpha1);
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
			if (thermal_shrinkage_flag)
			{
				bounds.set_parameter("THERMAL_EXPANSION1", alphat1);
				bounds.set_parameter("THERMAL_EXPANSION2", alphat2);
			}
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
			coh_mat.set_parameter("viscosity", viscosity);

			// XNCohesiveMaterial coh_mat;
			// coh_mat.set_parameter("SIGMA_CN", sigmac);
			// coh_mat.set_parameter("DELTA_C", deltac);
			// // coh_mat.set_parameter("delta_n", deltac);
			// // coh_mat.set_parameter("delta_t", deltac);
			// coh_mat.set_parameter("q", 1.0);	// Ratio of normal to tangential strengths

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
					coh_mat.set_parameter("SIGMA_C", sigma_c[m]);
					coh_mat.set_parameter("DELTA_C", delta_c[m]);
					coh_mat.set_parameter("viscosity", viscosity);
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



	// Assign BCs (FIXME: dOESN'T MAKE A GOOD SYSTEM IN 3D)
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
			BoundaryObject* boundary = problem->get_boundary();
			boundary->set_dirichlet_bcs_from_nodeset("right", 0, app_disp);
			boundary->set_dirichlet_bcs_from_nodeset("left", 0, 0.0);
			// boundary->set_dirichlet_bcs_from_nodeset("bottom", 1, 0.0);
			for(unsigned int d=0; d<mesh.dim(); ++d)
				boundary->set_dirichlet_bc(0, d, 0);

			// traction free on the left and right
			// id_type top_left_node = (x_elem+1) * y_elem;
			// BoundaryObject* boundary = mesh.get_boundary();
			// if (mesh.dim()>=2)
			// 	boundary->set_dirichlet_bc(top_left_node, 0, 0.0);
			// if (mesh.dim()==3)
			// 	boundary->set_dirichlet_bc(top_left_node, 2, 0.0);
			// for(unsigned int d=0; d<mesh.dim(); ++d)
			// 	boundary->set_dirichlet_bc(0, d, 0);
		}
	}



	// Add Sensitivity info
	{
		if (problem->getOptions()->getBoolOption("-SENSITIVITY"))
		{
			SensitivityRightLoad function;
			problem->addSensitivityFunction(&function);
			if (individual)
			{
				std::stringstream ss;
				for (id_type m=0; m<sigma_c.size(); ++m)
				{
					// Set the individual cohesive name
					ss.str("");
					ss << m;
					std::string coh_mat_name = "Epoxy-Carbon" + ss.str();
					problem->addMaterialSensitivityParameter(coh_mat_name, "sigma_c");
					// problem->addMaterialSensitivityParameter(coh_mat_name, "delta");
				}
			}
			else
			{
				problem->addMaterialSensitivityParameter("Epoxy-Carbon", "sigma");
				// problem->addMaterialSensitivityParameter("Epoxy-Carbon", "delta");
			}

			// // Shape sensitivity parameters
			// for (id_type i=0; i<mesh.n_inclusions(); ++i)
			// {
			// 	problem->addShapeSensitivityParameter(i, "X");
			// 	problem->addShapeSensitivityParameter(i, "Y");
			// 	problem->addShapeSensitivityParameter(i, "R");
			// }
		}
	}

	// Set up the problem
	PIGFEMPrint("\n\nSetting up the problem...");
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
	PIGFEMPrint("\n\nPROBLEM DETAILS\n---------------------------------------------------\n");
	PIGFEMPrint("Problem Type: " << problem_type_names[problem->get_type()] << std::endl);
	PIGFEMPrint("Number of Mesh Partitions: " << mesh.n_ranks() << std::endl);
	PIGFEMPrint("Number of Inclusions: " << mesh.n_inclusions() << std::endl);
	PIGFEMPrint("Number of Elements: " << problem->get_mesh()->n_global_elem() << std::endl);
	PIGFEMPrint("Number of Nodes: " << problem->get_mesh()->n_global_nodes()+problem->get_mesh()->n_global_enrich_nodes()  << " (" << problem->get_mesh()->n_global_enrich_nodes() << " enriched)" << std::endl);

	// Print some of the domain information
	PIGFEMPrint("The domain for this problem was :" << std::endl);
	for (id_type i=0; i<maxes.size(); ++i)
		PIGFEMPrint("\t" << i << ": (" << std::setprecision(16) << mins[i] << " - " << std::setprecision(16) << maxes[i] << ")\n" << std::setprecision(6));

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

}

























void add_circle(Mesh* mesh, Inclusion* circle, std::vector<double> vals)
{
	std::vector<double> center = {vals[0], vals[1]};
	circle->set_vec_parameter("center", center);
	circle->set_parameter("r", vals[2]);
	mesh->add_inclusion(circle);
}
void add_ellipse(Mesh* mesh, Inclusion* ellipse, std::vector<double> vals)
{
	std::vector<double> center = {vals[0], vals[1]};
	ellipse->set_vec_parameter("center", center);
	ellipse->set_parameter("a", vals[2]);
	ellipse->set_parameter("b", vals[3]);
	ellipse->set_parameter("alpha", vals[4]);
	mesh->add_inclusion(ellipse);
}
void add_sphere(Mesh* mesh, Inclusion* sphere, std::vector<double> vals)
{
	std::vector<double> center = {vals[0], vals[1], vals[2]};
	sphere->set_vec_parameter("center", center);
	sphere->set_parameter("r", vals[3]);
	mesh->add_inclusion(sphere);
}
void add_ellipsoid(Mesh* mesh, Inclusion* ellipsoid, std::vector<double> vals)
{
	std::vector<double> center = {vals[0], vals[1], vals[2]};
	ellipsoid->set_vec_parameter("center", center);
	ellipsoid->set_parameter("a", vals[3]);
	ellipsoid->set_parameter("b", vals[4]);
	ellipsoid->set_parameter("c", vals[5]);
	ellipsoid->set_parameter("alpha", vals[6]);
	ellipsoid->set_parameter("beta", vals[7]);
	ellipsoid->set_parameter("gamma", vals[8]);
	mesh->add_inclusion(ellipsoid);
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
	LinearElasticIsotropicMaterial inclusion;
	inclusion.set_parameter("E", E2);
	inclusion.set_parameter("nu", nu2);
	if (thermal_shrinkage_flag)
		inclusion.set_parameter("THERMAL_EXPANSION", alpha2);
	std::string mat_name = "Carbon";
	inclusion.set_name(mat_name);

	// Create the inclusion objects
	Circle_Inclusion circle;
	Ellipse_Inclusion ellipse;
	Sphere_Inclusion sphere;
	Ellipsoid_Inclusion ellipsoid;
	circle.set_material(&inclusion);
	ellipse.set_material(&inclusion);
	sphere.set_material(&inclusion);
	ellipsoid.set_material(&inclusion);

	// Read lines of comma-separated values until the end of the file
	id_type count = 0;
	std::stringstream ss;
	std::vector<std::pair<id_type, id_type> > inc_type_to_n_entries(5);
	inc_type_to_n_entries[0] = std::make_pair(5,7); // Plane
	inc_type_to_n_entries[1] = std::make_pair(3,5); // Circle
	inc_type_to_n_entries[2] = std::make_pair(5,7); // Ellipse
	inc_type_to_n_entries[3] = std::make_pair(4,6); // Sphere
	inc_type_to_n_entries[4] = std::make_pair(9,11); // Ellipsoid
	for (std::string line; std::getline(myfile, line); )
    {
        // Remove whitespace
        line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());

        // Parse into comma-separated values
        std::vector<std::string> strs = Utilities::splitString(line, ",");
        if (strs.size()==0)
        	break;
        id_type inc_type = std::stoi(strs[0]); // Type of inclusion (0=plane, 1=circle, 2=ellipse, 3=sphere, 4=ellipsoid)
        std::vector<double> vals(strs.size()-1);
        for (id_type i=1; i<strs.size(); ++i)
        	vals[i-1] = std::stod(strs[i]);


        // If the individual material properties are being specified, set the individual material name here
        if (vals.size()==inc_type_to_n_entries[inc_type].second)
        {
        	if (!individual && count!=0)	 // We've read in thing's already and they havent had sigma and delta specified...
        		err_message("Input file with multiple input formats!");
        	individual = true;
			ss.str("");
			ss << count;
			std::string mat_name = "Carbon" + ss.str();
			inclusion.set_name(mat_name); // Set the inclusion material name
		}

		// Check for a valid structure
		if ((vals.size()!=inc_type_to_n_entries[inc_type].first && !individual) ||
			(vals.size()!=inc_type_to_n_entries[inc_type].second && individual))
		{
			PIGFEMPrint("inc_type=" << inc_type << " n_entries=(" << inc_type_to_n_entries[inc_type].first << "," << inc_type_to_n_entries[inc_type].second << std::endl);
			err_message("Unknown structure of input csv file!");
		}

        // Depending on how many entries are on this line, we will add circles or spheres with or without uniques cohesive parameters
        double x, y, z, sig, del, mindim, maxdim;
        z = 0.0; sig = 0.0; del = 0.0;
        if (inc_type==1) // circle
        {
        	if (dimension == 3)
				err_message("Input file with multiple input formats!");
			dimension = 2;
			x = vals[0]; y = vals[1];
			add_circle(mesh, &circle, vals);
			mindim = vals[2];
			maxdim = vals[2];

			// If this is an individual properties simulation,
			if (vals.size()==3)
			{
				if (individual)
        			err_message("Input file with multiple input formats!");
        	}
        	else
        	{
        		sig = vals[3];
				del = vals[4];
        	}
        }
        else if (inc_type==2) // ellipse
        {
        	if (dimension == 3)
				err_message("Input file with multiple input formats!");
			dimension = 2;
			x = vals[0]; y = vals[1];
			add_ellipse(mesh, &ellipse, vals);
			mindim = std::min(vals[2], vals[3]);
			maxdim = std::max(vals[2], vals[3]);

			// If this is an individual properties simulation,
			if (vals.size()==5)
			{
				if (individual)
        			err_message("Input file with multiple input formats!");
        	}
        	else
        	{
        		sig = vals[5];
				del = vals[6];
        	}
        }
        else if (inc_type==3) // sphere
        {
        	if (dimension == 2)
				err_message("Input file with multiple input formats!");
			dimension = 3;
			x = vals[0]; y = vals[1]; z = vals[2];
			add_sphere(mesh, &sphere, vals);
			mindim = vals[3];
			maxdim = vals[3];

			// If this is an individual properties simulation,
			if (vals.size()==4)
			{
				if (individual)
        			err_message("Input file with multiple input formats!");
        	}
        	else
        	{
        		sig = vals[4];
				del = vals[5];
        	}
        }
        else if (inc_type==4) // ellipsoid
        {
        	if (dimension == 2)
				err_message("Input file with multiple input formats!");
			dimension = 3;
			x = vals[0]; y = vals[1]; z = vals[2];
			add_ellipsoid(mesh, &ellipsoid, vals);
			mindim = std::min(std::min(vals[3], vals[4]), vals[5]);
			maxdim = std::max(std::max(vals[3], vals[4]), vals[5]);

			// If this is an individual properties simulation,
			if (vals.size()==9)
			{
				if (individual)
        			err_message("Input file with multiple input formats!");
        	}
        	else
        	{
        		sig = vals[9];
				del = vals[10];
        	}
        }
        else
        	err_message("Unknown structure of input csv file!");

        // Set the maximum dimension values
        min_rad = std::min(min_rad, mindim);
        maxes[0] = std::max(maxes[0], x+maxdim);
        mins[0] = std::min(mins[0], x-maxdim);
        maxes[1] = std::max(maxes[1], y+maxdim);
        mins[1] = std::min(mins[1], y-maxdim);
        if (dimension == 3)
        {
        	maxes[2] = std::max(maxes[2], z+maxdim);
        	mins[2] = std::min(mins[2], z-maxdim);
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


