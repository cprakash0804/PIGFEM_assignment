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
#include "IGFEM_mesh.h"
#include "inclusion.h"
#include "plane_inclusion.h"
#include "ellipse_inclusion.h"
#include "ellipsoid_inclusion.h"
#include "elem.h"
#include "point1.h"
#include "edge2.h"
#include "tri3.h"
#include "quad4.h"
#include "tet4.h"
#include "hex8.h"
#include "material.h"
#include "material_lt.h"
#include "mpi.h"

using namespace std;

#define pi 3.141592653589799323
#define kappa1 0.03
#define kappa2 0.1
#define a 2.5
#define d_size 4
//#define z_size 5e-2
#define y_size 0.2
#define nelem 3
#define l_temp 0
#define r_temp 100
#define x_plane 2.11
#define y_plane 0.085
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

template <typename T>
void print_set(const std::set<T> set)
{
	if(set.size() >= 1)
	{
		std::cout << "{";
		for(auto it=set.begin(), end=set.end(); it!=end;)
		{
			cout << *it;
			if(++it != end)
				cout << ", ";
		}
		cout << "}";
	}
}

template <typename T>
void print_set(const std::set<std::pair<T, T> > set)
{
	if(set.size() >= 1)
	{
		std::cout << "{";
		for(auto it=set.begin(), end=set.end(); it!=end;)
		{
			cout << "{" << it->first << "," << it->second << "}";
			if(++it != end)
				cout << ", ";
		}
		cout << "}";
	}
}

std::vector<double> Analytical_temp(std::vector<double> gcoords)
{
	std::vector<double> temp;
	temp.push_back(((gcoords[0]+d_size)/(2*d_size))*(r_temp-l_temp) + l_temp);
	return temp;
}

std::vector<double> Analytical_IGFEM_temp(std::vector<double> gcoords)
{
	std::vector<double> temp;
	double s2 = (r_temp-l_temp)/(kappa2*(x_plane+d_size)/kappa1-(x_plane-d_size));
	double s1 = (kappa2/kappa1)*s2;
	if(gcoords[0] <= x_plane)
		temp.push_back(s1*(gcoords[0]+d_size)+l_temp);
	else
		temp.push_back(s2*(gcoords[0]-d_size)+r_temp);

	return temp;
}

std::vector<double> Analytical_IGFEM_temp1(std::vector<double> gcoords)
{
	std::vector<double> temp;
	double s2 = (r_temp-l_temp)/(kappa2*(y_plane)/kappa1-(y_plane-y_size));
	double s1 = (kappa2/kappa1)*s2;
	if(gcoords[1] <= y_plane)
		temp.push_back(s1*gcoords[1]+l_temp);
	else
		temp.push_back(s2*(gcoords[1]-y_size)+r_temp);

	return temp;
}


std::vector<double> Analytical_IGFEM_temp2(std::vector<double> gcoords)
{
	double d_size_diag = d_size*sqrt(2.0);
	double x_plane_diag = x_plane/sqrt(2.0);
	std::vector<double> temp;
	double s2 = (r_temp-l_temp)/(kappa2*(x_plane+d_size_diag)/kappa1-(x_plane-d_size_diag));
	double s1 = (kappa2/kappa1)*s2;
	// Some math to get the normal distance along the diagonal coordinate
	double dot = (gcoords[1]-gcoords[0])/sqrt(2.0);
	std::vector<double> diag_vec = {gcoords[0]+dot/sqrt(2.0), gcoords[1]-dot/sqrt(2.0)};
	double diag_dist = sqrt(pow(diag_vec[0],2) + pow(diag_vec[1],2));
	if(diag_vec[0] < 0)
		diag_dist = -1.0*diag_dist; 

	if(diag_dist <= x_plane_diag)
		temp.push_back(s1*(diag_dist+d_size_diag)+l_temp);
	else
		temp.push_back(s2*(diag_dist-d_size_diag)+r_temp);

	return temp;
}





int main (int argc, char* argv[])
{
/*
	// Create a mesh! start with 2D test
	Mesh mesh;
	mesh.init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	const clock_t g_start = clock();
	//mesh.generate_mesh(TET4, -d_size, d_size, nelem, -d_size, d_size, nelem, -z_size, z_size, nelem);
	mesh.generate_mesh(QUAD4, -d_size, d_size, nelem, -d_size, d_size, nelem);
	double gen_time = double(clock() - g_start)/CLOCKS_PER_SEC;

	// Partition the mesh. This has to be done before adding inclusions because I haven't implemented partitioning inclusions... Really should just call partition in generate_mesh
	const clock_t p_start = clock();
	mesh.partition_serial();
	double part_time = double(clock() - p_start)/CLOCKS_PER_SEC;

	// Call generate_node_elem. NOTE: This only generated the node_elem result for normal nodes. For enrichment nodes this is done in the mesh.add_enrichments function. This has to be done before calling the IGFEM routine though
	const clock_t ne_start = clock();
	mesh.generate_node_elem();
	double ne_time = double(clock() - ne_start)/CLOCKS_PER_SEC;

	// Assign Boundary Conditions post-partitioning
	mesh.set_dirichlet_bcs_from_nodeset("right", 0, r_temp);
	mesh.set_dirichlet_bcs_from_nodeset("left", 0, l_temp);

	// Add a material to the mesh post-partitioning
	LinearThermalMaterial material;
	material.set_parameter("kappa", kappa1);
	material.set_name("Epoxy");
	mesh.add_material(&material);
	mesh.set_material("Epoxy");


	// Now that Everything has been assigned, lets distribute the dofs
	const clock_t dof_start = clock();
	mesh.distribute_dofs_thermal();
	double dof_time = double(clock() - dof_start)/CLOCKS_PER_SEC;

	// Solve the system!
	double assem_time, solve_time;
	const clock_t s_start = clock();
	mesh.solve(assem_time, solve_time);
	double solve_acc_time = (double(clock() - s_start)/CLOCKS_PER_SEC) - assem_time - solve_time;

	// Compute L2-norm of the solution
	char buf[50];
	sprintf(buf, "FEM_thermal_%i.txt", nelem);
	mesh.print_solution( std::string(buf) );
	double L2 = mesh.Compute_L2_Norm(&Analytical_temp);

	// Collect the times on Process 0 for output
	double gen_times[size];
	double part_times[size];
	double ne_times[size];
	double dof_times[size];	
	double assem_times[size];
	double solve_times[size];
	double solve_acc_times[size];
	MPI_Gather(&gen_time, 1, MPI_DOUBLE, gen_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&part_time, 1, MPI_DOUBLE, part_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&ne_time, 1, MPI_DOUBLE, ne_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&dof_time, 1, MPI_DOUBLE, dof_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&assem_time, 1, MPI_DOUBLE, assem_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&solve_time, 1, MPI_DOUBLE, solve_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&solve_acc_time, 1, MPI_DOUBLE, solve_acc_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if(rank==0)
	{
		ofstream myfile;
		myfile.open("output3.txt", std::ofstream::app);
		myfile << "Solving " << mesh.n_dofs() << " DOFs with " << size << " processors." << endl;
		for(int i=0; i<size; ++i)
		{
			myfile << "Rank " << i << endl;
			myfile << "\tGENERATION TIME: " << gen_times[i] << " seconds\n";
			myfile << "\tPARTITIONING TIME: " << part_times[i] << " seconds\n";
			myfile << "\tNODE-ELEM GENERATION TIME: " << ne_times[i] << " seconds\n";
			myfile << "\tDOF DISTRIBUTION TIME: " << dof_times[i] << " seconds\n";
			myfile << "\tASSEMBLY TIME: " << assem_times[i] << " seconds\n";
			myfile << "\tSOLVE TIME: " << solve_times[i] << " seconds\n";
			myfile << "\tSOLVE ACCESORY TIME: " << solve_acc_times[i] << " seconds\n";
		}
		myfile << endl;
		myfile.close();
	}

	if(rank==0)
	{
		cout << "Number of Elements: " << mesh.n_elem() << endl;
		cout << "Number of Nodes: " << mesh.n_nodes() << endl;
		cout << "Number of dofs: " << mesh.n_global_dofs() << endl;
		cout << "L2 Norm = " << L2 << "\n\n";
	}

	mesh.Finalize();
*/
//===============================================================================================================================================================
// NORMAL MESH WORK ABOVE HERE





















// IGEFM MESH WORK BELOW HERE
//===============================================================================================================================================================

	// Create an IGFEM mesh!
	IGFEM_Mesh mesh;
	mesh.init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	const clock_t g_start = clock();
	//mesh.read_mesh("Tet_Heat2.inp");
	mesh.generate_mesh(QUAD4, -d_size, d_size, nelem, -d_size, d_size, nelem);
	double gen_time = double(clock() - g_start)/CLOCKS_PER_SEC;

	// Partition the mesh. This has to be done before adding inclusions because I haven't implemented partitioning inclusions... Really should just call partition in generate_mesh
	const clock_t p_start = clock();
	mesh.partition_serial();
	double part_time = double(clock() - p_start)/CLOCKS_PER_SEC;

	// Call generate_node_elem. NOTE: This only generated the node_elem result for normal nodes. For enrichment nodes this is done in the mesh.add_enrichments function. This has to be done before calling the IGFEM routine though
	const clock_t ne_start = clock();
	mesh.generate_node_elem();
	double ne_time = double(clock() - ne_start)/CLOCKS_PER_SEC;

	// Add a material to the mesh post-partitioning
	// NOTE: We don't really have to worry about adding maerials to each element like we did for a non-IGFEM mesh.
	//		 The element detection takes care of this by detectin that an element is part of the "matrix" that we define
	LinearThermalMaterial material;
	material.set_parameter("kappa", kappa1);
	material.set_name("Epoxy");
	mesh.add_material(&material);
	LinearThermalMaterial material2;
	material2.set_parameter("kappa", kappa2);
	material2.set_name("Glass");
	mesh.add_material(&material2); // Test adding inclusion with and without this line
	mesh.set_matrix_material("Epoxy");     // This is the first thing that is really inherent to IGFEM

	// Create the inclusion
	Plane_Inclusion plane;
	std::vector<double> vec = {x_plane, 0.0};
	plane.set_vec_parameter("Point", vec);
	vec = {1.0, 1.0};
	plane.set_vec_parameter("norm", vec);
	plane.set_material(&material2);
	mesh.add_inclusion(&plane);

	// Assign Boundary Conditions post-partitioning
	mesh.set_dirichlet_bcs_from_nodeset("right", 0, r_temp);
	mesh.set_dirichlet_bcs_from_nodeset("left", 0, l_temp);
	
	

	// This is where the real fun begins. I want to have one function that does all of these but I also want to be able to time all of them individually so I won't do that right now
	// Call nodal detection routine
	const clock_t an_start = clock();
	mesh.analyze_nodes();
	double an_time = double(clock() - an_start)/CLOCKS_PER_SEC;

	// Call the elemental detection routine
	const clock_t ae_start = clock();
	mesh.analyze_elements();
	double ae_time = double(clock() - ae_start)/CLOCKS_PER_SEC;

	// Call the add_enrichments function. This is the heavyweight. Does intersection finding, enrichment node generation, and integration elment generation
	const clock_t aen_start = clock();
	mesh.add_enrichments();
	double aen_time = double(clock() - aen_start)/CLOCKS_PER_SEC;


	// Now that Everything has been assigned, lets distribute the dofs
	const clock_t dof_start = clock();
	mesh.distribute_dofs_thermal();
	double dof_time = double(clock() - dof_start)/CLOCKS_PER_SEC;

	// Solve the system!
	double assem_time, solve_time;
	const clock_t s_start = clock();
	mesh.solve(assem_time, solve_time);
	double solve_acc_time = (double(clock() - s_start)/CLOCKS_PER_SEC) - assem_time - solve_time;


	char buf[50];
	sprintf(buf, "IGFEM_thermal_%i.txt", nelem);
	mesh.print_solution( std::string(buf) );

	// Compute L2-norm of the solution
	double L2 = mesh.Compute_L2_Norm(&Analytical_IGFEM_temp);


	// Collect the times on Process 0 for output
	double gen_times[size];
	double part_times[size];
	double ne_times[size];
	double an_times[size];
	double ae_times[size];
	double aen_times[size];
	double dof_times[size];	
	double assem_times[size];
	double solve_times[size];
	double solve_acc_times[size];
	MPI_Gather(&gen_time, 1, MPI_DOUBLE, gen_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&part_time, 1, MPI_DOUBLE, part_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&ne_time, 1, MPI_DOUBLE, ne_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&an_time, 1, MPI_DOUBLE, an_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&ae_time, 1, MPI_DOUBLE, ae_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&aen_time, 1, MPI_DOUBLE, aen_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&dof_time, 1, MPI_DOUBLE, dof_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&assem_time, 1, MPI_DOUBLE, assem_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&solve_time, 1, MPI_DOUBLE, solve_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&solve_acc_time, 1, MPI_DOUBLE, solve_acc_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if(rank==0)
	{
		ofstream myfile;
		myfile.open("output3.txt", std::ofstream::app);
		myfile << "Solving " << mesh.n_dofs() << " DOFs with " << size << " processors." << endl;
		for(int i=0; i<size; ++i)
		{
			myfile << "Rank " << i << endl;
			myfile << "\tGENERATION TIME: " << gen_times[i] << " seconds\n";
			myfile << "\tPARTITIONING TIME: " << part_times[i] << " seconds\n";
			myfile << "\tNODE-ELEM GENERATION TIME: " << ne_times[i] << " seconds\n";
			myfile << "\tNODE ANALYSIS TIME: " << an_times[i] << " seconds\n";
			myfile << "\tELEMENT ANALYSIS: " << ae_times[i] << " seconds\n";
			myfile << "\tENRICHMENT ADDITION TIME: " << aen_times[i] << " seconds\n";
			myfile << "\tDOF DISTRIBUTION TIME: " << dof_times[i] << " seconds\n";
			myfile << "\tASSEMBLY TIME: " << assem_times[i] << " seconds\n";
			myfile << "\tSOLVE TIME: " << solve_times[i] << " seconds\n";
			myfile << "\tSOLVE ACCESORY TIME: " << solve_acc_times[i] << " seconds\n";
		}
		myfile << endl;
		myfile.close();
	}

	if(rank==0)
	{
		cout << "Number of Elements: " << mesh.n_elem() << endl;
		cout << "Number of Nodes: " << mesh.n_nodes() << endl;
		cout << "Number of Enrichment Nodes: " << mesh.n_global_enrich_nodes() << endl;
		cout << "Number of dofs: " << mesh.n_global_dofs() << endl;
		cout << "L2 Norm = " << L2 << "\n\n";
	}

	mesh.Finalize();
}
