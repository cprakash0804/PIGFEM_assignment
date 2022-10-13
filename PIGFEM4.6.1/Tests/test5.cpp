#include <iostream>
#include <fstream>
#include <cstddef>
#include <sys/time.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <set>
#include <map>
#include "mesh.h"
#include "material_lei.h"   // I really shouldn't have to include this. Figure out why I do....
#include "mpi.h"

using namespace std;

#define E 69e9
#define nu 0.334
#define right_disp 0.1
#define x_dim 4.0
#define y_dim 4.0


template <typename T>
void print_vec(std::vector<T> vec)
{
	if(vec.size() >= 1)
	{
		std::cout << "{";
		for(unsigned int i=0; i<(vec.size()-1); ++i)
			cout << vec[i] << ", ";
		cout << vec[vec.size()-1] << "}";
	}
}


// Assumes plane strain
std::vector<double> TwoD_analytical_disp(std::vector<double> gcoords)
{
	std::vector<double> ret(2);
	double eps = right_disp/x_dim;
	ret[0] = eps*gcoords[0];
	ret[1] = (nu*eps/(nu-1.0))*gcoords[1];
	return ret;
}

int main (int argc, char* argv[])
{

	Mesh mesh;
	mesh.init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	const clock_t g_start = clock();
	mesh.generate_mesh(TRI3, 0.0, x_dim, 4, 0.0, y_dim, 4);
	double gen_time = double(clock() - g_start)/CLOCKS_PER_SEC;

	const clock_t p_start = clock();
	mesh.partition_serial();
	double part_time = double(clock() - p_start)/CLOCKS_PER_SEC;

	const clock_t ne_start = clock();
	mesh.generate_node_elem();
	double ne_time = double(clock() - ne_start)/CLOCKS_PER_SEC;

	// Add some elemsets post-partitioning
	std::set<unsigned int> e_set = {0, 1, 2, 3};
	mesh.add_elemset("bottom", e_set);
	e_set.clear();  e_set = {12, 13, 14, 15};
	mesh.add_elemset("top", e_set);

	// Add a material to the mesh post-partitioning
	LinearElasticIsotropicMaterial material;
	material.set_parameter("E", E);
	material.set_parameter("nu", nu);
	material.set_name("Aluminum");
	mesh.add_material(&material);

	// Set elemental materials
	mesh.set_material("Aluminum");

	// Assign Boundary Conditions post-partitioning
	mesh.set_dirichlet_bcs_from_nodeset("left", 0, 0); // Hold left
	mesh.set_dirichlet_bcs_from_nodeset("right", 0, right_disp); // Pull right
	mesh.set_dirichlet_bc(0, 1, 0); // stop y movement

	// Now that Everything has been assigned, lets distribute the dofs
	const clock_t dof_start = clock();
	mesh.distribute_dofs();
	double dof_time = double(clock() - dof_start)/CLOCKS_PER_SEC;

	// Solve the system!
	double assem_time, solve_time;
	const clock_t s_start = clock();
	mesh.solve(assem_time, solve_time);
	double solve_acc_time = (double(clock() - s_start)/CLOCKS_PER_SEC) - assem_time - solve_time;

	// Compute L2-norm of the solution
	double L2 = mesh.Compute_L2_Norm(&TwoD_analytical_disp);

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

	int rank_see = 0;
	if(rank==rank_see)
	{
		cout << "PARTITION " << rank_see << "\n--------------------------------------------------\n";
		
		cout << "Number of Elements: " << mesh.n_elem() << endl;
		cout << "Number of Nodes: " << mesh.n_nodes() << endl;

	/*
		// Output the elements of the mesh including the nodes that make them up
		cout << "Mesh:" << endl;
		for(unsigned int i=0; i<mesh.n_local_elem(); ++i)
			cout << *(mesh.get_elem_local(i));
	*/

	/*
		// Output nodal ownership and whether or not they are on a partition interface
		cout << "\n\nNodal ownership and interface information:\n";
		for(unsigned int i=0; i<mesh.n_nodes(); ++i)
		{
			unsigned int id = mesh.get_node_local(i)->get_id();
			cout << "Node " << id << "(" << mesh.get_node_owner_global(id) << "): ";
			if(!mesh.node_on_part_interface(id))
				cout << "Not on a partition interface!" << endl;
			else
			{
				cout << "Owned by " << mesh.n_parts_with_node(id) << " partitions. ";
				for(unsigned int j=0; j<mesh.n_parts_with_node(id); ++j)
					cout << mesh.get_part_int_from_node(id, j) << " ";
				cout << endl;
			}
		}
	*/

	/*
		// Output nodesets and element sets
		cout << "\n\nNodesets:\n";
		for(unsigned int n=0; n<mesh.n_nodesets(); ++n)
			cout << mesh.get_nodeset(n);
	*/

	/*
		cout << "\n\nElemsets:\n";
		for(unsigned int n=0; n<mesh.n_elemsets(); ++n)
			cout << mesh.get_elemset(n);
	*/

	/*
		// Output Material information
		cout << "\n\nMaterials:\n";
		for(unsigned int m=0; m<mesh.n_materials(); ++m)
		{
			Material* mat = mesh.get_material(m);
			cout << mat->get_name() << ": E=" << mat->get_parameter("E") << " K=" << mat->get_parameter("K") << endl;
		}
		cout << "Element Materials:\n";
		for(unsigned int e=0; e<mesh.n_local_elem(); ++e)
		{
			Elem* el = mesh.get_elem_local(e);
			unsigned int id = el->get_id();
			Material* mat = mesh.get_element_material_local(e);
			if(mat!=NULL)
				cout << "Element(" << id << "): Material(" << mat->get_name() << ")" << endl;
			else
				cout << "Element(" << id << "): Material(INVALID MATERIAL)" << endl;
		}
	*/

	/*
		// Output Boundary Condition information
		cout << "\n\nBoundary Conditions:\n";
		for(unsigned int n=0; n<mesh.n_local_nodes(); ++n)
		{
			unsigned int id = mesh.get_node_local(n)->get_id();
			for(unsigned int dof=0; dof<mesh.dim(); ++dof)
			{
				bool exists;
				double val = mesh.get_dirichlet_bc(id, dof, exists);
				if(exists)
					cout << "Node(" << id << ") dof(" << dof << ") val(" << val << ")" << endl;
			}
		}
	*/

	/*
		// Output Nodal degrees of freedom
		cout << "\n\nNodal Degrees of Freedom\n";
		cout << "Number of free degrees of freedom: " << mesh.n_global_free_dofs() << endl;
		cout << "Number of constrained degrees of freedom: " << mesh.n_global_const_dofs() << endl;
		cout << "Number of local free degrees of freedom: " << mesh.n_local_free_dofs() << endl;
		cout << "Number of local constrained degrees of freedom: " << mesh.n_local_const_dofs() << endl;
		for(unsigned int n=0; n<mesh.n_local_nodes(); ++n)
		{
			unsigned int id = mesh.get_node_local(n)->get_id();
			std::vector<unsigned int> dofs = mesh.get_nodal_global_dofs_local(n);
			cout << "Node(" << id << "): {";
			for(unsigned int d=0; d<dofs.size(); ++d)
			{
				cout << dofs[d];
				if(d!=(dofs.size()-1))
					cout << ", ";
			}
			cout << "}\n";
		}
	*/

	/*
		// Output the node_elem table (take 2)
		cout << "\n\nNode-Elem Connectivity Table:\n";
		for(unsigned int n=0; n<mesh.n_local_nodes(); ++n)
		{
			unsigned int id = mesh.get_node_local(n)->get_id();
			std::vector<unsigned int> ne = mesh.get_node_elem(id);
			cout << "Node(" << id << "): {";
			for(unsigned int i=0; i<ne.size(); ++i)
			{
				cout << ne[i];
				if(i!=(ne.size()-1))
					cout << ", ";
			}
			cout << "}\n";
		}
	*/
	
	/*
		// Output the _node_elem table
		cout << "Local Node_Elem table:" << endl;
		for(unsigned int n=0; n<mesh.n_local_nodes(); ++n)
		{
			unsigned int id = mesh.get_node_local(n)->get_id();
			std::vector<unsigned int> v = mesh.get_node_elem_local(n);
			if(v.size()!=0)
			cout << "Node(" << id << "): ";
			print_vec(v);
			cout << endl;
		}
		cout << "Remote Node_Elem table:" << endl;
		cout << "Table size: " << mesh.get_remote_node_elem_size() << endl;
		for(unsigned int i=0; i<mesh.n_nodes(); ++i)
		{
			unsigned int id = mesh.get_node_local(i)->get_id();
			for(int r=0; r<size; ++r)
			{
				std::vector<unsigned int> v = mesh.get_remote_node_elem(id, r);
				if(v.size()!=0)
				{
					cout << "Node(" << id << ") Rank(" << r << "): ";
					print_vec(v);
					cout << endl;
				}
			}
		}
	*/

	/*
		// Output the dofs for each node
		cout << "Nodal degrees of freedom:" << endl;
		for(unsigned int n=0; n<mesh.n_nodes(); ++n)
		{
			unsigned int id = mesh.get_node_local(n)->get_id();
			cout << "Node(" << id << "): ";
			print_vec(mesh.get_nodal_global_dofs_local(n));
			cout << endl;
		}
	*/

	/*
		// Output the solution
		cout << "\nSolution:\n\n";
		for(unsigned int n=0; n<mesh.n_local_nodes(); ++n)
		{
			unsigned int g_node = mesh.get_node_local(n)->get_id();
			for(unsigned int d=0; d<mesh.dim(); ++d)
			{
				std::cout << "Node(" << g_node << ") dof(" << d << "): " << mesh.get_solution_local(n, d) << std::endl;
			}
		}

		// Output the reactioc force
		cout << "\nReaction Force:\n\n";
		for(unsigned int n=0; n<mesh.n_local_nodes(); ++n)
		{
			unsigned int g_node = mesh.get_node_local(n)->get_id();
			for(unsigned int d=0; d<mesh.dim(); ++d)
			{
				std::cout << "Node(" << g_node << ") dof(" << d << "): " << mesh.get_force_local(n, d) << std::endl;
			}
		}
	*/
		cout << "\nL2 Norm = " << L2 << "\n\n";
	
	}
	mesh.print_solution("test5_out.txt");
	mesh.Finalize();
	
}
