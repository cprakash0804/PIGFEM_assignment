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
#define kappa1 1.15
#define kappa2 0.001
#define a 2.5
#define d_size 4
//#define z_size 5e-2
#define z_size 4
#define nelem 75
#define l_temp 0
#define r_temp 100
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



int main (int argc, char* argv[])
{

	// Create a mesh! start with 2D test
	Mesh mesh;
	mesh.init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	const clock_t r_start = clock();
	//mesh.read_mesh("Job-1.inp");
	mesh.generate_mesh(TRI3, -1.2, 1.6, 6, -1.2, 1.2, 5);
	double read_time = double(clock() - r_start)/CLOCKS_PER_SEC;

	const clock_t p_start = clock();
	mesh.partition_serial();
	double part_time = double(clock() - p_start)/CLOCKS_PER_SEC;

	const clock_t ne_start = clock();
	mesh.generate_node_elem();
	double ne_time = double(clock() - ne_start)/CLOCKS_PER_SEC;


	if(rank==1)
	{
		cout << "Number of Elements: " << mesh.n_elem() << endl;
		cout << "Number of Nodes: " << mesh.n_nodes() << endl;
		cout << "Read time: " << read_time << endl;
		cout << "Partition time: " << part_time << endl;
		cout << "Node-Elem time: " << ne_time << endl;


		// Output the elements of the mesh including the nodes that make them up
		cout << "Mesh:" << endl;
		for(unsigned int i=0; i<mesh.n_local_elem(); ++i)
			cout << *(mesh.get_elem_local(i));

		// Output nodesets and element sets
		cout << "\n\nNodesets:\n";
		for(Mesh::nodeset_iterator it=mesh.nodesets_begin(), end=mesh.nodesets_end(); it!=end; ++it)
		{
			cout << it->first << ": ";
			print_set(it->second);
			cout << endl;
		}

		cout << "\n\nElemsets:\n";
		for(Mesh::elemset_iterator it=mesh.elemsets_begin(), end=mesh.elemsets_end(); it!=end; ++it)
		{
			cout << it->first << ": ";
			print_set(it->second);
			cout << endl;
		}

		cout << "\n\nSidesets:\n";
		for(Mesh::sideset_iterator it=mesh.sidesets_begin(), end=mesh.sidesets_end(); it!=end; ++it)
		{
			cout << it->first << ": ";
			print_set(it->second);
			cout << endl;
		}

/*
		// Output nodesets and element sets
		cout << "\n\nNodesets:\n";
		for(unsigned int n=0; n<mesh.n_nodesets(); ++n)
		{
			print_set(mesh.get_nodeset(n));
		}
			


		cout << "\n\nElemsets:\n";
		for(unsigned int n=0; n<mesh.n_elemsets(); ++n)
			print_set(mesh.get_elemset(n));
*/

	}

	mesh.Finalize();
}
