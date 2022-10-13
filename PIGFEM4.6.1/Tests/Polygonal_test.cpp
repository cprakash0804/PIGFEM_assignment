#include <iostream>
#include <sys/time.h>
#include <vector>
#include "elem.h"
#include "Polygon.h"
#include "mpi.h"

using namespace std;


#define pi 3.141592653589799323
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
	Elem* pent = new Polygon(5);

	double coords[5][3] = {{1.0, 0.0},
        				   {1.0, 1.0},
        				   {0.3, 1.6},
        				   {-0.25, 0.5},
        				   {0.0, -1.25}};	
	std::vector<Node*> nodes(5);

	for (id_type n=0; n<5; ++n)
		nodes[n] = new Node(coords[n][0], coords[n][1], coords[n][2], n);

	pent->set_nodes(nodes);

	id_type nqp = pent->n_q_points();
	for (id_type qp=0; qp<nqp; ++qp)
	{
		// Get the quadrature point
		std::vector<double> rcoords;
		double W;
		pent->q_point(rcoords, W, qp);

		// Get the shape functions at this quad point
		std::vector<double> N;
		std::vector<std::vector<double> > dN;
		double J;
		pent->ShapeFunctions(rcoords, N, dN, J);
	}

	for (id_type n=0; n<5; ++n)
		delete nodes[n];

	delete pent;
}
