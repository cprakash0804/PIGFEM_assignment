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
#include "mpi.h"

using namespace std;

#define pi 3.141592653589799323
#define mu1 1.13e9
//#define mu2 1.13e9
#define mu2 31.15e9
#define nu1 0.33
//#define nu2 0.33
#define nu2 0.22

#define Force 50e6
#define x_disp 0.01

#define a 2.5
#define d_size 4
#define z_size 1


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

// From Goodier (1933) Plane strain
std::vector<double> TwoD_Analytical_Inclusion_Disp(std::vector<double> gcoords)
{
	double theta = atan2(gcoords[1], gcoords[0]);
	double r = sqrt(pow(gcoords[0], 2) + pow(gcoords[1], 2));
	double E1 = 2.0*mu1*(1.0+nu1);

	double A = a*a* (Force/(4.0*mu1)) * ((1.0-2.0*nu2)*mu1-(1.0-2.0*nu1)*mu2)/((1.0-2.0*nu2)*mu1+mu2);
	double B = a*a*a*a* (Force/(4.0*mu1)) * (mu1-mu2)/(mu1+(3.0-4.0*nu1)*mu2);
	double C = a*a* (Force/(2.0*mu1)) * (mu1-mu2)/(mu1+(3.0-4.0*nu1)*mu2);

	double u_r, u_t;
	if(r < a) // Inside the inclusion
	{
		double H = (1.0/(2.0*nu2*a*a)) * (Force/(4.0*mu2) - (mu1/mu2)*(3.0*B/(a*a*a*a)-C/(a*a)) - B/(a*a*a*a) - (C/(a*a))*(1.0-2.0*nu1) - (Force/(2.0*E1))*(1.0+nu1));
		double G = B/(a*a*a*a) + (C/(a*a))*(1.0-2.0*nu1) + (Force/(2.0*E1))*(1.0+nu1) + (2.0*nu2-3.0)*H*a*a;
		double F = (1.0/a)*(A/a + (-B/(a*a*a) + (2.0*C/a)*(1.0-nu1) - G*a - 2.0*nu2*H*a*a*a)*cos(2.0*theta) + (Force*a/(2.0*E1))*(1.0+nu1)*(1.0-2.0*nu1+cos(2.0*theta)));
		u_r = F*r + (G*r + 2.0*nu2*H*r*r*r)*cos(2.0*theta);
		u_t = -1.0*(G*r + (3.0-2.0*nu2)*H*r*r*r)*sin(2.0*theta);
	}
	else 	// Outside the inclusion
	{
		u_r = A/r + (-1.0*B/(r*r*r) + (2.0*C/r)*(1.0-nu1))*cos(2.0*theta) + (Force*r/(2.0*E1))*(1.0+nu1)*(1.0-2.0*nu1+cos(2.0*theta));
		u_t = -1.0*(B/(r*r*r) + (C/r)*(1.0-2.0*nu1))*sin(2.0*theta) - (Force*r/(2.0*E1))*(1.0+nu1)*sin(2.0*theta);
	}

	// Transform from cylindrical coordinates to cartesian
	std::vector<double> ret(2);
	ret[0] = u_r*cos(theta) - u_t*sin(theta);
	ret[1] = u_r*sin(theta) + u_t*cos(theta);
	return ret;
}

// From Goodier (1933)
std::vector<double> ThreeD_Analytical_Inclusion_Disp(std::vector<double> gcoords)
{
	double r = sqrt(pow(gcoords[0], 2) + pow(gcoords[1], 2) + pow(gcoords[2], 2));
	double r_th = sqrt(pow(gcoords[0], 2) + pow(gcoords[1], 2)); // Distance from z-axis parallel to x-y plane
	double phi = atan2(gcoords[1], gcoords[0]); // Angle of rotation about the z axis
	double theta = atan2(r_th, gcoords[2]);
	double E1 = 2.0*mu1*(1.0+nu1);

	double A = a*a*a*( (-1.0*Force/(8.0*mu1)) * ((mu1-mu2)/((7.0-5.0*nu1)*mu1+(8.0-10.0*nu1)*mu2)) * 
					   (((1.0-2.0*nu2)*(6.0-5.0*nu1)*2.0*mu1+(3.0+19.0*nu2-20.0*nu1*nu2)*mu2)/((1.0-2.0*nu2)*2.0*mu1+(1.0+nu2)*mu2)) +
					   (Force/(4.0*mu1)) * ((((1.0-nu1)*(1.0+nu2)/(1.0+nu1)-nu2)*mu2-(1.0-2.0*nu2)*mu1)/((1.0-2.0*nu2)*2.0*mu1+(1.0+nu2)*mu2)) );
	double B = a*a*a*a*a * (Force/(8.0*mu1)) * ((mu1-mu2)/((7.0-5.0*nu1)*mu1+(8.0-10.0*nu1)*mu2));
	double C = a*a*a * (Force/(8.0*mu1)) * ((5.0*(1.0-2.0*nu1)*(mu1-mu2))/((7.0-5.0*nu1)*mu1+(8.0-10.0*nu1)*mu2));

	double u_r, u_t;
	if(r < a) // Inside the inclusion
	{
		double ur_a = -1.0*A/(a*a) - 3.0*B/(a*a*a*a) + (((5.0-4.0*nu1)/(1.0-2.0*nu1))*(C/(a*a))-9.0*(B/(a*a*a*a)))*cos(2.0*theta) + (Force*a/(2.0*E1))*(1.0-nu1+(1.0+nu1)*cos(2.0*theta));
		double ut_a = -1.0*(2.0*C/(a*a) + 6.0*B/(a*a*a*a))*sin(2.0*theta) - (Force*a/(2.0*E1))*(1.0+nu1)*sin(2.0*theta);
		double sigrt_a = 2.0*mu1*(-1.0*(2.0*(1.0+nu1)/(1.0-2.0*nu1))*(C/(a*a*a))+24.0*B/(a*a*a*a*a))*sin(2.0*theta) - (Force/2.0)*sin(2.0*theta);
		double G = (-1.0/(6.0*nu2*a*a*a)) * (a*sigrt_a/(2.0*mu2)-ut_a)/sin(2.0*theta);
		double F = (1.0/3.0)*((-1.0/(2.0*mu2*sin(2.0*theta)))*sigrt_a - (7.0+2.0*nu2)*G*a*a);
		double H = (1.0/a)*(ur_a - F*a - 2.0*nu2*G*a*a*a - (3.0*F*a+6.0*nu2*G*a*a*a)*cos(2.0*theta));
		u_r = H*r + F*r + 2.0*nu2*G*r*r*r + (3.0*F*r+6.0*nu2*G*r*r*r)*cos(2.0*theta);
		u_t = -1.0*(3.0*F*r+(7.0-4.0*nu2)*G*r*r*r)*sin(2.0*theta);
	}
	else 	// Outside the inclusion
	{
		u_r = -1.0*A/(r*r) - 3.0*B/(r*r*r*r) + (((5.0-4.0*nu1)/(1.0-2.0*nu1))*(C/(r*r))-9.0*B/(r*r*r*r))*cos(2.0*theta) + (Force*r/(2.0*E1))*(1.0-nu1+(1.0+nu1)*cos(2.0*theta));
		u_t = -1.0*(2.0*C/(r*r) + 6.0*B/(r*r*r*r))*sin(2.0*theta) - (Force*r/(2.0*E1))*(1.0+nu1)*sin(2.0*theta);
	}

	// Transform from spherical coordinates to cartesian
	std::vector<double> ret(3);
	ret[0] = cos(phi)*(sin(theta)*u_r + cos(theta)*u_t);
	ret[1] = sin(phi)*(sin(theta)*u_r + cos(theta)*u_t);
	ret[2] = cos(theta)*u_r - sin(theta)*u_t;
	return ret;
}

std::vector<double> Analytical_x_tension(std::vector<double> gcoords)
{
	double eps_x = x_disp/(2.0*d_size);
	double eps_t = -nu1*eps_x;
	std::vector<double> u(gcoords.size());
	u[0] = eps_x*(gcoords[0]+d_size);

	if(gcoords.size()==2) // 2D
		u[1] = (-eps_t/(nu1-1.0))*(gcoords[1]+d_size); // Plane strain
	else if(gcoords.size()==3)
	{
		u[1] = (gcoords[1]+d_size)*eps_t;
		u[2] = (gcoords[2]+z_size)*eps_t;
	}
	else
		err_message("Unknown mesh dimension");

	return u;
}











#define nelem 50

int main (int argc, char* argv[])
{
	Mesh mesh(&argc, &argv);
	int rank = mesh.get_rank();
	int rank_see = 0;

	// Generate a mesh
	if(rank==rank_see)
		cout << "\nGenerating the mesh...";
	mesh.generate_mesh(TRI3, -d_size, d_size, nelem, -d_size, d_size, nelem);

	// Add a material to the mesh
	LinearElasticIsotropicMaterial material;
	material.set_parameter("mu", mu1);
	material.set_parameter("nu", nu1);
	material.set_name("Epoxy");
	material.set_plane_strain(true); // Goodier solution is for plane strain
	mesh.add_material(&material);
	mesh.set_material("Epoxy");

	// Assign BCs
	BoundaryObject* boundary = mesh.get_boundary();
	std::vector<std::string> sets;
	if(mesh.dim()==2)
		sets = {"left", "right", "top", "bottom"};
	else if(mesh.dim()==3)
		sets = {"left", "right", "front", "back", "top", "bottom"};
	else
		err_message("Unknown mesh dimension");

	for(unsigned int s=0; s<sets.size(); ++s)
	{
		std::set<unsigned int> curr_set = boundary->get_nodeset(sets[s]);
		for(auto it=curr_set.begin(), end=curr_set.end(); it!=end; ++it)
		{
			Node* node = mesh.get_node_global(*it);
			std::vector<double> gcoords;
			for(unsigned int d=0; d<mesh.dim(); ++d)
				gcoords.push_back( (*node)(d) );
			std::vector<double> disp;
			if(mesh.dim()==2)
				disp = TwoD_Analytical_Inclusion_Disp(gcoords);
			else if(mesh.dim()==3)
				disp = ThreeD_Analytical_Inclusion_Disp(gcoords);
			else
				err_message("Unknown mesh dimension");
			for(unsigned int d=0; d<mesh.dim(); ++d)
				boundary->set_dirichlet_bc(*it, d, disp[d]);
		}
	}

	// Create the inclusion
	LinearElasticIsotropicMaterial material2;
	material2.set_parameter("mu", mu2);
	material2.set_parameter("nu", nu2);
	material2.set_plane_strain(true);
	material2.set_name("Glass");
	if(mesh.dim()==2)
	{
		Ellipse_Inclusion ellipse;
		std::vector<double> center = {0.0, 0.0};
		ellipse.set_vec_parameter("center", center);
		ellipse.set_parameter("alpha", 0.0);
		ellipse.set_parameter("a", a);
		ellipse.set_parameter("b", a);
		ellipse.set_material(&material2);
		mesh.add_inclusion(&ellipse);
	}
	else if(mesh.dim()==3)
	{
		Ellipsoid_Inclusion ellipsoid;
		std::vector<double> center = {0, 0, 0};
		ellipsoid.set_vec_parameter("center", center);
		ellipsoid.set_parameter("alpha", 0.0);
		ellipsoid.set_parameter("beta", 0.0);
		ellipsoid.set_parameter("gamma", 0.0);
		ellipsoid.set_parameter("a", a);
		ellipsoid.set_parameter("b", a);
		ellipsoid.set_parameter("c", a);
		ellipsoid.set_material(&material2);
		mesh.add_inclusion(&ellipsoid);
	}
	else
		err_message("Unknown mesh dimension");

	// Define the problem
	ProblemLinearElasticity problem;
	problem.attach_mesh(&mesh);

	// Solve the problem!
	if(rank==rank_see)
		cout << "\nSolving the problem...";
	problem.solve_problem();

	// Compute the L2 Norm of the error
	if(rank==rank_see)
		cout << "\nComputing the L2 Norm of the error...";
	double L2;
	if(mesh.dim()==2)
		L2 = problem.Compute_L2_Norm(&TwoD_Analytical_Inclusion_Disp);
	else if(mesh.dim()==3)
		L2 = problem.Compute_L2_Norm(&ThreeD_Analytical_Inclusion_Disp);
	else
		err_message("Unknown mesh dimension");

	// Output result to a file
	char buf[50];
	sprintf(buf, "Output/Goodier_tri_%i.vtk", nelem);
	std::string filename( buf );
	problem.writeToVTK( filename );

	// Output more detailed results
	if(rank==rank_see)
	{
		cout << "\n\nPARTITION " << rank_see << "\n--------------------------------------------------\n";
		cout << "Number of Elements: " << problem.get_mesh()->n_global_elem() << endl;
		cout << "Number of Nodes: " << problem.get_mesh()->n_global_nodes() << endl;
		cout << "Number of Enrichment Nodes: " << problem.get_mesh()->n_global_enrich_nodes() << endl;
		cout << "Number of dofs: " << problem.get_dofs()->n_global_dofs() << endl;
		cout << "L2 Norm of the error: " << L2 << endl << endl;
	}


	problem.Finalize();
}
