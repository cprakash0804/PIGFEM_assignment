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
#include "edge2.h"
#include "tri3.h"
#include "quad4.h"
#include "tet4.h"
#include "hex8.h"
#include "material.h"
#include "material_lei.h"
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
#define a 2.5
#define d_size 4
//#define z_size 5e-2
#define z_size 4
#define nelem 200
#define x_disp 1e-1
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



void display_element_detection(Elem* elem, std::vector<int> detection)
{
	cout << "\n";
	
	// Dislpay the nodal detection
	cout << "Nodal Detection: "; print_vec(detection); cout << endl;

	// Display the info stored in the element
	cout << "\tElement Detector: " << elem->getElemDetector() << endl;
	cout << "\tIntersected Edges: "; print_vec(elem->getIntersectedEdges()); cout << endl;
	cout << "\tInclusion Number: " << elem->getInclusionNumber() << endl;
	std::vector<std::vector<unsigned int> > int_elems = elem->getIntegrationElemStructure();
	cout << "\tIntegration Elements: ";
	if(int_elems.size()!=0)
		print_vec(int_elems[0]);
	else
		cout << "{}";
	cout << endl;
	for(unsigned int e=1; e<int_elems.size(); ++e)
	{
		if(e!=0) cout << "\t                      ";
		print_vec(int_elems[e]); cout << endl;
	}
	std::vector<Material*> elem_mats = elem->getIntegrationElemMaterials();
	std::vector<std::string> elem_mat_names;
	for(unsigned int e=0; e<elem_mats.size(); ++e)
		elem_mat_names.push_back(elem_mats[e]->get_name());
	cout << "\tIntegration Element Materials: "; print_vec(elem_mat_names); cout << endl;
}


int main (int argc, char* argv[])
{
	/*
	* TEST INCLUSION FUNCTIONALITY
	*	1. Parameter setting/getting
	*	2. Domain Limits
	*	3. Surface normal
	*	4. Node detection
	*	5. Intersection finding
	*/

/*
	// Plane inclusion first
	// Parameter set/get
	Plane_Inclusion plane;
	std::vector<double> point = {1.0, 1.0, 1.0};
	std::vector<double> normal = {1.0, 1.0, 1.0};
	plane.set_vec_parameter("norm", normal);
	plane.set_vec_parameter("p", point);
	bool norm_good=true, point_good=true;
	std::vector<double> norm_test = plane.get_vec_parameter("normal");
	std::vector<double> pnt_test = plane.get_vec_parameter("point");
	if(normal.size() != norm_test.size())
		norm_good=false;
	else
	{
		for(unsigned int i=0; i<normal.size(); ++i)
		{
			if(normal[i] != norm_test[i])
				norm_good=false;
		}
	}
	if(point.size() != pnt_test.size())
		point_good=false;
	else
	{
		for(unsigned int i=0; i<point.size(); ++i)
		{
			if(point[i] != pnt_test[i])
				point_good=false;
		}
	}
	if(norm_good)
	{
		cout << "Plane Normal set/get was good!\n\tInput: "; print_vec(normal); cout << "\n\tOutput: "; print_vec(norm_test); cout << "\n";
	}
	else
	{
		cout << "Plane Normal set/get was not good!\n\tInput: "; print_vec(normal); cout << "\n\tOutput: "; print_vec(norm_test); cout << "\n";
	}
	if(point_good)
	{
		cout << "Plane Point set/get was good!\n\tInput: "; print_vec(point); cout << "\n\tOutput: "; print_vec(pnt_test); cout << "\n";
	}
	else
	{
		cout << "Plane Point set/get was not good!\n\tInput: "; print_vec(point); cout << "\n\tOutput: "; print_vec(pnt_test); cout << "\n";
	}
	// Domain limits
	std::vector<double> dl_plane = plane.get_domain_lims();
	cout << "Domain Limits of Plane with normal "; print_vec(normal); cout << "\n\tLimits = "; print_vec(dl_plane); cout << endl;
	// Surface normal
	std::vector<double> surf_norm_plane = plane.get_surface_normal(NULL);
	cout << "Surface Normal of Plane with normal "; print_vec(normal); cout << "\n\tNormal = "; print_vec(surf_norm_plane); cout << endl;
	// Nodal detection
	Node node0(0.0, 1.0, 1.0, 0);
	Node node1(2.0, 2.0, 2.0, 1);
	Node node2(point[0], point[1], point[2]+1e-14, 2);
	unsigned int detect0_plane = plane.detect_node(&node0);
	unsigned int detect1_plane = plane.detect_node(&node1);
	unsigned int detect2_plane = plane.detect_node(&node2);
	cout << "Nodal detection value for node " << node0 << "\tDetection = " << detect0_plane << endl;
	cout << "Nodal detection value for node " << node1 << "\tDetection = " << detect1_plane << endl;
	cout << "Nodal detection value for node " << node2 << "\tDetection = " << detect2_plane << endl;
	// Intersection finding
	std::vector<Node> nodes_plane = {node1, node0};
	std::vector<double> intersect_plane = plane.find_intersection(nodes_plane);
	cout << "Intersection point for node:\n" << "\t" <<nodes_plane[0] << "\t" << nodes_plane[1] <<  "\t"; print_vec(intersect_plane); cout << endl;



	// Ellipse inclusion nect
	// Parameter set/get
	cout << "\n\n\n\n";
	Ellipse_Inclusion ellipse;
	std::vector<double> center = {0.0, 0.0};
	double alpha = pi/2;
	double a = 2;
	double b = 1;
	ellipse.set_vec_parameter("center", center);
	ellipse.set_parameter("alpha", alpha);
	ellipse.set_parameter("a", a);
	ellipse.set_parameter("b", b);
	bool center_good=true, alpha_good=true, a_good=true, b_good=true;
	std::vector<double> center_test = ellipse.get_vec_parameter("center");
	double alpha_test = ellipse.get_parameter("alpha");
	double a_test = ellipse.get_parameter("a");
	double b_test = ellipse.get_parameter("b");
	if(center.size() != center_test.size())
		center_good=false;
	else
	{
		for(unsigned int i=0; i<center.size(); ++i)
		{
			if(center[i] != center_test[i])
				center_good=false;
		}
	}
	if(alpha!=alpha_test)
		alpha_good = false;
	if(a!=a_test)
		a_good = false;
	if(b!=b_test)
		b_good = false;
	if(center_good)
	{
		cout << "Ellipse Center set/get was good!\n\tInput: "; print_vec(center); cout << "\n\tOutput: "; print_vec(center_test); cout << "\n";
	}
	else
	{
		cout << "Ellipse Center set/get was not good!\n\tInput: "; print_vec(center); cout << "\n\tOutput: "; print_vec(center_test); cout << "\n";
	}
	if(alpha_good)
		cout << "Ellipse Alpha set/get was good!\n\tInput: " << alpha << "\n\tOutput: " << alpha_test << "\n";
	if(alpha_good)
		cout << "Ellipse A set/get was good!\n\tInput: " << a << "\n\tOutput: " << a_test << "\n";
	if(alpha_good)
		cout << "Ellipse B set/get was good!\n\tInput: " << b << "\n\tOutput: " << b_test << "\n";
	// Domain limits
	std::vector<double> dl_ellipse = ellipse.get_domain_lims();
	cout << "Domain Limits of Ellipse with parameters:\n\tCenter: "; print_vec(center); cout << "\n\talpha: " << alpha << " a: " << a << " b: " << b << "\n\tLimits = "; print_vec(dl_ellipse); cout << endl;
	// Surface normal
	Node surf_node(1.0, 0.0, 0.0, 0);
	std::vector<double> surf_norm_ellipse = ellipse.get_surface_normal(&surf_node);
	cout << "Surface Normal of Ellipse at " << surf_node << "\tNormal = "; print_vec(surf_norm_ellipse); cout << endl;
	// Nodal detection
	Node node3(2.0, 1.0, 0.0, 0);
	Node node4(2.0, 2.0, 0.0, 1);
	Node node5(center[0], center[1], 0.0, 2);
	unsigned int detect3_ellipse = ellipse.detect_node(&node3);
	unsigned int detect4_ellipse = ellipse.detect_node(&node4);
	unsigned int detect5_ellipse = ellipse.detect_node(&node5);
	cout << "Nodal detection value for node " << node3 << "\tDetection = " << detect3_ellipse << endl;
	cout << "Nodal detection value for node " << node4 << "\tDetection = " << detect4_ellipse << endl;
	cout << "Nodal detection value for node " << node5 << "\tDetection = " << detect5_ellipse << endl;
	// Intersection finding
	std::vector<Node> nodes_ellipse = {node4, node5};
	std::vector<double> intersect_ellipse = ellipse.find_intersection(nodes_ellipse);
	cout << "Intersection point for node:\n" << "\t" <<nodes_ellipse[0] << "\t" << nodes_ellipse[1] <<  "\t"; print_vec(intersect_ellipse); cout << endl;




	// Ellipsoid inclusion next
	// Parameter set/get
	cout << "\n\n\n\n";
	Ellipsoid_Inclusion ellipsoid;
	std::vector<double> center_ellipsoid = {0.0, 0.0, 0.0};
	ellipsoid.set_vec_parameter("center", center_ellipsoid);
	double alpha_ellipsoid=0.0, beta_ellipsoid=0.0, gamma_ellipsoid=pi/2, a_ellipsoid=3.0, b_ellipsoid=0.75, c_ellipsoid=1.0;
	ellipsoid.set_parameter("alpha", alpha_ellipsoid);
	ellipsoid.set_parameter("beta", beta_ellipsoid);
	ellipsoid.set_parameter("gamma", gamma_ellipsoid);
	ellipsoid.set_parameter("a", a_ellipsoid);
	ellipsoid.set_parameter("b", b_ellipsoid);
	ellipsoid.set_parameter("c", c_ellipsoid);
	bool center_good_ellipsoid=true;
	std::vector<double> center_test_ellipsoid = ellipsoid.get_vec_parameter("center");
	if(center_ellipsoid.size() != center_test_ellipsoid.size())
		center_good_ellipsoid=false;
	else
	{
		for(unsigned int i=0; i<center_ellipsoid.size(); ++i)
		{
			if(center_ellipsoid[i] != center_test_ellipsoid[i])
				center_good_ellipsoid=false;
		}
	}
	if(center_good_ellipsoid)
	{
		cout << "Ellipsoid Center set/get was good!\n\tInput: "; print_vec(center_ellipsoid); cout << "\n\tOutput: "; print_vec(center_test_ellipsoid); cout << "\n";
	}
	if(alpha_ellipsoid==ellipsoid.get_parameter("alpha"))
		cout << "Ellipsoid Alpha set/get was good!\n\tInput: " << alpha_ellipsoid << "\n\tOutput: " << ellipsoid.get_parameter("alpha") << "\n";
	if(beta_ellipsoid==ellipsoid.get_parameter("beta"))
		cout << "Ellipsoid beta set/get was good!\n\tInput: " << beta_ellipsoid << "\n\tOutput: " << ellipsoid.get_parameter("beta") << "\n";
	if(gamma_ellipsoid==ellipsoid.get_parameter("gamma"))
		cout << "Ellipsoid gamma set/get was good!\n\tInput: " << gamma_ellipsoid << "\n\tOutput: " << ellipsoid.get_parameter("gamma") << "\n";
	if(a_ellipsoid==ellipsoid.get_parameter("a"))
		cout << "Ellipsoid a set/get was good!\n\tInput: " << a_ellipsoid << "\n\tOutput: " << ellipsoid.get_parameter("a") << "\n";
	if(b_ellipsoid==ellipsoid.get_parameter("b"))
		cout << "Ellipsoid b set/get was good!\n\tInput: " << b_ellipsoid << "\n\tOutput: " << ellipsoid.get_parameter("b") << "\n";
	if(c_ellipsoid==ellipsoid.get_parameter("c"))
		cout << "Ellipsoid c set/get was good!\n\tInput: " << c_ellipsoid << "\n\tOutput: " << ellipsoid.get_parameter("c") << "\n";

	// Domain limits
	std::vector<double> dl_ellipsoid = ellipsoid.get_domain_lims();
	cout << "Domain Limits of Ellipsoid with parameters:\n\tCenter: "; print_vec(center_ellipsoid); cout << "\n\talpha: " << alpha_ellipsoid  << " beta: " << beta_ellipsoid << " gamma: " << gamma_ellipsoid << " a: " << a_ellipsoid << " b: " << b_ellipsoid << " c: " << c_ellipsoid << "\n\tLimits = "; print_vec(dl_ellipsoid); cout << endl;
	// Surface normal
	Node node6(1.0/sqrt(3), 1.0/sqrt(3), 1.0/sqrt(3), 0);
	std::vector<double> surf_norm_ellipsoid = ellipsoid.get_surface_normal(&node6);
	cout << "Surface Normal of Ellipsoid at " << node6 << "\tNormal = "; print_vec(surf_norm_ellipsoid); cout << endl;
	// Nodal detection
	Node node7(4.0, 3.0, 2.0, 1);
	Node node8(center_ellipsoid[0], center_ellipsoid[1], center_ellipsoid[2], 2);
	unsigned int detect6_ellipsoid = ellipsoid.detect_node(&node6);
	unsigned int detect7_ellipsoid = ellipsoid.detect_node(&node7);
	unsigned int detect8_ellipsoid = ellipsoid.detect_node(&node8);
	cout << "Nodal detection value for node " << node6 << "\tDetection = " << detect6_ellipsoid << endl;
	cout << "Nodal detection value for node " << node7 << "\tDetection = " << detect7_ellipsoid << endl;
	cout << "Nodal detection value for node " << node8 << "\tDetection = " << detect8_ellipsoid << endl;
	// Intersection finding
	std::vector<Node> nodes_ellipsoid = {node7, node8};
	std::vector<double> intersect_ellipsoid = ellipsoid.find_intersection(nodes_ellipsoid);
	cout << "Intersection point for node:\n" << "\t" <<nodes_ellipsoid[0] << "\t" << nodes_ellipsoid[1] <<  "\t"; print_vec(intersect_ellipsoid); cout << endl;
*/




	/*
	* TEST ELEMENT DETECTION FUNCTIONALITY
	*	1. test several (all?) configurations of node detection combinations for each element type
	*	2. Check the cut edge vector as well as the ntegration element structure
	*/

/*
	// Some things 'll use in every element
	Material* base_mat = new LinearElasticIsotropicMaterial;
	Material* mat0 = new LinearElasticIsotropicMaterial;
	Material* mat1 = new LinearElasticIsotropicMaterial;
	base_mat->set_name("Base Material");
	mat0->set_name("Material 0");
	mat1->set_name("Material 1");
	std::vector<Material*> inc_mats = {base_mat, mat0, mat1};
	Node node0(0.0, 0.0, 0.0, 0);
	Node node1(1.0, 0.0, 0.0, 1);
	Node node2(1.0, 1.0, 0.0, 2);
	Node node3(0.0, 1.0, 0.0, 3);
	Node node4(0.0, 0.0, 1.0, 4);
	Node node5(1.0, 0.0, 1.0, 5);
	Node node6(1.0, 1.0, 1.0, 6);
	Node node7(0.0, 1.0, 1.0, 7);
*/
/*
	// Start with an Edge2
	std::vector<Node*> edge2_nodes = {&node0, &node1};
	Edge2 edge(edge2_nodes, 0);
	cout << "\n\n\n" << edge; // Print element information
	std::vector<int> edge2_detection;
	edge2_detection = {-1, 0};
	edge.detection(edge2_detection, inc_mats);
	display_element_detection(&edge, edge2_detection);
	edge2_detection = {0, -1};
	edge.detection(edge2_detection, inc_mats);
	display_element_detection(&edge, edge2_detection);
	edge2_detection = {-1, -1};
	edge.detection(edge2_detection, inc_mats);
	display_element_detection(&edge, edge2_detection);
	edge2_detection = {0, 0};
	edge.detection(edge2_detection, inc_mats);
	display_element_detection(&edge, edge2_detection);

	// Now do a Tri3
	std::vector<Node*> tri3_nodes = {&node0, &node1, &node2};
	Tri3 tri(tri3_nodes, 0);
	cout << "\n\n\n" << tri; // Print element information
	std::vector<int> tri3_detection;
	tri3_detection = {-1, 0, 0};
	tri.detection(tri3_detection, inc_mats);			// One node separated
	display_element_detection(&tri, tri3_detection);
	tri3_detection = {0, -1, 0};
	tri.detection(tri3_detection, inc_mats);
	display_element_detection(&tri, tri3_detection);
	tri3_detection = {0, 0, 1};
	tri.detection(tri3_detection, inc_mats);
	display_element_detection(&tri, tri3_detection);
	tri3_detection = {-1, -1, -1};						// All in base
	tri.detection(tri3_detection, inc_mats);
	display_element_detection(&tri, tri3_detection);
	tri3_detection = {0, 0, 0};							// All in inclusion
	tri.detection(tri3_detection, inc_mats);
	display_element_detection(&tri, tri3_detection);
//	tri3_detection = {-1, 0, 1};						// 3 inclusions, should be error
//	tri.detection(tri3_detection, inc_mats);
//	display_element_detection(&tri, tri3_detection);


	// Now do a Quad4
	std::vector<Node*> quad4_nodes = {&node0, &node1, &node2, &node3};
	Quad4 quad(quad4_nodes, 0);
	cout << "\n\n\n" << quad; // Print element information
	std::vector<int> quad4_detection;
	quad4_detection = {-1, 0, 0, 0};					// One node separated
	quad.detection(quad4_detection, inc_mats);
	display_element_detection(&quad, quad4_detection);
	quad4_detection = {0, -1, 0, 0};
	quad.detection(quad4_detection, inc_mats);
	display_element_detection(&quad, quad4_detection);
	quad4_detection = {0, 0, 1, 0};
	quad.detection(quad4_detection, inc_mats);
	display_element_detection(&quad, quad4_detection);
	quad4_detection = {0, 0, 0, -1};
	quad.detection(quad4_detection, inc_mats);
	display_element_detection(&quad, quad4_detection);
	quad4_detection = {-1, -1, 0, 0};					// Two nodes separated
	quad.detection(quad4_detection, inc_mats);
	display_element_detection(&quad, quad4_detection);
	quad4_detection = {0, -1, -1, 0};
	quad.detection(quad4_detection, inc_mats);
	display_element_detection(&quad, quad4_detection);
	quad4_detection = {1, 1, 1, 1};						// All in one inclusion
	quad.detection(quad4_detection, inc_mats);
	display_element_detection(&quad, quad4_detection);
//	quad4_detection = {1, 0, 1, 0};						// Two opposite nodes separated, should be error
//	quad.detection(quad4_detection, inc_mats);
//	display_element_detection(&quad, quad4_detection);
//	quad4_detection = {-1, 0, 1, 0};					// Three inclusions, should be error
//	quad.detection(quad4_detection, inc_mats);
//	display_element_detection(&quad, quad4_detection);


	// Now do a Tet4
	std::vector<Node*> tet4_nodes = {&node0, &node1, &node2, &node5};
	Tet4 tet(tet4_nodes, 0);
	cout << "\n\n\n" << tet; // Print element information
	std::vector<int> tet4_detection;
	tet4_detection = {-1, 0, 0, 0};						// One node separated
	tet.detection(tet4_detection, inc_mats);
	display_element_detection(&tet, tet4_detection);
	tet4_detection = {0, 1, 0, 0};
	tet.detection(tet4_detection, inc_mats);
	display_element_detection(&tet, tet4_detection);
	tet4_detection = {0, 0, -1, 0};
	tet.detection(tet4_detection, inc_mats);
	display_element_detection(&tet, tet4_detection);
	tet4_detection = {0, 0, 0, 1};	
	tet.detection(tet4_detection, inc_mats);
	display_element_detection(&tet, tet4_detection);
	tet4_detection = {1, 1, 0, 0};						// Two nodes separated
	tet.detection(tet4_detection, inc_mats);
	display_element_detection(&tet, tet4_detection);
	tet4_detection = {-1, 1, 1, -1};	
	tet.detection(tet4_detection, inc_mats);
	display_element_detection(&tet, tet4_detection);
	tet4_detection = {0, -1, 0, -1};	
	tet.detection(tet4_detection, inc_mats);
	display_element_detection(&tet, tet4_detection);
	tet4_detection = {0, 0, 0, 0};						// All in one inclusion/base
	tet.detection(tet4_detection, inc_mats);
	display_element_detection(&tet, tet4_detection);
//	tet4_detection = {0, -1, 1, -1};					// Three inclusions, should be error
//	tet.detection(tet4_detection, inc_mats);
//	display_element_detection(&tet, tet4_detection);

*/




	/* ACTUALLY TEST THE MESH
	*	1. Generate a mesh of quad, tris, or tets (hex's aren't ready yet)
	*	2. Partition the mesh
	*	3. Add materials
	*	4. Add inclusions (Test this both with materials in the mesh already and also with inclusions with new materials
	*	5. Perform IGFEM process
	*		a. Node detection. Test with all interfaces as well as with nodes that lie on an inclusion interface
	*		b. Element Detection. Kinda already tested this I guess
	*		c. Enrichment node generation. Make sure I have the correct number and that they appear in the correct places.
	*		   Already tested intersection finding but its good to verify that it works here too
	*		d. Integration element generation. This is a big one, hasn't been tested at all yet
	*	6. Assemble and solve matrix
	*		a. Not quite sure how to test the assembly process... I don't have a structural IGFEM code to conpare it to
	*		b. Same with the solve process. Other than figuring out some way to plot the deformation, I can't think of anything.
	*		   Run a convergence test relative to some conforming mesh example?
	*		c. Goodier Solution comparison
	*/


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
	mesh.set_dirichlet_bcs_from_nodeset("right", 0, x_disp);
	mesh.set_dirichlet_bcs_from_nodeset("left", 0, 0.0);
	for(unsigned int d=0; d<mesh.dim(); ++d) // bottom left node
		mesh.set_dirichlet_bc(0, d, 0);
	if(mesh.dim()==3)
		mesh.set_dirichlet_bc(nelem*(nelem+1), 2, 0); // back left node

	// Add a material to the mesh post-partitioning
	// NOTE: We don't really have to worry about adding maerials to each element like we did for a non-IGFEM mesh.
	//		 The element detection takes care of this by detectin that an element is part of the "matrix" that we define
	LinearElasticIsotropicMaterial material;
	material.set_parameter("mu", mu1);
	material.set_parameter("nu", nu1);
	material.set_name("Epoxy");
	mesh.add_material(&material);
	mesh.set_material("Epoxy");     // This is the first thing that is really inherent to IGFEM


	// Now that Everything has been assigned, lets distribute the dofs
	const clock_t dof_start = clock();
	mesh.distribute_dofs_elasticity();
	double dof_time = double(clock() - dof_start)/CLOCKS_PER_SEC;

	// Solve the system!
	double assem_time, solve_time;
	const clock_t s_start = clock();
	mesh.solve(assem_time, solve_time);
	double solve_acc_time = (double(clock() - s_start)/CLOCKS_PER_SEC) - assem_time - solve_time;

	// Compute L2-norm of the solution
	char buf[50];
	sprintf(buf, "FEM_out_%i.txt", nelem);
	mesh.print_solution( std::string(buf) );
	double L2 = mesh.Compute_L2_Norm(&Analytical_x_tension);

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
*/
	


	// Create an IGFEM mesh!
	IGFEM_Mesh mesh;
	mesh.init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	const clock_t g_start = clock();
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
	// This must match the Goodier solution


	//mesh.set_dirichlet_bcs_from_nodeset("right", 0, x_disp);
	//mesh.set_dirichlet_bcs_from_nodeset("left", 0, 0.0);
	//for(unsigned int d=0; d<mesh.dim(); ++d)
	//	mesh.set_dirichlet_bc(0, d, 0);


	std::vector<std::string> sets;
	if(mesh.dim()==2)
		sets = {"left", "right", "top", "bottom"};
	else if(mesh.dim()==3)
		sets = {"left", "right", "front", "back", "top", "bottom"};
	else
		err_message("Unknown mesh dimension");

	for(unsigned int s=0; s<sets.size(); ++s)
	{
		std::set<unsigned int> curr_set = mesh.get_nodeset(sets[s]);
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
				mesh.set_dirichlet_bc(*it, d, disp[d]);
		}
	}


	// Add a material to the mesh post-partitioning
	// NOTE: We don't really have to worry about adding maerials to each element like we did for a non-IGFEM mesh.
	//		 The element detection takes care of this by detectin that an element is part of the "matrix" that we define
	LinearElasticIsotropicMaterial material;
	material.set_parameter("mu", mu1);
	material.set_parameter("nu", nu1);
	material.set_name("Epoxy");
	mesh.add_material(&material);
	mesh.set_matrix_material("Epoxy");     // This is the first thing that is really inherent to IGFEM
	LinearElasticIsotropicMaterial material2;
	material2.set_parameter("mu", mu2);
	material2.set_parameter("nu", nu2);
	material2.set_name("Glass");
	mesh.add_material(&material2); // Test adding inclusion with and without this line

	// Create the inclusion
	
	//Plane_Inclusion plane;
	//std::vector<double> vec = {1.0, 1.0};
	//plane.set_vec_parameter("Point", vec);
	//vec = {1.0, 1.0};
	//plane.set_vec_parameter("norm", vec);
	//plane.set_material(&material2);
	//mesh.add_inclusion(&plane);
	
	if(mesh.dim()==2)
	{
		Ellipse_Inclusion ellipse;
		std::vector<double> center = {0, 0};
		ellipse.set_vec_parameter("center", center);
		ellipse.set_parameter("alpha", 0.0);
		ellipse.set_parameter("a", a);
		ellipse.set_parameter("b", a);
		ellipse.set_material(&material2);
		mesh.add_inclusion(&ellipse);	// This will add the material to the mesh is line 596? is commented out
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
		mesh.add_inclusion(&ellipsoid);	// This will add the material to the mesh is line 596? is commented out
	}
	else
		err_message("Unknown mesh dimension");
	
	

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
	mesh.distribute_dofs_elasticity();
	double dof_time = double(clock() - dof_start)/CLOCKS_PER_SEC;

	// Solve the system!
	double assem_time, solve_time;
	const clock_t s_start = clock();
	mesh.solve(assem_time, solve_time);
	double solve_acc_time = (double(clock() - s_start)/CLOCKS_PER_SEC) - assem_time - solve_time;


	char buf[50];
	sprintf(buf, "IGFEM_out_%i.txt", nelem);
	mesh.print_solution( std::string(buf) );

	// Compute L2-norm of the solution
	//double L2 = mesh.Compute_L2_Norm(&Analytical_x_tension);

	double L2;
	if(mesh.dim()==2)
		L2 = mesh.Compute_L2_Norm(&TwoD_Analytical_Inclusion_Disp);
	else if(mesh.dim()==3)
		L2 = mesh.Compute_L2_Norm(&ThreeD_Analytical_Inclusion_Disp);
	else
		err_message("Unknown mesh dimension");


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



	// Print out things to see what they are
	int rank_see = 0;
	if(rank==rank_see)
	{
		//cout << "PARTITION " << rank_see << "\n--------------------------------------------------\n";
		cout << "Number of Elements: " << mesh.n_elem() << endl;
		cout << "Number of Nodes: " << mesh.n_nodes() << endl;
		cout << "Number of dofs: " << mesh.n_global_dofs() << endl;
/*
		// Print out the nodal detection result
		cout << "\n\nNodal Detection of the non-enrichment nodes:\n";
		print_vec(mesh.get_nodal_detection());


		// Print out all of the enrichment nodes
		cout << "\n\nEnrichment nodes in the local mesh:\n";
		for(unsigned int n=0; n<mesh.n_local_enrich_nodes(); ++n)
			cout << *mesh.get_enrich_node(n) << " Owner (" << mesh.get_enrich_node_owner_local(n) << ")\n";

		// Print out all integration elements
		cout  << "\n\nIntegration elements for intersected elements:\n";
		for(unsigned int e=0; e<mesh.n_local_elem(); ++e)
		{
			Elem* el = mesh.get_elem_local(e);
			if(el->getElemDetector()!=0)
			{
				cout << *el << "Integration elements for the above element:\n";
				for(unsigned int ie = 0; ie<el->n_integration_elem(); ++ie)
				{
					cout << el->get_int_elem_mat(ie)->get_name() << endl;
					cout << *el->get_integration_elem(ie);
				}
				cout << "\n";
			}
		}

		// Output the _node_elem table
		cout << "Local Node_Elem table:" << endl;
		for(unsigned int n=0; n<(mesh.n_local_nodes()+mesh.n_local_enrich_nodes()); ++n)
		{
			unsigned int id = mesh.get_node_local(n)->get_id();
			std::vector<unsigned int> v = mesh.get_node_elem_local(n);
			if(n< mesh.n_local_nodes())
				cout << "Node(" << id << "): ";
			else
				cout << "Enrichment Node(" << id << "): ";
			print_vec(v);
			cout << endl;
		}


		// Print out all of the nodal dofs
		// Output Nodal degrees of freedom
		cout << "\n\nNodal Degrees of Freedom\n";
		cout << "Number of free degrees of freedom: " << mesh.n_global_free_dofs() << endl;
		cout << "Number of constrained degrees of freedom: " << mesh.n_global_const_dofs() << endl;
		cout << "Number of local free degrees of freedom: " << mesh.n_local_free_dofs() << endl;
		cout << "Number of local constrained degrees of freedom: " << mesh.n_local_const_dofs() << endl;
		for(unsigned int n=0; n<(mesh.n_local_nodes()+mesh.n_local_enrich_nodes()); ++n)
		{
			unsigned int id = mesh.get_node_local(n)->get_id();
			if(n<mesh.n_local_nodes())
				cout << "Node(" << id << "): {";
			else
				cout << "Enrichment Node(" << id << "): {";
			
			std::vector<unsigned int> dofs = mesh.get_nodal_global_dofs_local(n);
			
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
		// Output the solution
		cout << "\nSolution:\n\n";
		for(unsigned int n=0; n<(mesh.n_local_nodes()+mesh.n_local_enrich_nodes()); ++n)
		{
			unsigned int g_node = mesh.get_node_local(n)->get_id();
			for(unsigned int d=0; d<mesh.dim(); ++d)
			{
				if(n<mesh.n_local_nodes())
					std::cout << "Node(";
				else
					std::cout << "Enrichment Node(";

				std::cout << g_node << ") dof(" << d << "): " << mesh.get_solution_local(n, d) << std::endl;
			}
		}

	
		// Output the reaction force
		cout << "\nReaction Force:\n\n";
		for(unsigned int n=0; n<(mesh.n_local_nodes()+mesh.n_local_enrich_nodes()); ++n)
		{
			unsigned int g_node = mesh.get_node_local(n)->get_id();
			for(unsigned int d=0; d<mesh.dim(); ++d)
			{
				if(n<mesh.n_local_nodes())
					std::cout << "Node(";
				else
					std::cout << "Enrichment Node(";

				std::cout << g_node << ") dof(" << d << "): " << mesh.get_force_local(n, d) << std::endl;
			}
		}
*/
		cout << "\nL2 Norm = " << L2 << "\n\n";

	}
	
	
	mesh.Finalize();
}
