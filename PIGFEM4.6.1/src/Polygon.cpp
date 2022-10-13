/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "Polygon.h"
#include <cmath>
#include <iostream>

id_type Polygon::side_nodes_map(id_type side, id_type node)
{
	return edge_nodes_map(side, node);
}

id_type Polygon::edge_nodes_map(id_type edge, id_type node)
{
	if(edge>=_nn)
		err_message("Please select a valid edge.");
	if(node>=2)
		err_message("Please select a valid node.");
	
	if (node == 0)
		return edge;
	else
	{
		if (edge < (_nn-1))
			return edge + 1;
		else
			return 0;
	}
}

Elem* Polygon::build_side(id_type side)
{
	if(side >= n_sides())
		err_message("Side number must be less than the number of sides.");

	Elem* el = build(EDGE2);
	std::vector<Node*> nodes;
	for(id_type n=0; n<2; ++n)
		nodes.push_back(_nodes[side_nodes_map(side,n)]);
	el->set_nodes(nodes);
	return el;
}


void Polygon::PolyTriangulate(std::vector<std::vector<double> >& p, std::vector<std::vector<id_type> >& Tri, std::vector<double> rcoords) const
{
	p.resize(_nn + 1);
	Tri.resize(_nn);

	for (id_type n=0; n<_nn; ++n)
	{
		p[n] = {cos(2.0 * PI * (n+1)/_nn), sin(2.0 * PI * (n+1)/_nn)};
		Tri[n] = {_nn, n, n+1};
	}
	p[_nn] = rcoords;
	Tri[_nn-1][2] = 0;
}


id_type Polygon::n_q_points(int order) const
{
	if(order < 0)
		err_message("Order of integration must be greater than or equal to 0.");
	
	Elem* tri = Elem::build(TRI3);
	id_type nqp = _nn * tri->n_q_points(order);
	delete tri;
	return nqp;
}
// Triangular quadrature is based off of the quadrature schemes presented in:
//	http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
// Note: all quadrature weights here are half of those presened in this paper.
void Polygon::q_point(std::vector<double>& coords, double& w, id_type qp, int order) const
{
	if(order < 0)
		err_message("Order of integration must be greater than or equal to 0.");
	if(qp<0 || qp>=n_q_points(order))
		err_message("Quadrature point selected must be less than the maximum number of quadrature points for the given order.");

	// Get the triangular quadrature point
	Elem* tri = Elem::build(TRI3);
	std::vector<double> tri_coords;
	double tri_w;
	id_type tri_nqp = tri->n_q_points(order);
	tri->q_point(tri_coords, tri_w, qp%tri_nqp, order);

	// Get the triangular shape function and gradients
	std::vector<double> tri_N = tri->compute_shape(tri_coords);
	std::vector<std::vector<double> > tri_dN = tri->compute_shape_grad(tri_coords);

	// Compute the subtriangles and get the right one
	id_type curr_tri = qp / tri_nqp;
	std::vector<std::vector<double> > p;
	std::vector<std::vector<id_type> > Tri;
	std::vector<double> origin = {0.0, 0.0};
	PolyTriangulate(p, Tri, origin); // Triangulate from the origin
	std::vector<std::vector<double> > pT(3);
	for (id_type n2=0; n2<3; ++n2)
		pT[n2] = p[Tri[curr_tri][n2]];

	// Compute the quadrature weight
	// Compute the mapping jacobian from the current triangle
	double J0[2][2] = {{0}};
	for (id_type i=0; i<2; ++i)
		for (id_type j=0; j<2; ++j)
			for(id_type n=0; n<3; ++n)
				J0[i][j] += pT[n][i] * tri_dN[n][j];
	double det_J0 = J0[0][0]*J0[1][1] - J0[0][1]*J0[1][0];
	w = det_J0 * tri_w;

	// Compute the quadrature point
	coords = {0.0, 0.0};
	for (id_type n=0; n<3; ++n)
	{
		coords[0] += tri_N[n] * pT[n][0];
		coords[1] += tri_N[n] * pT[n][1];
	}

	// Clean up memory
	delete tri;
}

std::vector<double> Polygon::compute_shape(const std::vector<double>& coords) const
{
	// Set up necessary structures
	std::vector<double> alpha(_nn);
	std::vector<double> A(_nn+1);

	// Compute the subtriangles
	std::vector<std::vector<double> > p;
	std::vector<std::vector<id_type> > Tri;
	PolyTriangulate(p, Tri, coords);

	// Compute the triangle sizes for each triangle
	for (id_type n=0; n<_nn; ++n)
	{
		std::vector<std::vector<double> > pT(3);
		for (id_type n2=0; n2<3; ++n2)
			pT[n2] = p[Tri[n][n2]];
		A[n+1] = 0.5*( (pT[1][0]*pT[2][1] - pT[1][1]*pT[2][0]) -
				 	   (pT[0][0]*pT[2][1] - pT[0][1]*pT[2][0]) + 
				 	   (pT[0][0]*pT[1][1] - pT[0][1]*pT[1][0]) );
	}
	A[0] = A[_nn];

	// Compute the alpha values
	double sum_alpha = 0.0;
	std::vector<double> sum_dalpha(2,0.0);
	for (id_type n=0; n<_nn; ++n)
	{
		alpha[n] = 1.0/(A[n] * A[n+1]);
		sum_alpha += alpha[n];
	}

	// Compute the shape functions
	for (id_type n=0; n<_nn; ++n)
		alpha[n] /= sum_alpha;
	return alpha;
}

std::vector<std::vector<double> > Polygon::compute_shape_grad(const std::vector<double>& coords) const
{
	// Fill in the shape functions
	std::vector<std::vector<double> > dNdxi(n_nodes());

	// Set up necessary structures
	std::vector<double> alpha(_nn);
	std::vector<std::vector<double> > dalpha(_nn);
	std::vector<double> A(_nn+1);
	std::vector<std::vector<double> > dA(_nn+1);
	for (id_type n=0; n<_nn; ++n){
		dalpha[n].resize(2);
		dA[n].resize(2);
	}
	dA[_nn].resize(2);


	// Compute the subtriangles
	std::vector<std::vector<double> > p;
	std::vector<std::vector<id_type> > Tri;
	PolyTriangulate(p, Tri, coords);


	// Compute the triangle sizes for each triangle
	for (id_type n=0; n<_nn; ++n)
	{
		std::vector<std::vector<double> > pT(3);
		for (id_type n2=0; n2<3; ++n2)
			pT[n2] = p[Tri[n][n2]];
		A[n+1] = 0.5*( (pT[1][0]*pT[2][1] - pT[1][1]*pT[2][0]) -
				 	   (pT[0][0]*pT[2][1] - pT[0][1]*pT[2][0]) + 
				 	   (pT[0][0]*pT[1][1] - pT[0][1]*pT[1][0]) );
		dA[n+1][0] = 0.5*(pT[2][1] - pT[1][1]);
		dA[n+1][1] = 0.5*(pT[1][0] - pT[2][0]);
	}
	A[0] = A[_nn];
	dA[0] = dA[_nn];

	// Compute the alpha values
	double sum_alpha = 0.0;
	std::vector<double> sum_dalpha(2,0.0);
	for (id_type n=0; n<_nn; ++n)
	{
		alpha[n] = 1.0/(A[n] * A[n+1]);
		dalpha[n][0] = -1.0*alpha[n] * (dA[n][0]/A[n] + dA[n+1][0]/A[n+1]);
		dalpha[n][1] = -1.0*alpha[n] * (dA[n][1]/A[n] + dA[n+1][1]/A[n+1]);
		sum_alpha += alpha[n];
		sum_dalpha[0] += dalpha[n][0];
		sum_dalpha[1] += dalpha[n][1];
	}

	// Compute the shape functions and the local gradients
	for (id_type n=0; n<n_nodes(); ++n)
	{
		double N = alpha[n] / sum_alpha;
		dNdxi[n].resize(2);
		dNdxi[n][0] = (dalpha[n][0] - N*sum_dalpha[0]) / sum_alpha;
		dNdxi[n][1] = (dalpha[n][1] - N*sum_dalpha[1]) / sum_alpha;
	}

	return dNdxi;
}

std::vector<DenseMatrix<double> > Polygon::compute_shape_grad_grad(const std::vector<double>& coords) const
{
	err_message("Shape function Laplacian only needs to be computed for parent elements and a Polygon should never be a parent element.");
}




bool Polygon::coords_inside(const std::vector<double>& rcoords)
{
	if (rcoords.size() != 2)
		err_message("Invalid coordinate size for coordinate inside check (POLYGON)!");

	// Figure out the angle and radius
	double theta = atan2(rcoords[1], rcoords[0]);
	double r = sqrt(rcoords[0]*rcoords[0] + rcoords[1]*rcoords[1]);

	double alpha = 2.0 * PI / _nn;
	double beta = fabs( std::fmod(theta, alpha) );
	double gamma = PI - beta - (PI-alpha)/2.0;
	double r_min = sin((PI-alpha)/2.0) / sin(gamma);

	if (r <= (r_min+_inside_tol))
		return true;
	else
		return false;
}



// IGFEM MATERIAL BELOW HERE
//==============================================================================================================================

void Polygon::detection(std::vector<int>& node_detection, std::vector<Material*>& mats, bool isCohesive)
{
	err_message("Element detection is not supported for n-gon elements.");
} // detection








void Polygon::refinement_nodes(std::vector<std::vector<id_type> >& refine_nodes)
{
	err_message("Refinement is not supported for n-gon elements.");
}
void Polygon::refinement_structure(std::vector<std::vector<id_type> >& structure, std::vector<elem_type>& types)
{
	err_message("Refinement is not supported for n-gon elements.");
}
