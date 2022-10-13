/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#include "CohQuad4.h"


id_type CohQuad4::n_q_points(int order) const
{
	if(order < 0)
		err_message("Order of integration must be greater than or equal to 0.");
	
	switch(order)
	{
		case 0:
		case 1:
			return 1;
		case 2:
		case 3:
			return 4;
		case 4:
		case 5:
			return 9;
		case 6:
		case 7:
			return 16;
		case 8:
		case 9:
			return 25;
		case 10:
		case 11:
			return 36;
		case 12:
		case 13:
			return 49;
		case 14:
		case 15:
			return 64;
		case 16:
		case 17:
			return 81;
		case 18:
		case 19:
			return 100;
		default:
			err_message("The selected quadrature order is not currently supported.");

	}
}
void CohQuad4::q_point(std::vector<double>& coords, double& w, id_type qp, int order) const
{
	if(order < 0)
		err_message("Order of integration must be greater than or equal to 0.");
	if(qp<0 || qp>=n_q_points(order))
		err_message("Quadrature point selected must be less than the maximum number of quadrature points for the given order.");
	
	coords.clear();
	CohesiveElem* edge = CohesiveElem::build(COHEDGE2);
	id_type nqp_edge = edge->n_q_points(order);
	
	std::vector<double> edge0_qp, edge1_qp;
	double edge0_w, edge1_w;
	edge->q_point(edge0_qp, edge0_w, (qp%nqp_edge), order);
	edge->q_point(edge1_qp, edge1_w, (qp/nqp_edge), order);
	delete edge;
	
	coords = {edge0_qp[0], edge1_qp[0]};
	w = edge0_w*edge1_w;
}




double CohQuad4::compute_shape(const std::vector<double>& coords, int n) const
{
	if(n<0 || n>3)
		err_message("Shape function number for an Quad4 must be between 0 and 3");
	double r = coords[0];
	double s = coords[1];
	if(n==0)
		return 0.25*(1.0-r)*(1.0-s);
	else if(n==1)
		return 0.25*(1.0+r)*(1.0-s);
	else if(n==2)
		return 0.25*(1.0+r)*(1.0+s);
	else if(n==3)
		return 0.25*(1.0-r)*(1.0+s);
	else
		err_message("Shape function selected must be less than 4 for a Quad4.");
}

std::vector<double> CohQuad4::compute_shape_grad(const std::vector<double>& coords, int n) const
{
	if(n<0 || n>3)
		err_message("Shape function number for an Quad4 must be between 0 and 3");
	
	double r = coords[0];
	double s = coords[1];
	std::vector<double> v;
	if(n==0)
		v = {-0.25*(1.0-s), -0.25*(1.0-r)};
	else if(n==1)
		v = {0.25*(1.0-s), -0.25*(1.0+r)};
	else if(n==2)
		v = {0.25*(1.0+s), 0.25*(1.0+r)};
	else if(n==3)
		v = {-0.25*(1.0+s), 0.25*(1.0-r)};
	else
		err_message("Shape function selected must be less than 4 for a Quad4.");
	return v;
}


// void CohQuad4::compute_rotation_general(const std::vector<double>& coords, const std::vector<std::vector<double> >& node_coords,
// 										DenseMatrix<double>& mat, double& dA) const
// {
// 	if(coords.size()!=2)
// 		err_message("Parent coordinates for a 2D element must only contain two values");

// 	// Initialize the array of shape function gradients wrt parent coords
// 	id_type n_nodes = node_coords.size();
// 	std::vector<std::vector<double> > parent_grad(n_nodes);
	
// 	for(id_type n=0; n<n_nodes; ++n)
// 		parent_grad[n] = compute_shape_grad(coords, n);
	
// 	// Compute dx/dr
// 	double dxdr[2][3] = {{0}};  // Initiaizes the whole array to 0
// 	for(int i=0; i<2; ++i)  // Loop over parent coords
// 	{
// 		for(int j=0; j<3; ++j)  // Loop over global coords
// 		{
// 			for(id_type n=0; n<n_nodes; ++n)
// 				dxdr[i][j] = dxdr[i][j] + parent_grad[n][i]*node_coords[n][j];
// 		}
// 	}

// 	// Cross product to form the normal vector
// 	std::vector<double> norm(3);
// 	norm[0] = dxdr[0][1]*dxdr[1][2] - dxdr[1][1]*dxdr[0][2];
// 	norm[1] = dxdr[1][0]*dxdr[0][2] - dxdr[0][0]*dxdr[1][2];
// 	norm[2] = dxdr[0][0]*dxdr[0][1] - dxdr[1][0]*dxdr[0][1];

// 	// Compute magnitude of each vec to normalize them
// 	double norm_mag = 0;
// 	double t1_mag = 0;
// 	double t2_mag = 0;
// 	for(int i=0; i<3; ++i)
// 	{
// 		norm_mag += pow(norm[i], 2);
// 		t1_mag += pow(dxdr[0][i], 2);
// 		t2_mag += pow(dxdr[1][i], 2);
// 	}
// 	norm_mag = sqrt(norm_mag);
// 	t1_mag = sqrt(t1_mag);
// 	t2_mag = sqrt(t2_mag);

// 	// Set the differential area (similar to jacobian for volumetric elements)
// 	dA = norm_mag; // ?????????????????????????

// 	// Store the roation matrix contributions
// 	mat.resize(3,3);
// 	mat(0,0) = norm[0]/norm_mag;
// 	mat(0,1) = norm[1]/norm_mag;
// 	mat(0,2) = norm[2]/norm_mag;
// 	mat(1,0) = dxdr[0][0]/t1_mag;
// 	mat(1,1) = dxdr[0][1]/t1_mag;
// 	mat(1,2) = dxdr[0][2]/t1_mag;
// 	mat(2,0) = dxdr[1][0]/t2_mag;
// 	mat(2,1) = dxdr[1][1]/t2_mag;
// 	mat(2,2) = dxdr[1][2]/t2_mag;
// }
