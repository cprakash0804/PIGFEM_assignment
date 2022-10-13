/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#include "CohEdge2.h"


id_type CohEdge2::n_q_points(int order) const
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
			return 2;
		case 4:
		case 5:
			return 3;
		case 6:
		case 7:
			return 4;
		case 8:
		case 9:
			return 5;
		case 10:
		case 11:
			return 6;
		case 12:
		case 13:
			return 7;
		case 14:
		case 15:
			return 8;
		case 16:
		case 17:
			return 9;
		case 18:
		case 19:
			return 10;
		default:
			err_message("The selected quadrature order is not currently supported.");

	}
}
void CohEdge2::q_point(std::vector<double>& coords, double& w, id_type qp, int order) const
{
	if(order < 0)
		err_message("Order of integration must be greater than or equal to 0.");
	if(qp<0 || qp>=n_q_points(order))
		err_message("Quadrature point selected must be less than the maximum number of quadrature points for the given order.");
	
	coords.clear();
	switch(order)
	{
		case 0:
		case 1:
		{
			coords = {0.0};
			w = 2.0;
			break;
		}
		case 2:
		case 3:
		{
			w = 1.0; // For all quadrature points in this case
			if(qp==0) {coords = {-0.57735026918962576450914878050196};}
			else if(qp==1) {coords = {0.57735026918962576450914878050196};}
			else err_message("Quadrature point for 3rd order quadrature of an Edge must be less than 2");
			break;
		}
		case 4:
		case 5:
		{
			if(qp==0) {coords = {-0.77459666924148337703585307995648}; w = 0.55555555555555555555555555555556;}
			else if(qp==1) {coords = {0.0}; w = 0.88888888888888888888888888888889;}
			else if(qp==2) {coords = {7.7459666924148337703585307995648}; w = 0.55555555555555555555555555555556;}
			else err_message("Quadrature point for 5th order quadrature of an Edge must be less than 3");
			break;
		}
		case 6:
		case 7:
		{
			if(qp==0) {coords = {-0.86113631159405257522394648889281}; w = 0.34785484513745385737306394922200;}
			else if(qp==1) {coords = {-0.33998104358485626480266575910324}; w = 0.65214515486254614262693605077800;}
			else if(qp==2) {coords = {0.33998104358485626480266575910324}; w = 0.65214515486254614262693605077800;}
			else if(qp==3) {coords = {0.86113631159405257522394648889281}; w = 0.34785484513745385737306394922200;}
			else err_message("Quadrature point for 7th order quadrature of an Edge must be less than 4");
			break;
		}
		case 8:
		case 9:
		{
			if(qp==0) {coords = {-0.9061798459386640}; w = 0.2369268850561891;}
			else if(qp==1) {coords = {-0.5384693101056831}; w = 0.4786286704993665;}
			else if(qp==2) {coords = {0.0}; w = 0.5688888888888889;}
			else if(qp==3) {coords = {0.5384693101056831}; w = 0.4786286704993665;}
			else if(qp==4) {coords = {0.9061798459386640}; w = 0.2369268850561891;}
			else err_message("Quadrature point for 9th order quadrature of an Edge must be less than 5");
			break;
		}
		case 10:
		case 11:
		{
			if(qp==0) {coords = {-0.9324695142031521}; w = 0.1713244923791704;}
			else if(qp==1) {coords = {-0.6612093864662645}; w = 0.3607615730481386;}
			else if(qp==2) {coords = {-0.2386191860831969}; w = 0.4679139345726910;}
			else if(qp==3) {coords = {0.2386191860831969}; w = 0.4679139345726910;}
			else if(qp==4) {coords = {0.6612093864662645}; w = 0.3607615730481386;}
			else if(qp==5) {coords = {0.9324695142031521}; w = 0.1713244923791704;}
			else err_message("Quadrature point for 11th order quadrature of an Edge must be less than 6");
			break;
		}
		case 12:
		case 13:
		{
			if(qp==0) {coords = {-0.9491079123427585}; w = 0.1294849661688697;}
			else if(qp==1) {coords = {-0.7415311855993945}; w = 0.2797053914892766;}
			else if(qp==2) {coords = {-0.4058451513773972}; w = 0.3818300505051189;}
			else if(qp==3) {coords = {0.0}; w = 0.4179591836734694;}
			else if(qp==4) {coords = {0.4058451513773972}; w = 0.3818300505051189;}
			else if(qp==5) {coords = {0.7415311855993945}; w = 0.2797053914892766;}
			else if(qp==6) {coords = {0.9491079123427585}; w = 0.1294849661688697;}
			else err_message("Quadrature point for 13th order quadrature of an Edge must be less than 7");
			break;
		}
		case 14:
		case 15:
		{
			if(qp==0) {coords = {-0.9602898564975363}; w = 0.1012285362903763;}
			else if(qp==1) {coords = {-0.7966664774136267}; w = 0.2223810344533745;}
			else if(qp==2) {coords = {-0.5255324099163290}; w = 0.3137066458778873;}
			else if(qp==3) {coords = {-0.1834346424956498}; w = 0.3626837833783620;}
			else if(qp==4) {coords = {0.1834346424956498}; w = 0.3626837833783620;}
			else if(qp==5) {coords = {0.5255324099163290}; w = 0.3137066458778873;}
			else if(qp==6) {coords = {0.7966664774136267}; w = 0.2223810344533745;}
			else if(qp==7) {coords = {0.9602898564975363}; w = 0.1012285362903763;}
			else err_message("Quadrature point for 15th order quadrature of an Edge must be less than 8");
			break;
		}
		case 16:
		case 17:
		{
			if(qp==0) {coords = {-0.9681602395076261}; w = 0.0812743883615744;}
			else if(qp==1) {coords = {-0.8360311073266358}; w = 0.1806481606948574;}
			else if(qp==2) {coords = {-0.6133714327005904}; w = 0.2606106964029354;}
			else if(qp==3) {coords = {-0.3242534234038089}; w = 0.3123470770400029;}
			else if(qp==4) {coords = {0.0}; w = 0.3302393550012598;}
			else if(qp==5) {coords = {0.3242534234038089}; w = 0.3123470770400029;}
			else if(qp==6) {coords = {0.6133714327005904}; w = 0.2606106964029354;}
			else if(qp==7) {coords = {0.8360311073266358}; w = 0.1806481606948574;}
			else if(qp==8) {coords = {0.9681602395076261}; w = 0.0812743883615744;}
			else err_message("Quadrature point for 17th order quadrature of an Edge must be less than 9");
			break;
		}
		case 18:
		case 19:
		{
			if(qp==0) {coords = {-0.9739065285171717}; w = 0.0666713443086881;}
			else if(qp==1) {coords = {-0.8650633666889845}; w = 0.1494513491505806;}
			else if(qp==2) {coords = {-0.6794095682990244}; w = 0.2190863625159820;}
			else if(qp==3) {coords = {-0.4333953941292472}; w = 0.2692667193099963;}
			else if(qp==4) {coords = {-0.1488743389816312}; w = 0.2955242247147529;}
			else if(qp==5) {coords = {0.1488743389816312}; w = 0.2955242247147529;}
			else if(qp==6) {coords = {0.4333953941292472}; w = 0.2692667193099963;}
			else if(qp==7) {coords = {0.6794095682990244}; w = 0.2190863625159820;}
			else if(qp==8) {coords = {0.8650633666889845}; w = 0.1494513491505806;}
			else if(qp==9) {coords = {0.9739065285171717}; w = 0.0666713443086881;}
			else err_message("Quadrature point for 19th order quadrature of an Edge must be less than 10");
			break;
		}
		default:
			err_message("The selected quadrature order is not currently supported.");
	}
}

double CohEdge2::compute_shape(const std::vector<double>& coords, int n) const
{
	if(n<0 || n>1)
		err_message("Shape function number for an Edge2 must be between 0 and 1");
	if(n==0)
		return 0.5*(1.0-coords[0]);
	else
		return 0.5*(1+coords[0]);
}

std::vector<double> CohEdge2::compute_shape_grad(const std::vector<double>& coords, int n) const
{
	if(n<0 || n>1)
		err_message("Shape function number for an Edge2 must be between 0 and 1");
	
	std::vector<double> v(1);
	if(n==0)
		v[0] = -0.5;
	else if(n==1)
		v[0] = 0.5;
	else
		err_message("Shape function number for an Edge2 must be between 0 and 1");
	return v;
}


// void CohEdge2::compute_rotation_general(const std::vector<double>& coords, const std::vector<std::vector<double> >& node_coords,
// 										DenseMatrix<double>& mat, double& dA) const
// {
// 	if(coords.size()!=1)
// 		err_message("Parent coordinates for a 1D element must only contain one value");
	
// 	// Compute dx/dr & dy/dr
// 	double dxdr = 0;
// 	double dydr = 0;
// 	for(id_type i=0; i<node_coords.size(); ++i) // Loop over the nodes
// 	{
// 		std::vector<double> parent_grad = compute_shape_grad(coords, i);
// 		dxdr = dxdr + parent_grad[0]*node_coords[i][0];
// 		dydr = dydr + parent_grad[0]*node_coords[i][1];
// 	}

// 	// Compute the tangent and normal vectors
// 	// Note, due to my desire for the normal vector to be first in the rotation matrix,
// 	//  the tangent vector needs to be reversed to maintain a right-handed coord system
// 	std::vector<double> t1 = {-1.0*dxdr, -1.0*dydr};
// 	std::vector<double> norm = {t1[1], -1.0*t1[0]};
// 	double mag = sqrt(pow(dxdr,2) + pow(dydr,2));

// 	// Set the differential area (similar to jacobian for volumetric elements)
// 	dA = mag; // ??????????????????????????????????

// 	mat.resize(2,2);
// 	// Normal vector
// 	mat(0,0) = norm[0]/mag;
// 	mat(0,1) = norm[1]/mag;
// 	// Tangent vector
// 	mat(1,0) = t1[0]/mag;
// 	mat(1,1) = t1[1]/mag;
// }


// Output helper function
std::vector<id_type> CohEdge2::plot_elem_ids()
{
	std::vector<id_type> v(n_nodes());
	v[0] = _nodes[0]->get_id();
	v[1] = _nodes[1]->get_id();
	v[2] = _nodes[3]->get_id();
	v[3] = _nodes[2]->get_id();
	return v;
}