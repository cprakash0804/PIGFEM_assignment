/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#include "CohPoint1.h"


id_type CohPoint1::n_q_points(int order) const
{
	return 1;
}
void CohPoint1::q_point(std::vector<double>& coords, double& w, id_type qp, int order) const
{
	coords.clear();
	coords.push_back(0.0);
	w = 1.0;
}

double CohPoint1::compute_shape(const std::vector<double>& coords, int n) const
{
	if(n!=0)
		err_message("Shape function number for an Edge2 must be 0");
	
	return 1.0;
}

std::vector<double> CohPoint1::compute_shape_grad(const std::vector<double>& coords, int n) const
{
	if(n!=0)
		err_message("Shape function number for an Edge2 must be 0");
	
	return std::vector<double>(1);
}

// void CohPoint1::compute_rotation_general(const std::vector<double>& coords, const std::vector<std::vector<double> >& node_coords,
// 										DenseMatrix<double>& mat, double& dA) const
// {
// 	mat.resize(1,1);
// 	mat(0,0) = 1.0;
// 	dA = 1; // ?????????????????????????????????
// }
