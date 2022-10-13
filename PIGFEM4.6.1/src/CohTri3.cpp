/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#include "CohTri3.h"


id_type CohTri3::n_q_points(int order) const
{
	if(order < 0)
		err_message("Order of integration must be greater than or equal to 0.");
	
	switch(order)
	{
		case 0:
		case 1:
			return 1;
		case 2:
			return 3;
		case 3:
		default:
			return 4;
		case 4:
			return 6;
		case 5:
			return 7;
		case 6:
			return 12;
		case 7:
			return 13;
		case 8:
			return 16;
	}
}
// Triangular quadrature is based off of the quadrature schemes presented in:
//	http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
// Note: all quadrature weights here are half of those presened in this paper.
void CohTri3::q_point(std::vector<double>& coords, double& w, id_type qp, int order) const
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
			if(qp==0) {coords = {0.333333333333333333333333333333333, 0.333333333333333333333333333333333}; w = 0.5;}
			else err_message("Quadrature point for 1st order quadrature of a Tri must be less than 1");
			break;
		}
		case 2:
		{
			if(qp==0) {coords = {0.166666666666666666666666666666667, 0.166666666666666666666666666666667}; w =  0.166666666666666666666666666666667;}
			else if(qp==1) {coords = {0.166666666666666666666666666666667, 0.666666666666666666666666666666667}; w = 0.166666666666666666666666666666667;}
			else if(qp==2) {coords = {0.666666666666666666666666666666667, 0.166666666666666666666666666666667}; w = 0.166666666666666666666666666666667;}
			else err_message("Quadrature point for 2nd order quadrature of a Tri must be less than 3");
			break;
		}
		default:    // I think I can put this not at the end...
		case 3:
		{
			/*
			if(qp==0) {coords = {1.5505102572168219018027159252941e-01L, 1.7855872826361642311703513337422e-01L}; w = 1.5902069087198858469718450103758e-01L;}
			else if(qp==1) {coords = {6.4494897427831780981972840747059e-01L, 7.5031110222608118177475598324603e-02L}; w = 9.0979309128011415302815498962418e-02L;}
			else if(qp==2) {coords = {1.5505102572168219018027159252941e-01L, 6.6639024601470138670269327409637e-01L}; w = 1.5902069087198858469718450103758e-01L;}
			else if(qp==3) {coords = {6.4494897427831780981972840747059e-01L, 2.8001991549907407200279599420481e-01L}; w = 9.0979309128011415302815498962418e-02L;}
			else err_message("Quadrature point for 3rd order quadrature of a Tri must be less than 4");
			break;
			*/
			if(qp==0) {coords = {0.333333333333333333333333333333333, 0.333333333333333333333333333333333}; w = -0.281250000000000000000000000000000;}
			else if(qp==1) {coords = {0.200000000000000000000000000000000, 0.200000000000000000000000000000000}; w = 0.260416666666666666666666666666665;}
			else if(qp==2) {coords = {0.200000000000000000000000000000000, 0.600000000000000000000000000000000}; w = 0.260416666666666666666666666666665;}
			else if(qp==3) {coords = {0.600000000000000000000000000000000, 0.200000000000000000000000000000000}; w = 0.260416666666666666666666666666665;}
			else err_message("Quadrature point for 3rd order quadrature of a Tri must be less than 4");
			break;
		}
		case 4:
		{
			if(qp==0) {coords = {0.44594849091597, 0.44594849091597}; w = 0.111690794839005;}
			else if(qp==1) {coords = {0.44594849091597, 0.10810301816807}; w = 0.111690794839005;}
			else if(qp==2) {coords = {0.10810301816807, 0.44594849091597}; w = 0.111690794839005;}
			else if(qp==3) {coords = {0.09157621350977, 0.09157621350977}; w = 0.054975871827660;}
			else if(qp==4) {coords = {0.09157621350977, 0.81684757298046}; w = 0.054975871827660;}
			else if(qp==5) {coords = {0.81684757298046, 0.09157621350977}; w = 0.054975871827660;}
			else err_message("Quadrature point for 4th order quadrature of a Tri must be less than 6");
				break;
		}
		case 5:
		{
			if(qp==0) {coords = {0.33333333333333, 0.33333333333333}; w = 0.112500000000000;}
			else if(qp==1) {coords = {0.47014206410511, 0.47014206410511}; w = 0.066197076394255;}
			else if(qp==2) {coords = {0.47014206410511, 0.05971587178977}; w = 0.066197076394255;}
			else if(qp==3) {coords = {0.05971587178977, 0.47014206410511}; w = 0.066197076394255;}
			else if(qp==4) {coords = {0.10128650732346, 0.10128650732346}; w = 0.062969590272415;}
			else if(qp==5) {coords = {0.10128650732346, 0.79742698535309}; w = 0.062969590272415;}
			else if(qp==6) {coords = {0.79742698535309, 0.10128650732346}; w = 0.062969590272415;}
			else err_message("Quadrature point for 5th order quadrature of a Tri must be less than 7");
			break;
		}
		case 6:
		{
			if(qp==0) {coords = {0.24928674517091, 0.24928674517091}; w = 0.11678627572638;}
			else if(qp==1) {coords = {0.24928674517091, 0.50142650965818}; w = 0.11678627572638;}
			else if(qp==2) {coords = {0.50142650965818, 0.24928674517091}; w = 0.11678627572638;}
			else if(qp==3) {coords = {0.06308901449150, 0.06308901449150}; w = 0.05084490637021;}
			else if(qp==4) {coords = {0.06308901449150, 0.87382197101700}; w = 0.05084490637021;}
			else if(qp==5) {coords = {0.87382197101700, 0.06308901449150}; w = 0.05084490637021;}
			else if(qp==6) {coords = {0.31035245103378, 0.63650249912140}; w = 0.08285107561837;}
			else if(qp==7) {coords = {0.63650249912140, 0.05314504984482}; w = 0.08285107561837;}
			else if(qp==8) {coords = {0.05314504984482, 0.31035245103378}; w = 0.08285107561837;}
			else if(qp==9) {coords = {0.63650249912140, 0.31035245103378}; w = 0.08285107561837;}
			else if(qp==10) {coords = {0.31035245103378, 0.05314504984482}; w = 0.08285107561837;}
			else if(qp==11) {coords = {0.05314504984482, 0.63650249912140}; w = 0.08285107561837;}
			else err_message("Quadrature point for 6th order quadrature of a Tri must be less than 12");
			break;
		}
		case 7:
		{
			if(qp==0) {coords = {0.33333333333333, 0.33333333333333}; w = -0.14957004446768;}
			else if(qp==1) {coords = {0.26034596607904, 0.26034596607904}; w = 0.17561525743321;}
			else if(qp==2) {coords = {0.26034596607904, 0.47930806784192}; w = 0.17561525743321;}
			else if(qp==3) {coords = {0.47930806784192, 0.26034596607904}; w = 0.17561525743321;}
			else if(qp==4) {coords = {0.06513010290222, 0.06513010290222}; w = 0.05334723560884;}
			else if(qp==5) {coords = {0.06513010290222, 0.86973979419557}; w = 0.05334723560884;}
			else if(qp==6) {coords = {0.86973979419557, 0.06513010290222}; w = 0.05334723560884;}
			else if(qp==7) {coords = {0.31286549600487, 0.63844418856981}; w = 0.07711376089026;}
			else if(qp==8) {coords = {0.63844418856981, 0.04869031542532}; w = 0.07711376089026;}
			else if(qp==9) {coords = {0.04869031542532, 0.31286549600487}; w = 0.07711376089026;}
			else if(qp==10) {coords = {0.63844418856981, 0.31286549600487}; w = 0.07711376089026;}
			else if(qp==11) {coords = {0.31286549600487, 0.04869031542532}; w = 0.07711376089026;}
			else if(qp==12) {coords = {0.04869031542532, 0.63844418856981}; w = 0.07711376089026;}
			else err_message("Quadrature point for 7th order quadrature of a Tri must be less than 13");
			break;
		}
		case 8:
		{
			if(qp==0) {coords = {0.33333333333333, 0.33333333333333}; w = 0.14431560767779;}
			else if(qp==1) {coords = {0.45929258829272, 0.45929258829272}; w = 0.09509163426728;}
			else if(qp==2) {coords = {0.45929258829272, 0.08141482341455}; w = 0.09509163426728;}
			else if(qp==3) {coords = {0.08141482341455, 0.45929258829272}; w = 0.09509163426728;}
			else if(qp==4) {coords = {0.17056930775176, 0.17056930775176}; w = 0.10321737053472;}
			else if(qp==5) {coords = {0.17056930775176, 0.65886138449648}; w = 0.10321737053472;}
			else if(qp==6) {coords = {0.65886138449648, 0.17056930775176}; w = 0.10321737053472;}
			else if(qp==7) {coords = {0.05054722831703, 0.05054722831703}; w = 0.03245849762320;}
			else if(qp==8) {coords = {0.05054722831703, 0.89890554336594}; w = 0.03245849762320;}
			else if(qp==9) {coords = {0.89890554336594, 0.05054722831703}; w = 0.03245849762320;}
			else if(qp==10) {coords = {0.26311282963464, 0.72849239295540}; w = 0.02723031417443;}
			else if(qp==11) {coords = {0.72849239295540, 0.00839477740996}; w = 0.02723031417443;}
			else if(qp==12) {coords = {0.00839477740996, 0.26311282963464}; w = 0.02723031417443;}
			else if(qp==13) {coords = {0.72849239295540, 0.26311282963464}; w = 0.02723031417443;}
			else if(qp==14) {coords = {0.26311282963464, 0.00839477740996}; w = 0.02723031417443;}
			else if(qp==15) {coords = {0.00839477740996, 0.72849239295540}; w = 0.02723031417443;}
			else err_message("Quadrature point for 8th order quadrature of a Tri must be less than 16");
			break;
		}
	}
}




double CohTri3::compute_shape(const std::vector<double>& coords, int n) const
{
	if(n<0 || n>3)
		err_message("Shape function number for an Tri3 must be between 0 and 3");
	double r = coords[0];
	double s = coords[1];
	if(n==0)
		return 1.0 - r - s;
	else if(n==1)
		return r;
	else if(n==2)
		return s;
	else
		err_message("Shape function selected must be less than 3 for a Tri3.");
}

std::vector<double> CohTri3::compute_shape_grad(const std::vector<double>& coords, int n) const
{
	if(n<0 || n>3)
		err_message("Shape function number for an Tri3 must be between 0 and 3");
	
	//double r = coords[0];
	//double s = coords[1];
	std::vector<double> v(2);
	if(n==0)
		v = {-1.0, -1.0};
	else if(n==1)
		v = {1.0, 0.0};
	else if(n==2)
		v = {0.0, 1.0};
	else
		err_message("Shape function selected must be less than 3 for a Tri3.");
	return v;
}


// void CohTri3::compute_rotation_general(const std::vector<double>& coords, const std::vector<std::vector<double> >& node_coords,
// 									   DenseMatrix<double>& mat, double& dA) const
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
// 				dxdr[i][j] += parent_grad[n][i]*node_coords[n][j];
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
