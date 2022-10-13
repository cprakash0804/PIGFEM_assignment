/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#include "tri3.h"
#include <iostream>

// ------------------------------------------------------------
const id_type Tri3::s_n_map[3][2] =
{
	{0, 1}, // Side 0
	{1, 2}, // Side 1
	{2, 0}  // Side 2
};

// Quad class static member initialization
const id_type Tri3::e_n_map[4][2] =
{
	{0, 1}, // Side 0
	{1, 2}, // Side 1
	{2, 0} // Side 2
};

id_type Tri3::side_nodes_map(id_type side, id_type node)
{
	if(side>=n_sides())
		err_message("Please select a valid side.");
	if(node>=2)
		err_message("Please select a valid node.");
	
	return s_n_map[side][node];
}

id_type Tri3::edge_nodes_map(id_type edge, id_type node)
{
	if(edge>=3)
		err_message("Please select a valid edge.");
	if(node>=2)
		err_message("Please select a valid node.");
	
	return e_n_map[edge][node];
}

Elem* Tri3::build_side(id_type side)
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


id_type Tri3::n_q_points(int order) const
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
void Tri3::q_point(std::vector<double>& coords, double& w, id_type qp, int order) const
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

std::vector<double> Tri3::compute_shape(const std::vector<double>& coords) const
{
	double r = coords[0];
	double s = coords[1];
	std::vector<double> N = {1.0 - r - s,
							 r,
							 s};
	return N;
}

std::vector<std::vector<double> > Tri3::compute_shape_grad(const std::vector<double>& coords) const
{
	std::vector<std::vector<double> > dN = {{-1.0, -1.0},
											{1.0, 0.0},
											{0.0, 1.0}};
	return dN;
}

std::vector<DenseMatrix<double> > Tri3::compute_shape_grad_grad(const std::vector<double>& coords) const
{
	DenseMatrix<double> ret(2,2);
	ret(0,0) = 0.0;	ret(0,1) = 0.0;
	ret(1,0) = 0.0;	ret(1,1) = 0.0;
	std::vector<DenseMatrix<double> > d2N(n_nodes());
	for (id_type n=0; n<n_nodes(); ++n)
		d2N[n] = ret;
	return d2N;
}


// Simpler inverse map for a triangle
// http://www4.ncsu.edu/~zhilin/TEACHING/MA587/chap9.pdf  (pg. 214)
std::vector<double> Tri3::inverse_map(std::vector<double> gcoords, bool& inside)
{
	std::vector<std::vector<double> > c = {{(*_nodes[0])(0), (*_nodes[0])(1)},
										   {(*_nodes[1])(0), (*_nodes[1])(1)},
										   {(*_nodes[2])(0), (*_nodes[2])(1)}};
	double A = 0.5*(c[0][0]*(c[1][1]-c[2][1]) + c[1][0]*(c[2][1]-c[0][1]) + c[2][0]*(c[0][1]-c[1][1]));
	if(A < 0.0)
		A = -1.0*A;
	
	std::vector<double> sol(2);
	sol[0] = (1.0/(2.0*A))*((c[2][1]-c[0][1])*(gcoords[0]-c[0][0]) - (c[2][0]-c[0][0])*(gcoords[1]-c[0][1]));
	sol[1] = (1.0/(2.0*A))*(-1.0*(c[1][1]-c[0][1])*(gcoords[0]-c[0][0]) + (c[1][0]-c[0][0])*(gcoords[1]-c[0][1]));

	inside = coords_inside(sol);
	return sol;
}


bool Tri3::coords_inside(const std::vector<double>& rcoords)
{
	if (rcoords.size() != 2)
		err_message("Invalid coordinate size for coordinate inside check (TRI3)!");

	if (rcoords[0] >= (0.0-_inside_tol) && rcoords[0] <= (1.0+_inside_tol))
	{
		if (rcoords[1] >= (0.0-_inside_tol) && rcoords[1] <= (1.0+_inside_tol))
		{
			if ((rcoords[0] + rcoords[1]) <= (1.0+_inside_tol))
				return true;
			else
				return false;
		}
		else
			return false;
	}
	else
		return false;
}




// IGFEM MATERIAL BELOW HERE
//==============================================================================================================================

void Tri3::detection(std::vector<int>& node_detection, std::vector<Material*>& mats, bool isCohesive)
{
	if(node_detection.size() != n_nodes())
		err_message("The number of nodes used for element detection must be the same as the number of nodes in the element.");

	// NOTE: the mats structure is a pointer to the material for every inclusion in the mesh, with the 0th elment being a pointer to the base/matrix material
	// For this reason, when we index the mats array with the nodal detection values (incusion indices) we need to add 1
	Material* matrix = mats[0];
	
	// Store many booleans that will be used many times
	int n0 = node_detection[0];
	int n1 = node_detection[1];
	int n2 = node_detection[2];

	//Storing whether the nodes belong to the same interface/inclusion as other nodes
	bool n0e1 = (n0 == n1); //bool n1e0 = n0e1;
	bool n0e2 = (n0 == n2); //bool n2e0 = n0e2;
	bool n1e2 = (n1 == n2); //bool n2e1 = n1e2;

	std::vector<int> vec = node_detection;
	std::sort(vec.begin(), vec.end());
	std::vector<int>::iterator it;
	it = std::unique(vec.begin(), vec.end());
	vec.resize( std::distance(vec.begin(), it) );
	int n_unique_nodes = vec.size();

	std::vector<int> set_counts;
	for(int i=0; i<n_unique_nodes; ++i)
		set_counts.push_back(std::count(node_detection.begin(), node_detection.end(), vec[i]));

	
	
	
	//All of the nodes are located in the same material, either the base material or an inclusion
	if(n_unique_nodes == 1)
	{
		_is_intersected = false;
		_cutEdge.clear();
		_enrichment_on_inclusion.clear(); // Should all be the same
		_int_elem_struct.clear(); // Hopefully this syntax works
		_coh_elem_struct.clear();
		_int_elem_mat.clear();
	}
	
	
	

	// Either one interface and the matrix or two interfaces
	// NOTE: We assume here that interfaces do not intersect here (they are separated by a region of matrix)
	else
	{
		_is_intersected = true;

		// All cases that are not cohesive
		if (!isCohesive)
		{
			if (n_unique_nodes == 2)
			{
				if (vec[0] < 0) // Only intersected by one inclusion
				{
					if(n1e2)		// node zero is separate
						{_cutEdge = {0, 2}; _int_elem_struct = {{0, 3, 4}, {1, 2, 4, 3}}; _int_elem_mat = {mats[n0+1], mats[n1+1]};}
					else if(n0e2)	// node one is separate
						{_cutEdge = {0, 1}; _int_elem_struct = {{1, 4, 3}, {0, 3, 4, 2}}; _int_elem_mat = {mats[n1+1], mats[n2+1]};}
					else if(n0e1)			// node two is separate
						{_cutEdge = {1, 2}; _int_elem_struct = {{2, 4, 3}, {0, 1, 3, 4}}; _int_elem_mat = {mats[n2+1], mats[n0+1]};}
					else
						err_message("TRIANGULAR ELEMENT HAS 1 NODE SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");

					_enrichment_on_inclusion.resize(_cutEdge.size());
					std::fill(_enrichment_on_inclusion.begin(), _enrichment_on_inclusion.end(), *std::max_element(vec.begin(), vec.end())); // All intersected by the same inclusion
				}

				// Intersected by 2 inclusions
				else
				{
					if(n1e2)		// node zero is separate
						{_cutEdge = {2,0,2,0}; _enrichment_on_inclusion = {n0,n0,n1,n1}; _int_elem_struct = {{0,4,3},{4,6,5,3},{1,2,5,6}}; _int_elem_mat = {mats[n0+1],matrix,mats[n1+1]};}
					else if(n0e2)	// node one is separate
						{_cutEdge = {0,1,0,1}; _enrichment_on_inclusion = {n1,n1,n2,n2}; _int_elem_struct = {{1,4,3},{4,6,5,3},{2,0,5,6}}; _int_elem_mat = {mats[n1+1],matrix,mats[n2+1]};}
					else if(n0e1)	// node two is separate
						{_cutEdge = {1,2,1,2}; _enrichment_on_inclusion = {n2,n2,n0,n0}; _int_elem_struct = {{2,4,3},{4,6,5,3},{0,1,5,6}}; _int_elem_mat = {mats[n2+1],matrix,mats[n0+1]};}
					else
						err_message("TRIANGULAR ELEMENT HAS 1 NODE SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
				}
			}


			// 3 unique nodal detection values
			else
			{
				// One of the nodes is in the matrix
				if (vec[0] < 0)
				{
					if (n0 < 0) // Node 0 is in the matrix
					{
						_cutEdge = {1, 0, 2, 1};
						_enrichment_on_inclusion = {n1,n1,n2,n2};
						_int_elem_mat = {mats[n1+1], matrix, mats[n2+1]};
						_int_elem_struct = {{1,3,4}, {0,4,3,6,5} ,{2,5,6}};
					}
					else if (n1 < 0) // Node 1 is in the matrix
					{
						_cutEdge = {2, 1, 0, 2};
						_enrichment_on_inclusion = {n2,n2,n0,n0};
						_int_elem_mat = {mats[n2+1], matrix, mats[n0+1]};
						_int_elem_struct = {{2,3,4}, {1,4,3,6,5} ,{0,5,6}};
					}
					else if (n2 < 0) // Node 2 is in the matrix
					{
						_cutEdge = {0, 2, 1, 0};
						_enrichment_on_inclusion = {n0,n0,n1,n1};
						_int_elem_mat = {mats[n0+1], matrix, mats[n1+1]};
						_int_elem_struct = {{0,3,4}, {2,4,3,6,5} ,{1,5,6}};
					}
					else
						err_message("Unknown manner of 2 inclusions cutting a Tri3.");
				}

				// All 3 nodes are in different inclusions
				else
				{
					_cutEdge = {0, 2, 2, 1, 1, 0};
					_enrichment_on_inclusion = {n0,n0,n2,n2,n1,n1};
					_int_elem_mat = {mats[n0+1], mats[n2+1], mats[n1+1], matrix};
					_int_elem_struct = {{0,3,4}, {2,5,6}, {1,7,8}, {8,7,6,5,4,3}};
				}
			}

			// Not cohesive So I can clear this
			_enrichment_nodes.resize(_cutEdge.size());
			_coh_elem_struct.clear();
		}




		// Cohesive
		else
		{
			if (n_unique_nodes == 2)
			{
				// At least one of the nodes is matrix
				if (vec[0] < 0)
				{
					if(n1e2)		// node zero is separate
					{
						_cutEdge = {2,0};
						_int_elem_mat = {mats[n0+1], mats[n1+1]};
						if(n0 < n1) {_int_elem_struct = {{0, 4, 3}, {1, 2, 5,6}}; _coh_elem_struct = {{3,4,5,6}};}
						else {_int_elem_struct = {{0, 6, 5}, {1, 2, 3, 4}}; _coh_elem_struct = {{4,3,6,5}};}
					}
					else if(n0e2)	// node one is separate
					{
						_cutEdge = {0, 1};
						_int_elem_mat = {mats[n1+1], mats[n0+1]};
						if(n1 < n0) {_int_elem_struct = {{1, 4, 3}, {2, 0, 5, 6}}; _coh_elem_struct = {{3,4,5,6}};}
						else {_int_elem_struct = {{1, 6, 5}, {2, 0, 3, 4}};  _coh_elem_struct = {{4,3,6,5}};}
					}
					else if(n0e1)			// node two is separate
					{
						_cutEdge = {1, 2};
						_int_elem_mat = {mats[n2+1], mats[n0+1]};
						if(n2 < n0) {_int_elem_struct = {{2, 4, 3}, {0, 1, 5, 6}}; _coh_elem_struct = {{3,4,5,6}};}
						else {_int_elem_struct = {{2, 6, 5}, {0, 1, 3, 4}}; _coh_elem_struct = {{4,3,6,5}};}
					}
					else
						err_message("TRIANGULAR ELEMENT HAS 1 NODE SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");

					_enrichment_on_inclusion.resize(_cutEdge.size());
					std::fill(_enrichment_on_inclusion.begin(), _enrichment_on_inclusion.end(), *std::max_element(vec.begin(), vec.end())); // All intersected by the same inclusion
				}

				// Intersected by 2 inclusion
				else
				{
					if(n1e2)		// node zero is separate
					{
						_cutEdge = {2,0,2,0};
						_enrichment_on_inclusion = {n0,n0,n1,n1};
						_int_elem_mat = {mats[n0+1], matrix, mats[n1+1]};
						_int_elem_struct = {{0,8,7}, {4,6,5,3}, {1,2,9,10}};
					}
					else if(n0e2)	// node one is separate
					{
						_cutEdge = {0,1,0,1};
						_enrichment_on_inclusion = {n1,n1,n2,n2};
						_int_elem_mat = {mats[n1+1], matrix, mats[n2+1]};
						_int_elem_struct = {{1,8,7}, {4,6,5,3},{2,0,9,10}};
					}
					else if(n0e1)			// node two is separate
					{
						_cutEdge = {1,2,1,2};
						_enrichment_on_inclusion = {n2,n2,n0,n0};
						_int_elem_mat = {mats[n2+1], matrix, mats[n0+1]};
						_int_elem_struct = {{2,8,7}, {4,6,5,3},{0,1,9,10}};
					}
					else
						err_message("TRIANGULAR ELEMENT HAS 1 NODE SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");

					// All cohesive structs are the same
					_coh_elem_struct = {{4,3,8,7}, {5,6,9,10}};
				}
			} // End 2 unique nodes



			// 3 unique nodal detection values
			else
			{
				// One of the nodes is in the matrix
				if (vec[0] < 0)
				{
					if (n0 < 0) // Node 0 is in the matrix
					{
						_cutEdge = {1, 0, 2, 1};
						_enrichment_on_inclusion = {n1,n1,n2,n2};
						_int_elem_mat = {mats[n1+1], matrix, mats[n2+1]};
						_int_elem_struct = {{1,7,8}, {0,4,3,6,5} ,{2,9,10}};
					}
					else if (n1 < 0) // Node 1 is in the matrix
					{
						_cutEdge = {2, 1, 0, 2};
						_enrichment_on_inclusion = {n2,n2,n0,n0};
						_int_elem_mat = {mats[n2+1], matrix, mats[n0+1]};
						_int_elem_struct = {{2,7,8}, {1,4,3,6,5} ,{0,9,10}};
					}
					else if (n2 < 0) // Node 2 is in the matrix
					{
						_cutEdge = {0, 2, 1, 0};
						_enrichment_on_inclusion = {n0,n0,n1,n1};
						_int_elem_mat = {mats[n0+1], matrix, mats[n1+1]};
						_int_elem_struct = {{0,7,8}, {2,4,3,6,5} ,{1,9,10}};
					}
					else
						err_message("Unknown manner of 2 inclusions cutting a Tri3.");

					// All the same
					_coh_elem_struct = {{3,4,7,8}, {5,6,9,10}};
				}

				// All 3 nodes are in different inclusions
				else
				{
					_cutEdge = {0, 2, 2, 1, 1, 0};
					_enrichment_on_inclusion = {n0,n0,n2,n2,n1,n1};
					_int_elem_mat = {mats[n0+1], mats[n2+1], mats[n1+1], matrix};
					_int_elem_struct = {{0,9,10}, {2,11,12}, {1,13,14}, {8,7,6,5,4,3}};
					_coh_elem_struct = {{3,4,9,10}, {5,6,11,12}, {7,8,13,14}};
				}
			} // End 3 unique nodes

			_enrichment_nodes.resize(_cutEdge.size() * 2);
			_coh_mat.resize(_coh_elem_struct.size());
		} // End is cohesive
	} // End is intersected

	
/*	
	// Either one interface and the matrix or two interfaces
	// NOTE: We assume here that interfaces do not intersect here (they are separated by a region of matrix)
	else
	{
		_is_intersected = true;

		// All cases that are not cohesive
		if (!isCohesive)
		{
			if (n_unique_nodes == 2)
			{
				if (vec[0] < 0) // Only intersected by one inclusion
				{
					if(n1e2)		// node zero is separate
						{_cutEdge = {0, 2}; _int_elem_struct = {{0, 3, 4}, {1, 2, 4, 3}}; _int_elem_mat = {mats[n0+1], mats[n1+1]};}
					else if(n0e2)	// node one is separate
						{_cutEdge = {0, 1}; _int_elem_struct = {{1, 4, 3}, {0, 3, 4, 2}}; _int_elem_mat = {mats[n1+1], mats[n0+1]};}
					else if(n0e1)			// node two is separate
						{_cutEdge = {1, 2}; _int_elem_struct = {{2, 4, 3}, {0, 1, 3, 4}}; _int_elem_mat = {mats[n2+1], mats[n0+1]};}
					else
						err_message("TRIANGULAR ELEMENT HAS 1 NODE SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");

					_enrichment_on_inclusion.resize(_cutEdge.size());
					std::fill(_enrichment_on_inclusion.begin(), _enrichment_on_inclusion.end(), *std::max_element(vec.begin(), vec.end())); // All intersected by the same inclusion
				}

				// Intersected by 2 inclusions
				else
				{
					if(n1e2)		// node zero is separate
						{_cutEdge = {2,0,2,0}; _enrichment_on_inclusion = {n0,n0,n1,n1}; _int_elem_struct = {{0,4,3},{4,6,5,3},{1,2,5,6}}; _int_elem_mat = {mats[n0+1],matrix,mats[n1+1]};}
					else if(n0e2)	// node one is separate
						{_cutEdge = {0,1,0,1}; _enrichment_on_inclusion = {n1,n1,n2,n2}; _int_elem_struct = {{1,4,3},{4,6,5,3},{2,0,5,6}}; _int_elem_mat = {mats[n1+1],matrix,mats[n2+1]};}
					else if(n0e1)	// node two is separate
						{_cutEdge = {1,2,1,2}; _enrichment_on_inclusion = {n2,n2,n0,n0}; _int_elem_struct = {{2,4,3},{4,6,5,3},{0,1,5,6}}; _int_elem_mat = {mats[n2+1],matrix,mats[n0+1]};}
					else
						err_message("TRIANGULAR ELEMENT HAS 1 NODE SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
				}
			}


			// 3 unique nodal detection values
			else
			{
				// One of the nodes is in the matrix
				if (vec[0] < 0)
				{
					if (n0 < 0) // Node 0 is in the matrix
					{
						_cutEdge = {1, 0, 2, 1};
						_enrichment_on_inclusion = {n1,n1,n2,n2};
						_int_elem_mat = {mats[n1+1], matrix, matrix, matrix, mats[n2+1]};
						_int_elem_struct = {{1,3,4}, {0,4,3}, {0,3,6}, {0,6,5} ,{2,5,6}};
					}
					else if (n1 < 0) // Node 1 is in the matrix
					{
						_cutEdge = {2, 1, 0, 2};
						_enrichment_on_inclusion = {n2,n2,n0,n0};
						_int_elem_mat = {mats[n2+1], matrix, matrix, matrix, mats[n0+1]};
						_int_elem_struct = {{2,3,4}, {1,4,3}, {1,3,6}, {1,6,5} ,{0,5,6}};
					}
					else if (n2 < 0) // Node 2 is in the matrix
					{
						_cutEdge = {0, 2, 1, 0};
						_enrichment_on_inclusion = {n0,n0,n1,n1};
						_int_elem_mat = {mats[n0+1], matrix, matrix, matrix, mats[n1+1]};
						_int_elem_struct = {{0,3,4}, {2,4,3}, {2,3,6}, {2,6,5} ,{1,5,6}};
					}
					else
						err_message("Unknown manner of 2 inclusions cutting a Tri3.");
				}

				// All 3 nodes are in different inclusions
				else
				{
					_cutEdge = {0, 2, 2, 1, 1, 0};
					_enrichment_on_inclusion = {n0,n0,n2,n2,n1,n1};
					_int_elem_mat = {mats[n0+1], mats[n2+1], mats[n1+1], matrix, matrix, matrix};
					_int_elem_struct = {{0,3,4}, {2,5,6}, {1,7,8}, {4,6,5}, {3,8,7}, {4,3,7,6}};
				}
			}

			// Not cohesive So I can clear this
			_enrichment_nodes.resize(_cutEdge.size());
			_coh_elem_struct.clear();
		}




		// Cohesive
		else
		{
			if (n_unique_nodes == 2)
			{
				// At least one of the nodes is matrix
				if (vec[0] < 0)
				{
					if(n1e2)		// node zero is separate
					{
						_cutEdge = {2,0};
						_int_elem_mat = {mats[n0+1], mats[n1+1]};
						if(n0 < n1) {_int_elem_struct = {{0, 4, 3}, {1, 2, 5,6}}; _coh_elem_struct = {{3,4,5,6}};}
						else {_int_elem_struct = {{0, 6, 5}, {1, 2, 3, 4}}; _coh_elem_struct = {{4,3,6,5}};}
					}
					else if(n0e2)	// node one is separate
					{
						_cutEdge = {0, 1};
						_int_elem_mat = {mats[n1+1], mats[n0+1]};
						if(n1 < n0) {_int_elem_struct = {{1, 4, 3}, {2, 0, 5, 6}}; _coh_elem_struct = {{3,4,5,6}};}
						else {_int_elem_struct = {{1, 6, 5}, {2, 0, 3, 4}};  _coh_elem_struct = {{4,3,6,5}};}
					}
					else if(n0e1)			// node two is separate
					{
						_cutEdge = {1, 2};
						_int_elem_mat = {mats[n2+1], mats[n0+1]};
						if(n2 < n0) {_int_elem_struct = {{2, 4, 3}, {0, 1, 5, 6}}; _coh_elem_struct = {{3,4,5,6}};}
						else {_int_elem_struct = {{2, 6, 5}, {0, 1, 3, 4}}; _coh_elem_struct = {{4,3,6,5}};}
					}
					else
						err_message("TRIANGULAR ELEMENT HAS 1 NODE SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");

					_enrichment_on_inclusion.resize(_cutEdge.size());
					std::fill(_enrichment_on_inclusion.begin(), _enrichment_on_inclusion.end(), *std::max_element(vec.begin(), vec.end())); // All intersected by the same inclusion
				}

				// Intersected by 2 inclusion
				else
				{
					if(n1e2)		// node zero is separate
					{
						_cutEdge = {2,0,2,0};
						_enrichment_on_inclusion = {n0,n0,n1,n1};
						_int_elem_mat = {mats[n0+1], matrix, mats[n1+1]};
						_int_elem_struct = {{0,8,7}, {4,6,5,3}, {1,2,9,10}};
					}
					else if(n0e2)	// node one is separate
					{
						_cutEdge = {0,1,0,1};
						_enrichment_on_inclusion = {n1,n1,n2,n2};
						_int_elem_mat = {mats[n1+1], matrix, mats[n2+1]};
						_int_elem_struct = {{1,8,7}, {4,6,5,3},{2,0,9,10}};
					}
					else if(n0e1)			// node two is separate
					{
						_cutEdge = {1,2,1,2};
						_enrichment_on_inclusion = {n2,n2,n0,n0};
						_int_elem_mat = {mats[n2+1], matrix, mats[n0+1]};
						_int_elem_struct = {{2,8,7}, {4,6,5,3},{0,1,9,10}};
					}
					else
						err_message("TRIANGULAR ELEMENT HAS 1 NODE SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");

					// All cohesive structs are the same
					_coh_elem_struct = {{4,3,8,7}, {5,6,9,10}};
				}
			} // End 2 unique nodes



			// 3 unique nodal detection values
			else
			{
				// One of the nodes is in the matrix
				if (vec[0] < 0)
				{
					if (n0 < 0) // Node 0 is in the matrix
					{
						_cutEdge = {1, 0, 2, 1};
						_enrichment_on_inclusion = {n1,n1,n2,n2};
						_int_elem_mat = {mats[n1+1], matrix, matrix, matrix, mats[n2+1]};
						_int_elem_struct = {{1,7,8}, {0,4,3}, {0,3,6}, {0,6,5} ,{2,9,10}};
					}
					else if (n1 < 0) // Node 1 is in the matrix
					{
						_cutEdge = {2, 1, 0, 2};
						_enrichment_on_inclusion = {n2,n2,n0,n0};
						_int_elem_mat = {mats[n2+1], matrix, matrix, matrix, mats[n0+1]};
						_int_elem_struct = {{2,7,8}, {1,4,3}, {1,3,6}, {1,6,5} ,{0,9,10}};
					}
					else if (n2 < 0) // Node 2 is in the matrix
					{
						_cutEdge = {0, 2, 1, 0};
						_enrichment_on_inclusion = {n0,n0,n1,n1};
						_int_elem_mat = {mats[n0+1], matrix, matrix, matrix, mats[n1+1]};
						_int_elem_struct = {{0,7,8}, {2,4,3}, {2,3,6}, {2,6,5} ,{1,9,10}};
					}
					else
						err_message("Unknown manner of 2 inclusions cutting a Tri3.");

					// All the same
					_coh_elem_struct = {{3,4,7,8}, {5,6,9,10}};
				}

				// All 3 nodes are in different inclusions
				else
				{
					_cutEdge = {0, 2, 2, 1, 1, 0};
					_enrichment_on_inclusion = {n0,n0,n2,n2,n1,n1};
					_int_elem_mat = {mats[n0+1], mats[n2+1], mats[n1+1], matrix, matrix, matrix};
					_int_elem_struct = {{0,9,10}, {2,11,12}, {1,13,14}, {4,6,5}, {3,8,7}, {4,3,7,6}};
					_coh_elem_struct = {{3,4,9,10}, {5,6,11,12}, {7,8,13,14}};
				}
			} // End 3 unique nodes

			_enrichment_nodes.resize(_cutEdge.size() * 2);
			_coh_mat.resize(_coh_elem_struct.size());
		} // End is cohesive
	} // End is intersected
*/

	// Set the boolean
	_detected = true;
} // detection








void Tri3::refinement_nodes(std::vector<std::vector<id_type> >& refine_nodes)
{
	std::vector<id_type> ids(n_nodes());
	for (id_type n=0; n<n_nodes(); ++n)
		ids[n] = _nodes[n]->get_id();
	refine_nodes = {{ids[0], ids[1]},
					{ids[1], ids[2]},
					{ids[2], ids[0]}};
}
void Tri3::refinement_structure(std::vector<std::vector<id_type> >& structure, std::vector<elem_type>& types)
{
	structure = {{0, 3, 5},
				 {3, 1, 4},
				 {5, 4, 2},
				 {3, 4, 5}};
	types = {TRI3, TRI3, TRI3, TRI3};
}
