/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#include "edge2.h"

const id_type Edge2::s_n_map[2][1] =
{
	{0}, // Side 0
	{1} // Side 1
};

const id_type Edge2::e_n_map[1][2] =
{
	{0, 1} // Side 0
};

id_type Edge2::side_nodes_map(id_type side, id_type node)
{
	if(side>=n_sides())
		err_message("Please select a valid side.");
	if(node>=1)
		err_message("Please select a valid node.");
	
	return s_n_map[side][node];
}

id_type Edge2::edge_nodes_map(id_type edge, id_type node)
{
	if(edge>=1)
		err_message("Please select a valid edge.");
	if(node>=2)
		err_message("Please select a valid node.");
	
	return e_n_map[edge][node];
}


Elem* Edge2::build_side(id_type side)
{
	if(side >= n_sides())
		err_message("Side number must be less than the number of sides.");

	Elem* el = build(POINT1);
	std::vector<Node*> nodes;
	for(id_type n=0; n<1; ++n)
		nodes.push_back(_nodes[side_nodes_map(side,n)]);
	el->set_nodes(nodes);
	return el;
}



id_type Edge2::n_q_points(int order) const
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
void Edge2::q_point(std::vector<double>& coords, double& w, id_type qp, int order) const
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

std::vector<double> Edge2::compute_shape(const std::vector<double>& coords) const
{
	double r = coords[0];
	std::vector<double> N = {0.5*(1.0-r),
							 0.5*(1.0+r)};
	return N;
}

std::vector<std::vector<double> > Edge2::compute_shape_grad(const std::vector<double>& coords) const
{
	std::vector<std::vector<double> > dN = {{-0.5},
											{0.5}};
	return dN;
}

std::vector<DenseMatrix<double> > Edge2::compute_shape_grad_grad(const std::vector<double>& coords) const
{
	DenseMatrix<double> ret(1,1);
	ret(0,0) = 0.0;
	std::vector<DenseMatrix<double> > d2N(n_nodes());
	for (id_type n=0; n<n_nodes(); ++n)
		d2N[n] = ret;
	return d2N;
}

bool Edge2::coords_inside(const std::vector<double>& rcoords)
{
	if (rcoords.size() != 1)
		err_message("Invalid coordinate size for coordinate inside check (EDGE2)!");

	if (rcoords[0] >= (-1.0-_inside_tol) && rcoords[0] <= (1.0+_inside_tol))
		return true;
	else
		return false;
}







// IGFEM MATERIAL BELOW HERE
//==============================================================================================================================

void Edge2::detection(std::vector<int>& node_detection, std::vector<Material*>& mats, bool isCohesive)
{
	if((node_detection.size()!=n_nodes()))
		err_message("The number of nodes used for element detection must be the same as the number of nodes in the element.");

	// NOTE: the mats structure is a pointer to the material for every inclusion in the mesh, with the 0th elment being a pointer to the base/matrix material
	// For this reason, when we index the mats array with the nodal detection values (incusion indices) we need to add 1
	Material* matrix = mats[0];

	// Store many booleans that will be used many times
	int n0 = node_detection[0];
	int n1 = node_detection[1];

	//Storing whether the nodes belong to the same interface/inclusion as other nodes
	//bool n0e1 = (n0 == n1); //bool n1e0 = n0e1;

	std::vector<int> vec = {n0, n1};
	std::sort(vec.begin(), vec.end());
	std::vector<int>::iterator it;
	it = std::unique(vec.begin(), vec.end());
	vec.resize( std::distance(vec.begin(), it) );
	int n_unique_nodes = vec.size();

	//All of the nodes are located in the same material, either the base material or an inclusion
	if(n_unique_nodes == 1)
	{
		_is_intersected = false;
		_cutEdge.clear();
		_enrichment_on_inclusion.clear();
		_int_elem_struct.clear();
		_coh_elem_struct.clear();
		_int_elem_mat.clear();
	}


	// 2 different nodal detection values
	// Could be 1 inclusion and the matrix or it could be two different inclusions
	// NOTE: We assume here that interfaces do not intersect here (they are separated by a region of matrix)
	else
	{
		_is_intersected = true;

		// Not cohesive
		if(!isCohesive)
		{
			// One of the nodes is in the matrix
			if (vec[0] < 0)
			{
				_cutEdge = {0};
				_enrichment_on_inclusion = {*std::max_element(vec.begin(), vec.end())};
				_int_elem_mat = {mats[n0+1], mats[n1+1]};
				_int_elem_struct = {{0,2}, {2,1}};
			}
			
			// Two different inclusions
			else
			{
				_cutEdge = {0, 0};
				_enrichment_on_inclusion = {n0, n1};
				_int_elem_mat = {mats[n0+1], matrix, mats[n1+1]};
				_int_elem_struct = {{0,2}, {2,3}, {3,1}};
			}

			_enrichment_nodes.resize(_cutEdge.size());
			_coh_elem_struct.clear();
		}

		// Cohesive
		else
		{
			// One of the nodes is in the matrix
			if (vec[0] < 0)
			{
				_cutEdge = {0};
				_enrichment_on_inclusion = {*std::max_element(vec.begin(), vec.end())};
				_int_elem_mat = {mats[n0+1], mats[n1+1]};
				if(n0 < n1)
					_int_elem_struct = {{0,2}, {3,1}};
				else
					_int_elem_struct = {{0,3}, {2,1}};

				_coh_elem_struct = {{2, 3}};
			}

			// Nodes lie in different inclusions
			else
			{
				_cutEdge = {0, 0};
				_enrichment_on_inclusion = {n0, n1};
				_int_elem_mat = {mats[n0+1], matrix, mats[n1+1]};
				_int_elem_struct = {{0,4}, {2,3}, {5,1}};
				_coh_elem_struct = {{2,4}, {3,5}};
			}

			_enrichment_nodes.resize(_cutEdge.size() * 2);
			_coh_mat.resize(_coh_elem_struct.size());
		}
	}



	_detected = true;

} // detection




void Edge2::refinement_nodes(std::vector<std::vector<id_type> >& refine_nodes)
{
	std::vector<id_type> ids(n_nodes());
	for (id_type n=0; n<n_nodes(); ++n)
		ids[n] = _nodes[n]->get_id();
	refine_nodes = {{ids[0], ids[1]}};
}
void Edge2::refinement_structure(std::vector<std::vector<id_type> >& structure, std::vector<elem_type>& types)
{
	structure = {{0, 2}, {2, 1}};
	types = {EDGE2, EDGE2};
}
