/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#include "prism6.h"
#include "tri3.h"
#include "edge2.h"


// ------------------------------------------------------------
// Prism6 class static member initializations
const int Prism6::s_n_map[5][4] =
{
	{0, 1, 2, -1},		// Side 0
	{0, 1, 4, 3},	// Side 1
	{1, 2, 5, 4},	// Side 2
	{2, 0, 3, 5},	// Side 3
	{3, 4, 5, -1}		// Side 4
};

const id_type Prism6::e_n_map[9][2] =
{
	{0, 1}, // Edge 0
	{1, 2}, // Edge 1
	{2, 0}, // Edge 2
	{0, 3}, // Edge 3
	{1, 4}, // Edge 4
	{2, 5}, // Edge 5
	{3, 4}, // Edge 6
	{4, 5}, // Edge 7
	{5, 3}	// Edge 8
};

id_type Prism6::side_nodes_map(id_type side, id_type node)
{
	if(side>=n_sides())
		err_message("Please select a valid side.");
	if(node>=4)
		err_message("Please select a valid node.");
	
	return s_n_map[side][node];
}

id_type Prism6::edge_nodes_map(id_type edge, id_type node)
{
	if(edge>=n_edges())
		err_message("Please select a valid edge.");
	if(node>=2)
		err_message("Please select a valid node.");
	
	return e_n_map[edge][node];
}

id_type Prism6::n_nodes_on_side(id_type side)
{
	if(side==0 || side==4)
		return 3;
	else
		return 4;
}

Elem* Prism6::build_side(id_type side)
{
	if(side >= n_sides())
		err_message("Side number must be less than the number of sides.");

	Elem* el;
	std::vector<Node*> nodes;
	switch(side)
	{
		case 1:
		case 2:
		case 3:
			el = build(QUAD4);
			for(id_type n=0; n<4; ++n)
				nodes.push_back(_nodes[side_nodes_map(side,n)]);
			el->set_nodes(nodes);
			return el;
		case 0:
		case 4:
			el = build(TRI3);
			for(id_type n=0; n<3; ++n)
				nodes.push_back(_nodes[side_nodes_map(side,n)]);
			el->set_nodes(nodes);
			return el;
		default:
			err_message("Invalid side (Prism6)");
	}
}

// Quadrature rules for the prism are a composite of the rules for a triangle and for an edge
id_type Prism6::n_q_points(int order) const
{
	if(order < 0)
		err_message("Order of integration must be greater than or equal to 0.");
	
	Elem* tri = Elem::build(TRI3);
	Elem* edge = Elem::build(EDGE2);
	int nqp_tri = tri->n_q_points(order);
	int nqp_edge = edge->n_q_points(order);
	delete tri;
	delete edge;
	return nqp_tri*nqp_edge;
}

void Prism6::q_point(std::vector<double>& coords, double& w, id_type qp, int order) const
{
	if(order < 0)
		err_message("Order of integration must be greater than or equal to 0.");
	if(qp<0 || qp>=n_q_points(order))
		err_message("Quadrature point selected must be less than the maximum number of quadrature points for the given order.");
	
	coords.clear();
	Elem* tri = Elem::build(TRI3);
	Elem* edge = Elem::build(EDGE2);
	id_type nqp_tri = tri->n_q_points(order);
	
	std::vector<double> tri_qp, edge_qp;
	double tri_w, edge_w;
	tri->q_point(tri_qp, tri_w, (qp%nqp_tri), order);
	edge->q_point(edge_qp, edge_w, (qp/nqp_tri), order);
	delete tri;
	delete edge;

	coords = {tri_qp[0], tri_qp[1], edge_qp[0]};
	w = tri_w*edge_w;
}


std::vector<double> Prism6::compute_shape(const std::vector<double>& coords) const
{
	double r = coords[0];
	double s = coords[1];
	double t = coords[2];
	std::vector<double> N = {(1.-r-s) * .5*(1.-t),
							 r * 0.5*(1.-t),
							 s * 0.5*(1.-t),
							 (1.-r-s) * 0.5*(1.+t),
							 r * 0.5*(1.+t),
							 s * 0.5*(1.+t)};
	return N;
}

std::vector<std::vector<double> > Prism6::compute_shape_grad(const std::vector<double>& coords) const
{
	double r = coords[0];
	double s = coords[1];
	double t = coords[2];
	std::vector<std::vector<double> > dN = {{-0.5*(1.0-t), -0.5*(1.0-t), -0.5*(1.0-r-s)},
											{0.5*(1.0-t), 0.0, -0.5*r},
											{0.0, 0.5*(1.0-t), -0.5*s},
											{-0.5*(1.0+t), -0.5*(1.0+t), 0.5*(1.0-r-s)},
											{0.5*(1.0+t), 0.0, 0.5*r},
											{0.0, 0.5*(1.0+t), 0.5*s}};
	return dN;
}

std::vector<DenseMatrix<double> > Prism6::compute_shape_grad_grad(const std::vector<double>& coords) const
{
	std::vector<DenseMatrix<double> > d2N(n_nodes());
	// Node 0
		d2N[0](0,0) = 0.0;		d2N[0](0,1) = 0.0;		d2N[0](0,2) = 0.5;
		d2N[0](1,0) = 0.0;		d2N[0](1,1) = 0.0;		d2N[0](1,2) = 0.5;
		d2N[0](2,0) = 0.5;		d2N[0](2,1) = 0.5;		d2N[0](2,2) = 0.0;
	// Node 1
		d2N[1](0,0) = 0.0;		d2N[1](0,1) = 0.0;		d2N[1](0,2) = -0.5;
		d2N[1](1,0) = 0.0;		d2N[1](1,1) = 0.0;		d2N[1](1,2) = 0.0;
		d2N[1](2,0) = -0.5;		d2N[1](2,1) = 0.0;		d2N[1](2,2) = 0.0;
	// Node 2
		d2N[2](0,0) = 0.0;		d2N[2](0,1) = 0.0;		d2N[2](0,2) = 0.0;
		d2N[2](1,0) = 0.0;		d2N[2](1,1) = 0.0;		d2N[2](1,2) = -0.5;
		d2N[2](2,0) = 0.0;		d2N[2](2,1) = -0.5;		d2N[2](2,2) = 0.0;
	// Node 3
		d2N[3](0,0) = 0.0;		d2N[3](0,1) = 0.0;		d2N[3](0,2) = -0.5;
		d2N[3](1,0) = 0.0;		d2N[3](1,1) = 0.0;		d2N[3](1,2) = -0.5;
		d2N[3](2,0) = -0.5;		d2N[3](2,1) = -0.5;		d2N[3](2,2) = 0.0;
	// Node 4
		d2N[4](0,0) = 0.0;		d2N[4](0,1) = 0.0;		d2N[4](0,2) = 0.5;
		d2N[4](1,0) = 0.0;		d2N[4](1,1) = 0.0;		d2N[4](1,2) = 0.0;
		d2N[4](2,0) = 0.5;		d2N[4](2,1) = 0.0;		d2N[4](2,2) = 0.0;
	// Node 5
		d2N[5](0,0) = 0.0;		d2N[5](0,1) = 0.0;		d2N[5](0,2) = 0.0;
		d2N[5](1,0) = 0.0;		d2N[5](1,1) = 0.0;		d2N[5](1,2) = 0.5;
		d2N[5](2,0) = 0.0;		d2N[5](2,1) = 0.5;		d2N[5](2,2) = 0.0;
	return d2N;
}

bool Prism6::coords_inside(const std::vector<double>& rcoords)
{
	if (rcoords.size() != 3)
		err_message("Invalid coordinate size for coordinate inside check (PRISM6)!");

	if (rcoords[0] >= (0.0-_inside_tol) && rcoords[0] <= (1.0+_inside_tol))
	{
		if (rcoords[1] >= (0.0-_inside_tol) && rcoords[1] <= (1.0+_inside_tol))
		{
			if ((rcoords[0] + rcoords[1]) <= (1.0+_inside_tol))
			{
				if (rcoords[2] >= (-1.0-_inside_tol) && rcoords[2] <= (1.0+_inside_tol))
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
	else
		return false;
}





void Prism6::detection(std::vector<int>& node_detection, std::vector<Material*>& mats, bool isCohesive)
{
	err_message("Element detection is not presently supported for Prism6 elements.");
/*
	if(node_detection.size() != n_nodes())
		err_message("The number of nodes used for element detection must be the same as the number of nodes in the element.");

	// NOTE: the mats structure is a pointer to the material for every inclusion in the mesh, with the 0th element being a pointer to the base/matrix material
	// For this reason, when we index the mats array with the nodal detection values (incusion indices) we need to add 1
	id_type n_inclusions = mats.size();

	// Store the nodal detector values
	int n0 = node_detection[0];
	int n1 = node_detection[1];
	int n2 = node_detection[2];
	int n3 = node_detection[3];
	int n4 = node_detection[4];
	int n5 = node_detection[5];

	// Store many booleans that will be used many times

	//Storing whether the nodes belong to the same interface/inclusion as other nodes
	bool n0e1 = (n0 == n1); bool n1e0 = n0e1;
	bool n0e2 = (n0 == n2); bool n2e0 = n0e2;
	bool n0e3 = (n0 == n3); bool n3e0 = n0e3;
	bool n0e4 = (n0 == n4); bool n4e0 = n0e4;
	bool n0e5 = (n0 == n5); bool n5e0 = n0e5;

	bool n1e2 = (n1 == n2); bool n2e1 = n1e2;
	bool n1e3 = (n1 == n3); bool n3e1 = n1e3;
	bool n1e4 = (n1 == n4); bool n4e1 = n1e4;
	bool n1e5 = (n1 == n5); bool n5e1 = n1e5;

	bool n2e3 = (n2 == n3); bool n3e2 = n2e3;
	bool n2e4 = (n2 == n4); bool n4e2 = n2e4;
	bool n2e5 = (n2 == n5); bool n5e2 = n2e5;

	bool n3e4 = (n3 == n4); bool n4e3 = n3e4;
	bool n3e5 = (n3 == n5); bool n5e3 = n3e5;

	bool n4e5 = (n4 == n5); bool n5e4 = n4e5;

	std::vector<int> vec = node_detection;
	std::sort(vec.begin(), vec.end());
	std::vector<int>::iterator it;
	it = std::unique(vec.begin(), vec.end());
	vec.resize( std::distance(vec.begin(), it) );
	int n_unique_nodes = vec.size();

	//Check if an element is intersected by more than one interface
	// If so, the element needs to be refined
	if(n_unique_nodes > 2)
	{
		//Do something about refining the mesh here...
		err_message("An element is cut by more than one interface. This is not currently supported.");
	}

	std::vector<int> set_counts;
	for(int i=0; i<n_unique_nodes; ++i)
		set_counts.push_back(std::count(node_detection.begin(), node_detection.end(), vec[i]));


	
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




	else // The interface cuts the element somehow
	{
		_is_intersected = true;

		// Not cohesive
		if (!isCohesive)
		{
			// 2 unique nodal detection values
			if (n_unique_nodes == 2)
			{
				// One of the materials is the matrix (Only have one inclusions cutting this element)
				if (vec[0] < 0)
				{
					_enrichment_on_inclusion.resize(_cutEdge.size());
					std::fill(_enrichment_on_inclusion.begin(), _enrichment_on_inclusion.end(), *std::max_element(vec.begin(), vec.end())); // All intersected by the same inclusion
				}

				// 2 Different inclusions cut this element somehow covering all nodes
				else
				{

				}
			} // End 2 unique nodes



			// 3 unique nodal detection values
			else if (n_unique_nodes == 3)
			{
				// One of the materials is the matrix
				if (vec[0] < 0)
				{}

				// Otherwise 3 distinct inclusions cut this element covering all nodes
				else
				{}
			} // End 3 unique nodes



			// 4 unique nodal detection values
			else if (n_unique_nodes == 4)
			{
				// One of the materials is the matrix
				if (vec[0] < 0)
				{}

				// Otherwise 4 distinct inclusions cut this element covering all nodes
				else
				{}
			} // End 4 unique nodes



			// 5 unique nodal detection values
			else if (n_unique_nodes == 5)
			{
				// One of the materials is the matrix
				if (vec[0] < 0)
				{}

				// Otherwise 5 distinct inclusions cut this element covering all nodes
				else
				{}
			} // End 5 unique nodes



			// 6 unique nodal detection values
			else if (n_unique_nodes == 6)
			{
				// One of the materials is the matrix
				if (vec[0] < 0)
				{}

				// Otherwise 6 distinct inclusions cut this element covering all nodes
				else
				{}
			} // End 6 unique nodes

			_coh_elem_struct.clear();
			_enrichment_nodes.resize(_cutEdge.size());
		} // End is not cohesive




		// Cohesive
		else
		{
			_enrichment_nodes.resize(_cutEdge.size() * 2);
			_coh_mat.resize(_coh_elem_struct.size());
		} // End is cohesive
	} // End is intersected

	_detected = true;
*/
} // detection




void Prism6::refinement_nodes(std::vector<std::vector<id_type> >& refine_nodes)
{
	std::vector<id_type> ids(n_nodes());
	for (id_type n=0; n<n_nodes(); ++n)
		ids[n] = _nodes[n]->get_id();
	refine_nodes = {{ids[0], ids[3]},
					{ids[1], ids[4]},
					{ids[2], ids[5]},
					{ids[0], ids[1]},
					{ids[1], ids[2]},
					{ids[2], ids[0]},
					{ids[0], ids[1], ids[4], ids[3]},
					{ids[1], ids[2], ids[5], ids[4]},
					{ids[2], ids[0], ids[3], ids[5]},
					{ids[3], ids[4]},
					{ids[4], ids[5]},
					{ids[5], ids[3]}};
}
void Prism6::refinement_structure(std::vector<std::vector<id_type> >& structure, std::vector<elem_type>& types)
{
	structure = {{0, 9, 11, 6, 12, 14},
				 {1, 10, 9, 7, 13, 12},
				 {2, 11, 10, 8, 14, 13},
				 {9, 10, 11, 12, 13, 14},
				 {6, 12, 14, 3, 15, 17},
				 {7, 13, 12, 4, 16, 15},
				 {8, 14, 13, 5, 17, 16},
				 {12, 13, 14, 15, 16, 17}};
	types = {PRISM6, PRISM6, PRISM6, PRISM6, PRISM6, PRISM6, PRISM6, PRISM6};
}
