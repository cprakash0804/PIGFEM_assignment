/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#include "hex8.h"
#include "edge2.h"


// ------------------------------------------------------------
// Hex8 class static member initializations
const id_type Hex8::s_n_map[6][4] =
{
	{0, 3, 2, 1}, // Side 0
	{0, 1, 5, 4}, // Side 1
	{1, 2, 6, 5}, // Side 2
	{2, 3, 7, 6}, // Side 3
	{3, 0, 4, 7}, // Side 4
	{4, 5, 6, 7}  // Side 5
};

const id_type Hex8::e_n_map[12][2] =
{
	{0, 1}, // Edge 0
	{1, 2}, // Edge 1
	{2, 3}, // Edge 2
	{0, 3}, // Edge 3
	{0, 4}, // Edge 4
	{1, 5}, // Edge 5
	{2, 6}, // Edge 6
	{3, 7}, // Edge 7
	{4, 5}, // Edge 8
	{5, 6}, // Edge 9
	{6, 7}, // Edge 10
	{4, 7}  // Edge 11
};

id_type Hex8::side_nodes_map(id_type side, id_type node)
{
	if(side>=n_sides())
		err_message("Please select a valid side.");
	if(node>=4)
		err_message("Please select a valid node.");
	
	return s_n_map[side][node];
}

id_type Hex8::edge_nodes_map(id_type edge, id_type node)
{
	if(edge>=n_edges())
		err_message("Please select a valid edge.");
	if(node>=2)
		err_message("Please select a valid node.");
	
	return e_n_map[edge][node];
}

Elem* Hex8::build_side(id_type side)
{
	if(side >= n_sides())
		err_message("Side number must be less than the number of sides.");

	Elem* el = build(QUAD4);
	std::vector<Node*> nodes;
	for(id_type n=0; n<4; ++n)
		nodes.push_back(_nodes[side_nodes_map(side,n)]);
	el->set_nodes(nodes);
	return el;
}

  
id_type Hex8::n_q_points(int order) const
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
			return 8;
		case 4:
		case 5:
			return 27;
		case 6:
		case 7:
			return 64;
		case 8:
		case 9:
			return 125;
		case 10:
		case 11:
			return 216;
		case 12:
		case 13:
			return 343;
		case 14:
		case 15:
			return 512;
		case 16:
		case 17:
			return 729;
		case 18:
		case 19:
			return 1000;
		default:
			err_message("The selected quadrature order is not currently supported.");

	}
}
void Hex8::q_point(std::vector<double>& coords, double& w, id_type qp, int order) const
{
	if(order < 0)
		err_message("Order of integration must be greater than or equal to 0.");
	if(qp<0 || qp>=n_q_points(order))
		err_message("Quadrature point selected must be less than the maximum number of quadrature points for the given order.");
	
	coords.clear();
	Elem* edge = Elem::build(EDGE2);
	id_type nqp_edge = edge->n_q_points(order);
	
	std::vector<double> edge0_qp, edge1_qp, edge2_qp;
	double edge0_w, edge1_w, edge2_w;
	edge->q_point(edge0_qp, edge0_w, (qp%nqp_edge), order);
	edge->q_point(edge1_qp, edge1_w, ((qp/nqp_edge)%nqp_edge), order);
	edge->q_point(edge2_qp, edge2_w, (qp/(nqp_edge*nqp_edge)), order);
	delete edge;
	
	coords = {edge0_qp[0], edge1_qp[0], edge2_qp[0]};
	w = edge0_w*edge1_w*edge2_w;
}


std::vector<double> Hex8::compute_shape(const std::vector<double>& coords) const
{
	double r = coords[0];
	double s = coords[1];
	double t = coords[2];
	std::vector<double> N = {0.125*(1.0-r)*(1.0-s)*(1.0-t),
							 0.125*(1.0+r)*(1.0-s)*(1.0-t),
							 0.125*(1.0+r)*(1.0+s)*(1.0-t),
							 0.125*(1.0-r)*(1.0+s)*(1.0-t),
							 0.125*(1.0-r)*(1.0-s)*(1.0+t),
							 0.125*(1.0+r)*(1.0-s)*(1.0+t),
							 0.125*(1.0+r)*(1.0+s)*(1.0+t),
							 0.125*(1.0-r)*(1.0+s)*(1.0+t)};
	return N;
}

std::vector<std::vector<double> > Hex8::compute_shape_grad(const std::vector<double>& coords) const
{
	double r = coords[0];
	double s = coords[1];
	double t = coords[2];
	std::vector<std::vector<double> > dN = {{-0.125*(1.0-s)*(1-t), -0.125*(1.0-r)*(1-t), -0.125*(1.0-r)*(1.0-s)},
											{0.125*(1.0-s)*(1-t), -0.125*(1.0+r)*(1-t), -0.125*(1.0+r)*(1.0-s)},
											{0.125*(1.0+s)*(1-t), 0.125*(1.0+r)*(1-t), -0.125*(1.0+r)*(1.0+s)},
											{-0.125*(1.0+s)*(1-t), 0.125*(1.0-r)*(1-t), -0.125*(1.0-r)*(1.0+s)},
											{-0.125*(1.0-s)*(1+t), -0.125*(1.0-r)*(1+t), 0.125*(1.0-r)*(1.0-s)},
											{0.125*(1.0-s)*(1+t), -0.125*(1.0+r)*(1+t), 0.125*(1.0+r)*(1.0-s)},
											{0.125*(1.0+s)*(1+t), 0.125*(1.0+r)*(1+t), 0.125*(1.0+r)*(1.0+s)},
											{-0.125*(1.0+s)*(1+t), 0.125*(1.0-r)*(1+t), 0.125*(1.0-r)*(1.0+s)}};
	return dN;
}

std::vector<DenseMatrix<double> > Hex8::compute_shape_grad_grad(const std::vector<double>& coords) const
{
	double r = coords[0];
	double s = coords[1];
	double t = coords[2];
	std::vector<DenseMatrix<double> > d2N(n_nodes());
	// Node 0
		d2N[0](0,0) = 0.0;				d2N[0](0,1) = 0.125*(1.0-t);	d2N[0](0,2) = 0.125*(1.0-s);
		d2N[0](1,0) = 0.125*(1.0-t);	d2N[0](1,1) = 0.0;				d2N[0](1,2) = 0.125*(1.0-r);
		d2N[0](2,0) = 0.125*(1.0-s);	d2N[0](2,1) = 0.125*(1.0-r);	d2N[0](2,2) = 0.0;
	// Node 1
		d2N[1](0,0) = 0.0;				d2N[1](0,1) = -0.125*(1.0-t);	d2N[1](0,2) = -0.125*(1.0-s);
		d2N[1](1,0) = -0.125*(1.0-t);	d2N[1](1,1) = 0.0;				d2N[1](1,2) = 0.125*(1.0+r);
		d2N[1](2,0) = -0.125*(1.0-s);	d2N[1](2,1) = 0.125*(1.0+r);	d2N[1](2,2) = 0.0;
	// Node 2
		d2N[2](0,0) = 0.0;				d2N[2](0,1) = 0.125*(1.0-t);	d2N[2](0,2) = -0.125*(1.0+s);
		d2N[2](1,0) = 0.125*(1.0-t);	d2N[2](1,1) = 0.0;				d2N[2](1,2) = -0.125*(1.0+r);
		d2N[2](2,0) = -0.125*(1.0+s);	d2N[2](2,1) = -0.125*(1.0+r);	d2N[2](2,2) = 0.0;
	// Node 3
		d2N[3](0,0) = 0.0;				d2N[3](0,1) = -0.125*(1.0-t);	d2N[3](0,2) = 0.125*(1.0+s);
		d2N[3](1,0) = -0.125*(1.0-t);	d2N[3](1,1) = 0.0;				d2N[3](1,2) = -0.125*(1.0-r);
		d2N[3](2,0) = 0.125*(1.0+s);	d2N[3](2,1) = -0.125*(1.0-r);	d2N[3](2,2) = 0.0;
	// Node 4
		d2N[4](0,0) = 0.0;				d2N[4](0,1) = 0.125*(1.0+t);	d2N[4](0,2) = -0.125*(1.0-s);
		d2N[4](1,0) = 0.125*(1.0+t);	d2N[4](1,1) = 0.0;				d2N[4](1,2) = -0.125*(1.0-r);
		d2N[4](2,0) = -0.125*(1.0-s);	d2N[4](2,1) = -0.125*(1.0-r);	d2N[4](2,2) = 0.0;
	// Node 5
		d2N[5](0,0) = 0.0;				d2N[5](0,1) = -0.125*(1.0+t);	d2N[5](0,2) = 0.125*(1.0-s);
		d2N[5](1,0) = -0.125*(1.0+t);	d2N[5](1,1) = 0.0;				d2N[5](1,2) = -0.125*(1.0+r);
		d2N[5](2,0) = 0.125*(1.0-s);	d2N[5](2,1) = -0.125*(1.0+r);	d2N[5](2,2) = 0.0;
	// Node 6
		d2N[6](0,0) = 0.0;				d2N[6](0,1) = 0.125*(1.0+t);	d2N[6](0,2) = 0.125*(1.0+s);
		d2N[6](1,0) = 0.125*(1.0+t);	d2N[6](1,1) = 0.0;				d2N[6](1,2) = 0.125*(1.0+r);
		d2N[6](2,0) = 0.125*(1.0+s);	d2N[6](2,1) = 0.125*(1.0+r);	d2N[6](2,2) = 0.0;
	// Node 7
		d2N[7](0,0) = 0.0;				d2N[7](0,1) = -0.125*(1.0+t);	d2N[7](0,2) = -0.125*(1.0+s);
		d2N[7](1,0) = -0.125*(1.0+t);	d2N[7](1,1) = 0.0;				d2N[7](1,2) = 0.125*(1.0-r);
		d2N[7](2,0) = -0.125*(1.0+s);	d2N[7](2,1) = 0.125*(1.0-r);	d2N[7](2,2) = 0.0;
	return d2N;
}


bool Hex8::coords_inside(const std::vector<double>& rcoords)
{
	if (rcoords.size() != 3)
		err_message("Invalid coordinate size for coordinate inside check (HEX8)!");

	if (rcoords[0] >= (-1.0-_inside_tol) && rcoords[0] <= (1.0+_inside_tol))
	{
		if (rcoords[1] >= (-1.0-_inside_tol) && rcoords[1] <= (1.0+_inside_tol))
		{
			if (rcoords[2] >= (-1.0-_inside_tol) && rcoords[2] <= (1.0+_inside_tol))
			{
				return true;
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





void Hex8::detection(std::vector<int>& node_detection, std::vector<Material*>& mats, bool isCohesive)
{
	err_message("Hexahedral element detection is not currently supported");
/*
	if(node_detection.size() != n_nodes())
		err_message("The number of nodes used for element detection must be the same as the number of nodes in the element.");

	// NOTE: the mats structure is a pointer to the material for every inclusion in the mesh, with the 0th elment being a pointer to the base/matrix material
	// For this reason, when we index the mats array with the nodal detection values (incusion indices) we need to add 1
	//id_type n_inclusions = mats.size();

	// Store the nodal detector values
	int n0 = node_detection[0];
	int n1 = node_detection[1];
	int n2 = node_detection[2];
	int n3 = node_detection[3];
	int n4 = node_detection[4];
	int n5 = node_detection[5];
	int n6 = node_detection[6];
	int n7 = node_detection[7];

	// Store many booleans that will be used many times

	//Storing whether the nodes belong to the same interface/inclusion as other nodes
	bool n0e1 = (n0 == n1); //bool n1e0 = n0e1;
	bool n0e2 = (n0 == n2); //bool n2e0 = n0e2;
	bool n0e3 = (n0 == n3); //bool n3e0 = n0e3;
	bool n0e4 = (n0 == n4); //bool n4e0 = n0e4;
	bool n0e5 = (n0 == n5); //bool n5e0 = n0e5;
	bool n0e6 = (n0 == n6); //bool n6e0 = n0e6;
	bool n0e7 = (n0 == n7); //bool n7e0 = n0e7;

	bool n1e2 = (n1 == n2); //bool n2e1 = n1e2;
	bool n1e3 = (n1 == n3); //bool n3e1 = n1e3;
	bool n1e4 = (n1 == n4); //bool n4e1 = n1e4;
	bool n1e5 = (n1 == n5); //bool n5e1 = n1e5;
	bool n1e6 = (n1 == n6); //bool n6e1 = n1e6;
	bool n1e7 = (n1 == n7); //bool n7e1 = n1e7;

	bool n2e3 = (n2 == n3); //bool n3e2 = n2e3;
	bool n2e4 = (n2 == n4); //bool n4e2 = n2e4;
	bool n2e5 = (n2 == n5); //bool n5e2 = n2e5;
	bool n2e6 = (n2 == n6); //bool n6e2 = n2e6;
	bool n2e7 = (n2 == n7); //bool n7e2 = n2e7;

	bool n3e4 = (n3 == n4); //bool n4e3 = n3e4;
	bool n3e5 = (n3 == n5); //bool n5e3 = n3e5;
	bool n3e6 = (n3 == n6); //bool n6e3 = n3e6;
	bool n3e7 = (n3 == n7); //bool n7e3 = n3e7;

	bool n4e5 = (n4 == n5); //bool n5e4 = n4e5;
	bool n4e6 = (n4 == n6); //bool n6e4 = n4e6;
	bool n4e7 = (n4 == n7); //bool n7e4 = n4e7;

	bool n5e6 = (n5 == n6); //bool n6e5 = n5e6;
	bool n5e7 = (n5 == n7); //bool n7e5 = n5e7;

	bool n6e7 = (n6 == n7); //bool n7e6 = n6e7;

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
					// 3-edge cut
					// One corner is cut off so the three edges connecting it must be cut
					if((set_counts[0]==1) || (set_counts[0]==7))
					{
						if(!n0e1 && n1e2) // Node 0 is seperate
							_cutEdge = {0, 3, 4};
						else if(!n0e1 && n0e2) // Node 1 is seperate
							_cutEdge = {0, 1, 5};
						else if(!n0e2 && n0e1) // Node 2 is seperate
							_cutEdge = {1, 2, 6};
						else if(!n0e3 && n0e1) // Node 3 is seperate
							_cutEdge = {2, 3, 7};
						else if(!n0e4 && n0e1) // Node 4 is seperate
							_cutEdge = {4, 8, 11};
						else if(!n0e5 && n0e1) // Node 5 is seperate
							_cutEdge = {5, 8, 9};
						else if(!n0e6 && n0e1) // Node 6 is seperate
							_cutEdge = {6, 9, 10};
						else if(!n0e7 && n0e1) // Node 7 is seperate
							_cutEdge = {7, 10, 11};
						else
						{
							//This should hopefully never happen
							//If it does it probably means the mesh needs to be refined somehow...
							err_message("HEXAHEDRAL ELEMENT HAS 1 NODE SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
						}
					}

					// 4-edge cut, type 1
					// A triangular prism is cut by cutting off an edge effectively
					else if((set_counts[0]==2) || (set_counts[0]==6))
					{
						if(n0e1 && !n0e2 && n2e3 && n2e4) //edge 0
							_cutEdge = {1, 3, 4, 5};
						else if(n1e2 && !n0e1 && n0e3 && n0e4) //edge 1
							_cutEdge = {0, 2, 5, 6};
						else if(n2e3 && !n0e2 && n0e1 && n0e4) //edge 2
							_cutEdge = {1, 3, 6, 7};
						else if(n0e3 && !n0e1 && n1e2 && n1e4) //edge 3
							_cutEdge = {0, 2, 4, 7};
						else if(n0e4 && !n0e1 && n1e2 && n1e3) //edge 4
							_cutEdge = {0, 3, 8, 11};
						else if(n1e5 && !n0e1 && n0e2 && n0e3) //edge 5
							_cutEdge = {0, 1, 8, 9};
						else if(n2e6 && !n0e2 && n0e1 && n0e3) //edge 6
							_cutEdge = {1, 2, 9, 10};
						else if(n3e7 && !n0e3 && n0e1 && n0e2) //edge 7
							_cutEdge = {2, 3, 10, 11};
						else if(n4e5 && !n0e4 && n0e1 && n0e2) //edge 8
							_cutEdge = {4, 5, 9, 11};
						else if(n5e6 && !n0e5 && n0e1 && n0e2) //edge 9
							_cutEdge = {5, 6, 8, 10};
						else if(n6e7 && !n0e6 && n0e1 && n0e2) //edge 10
							_cutEdge = {6, 7, 9, 11};
						else if(n4e7 && !n0e4 && n0e1 && n0e2) //edge 11
							_cutEdge = {4, 7, 8, 10};
						else
						{
							//This should hopefully never happen
							//If it does it probably means the mesh needs to be refined somehow...
							err_message("HEXAHEDRAL ELEMENT HAS 2 NODES SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
						}
					}

					// 5-edge cut
					// The nodes are split into sets of 3 neighboring nodes and the remaining 5 nodes by an interface that cuts 5 edges
					// 24 cases...
					// The odrer that the edges are put into the _cutEdge vector is important here. It is the order in which the cutting surface will be separated into triangles
					// Each successive triad of sides represents one of the trangles that the surface will be composed of
					else if((set_counts[0]==3) ||(set_counts[0]==5))
					{
						if(n0e1 && n0e3 && !n0e2 && n2e4 && n2e5 && n2e6 && n2e7) // {0, 1, 3}/{2, 4, 5, 6, 7} set (2-4 diagonal, moved down)
							_cutEdge = {4, 5, 7, 1, 2};
						else if(n0e1 && n0e4 && !n0e2 && n2e3 && n2e5 && n2e6 && n2e7) // {0, 1, 4}/{2, 3, 5, 6, 7} set
							_cutEdge = {3, 1, 11, 5, 8};
						else if(n0e3 && n0e4 && !n0e1 && n1e2 && n1e5 && n1e6 && n1e7) // {0, 3, 4}/{1, 2, 5, 6, 7} set
							_cutEdge = {0, 2, 8, 7, 11};
						else if(n0e1 && n0e2 && !n1e3 && n3e4 && n3e5 && n3e6 && n3e7) // {1, 0, 2}/{3, 4, 5, 6, 7} set (3-5 diagonal, moved down)
							_cutEdge = {5, 4, 6, 3, 2};
						else if(n0e1 && n0e5 && !n1e2 && n2e3 && n2e4 && n2e6 && n2e7) // {1, 0, 5}/{2, 3, 4, 6, 7} set
							_cutEdge = {1, 3, 9, 4, 8};
						else if(n1e2 && n1e5 && !n0e1 && n0e3 && n0e4 && n0e6 && n0e7) // {1, 2, 5}/{0, 3, 4, 6, 7} set
							_cutEdge = {0, 2, 8, 6, 9};
						else if(n1e2 && n1e3 && !n0e1 && n0e4 && n0e5 && n0e6 && n0e7) // {2, 1, 3}/{0, 4, 5, 6, 7} set (0-6 diagonal, moved down)
							_cutEdge = {6, 5, 7, 0, 3};
						else if(n1e2 && n1e6 && !n0e1 && n0e3 && n0e4 && n0e5 && n0e7) // {2, 1, 6}/{0, 3, 4, 5, 7} set
							_cutEdge = {2, 0, 10, 5, 9};
						else if(n2e3 && n2e6 && !n0e2 && n0e1 && n0e4 && n0e5 && n0e7) // {2, 3, 6}/{0, 1, 4, 5, 7} set
							_cutEdge = {1, 3, 9, 7, 10};
						else if(n0e2 && n0e3 && !n0e1 && n1e4 && n1e5 && n1e6 && n1e7) // {3, 0, 2}/{1, 4, 5, 6, 7} set (1-7 diagonal, moved down)
							_cutEdge = {7, 4, 6, 0, 1};
						else if(n0e3 && n0e7 && !n0e1 && n1e2 && n1e4 && n1e5 && n1e6) // {3, 0, 7}/{1, 2, 4, 5, 6} set
							_cutEdge = {2, 0, 10, 4, 11};
						else if(n2e3 && n2e7 && !n0e2 && n0e1 && n0e4 && n0e5 && n0e6) // {3, 2, 7}/{0, 1, 4, 5, 6} set
							_cutEdge = {3, 1, 11, 6, 10};
						else if(n4e5 && n4e7 && !n0e4 && n0e1 && n0e2 && n0e3 && n0e6) // {4, 5, 7}/{0, 1, 2, 3, 6} set (0-6 diagonal, moved up)
							_cutEdge = {4, 5, 7, 9, 10};
						else if(n0e4 && n0e5 && !n0e1 && n1e2 && n1e3 && n1e6 && n1e7) // {4, 0, 5}/{1, 2, 3, 6, 7} set
							_cutEdge = {11, 9, 3, 5, 0};
						else if(n0e4 && n0e7 && !n0e1 && n1e2 && n1e3 && n1e5 && n1e6) // {4, 0, 7}/{1, 2, 3, 5, 6} set
							_cutEdge = {8, 0, 10, 3, 7};
						else if(n4e5 && n4e6 && !n0e4 && n0e1 && n0e2 && n0e3 && n0e7) // {5, 4, 6}/{0, 1, 2, 3, 7} set (1-7 diagonal, moved up)
							_cutEdge = {5, 4, 6, 11, 10};
						else if(n1e5 && n1e4 && !n0e1 && n0e2 && n0e3 && n0e6 && n0e7) // {5, 1, 4}/{0, 2, 3, 6, 7} set
							_cutEdge = {9, 1, 11, 0, 4};
						else if(n1e5 && n1e6 && !n0e1 && n0e2 && n0e3 && n0e4 && n0e7) // {5, 1, 6}/{0, 2, 3, 4, 7} set
							_cutEdge = {8, 0, 10, 1, 6};
						else if(n5e6 && n5e7 && !n0e5 && n0e1 && n0e2 && n0e3 && n0e4) // {6, 5, 7}/{0, 1, 2, 3, 4} set (2-4 diagonal, moved up)
							_cutEdge = {6, 5, 7, 8, 11};
						else if(n2e6 && n2e5 && !n0e2 && n0e1 && n0e3 && n0e4 && n0e7) // {6, 2, 5}/{0, 1, 3, 4, 7} set
							_cutEdge = {10, 2, 8, 1, 5};
						else if(n2e6 && n2e7 && !n0e2 && n0e1 && n0e3 && n0e4 && n0e5) // {6, 2, 7}/{0, 1, 3, 4, 5} set
							_cutEdge = {9, 1, 11, 2, 7};
						else if(n3e4 && n3e7 && !n0e3 && n0e1 && n0e2 && n0e5 && n0e6) // {7, 3, 4}/{0, 1, 2, 5, 6} set
							_cutEdge = {10, 8, 2, 4, 3};
						else if(n3e6 && n3e7 && !n0e3 && n0e1 && n0e2 && n0e4 && n0e5) // {7, 3, 6}/{0, 1, 2, 4, 5} set
							_cutEdge = {11, 3, 9, 2, 6};
						else if(n4e6 && n4e7 && !n0e4 && n0e1 && n0e2 && n0e3 && n0e5) // {7, 4, 6}/{0, 1, 2, 3, 5} set (3-5 diagonal, moved up)
							_cutEdge = {7, 4, 6, 8, 9};
						else
						{
							//This should hopefully never happen
							//If it does it probably means the mesh needs to be refined somehow...
							err_message("HEXAHEDRAL ELEMENT HAS 3 NODES SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
						}
					}

					// Either 4-edge cut type 2 or a six edge cut
					else if(set_counts[0] == 4)
					{
						//First three cases cover if one entire face is cut away
						if (n0e1 && n0e2 && n0e3 && !n0e4 && n4e5 && n4e6 && n4e7) // {0, 1, 2, 3}/{5, 6, 7, 8} face
							{_cutEdge = {4, 5, 6, 7};}
						else if(n0e3 && n3e4 && n4e7 && !n0e1 && n1e2 && n2e5 && n5e6) // {0, 3, 4, 7}/{1, 2, 5, 6} face
							{_cutEdge = {0, 2, 8, 10};}
						else if(n0e1 && n1e4 && n4e5 && !n0e2 && n2e3 && n3e6 && n6e7) // {0, 1, 4, 5}/{2, 3, 6, 7} face
							{_cutEdge = {1, 3, 9, 11};}
						
						// Next 4 cases cover a diagonal inclusion face cutting away 4 neighboring nodes in a pyramid shape
						else if(n0e1 && n1e3 && n3e4 && !n0e2 && n2e5 && n5e6 && n6e7) // {0, 1, 3, 4}/{6, 2, 5, 7} set
							{_cutEdge = {1, 2, 5, 7, 8, 11};}
						else if(n0e1 && n1e2 && n2e5 && !n0e3 && n3e4 && n4e6 && n6e7) // {1, 0, 2, 5}/{7, 3, 4, 6} set
							{_cutEdge = {3, 2, 4, 6, 8, 9};}
						else if(n0e4 && n4e5 && n5e7 && !n0e1 && n1e2 && n2e3 && n3e6) // {4, 0, 5, 7}/{2, 1, 3, 6} set
							{_cutEdge = {3, 0, 7, 5, 10, 9};}
						else if(n0e2 && n2e3 && n3e7 && !n0e1 && n1e4 && n4e5 && n5e6) // {3, 0, 2, 7}/{5, 1, 4, 6} set
							{_cutEdge = {0, 1, 4, 6, 11, 10};}
						else
						{
							//This should hopefully never happen
							//If it does it probably means the mesh needs to be refined somehow...
							err_message("HEXAHEDRAL ELEMENT HAS 4 NODES SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
						}
					}

					// Final Catch-all case
					else
					{
						//This REALLY should never happen
						//No idea what went wrong to get here. The sorting of the set_counts vector should never allow us to get here
						err_message("Hex8 detection. We'll never get here!");
					}

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



			// 7 unique nodal detection values
			else if (n_unique_nodes == 7)
			{
				// One of the materials is the matrix
				if (vec[0] < 0)
				{}

				// Otherwise 7 distinct inclusions cut this element covering all nodes
				else
				{}
			} // End 7 unique nodes



			// 8 unique nodal detection values
			else if (n_unique_nodes == 8)
			{
				// One of the materials is the matrix
				if (vec[0] < 0)
				{}

				// Otherwise 8 distinct inclusions cut this element covering all nodes
				else
				{}
			} // End 8 unique nodes

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



void Hex8::refinement_nodes(std::vector<std::vector<id_type> >& refine_nodes)
{
	std::vector<id_type> ids(n_nodes());
	for (id_type n=0; n<n_nodes(); ++n)
		ids[n] = _nodes[n]->get_id();
	refine_nodes = {{ids[0], ids[4]},
					{ids[1], ids[5]},
					{ids[2], ids[6]},
					{ids[3], ids[7]},
					{ids[0], ids[1]},
					{ids[1], ids[2]},
					{ids[2], ids[3]},
					{ids[3], ids[0]},
					{ids[0], ids[1], ids[2], ids[3]},
					{ids[0], ids[1], ids[5], ids[4]},
					{ids[1], ids[2], ids[6], ids[5]},
					{ids[2], ids[3], ids[7], ids[6]},
					{ids[3], ids[0], ids[4], ids[7]},
					{ids[0], ids[1], ids[2], ids[3], ids[4], ids[5], ids[6], ids[7]},
					{ids[4], ids[5]},
					{ids[5], ids[6]},
					{ids[6], ids[7]},
					{ids[7], ids[4]},
					{ids[4], ids[5], ids[6], ids[7]}};
}
void Hex8::refinement_structure(std::vector<std::vector<id_type> >& structure, std::vector<elem_type>& types)
{
	structure = {{0, 12, 16, 15, 8, 17, 21, 20},
				 {12, 1, 13, 16, 17, 9, 18, 21},
				 {16, 13, 2, 14, 21, 18, 10, 19},
				 {15, 16, 14, 3, 20, 21, 19, 11},
				 {8, 17, 21, 20, 4, 22, 26, 25},
				 {17, 9, 18, 21, 22, 5, 23, 26},
				 {21, 18, 10, 19, 26, 23, 6, 24},
				 {20, 21, 19, 11, 25, 26, 24, 7}};
	types = {HEX8, HEX8, HEX8, HEX8, HEX8, HEX8, HEX8, HEX8};
}
