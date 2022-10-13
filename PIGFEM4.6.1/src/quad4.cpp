/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#include "quad4.h"
#include <iostream>

// ------------------------------------------------------------
// Quad class static member initialization
const id_type Quad4::s_n_map[4][2] =
{
	{0, 1}, // Side 0
	{1, 2}, // Side 1
	{2, 3}, // Side 2
	{3, 0}  // Side 3
};

const id_type Quad4::e_n_map[4][2] =
{
	{0, 1}, // Side 0
	{1, 2}, // Side 1
	{2, 3}, // Side 2
	{3, 0}  // Side 3
};

id_type Quad4::side_nodes_map(id_type side, id_type node)
{
	if(side>=n_sides())
		err_message("Please select a valid side.");
	if(node>=2)
		err_message("Please select a valid node.");
	
	return s_n_map[side][node];
}

id_type Quad4::edge_nodes_map(id_type edge, id_type node)
{
	if(edge>=4)
		err_message("Please select a valid edge.");
	if(node>=2)
		err_message("Please select a valid node.");
	
	return e_n_map[edge][node];
}

Elem* Quad4::build_side(id_type side)
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

// Quadrature rules for the quad are a composite of the rules for an edge
id_type Quad4::n_q_points(int order) const
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

void Quad4::q_point(std::vector<double>& coords, double& w, id_type qp, int order) const
{
	if(order < 0)
		err_message("Order of integration must be greater than or equal to 0.");
	if(qp<0 || qp>=n_q_points(order))
		err_message("Quadrature point selected must be less than the maximum number of quadrature points for the given order.");
	
	coords.clear();
	Elem* edge = Elem::build(EDGE2);
	id_type nqp_edge = edge->n_q_points(order);
	
	std::vector<double> edge0_qp, edge1_qp;
	double edge0_w, edge1_w;
	edge->q_point(edge0_qp, edge0_w, (qp%nqp_edge), order);
	edge->q_point(edge1_qp, edge1_w, (qp/nqp_edge), order);
	delete edge;
	
	coords = {edge0_qp[0], edge1_qp[0]};
	w = edge0_w*edge1_w;
}

std::vector<double> Quad4::compute_shape(const std::vector<double>& coords) const
{
	double r = coords[0];
	double s = coords[1];
	std::vector<double> N = {0.25*(1.0-r)*(1.0-s),
							0.25*(1.0+r)*(1.0-s),
							0.25*(1.0+r)*(1.0+s),
							0.25*(1.0-r)*(1.0+s)};
	return N;
}

std::vector<std::vector<double> > Quad4::compute_shape_grad(const std::vector<double>& coords) const
{
	double r = coords[0];
	double s = coords[1];
	std::vector<std::vector<double> > dN = {{-0.25*(1.0-s), -0.25*(1.0-r)},
											{0.25*(1.0-s), -0.25*(1.0+r)},
											{0.25*(1.0+s), 0.25*(1.0+r)},
											{-0.25*(1.0+s), 0.25*(1.0-r)}};
	return dN;
}

std::vector<DenseMatrix<double> > Quad4::compute_shape_grad_grad(const std::vector<double>& coords) const
{
	std::vector<DenseMatrix<double> > d2N(n_nodes());
	for (id_type n=0; n<n_nodes(); ++n)
	{
		d2N[n].resize(2,2);
		d2N[n](0,1) = 0.25 * pow(-1, n);
		d2N[n](1,0) = 0.25 * pow(-1, n);
	}
	return d2N;
}


bool Quad4::coords_inside(const std::vector<double>& rcoords)
{
	if (rcoords.size() != 2)
		err_message("Invalid coordinate size for coordinate inside check (QUAD4)!");

	if (rcoords[0] >= (-1.0-_inside_tol) && rcoords[0] <= (1.0+_inside_tol))
	{
		if (rcoords[1] >= (-1.0-_inside_tol) && rcoords[1] <= (1.0+_inside_tol))
		{
			return true;
		}
		else
			return false;
	}
	else
		return false;
}




// IGFEM MATERIAL BELOW HERE
//==============================================================================================================================

void Quad4::detection(std::vector<int>& node_detection, std::vector<Material*>& mats, bool isCohesive)
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
	int n3 = node_detection[3];

	//Storing whether the nodes belong to the same interface/inclusion as other nodes
	bool n0e1 = (n0 == n1); //bool n1e0 = n0e1;
	bool n0e2 = (n0 == n2); //bool n2e0 = n0e2;
	bool n0e3 = (n0 == n3); //bool n3e0 = n0e3;
	bool n1e2 = (n1 == n2); //bool n2e1 = n1e2;
	bool n1e3 = (n1 == n3); //bool n3e1 = n1e3;
	bool n2e3 = (n2 == n3); //bool n3e2 = n2e3;

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
		_enrichment_on_inclusion.clear();
		_int_elem_struct.clear();
		_coh_elem_struct.clear();
		_int_elem_mat.clear();
	}



	// ==============================================================================
	// CASES WRITTEN USING POLYGONAL INTEGRATION ELEMENTS
	// ==============================================================================


	// This element is intersected somehow
	else
	{
		_is_intersected = true;

		// Not cohesive cases
		if (!isCohesive)
		{
			// 2 unique nodal detection values
			if (n_unique_nodes == 2)
			{
				// One of the materials is the matrix (Only have one inclusions cutting this element)
				if (vec[0] < 0)
				{
					// Three nodes are located in one material and one in another
					if((set_counts[0]==1) || (set_counts[0]==3))
					{
						if(n1e2 && n1e3)		// node zero is separate
							{_cutEdge = {3,0}; _int_elem_struct = {{0,5,4},{1,2,3,4,5}}; _int_elem_mat = {mats[n0+1],mats[n2+1]};}
						else if(n0e2 && n0e3)	// node one is separate
							{_cutEdge = {0,1}; _int_elem_struct = {{1,5,4},{2,3,0,4,5}}; _int_elem_mat = {mats[n1+1],mats[n3+1]};}
						else if(n0e1 && n0e3)	// node two is separate
							{_cutEdge = {1,2}; _int_elem_struct = {{2,5,4},{3,0,1,4,5}}; _int_elem_mat = {mats[n2+1],mats[n0+1]};}
						else if(n0e1 && n1e2)		// node three is separate
							{_cutEdge = {2,3}; _int_elem_struct = {{3,5,4},{0,1,2,4,5}}; _int_elem_mat = {mats[n3+1],mats[n1+1]};}
						else
							err_message("QUADRILATERAL ELEMENT HAS 1 NODE SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
					}
					
					// Two cases:
					//	1. Two adacent nodes are in one material and the other two nodes are in another, this splitting this into 2 Quad4's-direction
					//	2. Two opposite nodes are in the same material and the other two nodes are in the same material. This means the interface is highly curved and an error should be thrown (will later lead to mesh refinement)
					else if(set_counts[0]==2)
					{
						if(n0e1)			// {0,1}/{2,3}
							{_cutEdge = {1, 3}; _int_elem_struct = {{0, 1, 4, 5}, {2, 3, 5, 4}}; _int_elem_mat = {mats[n0+1], mats[n2+1]};}
						else if(n1e2)		// {0,3}/{1,2}
							{_cutEdge = {2, 0}; _int_elem_struct = {{1, 2, 4, 5}, {3, 0, 5, 4}}; _int_elem_mat = {mats[n1+1], mats[n3+1]};}
						else				// Interface is highly curved
							err_message("QUADRILATERAL ELEMENT HAS 2 NODES SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
					}
					
					else
						err_message("Quad4 element detector: We'll never get here! (Hopefully)");

					_enrichment_on_inclusion.resize(_cutEdge.size());
					std::fill(_enrichment_on_inclusion.begin(), _enrichment_on_inclusion.end(), *std::max_element(vec.begin(), vec.end())); // All intersected by the same inclusion
				}



				// Otherwise there are two inclusions cutting this element
				else
				{
					// Three nodes are located in one material and one in another
					if((set_counts[0]==1) || (set_counts[0]==3))
					{
						if(n1e2 && n1e3)		// node zero is separate
							{_cutEdge = {0,3,3,0}; _enrichment_on_inclusion = {n0,n0,n2,n2}; _int_elem_struct = {{0,4,5},{5,4,7,6},{1,2,3,6,7}}; _int_elem_mat = {mats[n0+1],matrix,mats[n2+1]};}
						else if(n0e2 && n0e3)	// node one is separate
							{_cutEdge = {1,0,0,1}; _enrichment_on_inclusion = {n1,n1,n3,n3}; _int_elem_struct = {{1,4,5},{5,4,7,6},{2,3,0,6,7}}; _int_elem_mat = {mats[n1+1],matrix,mats[n3+1]};}
						else if(n0e1 && n0e3)	// node two is separate
							{_cutEdge = {2,1,1,2}; _enrichment_on_inclusion = {n2,n2,n0,n0}; _int_elem_struct = {{2,4,5},{5,4,7,6},{3,0,1,6,7}}; _int_elem_mat = {mats[n2+1],matrix,mats[n0+1]};}
						else if(n0e1 && n1e2)		// node three is separate
							{_cutEdge = {3,2,2,3}; _enrichment_on_inclusion = {n3,n3,n1,n1}; _int_elem_struct = {{3,4,5},{5,4,7,6},{0,1,2,6,7}}; _int_elem_mat = {mats[n3+1],matrix,mats[n1+1]};}
						else
							err_message("QUADRILATERAL ELEMENT HAS 1 NODE SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
					}
					
					// Two cases:
					//	1. Two adacent nodes are in one material and the other two nodes are in another, this splitting this into 2 Quad4's-direction
					//	2. Two opposite nodes are in the same material and the other two nodes are in the same material. This means the interface is highly curved and an error should be thrown (will later lead to mesh refinement)
					else if(set_counts[0]==2)
					{
						if(n0e1)			// {0,1}/{2,3}
							{_cutEdge = {1,3,3,1}; _enrichment_on_inclusion = {n0,n0,n2,n2}; _int_elem_struct = {{0,1,4,5},{5,4,7,6},{2,3,6,7}}; _int_elem_mat = {mats[n0+1],matrix,mats[n2+1]};}
						else if(n1e2)		// {0,3}/{1,2}
							{_cutEdge = {2,0,0,2}; _enrichment_on_inclusion = {n1,n1,n3,n3}; _int_elem_struct = {{1,2,4,5},{5,4,7,6},{3,0,6,7}}; _int_elem_mat = {mats[n1+1],matrix,mats[n3+1]};}
						else				// Interface is highly curved
							err_message("QUADRILATERAL ELEMENT HAS 2 NODES SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
					}
					
					else
						err_message("Quad4 element detector: We'll never get here! (Hopefully)");


					_enrichment_nodes.resize(4);
				}
			} // End 2 unique nodes



			// 3 unique nodal detection values
			else if (n_unique_nodes == 3)
			{
				// One of the materials is the matrix (2 options)
				if (vec[0] < 0)
				{
					// Case where 2 of the nodes are in the matrix and the remaining 2 nodes are in different inclusions
					if (set_counts[0] == 2)
					{
						if (n0e1)		// Nodes 0 and 1 are in the matrix
							{_cutEdge = {3,2,2,1}; _enrichment_on_inclusion = {n3,n3,n2,n2}; _int_elem_struct = {{3,4,5},{2,6,7},{0,1,7,6,5,4}}; _int_elem_mat = {mats[n3+1],mats[n2+1],matrix};}
						else if (n1e2)	// Nodes 1 and 2 are in the matrix
							{_cutEdge = {0,3,3,2}; _enrichment_on_inclusion = {n0,n0,n3,n3}; _int_elem_struct = {{0,4,5},{3,6,7},{1,2,7,6,5,4}}; _int_elem_mat = {mats[n0+1],mats[n3+1],matrix};}
						else if (n2e3)	// Nodes 2 and 3 are in the matrix
							{_cutEdge = {1,0,0,3}; _enrichment_on_inclusion = {n1,n1,n0,n0}; _int_elem_struct = {{1,4,5},{0,6,7},{2,3,7,6,5,4}}; _int_elem_mat = {mats[n1+1],mats[n0+1],matrix};}
						else if (n0e3)	// Nodes 3 and 0 are in the matrix
							{_cutEdge = {2,1,1,0}; _enrichment_on_inclusion = {n2,n2,n1,n1}; _int_elem_struct = {{2,4,5},{1,6,7},{3,0,7,6,5,4}}; _int_elem_mat = {mats[n2+1],mats[n1+1],matrix};}
						else if (n0e2)	// Nodes 0 and 2 are in the matrix
							{_cutEdge = {1,0,3,2}; _enrichment_on_inclusion = {n1,n1,n3,n3}; _int_elem_struct = {{1,4,5},{3,6,7},{2,7,6,0,5,4}}; _int_elem_mat = {mats[n1+1],mats[n3+1],matrix};}
						else if (n1e3)	// Nodes 1 and 3 are in the matrix
							{_cutEdge = {0,3,2,1}; _enrichment_on_inclusion = {n0,n0,n2,n2}; _int_elem_struct = {{0,4,5},{2,6,7},{1,7,6,3,5,4}}; _int_elem_mat = {mats[n0+1],mats[n2+1],matrix};}
						else
							err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
					}

					// Otherwise, the matrix is only on one of the nodes and two of the other nodes are covered by the same inclusion
					else
					{
						if (n0e1)
						{
							if (n2<0)
								{_cutEdge = {1,3,3,2}; _enrichment_on_inclusion = {n0,n0,n3,n3}; _int_elem_struct = {{0,1,4,5},{3,6,7},{2,7,6,5,4}}; _int_elem_mat = {mats[n0+1],mats[n3+1],matrix};}
							else if (n3<0)
								{_cutEdge = {1,3,2,1}; _enrichment_on_inclusion = {n0,n0,n2,n2}; _int_elem_struct = {{0,1,4,5},{2,6,7},{3,5,4,7,6}}; _int_elem_mat = {mats[n0+1],mats[n2+1],matrix};}
							else err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
						}
						else if (n1e2)
						{
							if (n3<0)
								{_cutEdge = {2,0,0,3}; _enrichment_on_inclusion = {n1,n1,n0,n0}; _int_elem_struct = {{1,2,4,5},{0,6,7},{3,7,6,5,4}}; _int_elem_mat = {mats[n1+1],mats[n0+1],matrix};}
							else if (n0<0)
								{_cutEdge = {2,0,3,2}; _enrichment_on_inclusion = {n1,n1,n3,n3}; _int_elem_struct = {{1,2,4,5},{3,6,7},{0,5,4,7,6}}; _int_elem_mat = {mats[n1+1],mats[n3+1],matrix};}
							else err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
						}
						else if (n2e3)
						{
							if (n0<0)
								{_cutEdge = {3,1,1,0}; _enrichment_on_inclusion = {n2,n2,n1,n1}; _int_elem_struct = {{2,3,4,5},{1,6,7},{0,7,6,5,4}}; _int_elem_mat = {mats[n2+1],mats[n1+1],matrix};}
							else if (n1<0)
								{_cutEdge = {3,1,0,3}; _enrichment_on_inclusion = {n2,n2,n0,n0}; _int_elem_struct = {{2,3,4,5},{0,6,7},{1,5,4,7,6}}; _int_elem_mat = {mats[n2+1],mats[n0+1],matrix};}
							else err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
						}
						else if (n0e3)
						{
							if (n1<0)
								{_cutEdge = {0,2,2,1}; _enrichment_on_inclusion = {n3,n3,n2,n2}; _int_elem_struct = {{3,0,4,5},{2,6,7},{1,7,6,5,4}}; _int_elem_mat = {mats[n3+1],mats[n2+1],matrix};}
							else if (n2<0)
								{_cutEdge = {0,2,1,0}; _enrichment_on_inclusion = {n3,n3,n1,n1}; _int_elem_struct = {{3,0,4,5},{1,6,7},{2,5,4,7,6}}; _int_elem_mat = {mats[n3+1],mats[n1+1],matrix};}
							else err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
						}
						else
							err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
					}
				} // End contains matrix

				// Cut by 3 interface (one of them covers 2 of the nodes)
				else
				{
					if (n0e1)
						{_cutEdge = {1,3,3,2,2,1}; _enrichment_on_inclusion = {n0,n0,n3,n3,n2,n2}; _int_elem_struct = {{0,1,4,5},{3,6,7},{2,8,9},{9,8,7,6,5,4}}; _int_elem_mat = {mats[n0+1],mats[n3+1],mats[n2+1],matrix};}
					else if (n1e2)
						{_cutEdge = {2,0,0,3,3,2}; _enrichment_on_inclusion = {n1,n1,n0,n0,n3,n3}; _int_elem_struct = {{1,2,4,5},{0,6,7},{3,8,9},{9,8,7,6,5,4}}; _int_elem_mat = {mats[n1+1],mats[n0+1],mats[n3+1],matrix};}
					else if (n2e3)
						{_cutEdge = {3,1,1,0,0,3}; _enrichment_on_inclusion = {n2,n2,n1,n1,n0,n0}; _int_elem_struct = {{2,3,4,5},{1,6,7},{0,8,9},{9,8,7,6,5,4}}; _int_elem_mat = {mats[n2+1],mats[n1+1],mats[n0+1],matrix};}
					else if (n0e3)
						{_cutEdge = {0,2,2,1,1,0}; _enrichment_on_inclusion = {n3,n3,n2,n2,n1,n1}; _int_elem_struct = {{3,0,4,5},{2,6,7},{1,8,9},{9,8,7,6,5,4}}; _int_elem_mat = {mats[n3+1],mats[n2+1],mats[n1+1],matrix};}
					else
						err_message("Unknown manner of 3 inclusions intersecting an element.");
				}
			} // End 3 unique nodes



			// 4 unique nodal detection values
			else
			{
				// One of the nodal detection points is in the matrix which means we have 3 inclusions and the matrix here
				if (vec[0] < 0)
				{
					if (n0 < 0)			// Node 0 is in the matrix
						{_cutEdge = {3,2,2,1,1,0}; _enrichment_on_inclusion = {n3,n3,n2,n2,n1,n1}; _int_elem_struct = {{3,4,5},{2,6,7},{1,8,9},{0,9,8,7,6,5,4}}; _int_elem_mat = {mats[n3+1],mats[n2+1],mats[n1+1],matrix};}
					else if (n1 < 0)	// Node 1 is in the matrix
						{_cutEdge = {0,3,3,2,2,1}; _enrichment_on_inclusion = {n0,n0,n3,n3,n2,n2}; _int_elem_struct = {{0,4,5},{3,6,7},{2,8,9},{1,9,8,7,6,5,4}}; _int_elem_mat = {mats[n0+1],mats[n3+1],mats[n2+1],matrix};}
					else if (n2 < 0)	// Node 2 is in the matrix
						{_cutEdge = {1,0,0,3,3,2}; _enrichment_on_inclusion = {n1,n1,n0,n0,n3,n3}; _int_elem_struct = {{1,4,5},{0,6,7},{3,8,9},{2,9,8,7,6,5,4}}; _int_elem_mat = {mats[n1+1],mats[n0+1],mats[n3+1],matrix};}
					else if (n3 < 0)	// Node 3 is in the matrix
						{_cutEdge = {2,1,1,0,0,3}; _enrichment_on_inclusion = {n2,n2,n1,n1,n0,n0}; _int_elem_struct = {{2,4,5},{1,6,7},{0,8,9},{3,9,8,7,6,5,4}}; _int_elem_mat = {mats[n2+1],mats[n1+1],mats[n0+1],matrix};}
					else
						err_message("Unknown manner of 3 inclusions and the matrix intersecting the element.");
				}

				// Otherwise 4 inclusions cut this element
				else
				{
					_cutEdge = {3,2,2,1,1,0,0,3};
					_enrichment_on_inclusion = {n3,n3,n2,n2,n1,n1,n0,n0};
					_int_elem_struct = {{3,4,5},{2,6,7},{1,8,9},{0,10,11},{11,10,9,8,7,6,5,4}};
					_int_elem_mat = {mats[n3+1],mats[n2+1],mats[n1+1],mats[n0+1],matrix};
				}
			}

			// Not cohesive so we can clear this
			_enrichment_nodes.resize(_cutEdge.size());
			_coh_elem_struct.clear();
		}




		// Cohesive
		else
		{
			// 2 unique nodal detection values
			if (n_unique_nodes == 2)
			{
				// One of the materials is the matrix (Only have one inclusions cutting this element)
				if (vec[0] < 0)
				{
					// Three nodes are located in one material and one in another
					if((set_counts[0]==1) || (set_counts[0]==3))
					{
						if(n1e2 && n1e3)		// node zero is separate
						{
							_cutEdge = {3,0};
							_int_elem_mat = {mats[n0+1], mats[n2+1]};
							if(n0 < n1) {_int_elem_struct = {{0,5,4},{1,2,3,6,7}}; _coh_elem_struct = {{4,5,6,7}};}
							else {_int_elem_struct = {{0,7,6},{1,2,3,4,5}}; _coh_elem_struct = {{5,4,7,6}};}
						}
						else if(n0e2 && n0e3)	// node one is separate
						{
							_cutEdge = {0,1};
							_int_elem_mat = {mats[n1+1], mats[n3+1]}; 
							if(n1 < n0) {_int_elem_struct = {{1,5,4},{2,3,0,6,7}}; _coh_elem_struct = {{4,5,6,7}};}
							else {_int_elem_struct = {{1,7,6},{2,3,0,4,5}}; _coh_elem_struct = {{5,4,7,6}};}
						}
						else if(n0e1 && n0e3)	// node two is separate
						{
							_cutEdge = {1,2};
							_int_elem_mat = {mats[n2+1],mats[n0+1]};
							if(n2 < n0) {_int_elem_struct = {{2,5,4},{3,0,1,6,7}}; _coh_elem_struct = {{4,5,6,7}};}
							else {_int_elem_struct = {{2,7,6},{3,0,1,4,5}}; _coh_elem_struct = {{5,4,7,6}};}
						}
						else if(n0e1 && n1e2)		// node three is separate
						{
							_cutEdge = {2,3};
							_int_elem_mat = {mats[n3+1],mats[n1+1]};
							if(n3 < n0) {_int_elem_struct = {{3,5,4},{0,1,2,6,7}}; _coh_elem_struct = {{4,5,6,7}};}
							else {_int_elem_struct = {{3,7,6},{0,1,2,4,5}}; _coh_elem_struct = {{5,4,7,6}};}
						}
						else
							err_message("QUADRILATERAL ELEMENT HAS 1 NODE SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
					}
					
					// Two cases:
					//	1. Two adacent nodes are in one material and the other two nodes are in another, this splitting this into 2 Quad4's-direction
					//	2. Two opposite nodes are in the same material and the other two nodes are in the same material. This means the interface is highly curved and an error should be thrown (will later lead to mesh refinement)
					//		- This could also be handled similar to how 2 inclusions is handled
					else if(set_counts[0]==2)
					{
						if(n0e1)			// {0,1}/{2,3}
						{
							_cutEdge = {1, 3};
							_int_elem_mat = {mats[n0+1], mats[n2+1]};
							if(n0 < n2) {_int_elem_struct = {{0,1,4,5},{2,3,7,6}}; _coh_elem_struct = {{5,4,7,6}};}
							else {_int_elem_struct = {{0,1,6,7},{2,3,5,4}}; _coh_elem_struct = {{4,5,6,7}};}
						}
						else if(n1e2)		// {0,3}/{1,2}
						{
							_cutEdge = {2, 0};
							_int_elem_mat = {mats[n0+1], mats[n1+1]};
							if(n0 < n1) {_int_elem_struct = {{3,0,5,4},{1,2,6,7}}; _coh_elem_struct = {{4,5,6,7}};}
							else {_int_elem_struct = {{3,0,7,6},{1,2,4,5}}; _coh_elem_struct = {{5,4,7,6}};}
						}
						// FIXME: I would like to add 2 extra cases here
						else				// Interface is highly curved
							err_message("QUADRILATERAL ELEMENT HAS 2 NODES SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
					}
				
					else
						err_message("Quad4 element detector: We'll never get here! (Hopefully)");

					_enrichment_on_inclusion.resize(_cutEdge.size());
					std::fill(_enrichment_on_inclusion.begin(), _enrichment_on_inclusion.end(), *std::max_element(vec.begin(), vec.end())); // All intersected by the same inclusion
				}




				// Otherwise, there are two inclusions cutting this element
				else
				{
					// Three nodes are located in one material and one in another
					if((set_counts[0]==1) || (set_counts[0]==3))
					{
						if(n1e2 && n1e3)		// node zero is separate
							{_cutEdge = {0,3,3,0}; _enrichment_on_inclusion = {n0,n0,n2,n2}; _int_elem_struct = {{0,8,9},{5,4,7,6},{1,2,3,10,11}}; _int_elem_mat = {mats[n0+1],matrix,mats[n2+1]};}
						else if(n0e2 && n0e3)	// node one is separate
							{_cutEdge = {1,0,0,1}; _enrichment_on_inclusion = {n1,n1,n3,n3}; _int_elem_struct = {{1,8,9},{5,4,7,6},{2,3,0,10,11}}; _int_elem_mat = {mats[n1+1],matrix,mats[n3+1]};}
						else if(n0e1 && n0e3)	// node two is separate
							{_cutEdge = {2,1,1,2}; _enrichment_on_inclusion = {n2,n2,n0,n0}; _int_elem_struct = {{2,8,9},{5,4,7,6},{3,0,1,10,11}}; _int_elem_mat = {mats[n2+1],matrix,mats[n0+1]};}
						else if(n0e1 && n1e2)		// node three is separate
							{_cutEdge = {3,2,2,3}; _enrichment_on_inclusion = {n3,n3,n1,n1}; _int_elem_struct = {{3,8,9},{5,4,7,6},{0,1,2,10,11}}; _int_elem_mat = {mats[n3+1],matrix,mats[n1+1]};}
						else
							err_message("QUADRILATERAL ELEMENT HAS 1 NODE SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
						// Cohesive structure is the same for all of these cases
						_coh_elem_struct = {{4,5,8,9}, {6,7,10,11}};
					}
					
					// Two cases:
					//	1. Two adacent nodes are in one material and the other two nodes are in another, this splitting this into 2 Quad4's-direction
					//	2. Two opposite nodes are in the same material and the other two nodes are in the same material. This means the interface is highly curved and an error should be thrown (will later lead to mesh refinement)
					else if(set_counts[0]==2)
					{
						if(n0e1)			// {0,1}/{2,3}
							{_cutEdge = {1,3,3,1}; _enrichment_on_inclusion = {n0,n0,n2,n2}; _int_elem_struct = {{0,1,8,9},{5,4,7,6},{2,3,10,11}}; _int_elem_mat = {mats[n0+1],matrix,mats[n2+1]};}
						else if(n1e2)		// {0,3}/{1,2}
							{_cutEdge = {2,0,0,2}; _enrichment_on_inclusion = {n1,n1,n3,n3}; _int_elem_struct = {{1,2,8,9},{5,4,7,6},{3,0,10,11}}; _int_elem_mat = {mats[n1+1],matrix,mats[n3+1]};}
						else				// Interface is highly curved
							err_message("QUADRILATERAL ELEMENT HAS 2 NODES SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
						// Cohesive structure is the same for all of these cases
						_coh_elem_struct = {{4,5,8,9}, {6,7,10,11}};
					}
					
					else
						err_message("Quad4 element detector: We'll never get here! (Hopefully)");
				}
			} // End 2 unique values



			// 3 unique nodal detection values
			else if (n_unique_nodes == 3)
			{
				// One of the materials is the matrix (2 options)
				if (vec[0] < 0)
				{
					// Case where 2 of the nodes are in the matrix and the remaining 2 nodes are in different inclusions
					if (set_counts[0] == 2)
					{
						if (n0e1)		// Nodes 0 and 1 are in the matrix
							{_cutEdge = {3,2,2,1}; _enrichment_on_inclusion = {n3,n3,n2,n2}; _int_elem_struct = {{3,8,9},{2,10,11},{0,1,7,6,5,4}}; _int_elem_mat = {mats[n3+1],mats[n2+1],matrix};}
						else if (n1e2)	// Nodes 1 and 2 are in the matrix
							{_cutEdge = {0,3,3,2}; _enrichment_on_inclusion = {n0,n0,n3,n3}; _int_elem_struct = {{0,8,9},{3,10,11},{1,2,7,6,5,4}}; _int_elem_mat = {mats[n0+1],mats[n3+1],matrix};}
						else if (n2e3)	// Nodes 2 and 3 are in the matrix
							{_cutEdge = {1,0,0,3}; _enrichment_on_inclusion = {n1,n1,n0,n0}; _int_elem_struct = {{1,8,9},{0,10,11},{2,3,7,6,5,4}}; _int_elem_mat = {mats[n1+1],mats[n0+1],matrix};}
						else if (n0e3)	// Nodes 3 and 0 are in the matrix
							{_cutEdge = {2,1,1,0}; _enrichment_on_inclusion = {n2,n2,n1,n1}; _int_elem_struct = {{2,8,9},{1,10,11},{3,0,7,6,5,4}}; _int_elem_mat = {mats[n2+1],mats[n1+1],matrix};}
						else if (n0e2)	// Nodes 0 and 2 are in the matrix
							{_cutEdge = {1,0,3,2}; _enrichment_on_inclusion = {n1,n1,n3,n3}; _int_elem_struct = {{1,8,9},{3,10,11},{2,7,6,0,5,4}}; _int_elem_mat = {mats[n1+1],mats[n3+1],matrix};}
						else if (n1e3)	// Nodes 1 and 3 are in the matrix
							{_cutEdge = {0,3,2,1}; _enrichment_on_inclusion = {n0,n0,n2,n2}; _int_elem_struct = {{0,8,9},{2,10,11},{1,7,6,3,5,4}}; _int_elem_mat = {mats[n0+1],mats[n2+1],matrix};}
						else
							err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
					}

					// Otherwise, the matrix is only on one of the nodes and two of the other nodes are covered by the same inclusion
					else
					{
						if (n0e1)
						{
							if (n2<0)
								{_cutEdge = {1,3,3,2}; _enrichment_on_inclusion = {n0,n0,n3,n3}; _int_elem_struct = {{0,1,8,9},{3,10,11},{2,7,6,5,4}}; _int_elem_mat = {mats[n0+1],mats[n3+1],matrix};}
							else if (n3<0)
								{_cutEdge = {1,3,2,1}; _enrichment_on_inclusion = {n0,n0,n2,n2}; _int_elem_struct = {{0,1,8,9},{2,10,11},{3,5,4,7,6}}; _int_elem_mat = {mats[n0+1],mats[n2+1],matrix};}
							else err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
						}
						else if (n1e2)
						{
							if (n3<0)
								{_cutEdge = {2,0,0,3}; _enrichment_on_inclusion = {n1,n1,n0,n0}; _int_elem_struct = {{1,2,8,9},{0,10,11},{3,7,6,5,4}}; _int_elem_mat = {mats[n1+1],mats[n0+1],matrix};}
							else if (n0<0)
								{_cutEdge = {2,0,3,2}; _enrichment_on_inclusion = {n1,n1,n3,n3}; _int_elem_struct = {{1,2,8,9},{3,10,11},{0,5,4,7,6}}; _int_elem_mat = {mats[n1+1],mats[n3+1],matrix};}
							else err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
						}
						else if (n2e3)
						{
							if (n0<0)
								{_cutEdge = {3,1,1,0}; _enrichment_on_inclusion = {n2,n2,n1,n1}; _int_elem_struct = {{2,3,8,9},{1,10,11},{0,7,6,5,4}}; _int_elem_mat = {mats[n2+1],mats[n1+1],matrix};}
							else if (n1<0)
								{_cutEdge = {3,1,0,3}; _enrichment_on_inclusion = {n2,n2,n0,n0}; _int_elem_struct = {{2,3,8,9},{0,10,11},{1,5,4,7,6}}; _int_elem_mat = {mats[n2+1],mats[n0+1],matrix};}
							else err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
						}
						else if (n0e3)
						{
							if (n1<0)
								{_cutEdge = {0,2,2,1}; _enrichment_on_inclusion = {n3,n3,n2,n2}; _int_elem_struct = {{3,0,8,9},{2,10,11},{1,7,6,5,4}}; _int_elem_mat = {mats[n3+1],mats[n2+1],matrix};}
							else if (n2<0)
								{_cutEdge = {0,2,1,0}; _enrichment_on_inclusion = {n3,n3,n1,n1}; _int_elem_struct = {{3,0,8,9},{1,10,11},{2,5,4,7,6}}; _int_elem_mat = {mats[n3+1],mats[n1+1],matrix};}
							else err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
						}
						else
							err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
					}

					// All cohesive structures are the same here
					_coh_elem_struct = {{4,5,8,9}, {6,7,10,11}};
				}


				// Cut by 3 interface (one of them covers 2 of the nodes)
				else
				{
					if (n0e1)
						{_cutEdge = {1,3,3,2,2,1}; _enrichment_on_inclusion = {n0,n0,n3,n3,n2,n2}; _int_elem_struct = {{0,1,10,11},{3,12,13},{2,14,15},{9,8,7,6,5,4}}; _int_elem_mat = {mats[n0+1],mats[n3+1],mats[n2+1],matrix};}
					else if (n1e2)
						{_cutEdge = {2,0,0,3,3,2}; _enrichment_on_inclusion = {n1,n1,n0,n0,n3,n3}; _int_elem_struct = {{1,2,10,11},{0,12,13},{3,14,15},{9,8,7,6,5,4}}; _int_elem_mat = {mats[n1+1],mats[n0+1],mats[n3+1],matrix};}
					else if (n2e3)
						{_cutEdge = {3,1,1,0,0,3}; _enrichment_on_inclusion = {n2,n2,n1,n1,n0,n0}; _int_elem_struct = {{2,3,10,11},{1,12,13},{0,14,15},{9,8,7,6,5,4}}; _int_elem_mat = {mats[n2+1],mats[n1+1],mats[n0+1],matrix};}
					else if (n0e3)
						{_cutEdge = {0,2,2,1,1,0}; _enrichment_on_inclusion = {n3,n3,n2,n2,n1,n1}; _int_elem_struct = {{3,0,10,11},{2,12,13},{1,14,15},{9,8,7,6,5,4}}; _int_elem_mat = {mats[n3+1],mats[n2+1],mats[n1+1],matrix};}
					else
						err_message("Unknown manner of 3 inclusions intersecting an element.");

					_coh_elem_struct = {{4,5,10,11},{6,7,12,13},{8,9,14,15}};
				}
			} // End 3 unique nodes



			// 4 unique nodal detection values
			else
			{
				// One of the nodal detection points is in the matrix which means we have 3 inclusions and the matrix here
				if (vec[0] < 0)
				{
					if (n0 < 0)			// Node 0 is in the matrix
						{_cutEdge = {3,2,2,1,1,0}; _enrichment_on_inclusion = {n3,n3,n2,n2,n1,n1}; _int_elem_struct = {{3,10,11},{2,12,13},{1,14,15},{0,9,8,7,6,5,4}}; _int_elem_mat = {mats[n3+1],mats[n2+1],mats[n1+1],matrix};}
					else if (n1 < 0)	// Node 1 is in the matrix
						{_cutEdge = {0,3,3,2,2,1}; _enrichment_on_inclusion = {n0,n0,n3,n3,n2,n2}; _int_elem_struct = {{0,10,11},{3,12,13},{2,14,15},{1,9,8,7,6,5,4}}; _int_elem_mat = {mats[n0+1],mats[n3+1],mats[n2+1],matrix};}
					else if (n2 < 0)	// Node 2 is in the matrix
						{_cutEdge = {1,0,0,3,3,2}; _enrichment_on_inclusion = {n1,n1,n0,n0,n3,n3}; _int_elem_struct = {{1,10,11},{0,12,13},{3,14,15},{2,9,8,7,6,5,4}}; _int_elem_mat = {mats[n1+1],mats[n0+1],mats[n3+1],matrix};}
					else if (n3 < 0)	// Node 3 is in the matrix
						{_cutEdge = {2,1,1,0,0,3}; _enrichment_on_inclusion = {n2,n2,n1,n1,n0,n0}; _int_elem_struct = {{2,10,11},{1,12,13},{0,14,15},{3,9,8,7,6,5,4}}; _int_elem_mat = {mats[n2+1],mats[n1+1],mats[n0+1],matrix};}
					else
						err_message("Unknown manner of 3 inclusions and the matrix intersecting the element.");

					_coh_elem_struct = {{4,5,10,11},{6,7,12,13},{8,9,14,15}};
				}

				// Otherwise 4 inclusions cut this element
				else
				{
					_cutEdge = {3,2,2,1,1,0,0,3};
					_enrichment_on_inclusion = {n3,n3,n2,n2,n1,n1,n0,n0};
					_int_elem_struct = {{3,12,13},{2,14,15},{1,16,17},{0,18,19},{11,10,9,8,7,6,5,4}};
					_int_elem_mat = {mats[n3+1],mats[n2+1],mats[n1+1],mats[n0+1],matrix};
				}
			} // End 4 unique nodes

			_enrichment_nodes.resize(_cutEdge.size() * 2);
			_coh_mat.resize(_coh_elem_struct.size());
		} // End is cohesive
	} // End is intersected













	// ==============================================================================
	// CASES WRITTEN NOT USING POLYGONAL INTEGRATION ELEMENTS
	// ==============================================================================

/*

	// This element is intersected somehow
	else
	{
		_is_intersected = true;

		// Not cohesive cases
		if (!isCohesive)
		{
			// 2 unique nodal detection values
			if (n_unique_nodes == 2)
			{
				// One of the materials is the matrix (Only have one inclusions cutting this element)
				if (vec[0] < 0)
				{
					// Three nodes are located in one material and one in another
					if((set_counts[0]==1) || (set_counts[0]==3))
					{
						if(n1e2 && n1e3)		// node zero is separate
							{_cutEdge = {3, 0}; _int_elem_struct = {{0, 5, 4}, {1, 2, 5}, {2, 3, 4, 5}}; _int_elem_mat = {mats[n0+1], mats[n1+1], mats[n2+1]};}
						else if(n0e2 && n0e3)	// node one is separate
							{_cutEdge = {0, 1}; _int_elem_struct = {{1, 5, 4}, {2, 3, 5}, {3, 0, 4, 5}}; _int_elem_mat = {mats[n1+1], mats[n2+1], mats[n3+1]};}
						else if(n0e1 && n0e3)	// node two is separate
							{_cutEdge = {1, 2}; _int_elem_struct = {{2, 5, 4}, {3, 0, 5}, {0, 1, 4, 5}}; _int_elem_mat = {mats[n2+1], mats[n3+1], mats[n0+1]};}
						else if(n0e1 && n1e2)		// node three is separate
							{_cutEdge = {2, 3}; _int_elem_struct = {{3, 5, 4}, {0, 1, 5}, {1, 2, 4, 5}}; _int_elem_mat = {mats[n3+1], mats[n0+1], mats[n1+1]};}
						else
							err_message("QUADRILATERAL ELEMENT HAS 1 NODE SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
					}
					
					// Two cases:
					//	1. Two adacent nodes are in one material and the other two nodes are in another, this splitting this into 2 Quad4's-direction
					//	2. Two opposite nodes are in the same material and the other two nodes are in the same material. This means the interface is highly curved and an error should be thrown (will later lead to mesh refinement)
					else if(set_counts[0]==2)
					{
						if(n0e1)			// {0,1}/{2,3}
							{_cutEdge = {1, 3}; _int_elem_struct = {{0, 1, 4, 5}, {2, 3, 5, 4}}; _int_elem_mat = {mats[n0+1], mats[n2+1]};}
						else if(n1e2)		// {0,3}/{1,2}
							{_cutEdge = {2, 0}; _int_elem_struct = {{1, 2, 4, 5}, {3, 0, 5, 4}}; _int_elem_mat = {mats[n1+1], mats[n3+1]};}
						else				// Interface is highly curved
							err_message("QUADRILATERAL ELEMENT HAS 2 NODES SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
					}
					
					else
						err_message("Quad4 element detector: We'll never get here! (Hopefully)");

					_enrichment_on_inclusion.resize(_cutEdge.size());
					std::fill(_enrichment_on_inclusion.begin(), _enrichment_on_inclusion.end(), *std::max_element(vec.begin(), vec.end())); // All intersected by the same inclusion
				}



				// Otherwise there are two inclusions cutting this element
				else
				{
					// Three nodes are located in one material and one in another
					if((set_counts[0]==1) || (set_counts[0]==3))
					{
						if(n1e2 && n1e3)		// node zero is separate
							{_cutEdge = {0,3,3,0}; _enrichment_on_inclusion = {n0,n0,n2,n2}; _int_elem_struct = {{0,4,5},{5,4,7,6},{1,2,7},{2,6,7},{2,3,6}}; _int_elem_mat = {mats[n0+1],matrix,mats[n2+2],mats[n2+1],mats[n2+1]};}
						else if(n0e2 && n0e3)	// node one is separate
							{_cutEdge = {1,0,0,1}; _enrichment_on_inclusion = {n1,n1,n3,n3}; _int_elem_struct = {{1,4,5},{5,4,7,6},{2,3,7},{3,6,7},{3,0,6}}; _int_elem_mat = {mats[n1+1],matrix,mats[n3+2],mats[n3+1],mats[n3+1]};}
						else if(n0e1 && n0e3)	// node two is separate
							{_cutEdge = {2,1,1,2}; _enrichment_on_inclusion = {n2,n2,n0,n0}; _int_elem_struct = {{2,4,5},{5,4,7,6},{3,0,7},{0,6,7},{0,1,6}}; _int_elem_mat = {mats[n2+1],matrix,mats[n0+2],mats[n0+1],mats[n0+1]};}
						else if(n0e1 && n1e2)		// node three is separate
							{_cutEdge = {3,2,2,3}; _enrichment_on_inclusion = {n3,n3,n1,n1}; _int_elem_struct = {{3,4,5},{5,4,7,6},{0,1,7},{1,6,7},{1,2,6}}; _int_elem_mat = {mats[n3+1],matrix,mats[n1+2],mats[n1+1],mats[n1+1]};}
						else
							err_message("QUADRILATERAL ELEMENT HAS 1 NODE SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
					}
					
					// Two cases:
					//	1. Two adacent nodes are in one material and the other two nodes are in another, this splitting this into 2 Quad4's-direction
					//	2. Two opposite nodes are in the same material and the other two nodes are in the same material. This means the interface is highly curved and an error should be thrown (will later lead to mesh refinement)
					else if(set_counts[0]==2)
					{
						if(n0e1)			// {0,1}/{2,3}
							{_cutEdge = {1,3,3,1}; _enrichment_on_inclusion = {n0,n0,n2,n2}; _int_elem_struct = {{0,1,4,5},{5,4,7,6},{2,3,6,7}}; _int_elem_mat = {mats[n0+1],matrix,mats[n2+1]};}
						else if(n1e2)		// {0,3}/{1,2}
							{_cutEdge = {2,0,0,2}; _enrichment_on_inclusion = {n1,n1,n3,n3}; _int_elem_struct = {{1,2,4,5},{5,4,7,6},{3,0,6,7}}; _int_elem_mat = {mats[n1+1],matrix,mats[n3+1]};}
						else				// Interface is highly curved
							err_message("QUADRILATERAL ELEMENT HAS 2 NODES SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
					}
					
					else
						err_message("Quad4 element detector: We'll never get here! (Hopefully)");


					_enrichment_nodes.resize(4);
				}
			} // End 2 unique nodes



			// 3 unique nodal detection values
			else if (n_unique_nodes == 3)
			{
				// One of the materials is the matrix (2 options)
				if (vec[0] < 0)
				{
					// Case where 2 of the nodes are in the matrix and the remaining 2 nodes are in different inclusions
					if (set_counts[0] == 2)
					{
						if (n0e1)		// Nodes 0 and 1 are in the matrix
							{_cutEdge = {3,2,2,1}; _enrichment_on_inclusion = {n3,n3,n2,n2}; _int_elem_struct = {{3,4,5},{2,6,7},{0,1,7,4},{5,4,7,6}}; _int_elem_mat = {mats[n3+1],mats[n2+1],matrix,matrix};}
						else if (n1e2)	// Nodes 1 and 2 are in the matrix
							{_cutEdge = {0,3,3,2}; _enrichment_on_inclusion = {n0,n0,n3,n3}; _int_elem_struct = {{0,4,5},{3,6,7},{1,2,7,4},{5,4,7,6}}; _int_elem_mat = {mats[n0+1],mats[n3+1],matrix,matrix};}
						else if (n2e3)	// Nodes 2 and 3 are in the matrix
							{_cutEdge = {1,0,0,3}; _enrichment_on_inclusion = {n1,n1,n0,n0}; _int_elem_struct = {{1,4,5},{0,6,7},{2,3,7,4},{5,4,7,6}}; _int_elem_mat = {mats[n1+1],mats[n0+1],matrix,matrix};}
						else if (n0e3)	// Nodes 3 and 0 are in the matrix
							{_cutEdge = {2,1,1,0}; _enrichment_on_inclusion = {n2,n2,n1,n1}; _int_elem_struct = {{2,4,5},{1,6,7},{3,0,7,4},{5,4,7,6}}; _int_elem_mat = {mats[n2+1],mats[n1+1],matrix,matrix};}
						else if (n1e3)	// Nodes 1 and 3 are in the matrix
							{_cutEdge = {0,3,2,1}; _enrichment_on_inclusion = {n0,n0,n2,n2}; _int_elem_struct = {{0,4,5},{2,6,7},{3,5,6},{5,4,7,6},{1,7,4}}; _int_elem_mat = {mats[n0+1],mats[n2+1],matrix,matrix,matrix};}
						else if (n0e2)	// Nodes 0 and 2 are in the matrix
							{_cutEdge = {1,0,3,2}; _enrichment_on_inclusion = {n1,n1,n3,n3}; _int_elem_struct = {{1,4,5},{3,6,7},{0,5,6},{5,4,7,6},{2,7,4}}; _int_elem_mat = {mats[n1+1],mats[n3+1],matrix,matrix,matrix};}
						else
							err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
					}

					// Otherwise, the matrix is only on one of the nodes and two of the other nodes are covered by the same inclusion
					else
					{
						if (n0e1)
						{
							if (n2<0)
								{_cutEdge = {1,3,3,2}; _enrichment_on_inclusion = {n0,n0,n3,n3}; _int_elem_struct = {{0,1,4,5},{3,6,7},{2,7,4},{5,4,7,6}}; _int_elem_mat = {mats[n0+1],mats[n3+1],matrix,matrix};}
							else if (n3<0)
								{_cutEdge = {1,3,2,1}; _enrichment_on_inclusion = {n0,n0,n2,n2}; _int_elem_struct = {{0,1,4,5},{2,6,7},{3,5,6},{5,4,7,6}}; _int_elem_mat = {mats[n0+1],mats[n2+1],matrix,matrix};}
							else err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
						}
						else if (n1e2)
						{
							if (n3<0)
								{_cutEdge = {2,0,0,3}; _enrichment_on_inclusion = {n1,n1,n0,n0}; _int_elem_struct = {{1,2,4,5},{0,6,7},{3,7,4},{5,4,7,6}}; _int_elem_mat = {mats[n1+1],mats[n0+1],matrix,matrix};}
							else if (n0<0)
								{_cutEdge = {2,0,3,2}; _enrichment_on_inclusion = {n1,n1,n3,n3}; _int_elem_struct = {{1,2,4,5},{3,6,7},{0,5,6},{5,4,7,6}}; _int_elem_mat = {mats[n1+1],mats[n3+1],matrix,matrix};}
							else err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
						}
						else if (n2e3)
						{
							if (n0<0)
								{_cutEdge = {3,1,1,0}; _enrichment_on_inclusion = {n2,n2,n1,n1}; _int_elem_struct = {{2,3,4,5},{1,6,7},{0,7,4},{5,4,7,6}}; _int_elem_mat = {mats[n2+1],mats[n1+1],matrix,matrix};}
							else if (n1<0)
								{_cutEdge = {3,1,0,3}; _enrichment_on_inclusion = {n2,n2,n0,n0}; _int_elem_struct = {{2,3,4,5},{0,6,7},{1,5,6},{5,4,7,6}}; _int_elem_mat = {mats[n2+1],mats[n0+1],matrix,matrix};}
							else err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
						}
						else if (n0e3)
						{
							if (n1<0)
								{_cutEdge = {0,2,2,1}; _enrichment_on_inclusion = {n3,n3,n2,n2}; _int_elem_struct = {{3,0,4,5},{2,6,7},{1,7,4},{5,4,7,6}}; _int_elem_mat = {mats[n3+1],mats[n2+1],matrix,matrix};}
							else if (n2<0)
								{_cutEdge = {0,2,1,0}; _enrichment_on_inclusion = {n3,n3,n1,n1}; _int_elem_struct = {{3,0,4,5},{1,6,7},{2,5,6},{5,4,7,6}}; _int_elem_mat = {mats[n3+1],mats[n1+1],matrix,matrix};}
							else err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
						}
						else
							err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
					}
				} // End contains matrix

				// Cut by 3 interface (one of them covers 2 of the nodes)
				else
				{
					if (n0e1)
						{_cutEdge = {1,3,3,2,2,1}; _enrichment_on_inclusion = {n0,n0,n3,n3,n2,n2}; _int_elem_struct = {{0,1,4,5},{3,6,7},{2,8,9},{5,7,6},{4,9,8},{5,4,8,7}}; _int_elem_mat = {mats[n0+1],mats[n3+1],mats[n2+1],matrix,matrix,matrix};}
					else if (n1e2)
						{_cutEdge = {2,0,0,3,3,2}; _enrichment_on_inclusion = {n1,n1,n0,n0,n3,n3}; _int_elem_struct = {{1,2,4,5},{0,6,7},{3,8,9},{5,7,6},{4,9,8},{5,4,8,7}}; _int_elem_mat = {mats[n1+1],mats[n0+1],mats[n3+1],matrix,matrix,matrix};}
					else if (n2e3)
						{_cutEdge = {3,1,1,0,0,3}; _enrichment_on_inclusion = {n2,n2,n1,n1,n0,n0}; _int_elem_struct = {{2,3,4,5},{1,6,7},{0,8,9},{5,7,6},{4,9,8},{5,4,8,7}}; _int_elem_mat = {mats[n2+1],mats[n1+1],mats[n0+1],matrix,matrix,matrix};}
					else if (n0e3)
						{_cutEdge = {0,2,2,1,1,0}; _enrichment_on_inclusion = {n3,n3,n2,n2,n1,n1}; _int_elem_struct = {{3,0,4,5},{2,6,7},{1,8,9},{5,7,6},{4,9,8},{5,4,8,7}}; _int_elem_mat = {mats[n3+1],mats[n2+1],mats[n1+1],matrix,matrix,matrix};}
					else
						err_message("Unknown manner of 3 inclusions intersecting an element.");
				}
			} // End 3 unique nodes



			// 4 unique nodal detection values
			else
			{
				// One of the nodal detection points is in th ematrix which means we have 3 inclusions and the matrix here
				if (vec[0] < 0)
				{
					if (n0 < 0)			// Node 0 is in the matrix
						{_cutEdge = {3,2,2,1,1,0}; _enrichment_on_inclusion = {n3,n3,n2,n2,n1,n1}; _int_elem_struct = {{3,4,5},{2,6,7},{1,8,9},{0,9,8,7},{0,7,6},{0,6,5,4}}; _int_elem_mat = {mats[n3+1],mats[n2+1],mats[n1+1],matrix,matrix,matrix};}
					else if (n1 < 0)	// Node 1 is in the matrix
						{_cutEdge = {0,3,3,2,2,1}; _enrichment_on_inclusion = {n0,n0,n3,n3,n2,n2}; _int_elem_struct = {{0,4,5},{3,6,7},{2,8,9},{1,9,8,7},{1,7,6},{1,6,5,4}}; _int_elem_mat = {mats[n0+1],mats[n3+1],mats[n2+1],matrix,matrix,matrix};}
					else if (n2 < 0)	// Node 2 is in the matrix
						{_cutEdge = {1,0,0,3,3,2}; _enrichment_on_inclusion = {n1,n1,n0,n0,n3,n3}; _int_elem_struct = {{1,4,5},{0,6,7},{3,8,9},{2,9,8,7},{2,7,6},{2,6,5,4}}; _int_elem_mat = {mats[n1+1],mats[n0+1],mats[n3+1],matrix,matrix,matrix};}
					else if (n3 < 0)	// Node 3 is in the matrix
						{_cutEdge = {2,1,1,0,0,3}; _enrichment_on_inclusion = {n2,n2,n1,n1,n0,n0}; _int_elem_struct = {{2,4,5},{1,6,7},{0,8,9},{3,9,8,7},{3,7,6},{3,6,5,4}}; _int_elem_mat = {mats[n2+1],mats[n1+1],mats[n0+1],matrix,matrix,matrix};}
					else
						err_message("Unknown manner of 3 inclusions and the matrix intersecting the element.");
				}

				// Otherwise 4 inclusions cut this element
				else
				{
					_cutEdge = {3,2,2,1,1,0,0,3};
					_enrichment_on_inclusion = {n3,n3,n2,n2,n1,n1,n0,n0};
					_int_elem_struct = {{3,4,5},{2,6,7},{1,8,9},{0,10,11},{11,10,9,8},{11,8,7,4},{7,6,5,4}};
					_int_elem_mat = {mats[n3+1],mats[n2+1],mats[n1+1],mats[n0+1],matrix,matrix,matrix};
				}
			}

			// Not cohesive so we can clear this
			_enrichment_nodes.resize(_cutEdge.size());
			_coh_elem_struct.clear();
		}




		// Cohesive
		else
		{
			// 2 unique nodal detection values
			if (n_unique_nodes == 2)
			{
				// One of the materials is the matrix (Only have one inclusions cutting this element)
				if (vec[0] < 0)
				{
					// Three nodes are located in one material and one in another
					if((set_counts[0]==1) || (set_counts[0]==3))
					{
						if(n1e2 && n1e3)		// node zero is separate
						{
							_cutEdge = {3,0};
							_int_elem_mat = {mats[n0+1], mats[n1+1], mats[n2+1]};
							if(n0 < n1) {_int_elem_struct = {{0,5,4},{1,2,7},{2,3,6,7}}; _coh_elem_struct = {{4,5,6,7}};}
							else {_int_elem_struct = {{0,7,6},{1,2,5},{2,3,4,5}}; _coh_elem_struct = {{5,4,7,6}};}
						}
						else if(n0e2 && n0e3)	// node one is separate
						{
							_cutEdge = {0,1};
							_int_elem_mat = {mats[n1+1], mats[n2+1], mats[n3+1]};
							if(n1 < n2) {_int_elem_struct = {{1,5,4},{2,3,7},{3,0,6,7}}; _coh_elem_struct = {{4,5,6,7}};}
							else {_int_elem_struct = {{1,7,6},{2,3,5},{3,0,4,5}}; _coh_elem_struct = {{5,4,7,6}};}
						}
						else if(n0e1 && n0e3)	// node two is separate
						{
							_cutEdge = {1,2};
							_int_elem_mat = {mats[n2+1], mats[n3+1], mats[n0+1]};
							if(n2 < n3) {_int_elem_struct = {{2,5,4},{3,0,7},{0,1,6,7}}; _coh_elem_struct = {{4,5,6,7}};}
							else {_int_elem_struct = {{2,7,6},{3,0,5},{0,1,4,5}}; _coh_elem_struct = {{5,4,7,6}};}
						}
						else if(n0e1 && n1e2)		// node three is separate
						{
							_cutEdge = {2,3};
							_int_elem_mat = {mats[n3+1], mats[n0+1], mats[n1+1]};
							if(n3 < n0) {_int_elem_struct = {{3,5,4},{0,1,7},{1,2,6,7}}; _coh_elem_struct = {{4,5,6,7}};}
							else {_int_elem_struct = {{3,7,6},{0,1,5},{1,2,4,5}}; _coh_elem_struct = {{5,4,7,6}};}
						}
						else
							err_message("QUADRILATERAL ELEMENT HAS 1 NODE SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
					}
					
					// Two cases:
					//	1. Two adacent nodes are in one material and the other two nodes are in another, this splitting this into 2 Quad4's-direction
					//	2. Two opposite nodes are in the same material and the other two nodes are in the same material. This means the interface is highly curved and an error should be thrown (will later lead to mesh refinement)
					//		- This could also be handled similar to how 2 inclusions is handled
					else if(set_counts[0]==2)
					{
						if(n0e1)			// {0,1}/{2,3}
						{
							_cutEdge = {1, 3};
							_int_elem_mat = {mats[n0+1], mats[n2+1]};
							if(n0 < n2) {_int_elem_struct = {{0,1,4,5},{2,3,7,6}}; _coh_elem_struct = {{5,4,7,6}};}
							else {_int_elem_struct = {{0,1,6,7},{2,3,5,4}}; _coh_elem_struct = {{4,5,6,7}};}
						}
						else if(n1e2)		// {0,3}/{1,2}
						{
							_cutEdge = {2, 0};
							_int_elem_mat = {mats[n1+1], mats[n3+1]};
							if(n0 < n1) {_int_elem_struct = {{3,0,5,4},{1,2,6,7}}; _coh_elem_struct = {{4,5,6,7}};}
							else {_int_elem_struct = {{1,2,4,5},{3,0,7,6}}; _coh_elem_struct = {{5,4,7,6}};}
						}
						else				// Interface is highly curved
							err_message("QUADRILATERAL ELEMENT HAS 2 NODES SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
					}
				
					else
						err_message("Quad4 element detector: We'll never get here! (Hopefully)");

					_enrichment_on_inclusion.resize(_cutEdge.size());
					std::fill(_enrichment_on_inclusion.begin(), _enrichment_on_inclusion.end(), *std::max_element(vec.begin(), vec.end())); // All intersected by the same inclusion
				}




				// Otherwise, there are two inclusions cutting this element
				else
				{
					// Three nodes are located in one material and one in another
					if((set_counts[0]==1) || (set_counts[0]==3))
					{
						if(n1e2 && n1e3)		// node zero is separate
							{_cutEdge = {0,3,3,0}; _enrichment_on_inclusion = {n0,n0,n2,n2}; _int_elem_struct = {{0,8,9},{5,4,7,6},{1,2,11},{2,10,11},{2,3,10}}; _int_elem_mat = {mats[n0+1],matrix,mats[n2+2],mats[n2+1],mats[n2+1]};}
						else if(n0e2 && n0e3)	// node one is separate
							{_cutEdge = {1,0,0,1}; _enrichment_on_inclusion = {n1,n1,n3,n3}; _int_elem_struct = {{1,8,9},{5,4,7,6},{2,3,11},{3,10,11},{3,0,10}}; _int_elem_mat = {mats[n1+1],matrix,mats[n3+2],mats[n3+1],mats[n3+1]};}
						else if(n0e1 && n0e3)	// node two is separate
							{_cutEdge = {2,1,1,2}; _enrichment_on_inclusion = {n2,n2,n0,n0}; _int_elem_struct = {{2,8,9},{5,4,7,6},{3,0,11},{0,10,11},{0,1,10}}; _int_elem_mat = {mats[n2+1],matrix,mats[n0+2],mats[n0+1],mats[n0+1]};}
						else if(n0e1 && n1e2)		// node three is separate
							{_cutEdge = {3,2,2,3}; _enrichment_on_inclusion = {n3,n3,n1,n1}; _int_elem_struct = {{3,8,9},{5,4,7,6},{0,1,11},{1,10,11},{1,2,10}}; _int_elem_mat = {mats[n3+1],matrix,mats[n1+2],mats[n1+1],mats[n1+1]};}
						else
							err_message("QUADRILATERAL ELEMENT HAS 1 NODE SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
						// Cohesive structure is the same for all of these cases
						_coh_elem_struct = {{4,5,8,9}, {6,7,10,11}};
					}
					
					// Two cases:
					//	1. Two adacent nodes are in one material and the other two nodes are in another, this splitting this into 2 Quad4's-direction
					//	2. Two opposite nodes are in the same material and the other two nodes are in the same material. This means the interface is highly curved and an error should be thrown (will later lead to mesh refinement)
					else if(set_counts[0]==2)
					{
						if(n0e1)			// {0,1}/{2,3}
							{_cutEdge = {1,3,3,1}; _enrichment_on_inclusion = {n0,n0,n2,n2}; _int_elem_struct = {{0,1,8,9},{5,4,7,6},{2,3,10,11}}; _int_elem_mat = {mats[n0+1],matrix,mats[n2+1]};}
						else if(n1e2)		// {0,3}/{1,2}
							{_cutEdge = {2,0,0,2}; _enrichment_on_inclusion = {n1,n1,n3,n3}; _int_elem_struct = {{1,2,8,9},{5,4,7,6},{3,0,10,11}}; _int_elem_mat = {mats[n1+1],matrix,mats[n3+1]};}
						else				// Interface is highly curved
							err_message("QUADRILATERAL ELEMENT HAS 2 NODES SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
						// Cohesive structure is the same for all of these cases
						_coh_elem_struct = {{4,5,8,9}, {6,7,10,11}};
					}
					
					else
						err_message("Quad4 element detector: We'll never get here! (Hopefully)");
				}
			} // End 2 unique values



			// 3 unique nodal detection values
			else if (n_unique_nodes == 3)
			{
				// One of the materials is the matrix (2 options)
				if (vec[0] < 0)
				{
					// Case where 2 of the nodes are in the matrix and the remaining 2 nodes are in different inclusions
					if (set_counts[0] == 2)
					{
						if (n0e1)		// Nodes 0 and 1 are in the matrix
							{_cutEdge = {3,2,2,1}; _enrichment_on_inclusion = {n3,n3,n2,n2}; _int_elem_struct = {{3,8,9},{2,10,11},{0,1,7,4},{5,4,7,6}}; _int_elem_mat = {mats[n3+1],mats[n2+1],matrix,matrix};}
						else if (n1e2)	// Nodes 1 and 2 are in the matrix
							{_cutEdge = {0,3,3,2}; _enrichment_on_inclusion = {n0,n0,n3,n3}; _int_elem_struct = {{0,8,9},{3,10,11},{1,2,7,4},{5,4,7,6}}; _int_elem_mat = {mats[n0+1],mats[n3+1],matrix,matrix};}
						else if (n2e3)	// Nodes 2 and 3 are in the matrix
							{_cutEdge = {1,0,0,3}; _enrichment_on_inclusion = {n1,n1,n0,n0}; _int_elem_struct = {{1,8,9},{0,10,11},{2,3,7,4},{5,4,7,6}}; _int_elem_mat = {mats[n1+1],mats[n0+1],matrix,matrix};}
						else if (n0e3)	// Nodes 3 and 0 are in the matrix
							{_cutEdge = {2,1,1,0}; _enrichment_on_inclusion = {n2,n2,n1,n1}; _int_elem_struct = {{2,8,9},{1,10,11},{3,0,7,4},{5,4,7,6}}; _int_elem_mat = {mats[n2+1],mats[n1+1],matrix,matrix};}
						else if (n1e3)	// Nodes 0 and 2 are in the matrix
							{_cutEdge = {0,3,2,1}; _enrichment_on_inclusion = {n0,n0,n2,n2}; _int_elem_struct = {{0,8,9},{2,10,11},{3,5,6},{5,4,7,6},{1,7,4}}; _int_elem_mat = {mats[n0+1],mats[n2+1],matrix,matrix,matrix};}
						else if (n0e2)	// Nodes 1 and 3 are in the matrix
							{_cutEdge = {1,0,3,2}; _enrichment_on_inclusion = {n1,n1,n3,n3}; _int_elem_struct = {{1,8,9},{3,10,11},{0,5,6},{5,4,7,6},{2,7,4}}; _int_elem_mat = {mats[n3+1],mats[n2+1],matrix,matrix,matrix};}
						else
							err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
					}

					// Otherwise, the matrix is only on one of the nodes and two of the other nodes are covered by the same inclusion
					else
					{
						if (n0e1)
						{
							if (n2<0)
								{_cutEdge = {1,3,3,2}; _enrichment_on_inclusion = {n0,n0,n3,n3}; _int_elem_struct = {{0,1,8,9},{3,10,11},{2,7,4},{5,4,7,6}}; _int_elem_mat = {mats[n0+1],mats[n3+1],matrix,matrix};}
							else if (n3<0)
								{_cutEdge = {1,3,2,1}; _enrichment_on_inclusion = {n0,n0,n2,n2}; _int_elem_struct = {{0,1,8,9},{2,10,11},{3,5,6},{5,4,7,6}}; _int_elem_mat = {mats[n0+1],mats[n2+1],matrix,matrix};}
							else err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
						}
						else if (n1e2)
						{
							if (n3<0)
								{_cutEdge = {2,0,0,3}; _enrichment_on_inclusion = {n1,n1,n0,n0}; _int_elem_struct = {{1,2,8,9},{0,10,11},{3,7,4},{5,4,7,6}}; _int_elem_mat = {mats[n1+1],mats[n0+1],matrix,matrix};}
							else if (n0<0)
								{_cutEdge = {2,0,3,2}; _enrichment_on_inclusion = {n1,n1,n3,n3}; _int_elem_struct = {{1,2,8,9},{3,10,11},{0,5,6},{5,4,7,6}}; _int_elem_mat = {mats[n1+1],mats[n3+1],matrix,matrix};}
							else err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
						}
						else if (n2e3)
						{
							if (n0<0)
								{_cutEdge = {3,1,1,0}; _enrichment_on_inclusion = {n2,n2,n1,n1}; _int_elem_struct = {{2,3,8,9},{1,10,11},{0,7,4},{5,4,7,6}}; _int_elem_mat = {mats[n2+1],mats[n1+1],matrix,matrix};}
							else if (n1<0)
								{_cutEdge = {3,1,0,3}; _enrichment_on_inclusion = {n2,n2,n0,n0}; _int_elem_struct = {{2,3,8,9},{0,10,11},{1,5,6},{5,4,7,6}}; _int_elem_mat = {mats[n2+1],mats[n0+1],matrix,matrix};}
							else err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
						}
						else if (n0e3)
						{
							if (n1<0)
								{_cutEdge = {0,2,2,1}; _enrichment_on_inclusion = {n3,n3,n2,n2}; _int_elem_struct = {{3,0,8,9},{2,10,11},{1,7,4},{5,4,7,6}}; _int_elem_mat = {mats[n3+1],mats[n2+1],matrix,matrix};}
							else if (n2<0)
								{_cutEdge = {0,2,1,0}; _enrichment_on_inclusion = {n3,n3,n1,n1}; _int_elem_struct = {{3,0,8,9},{1,10,11},{2,5,6},{5,4,7,6}}; _int_elem_mat = {mats[n3+1],mats[n1+1],matrix,matrix};}
							else err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
						}
						else
							err_message("Unkown manner of 2 inclusions and the matrix intersecting an element.");
					}

					// All cohesive structures are the same here
					_coh_elem_struct = {{4,5,8,9}, {6,7,10,11}};
				}


				// Cut by 3 interface (one of them covers 2 of the nodes)
				else
				{
					if (n0e1)
						{_cutEdge = {1,3,3,2,2,1}; _enrichment_on_inclusion = {n0,n0,n3,n3,n2,n2}; _int_elem_struct = {{0,1,10,11},{3,12,13},{2,14,15},{5,7,6},{4,9,8},{5,4,8,7}}; _int_elem_mat = {mats[n0+1],mats[n3+1],mats[n2+1],matrix,matrix,matrix};}
					else if (n1e2)
						{_cutEdge = {2,0,0,3,3,2}; _enrichment_on_inclusion = {n1,n1,n0,n0,n3,n3}; _int_elem_struct = {{1,2,10,11},{0,12,13},{3,14,15},{5,7,6},{4,9,8},{5,4,8,7}}; _int_elem_mat = {mats[n1+1],mats[n0+1],mats[n3+1],matrix,matrix,matrix};}
					else if (n2e3)
						{_cutEdge = {3,1,1,0,0,3}; _enrichment_on_inclusion = {n2,n2,n1,n1,n0,n0}; _int_elem_struct = {{2,3,10,11},{1,12,13},{0,14,15},{5,7,6},{4,9,8},{5,4,8,7}}; _int_elem_mat = {mats[n2+1],mats[n1+1],mats[n0+1],matrix,matrix,matrix};}
					else if (n0e3)
						{_cutEdge = {0,2,2,1,1,0}; _enrichment_on_inclusion = {n3,n3,n2,n2,n1,n1}; _int_elem_struct = {{3,0,10,11},{2,12,13},{1,14,15},{5,7,6},{4,9,8},{5,4,8,7}}; _int_elem_mat = {mats[n3+1],mats[n2+1],mats[n1+1],matrix,matrix,matrix};}
					else
						err_message("Unknown manner of 3 inclusions intersecting an element.");

					_coh_elem_struct = {{4,5,10,11},{6,7,12,13},{8,9,14,15}};
				}
			} // End 3 unique nodes



			// 4 unique nodal detection values
			else
			{
				// One of the nodal detection points is in th ematrix which means we have 3 inclusions and the matrix here
				if (vec[0] < 0)
				{
					if (n0 < 0)			// Node 0 is in the matrix
						{_cutEdge = {3,2,2,1,1,0}; _enrichment_on_inclusion = {n3,n3,n2,n2,n1,n1}; _int_elem_struct = {{3,10,11},{2,12,13},{1,14,15},{0,9,8,7},{0,7,6},{0,6,5,4}}; _int_elem_mat = {mats[n3+1],mats[n2+1],mats[n1+1],matrix,matrix,matrix};}
					else if (n1 < 0)	// Node 1 is in the matrix
						{_cutEdge = {0,3,3,2,2,1}; _enrichment_on_inclusion = {n0,n0,n3,n3,n2,n2}; _int_elem_struct = {{0,10,11},{3,12,13},{2,14,15},{1,9,8,7},{1,7,6},{1,6,5,4}}; _int_elem_mat = {mats[n0+1],mats[n3+1],mats[n2+1],matrix,matrix,matrix};}
					else if (n2 < 0)	// Node 2 is in the matrix
						{_cutEdge = {1,0,0,3,3,2}; _enrichment_on_inclusion = {n1,n1,n0,n0,n3,n3}; _int_elem_struct = {{1,10,11},{0,12,13},{3,14,15},{2,9,8,7},{2,7,6},{2,6,5,4}}; _int_elem_mat = {mats[n1+1],mats[n0+1],mats[n3+1],matrix,matrix,matrix};}
					else if (n3 < 0)	// Node 3 is in the matrix
						{_cutEdge = {2,1,1,0,0,3}; _enrichment_on_inclusion = {n2,n2,n1,n1,n0,n0}; _int_elem_struct = {{2,10,11},{1,12,13},{0,14,15},{3,9,8,7},{3,7,6},{3,6,5,4}}; _int_elem_mat = {mats[n2+1],mats[n1+1],mats[n0+1],matrix,matrix,matrix};}
					else
						err_message("Unknown manner of 3 inclusions and the matrix intersecting the element.");

					_coh_elem_struct = {{4,5,10,11},{6,7,12,13},{8,9,14,15}};
				}

				// Otherwise 4 inclusions cut this element
				else
				{
					_cutEdge = {3,2,2,1,1,0,0,3};
					_enrichment_on_inclusion = {n3,n3,n2,n2,n1,n1,n0,n0};
					_int_elem_struct = {{3,12,13},{2,14,15},{1,16,17},{0,18,19},{11,10,9,8},{11,8,7,4},{7,6,5,4}};
					_int_elem_mat = {mats[n3+1],mats[n2+1],mats[n1+1],mats[n0+1],matrix,matrix,matrix};
				}
			} // End 4 unique nodes

			_enrichment_nodes.resize(_cutEdge.size() * 2);
			_coh_mat.resize(_coh_elem_struct.size());
		} // End is cohesive
	} // End is intersected
*/

	_detected = true;
	
} // detection








void Quad4::refinement_nodes(std::vector<std::vector<id_type> >& refine_nodes)
{
	std::vector<id_type> ids(n_nodes());
	for (id_type n=0; n<n_nodes(); ++n)
		ids[n] = _nodes[n]->get_id();
	refine_nodes = {{ids[0], ids[1]},
					{ids[1], ids[2]},
					{ids[2], ids[3]},
					{ids[3], ids[0]},
					{ids[0], ids[1], ids[2], ids[3]}};
}
void Quad4::refinement_structure(std::vector<std::vector<id_type> >& structure, std::vector<elem_type>& types)
{
	structure = {{0, 4, 8, 7},
				 {4, 1, 5, 8},
				 {8, 5, 2, 6},
				 {7, 8, 6, 3}};
	types = {QUAD4, QUAD4, QUAD4, QUAD4};
}
