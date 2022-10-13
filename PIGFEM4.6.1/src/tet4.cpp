/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#include "tet4.h"

// ------------------------------------------------------------
// Tet4 class static member initializations
const id_type Tet4::s_n_map[4][3] =
{
	{0, 2, 1}, // Side 0
	{0, 1, 3}, // Side 1
	{1, 2, 3}, // Side 2
	{2, 0, 3}  // Side 3
};

const id_type Tet4::e_n_map[6][2] =
{
	{0, 1}, // Edge 0
	{1, 2}, // Edge 1
	{0, 2}, // Edge 2
	{0, 3}, // Edge 3
	{1, 3}, // Edge 4
	{2, 3}  // Edge 5
};

id_type Tet4::side_nodes_map(id_type side, id_type node)
{
	if(side>=n_sides())
		err_message("Please select a valid side.");
	if(node>=3)
		err_message("Please select a valid node.");
	
	return s_n_map[side][node];
}

id_type Tet4::edge_nodes_map(id_type edge, id_type node)
{
	if(edge>=n_edges())
		err_message("Please select a valid edge.");
	if(node>=2)
		err_message("Please select a valid node.");
	
	return e_n_map[edge][node];
}

Elem* Tet4::build_side(id_type side)
{
	if(side >= n_sides())
		err_message("Side number must be less than the number of sides.");

	Elem* el = build(TRI3);
	std::vector<Node*> nodes;
	for(id_type n=0; n<3; ++n)
		nodes.push_back(_nodes[side_nodes_map(side,n)]);
	el->set_nodes(nodes);
	return el;
}


id_type Tet4::n_q_points(int order) const
{
	if(order < 0)
		err_message("Order of integration must be greater than or equal to 0.");
	
	switch(order)
	{
		case 0:
		case 1:
			return 1;
		case 2:
			return 4;
		default:
		case 3:
			return 5;
		case 4:
			return 11;
		case 5:
			return 15;
	}
}
void Tet4::q_point(std::vector<double>& coords, double& w, id_type qp, int order) const
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
			if(qp==0) {coords = {0.25, 0.25, 0.25}; w = 1.0;}
			else err_message("Quadrature point for 1st order quadrature of a Tet must be less than 1");
			break;
		}
		case 2:
		{
			if(qp==0) {coords = {0.138196601125011, 0.138196601125011, 0.138196601125011}; w = 0.25;}
			else if(qp==1) {coords = {0.585410196624969, 0.138196601125011, 0.138196601125011}; w = 0.25;}
			else if(qp==2) {coords = {0.138196601125011, 0.585410196624969, 0.138196601125011}; w = 0.25;}
			else if(qp==3) {coords = {0.138196601125011, 0.138196601125011, 0.585410196624969}; w = 0.25;}
			else err_message("Quadrature point for 2nd order quadrature of a Tet must be less than 4");
			break;
		}
		default:
		case 3:
		{
			if(qp==0) {coords = {0.25, 0.25, 0.25}; w = -0.8;}
			else if(qp==1) {coords = {0.5, 0.166666666666667, 0.166666666666667}; w = 0.45;}
			else if(qp==2) {coords = {0.166666666666667, 0.166666666666667, 0.166666666666667}; w = 0.45;}
			else if(qp==3) {coords = {0.166666666666667, 0.166666666666667, 0.5}; w = 0.45;}
			else if(qp==4) {coords = {0.166666666666667, 0.5, 0.166666666666667}; w = 0.45;}
			else err_message("Quadrature point for 3rd order quadrature of a Tet must be less than 5");
			break;
		}
		case 4:
		{
			if(qp==0) {coords = {0.25, 0.25, 0.25}; w = -0.013155555555556;}
			else if(qp==1) {coords = {0.071428571428571, 0.071428571428571, 0.071428571428571}; w = 0.007622222222222;}
			else if(qp==2) {coords = {0.785714285714286, 0.071428571428571, 0.071428571428571}; w = 0.007622222222222;}
			else if(qp==3) {coords = {0.071428571428571, 0.785714285714286, 0.071428571428571}; w = 0.007622222222222;}
			else if(qp==4) {coords = {0.071428571428571, 0.071428571428571, 0.785714285714286}; w = 0.007622222222222;}
			else if(qp==5) {coords = {0.399403576166799, 0.100596423833201, 0.100596423833201}; w = 0.024888888888889;}
			else if(qp==6) {coords = {0.100596423833201, 0.399403576166799, 0.100596423833201}; w = 0.024888888888889;}
			else if(qp==7) {coords = {0.399403576166799, 0.399403576166799, 0.100596423833201}; w = 0.024888888888889;}
			else if(qp==8) {coords = {0.100596423833201, 0.100596423833201, 0.399403576166799}; w = 0.024888888888889;}
			else if(qp==9) {coords = {0.399403576166799, 0.100596423833201, 0.399403576166799}; w = 0.024888888888889;}
			else if(qp==10) {coords = {0.100596423833201, 0.399403576166799, 0.399403576166799}; w = 0.024888888888889;}
			else err_message("Quadrature point for 4th order quadrature of a Tet must be less than 11");
			break;
		}
		case 5:
		{
			if(qp==0) {coords = {0.25, 0.25, 0.25}; w = 0.030283678097089;}
			else if(qp==1) {coords = {0.333333333333333, 0.333333333333333, 0.333333333333333}; w = 0.006026785714286;}
			else if(qp==2) {coords = {0.0, 0.333333333333333, 0.333333333333333}; w = 0.006026785714286;}
			else if(qp==3) {coords = {0.333333333333333, 0.0, 0.333333333333333}; w = 0.006026785714286;}
			else if(qp==4) {coords = {0.333333333333333, 0.333333333333333, 0.0}; w = 0.006026785714286;}
			else if(qp==5) {coords = {0.090909090909091, 0.090909090909091, 0.090909090909091}; w = 0.011645249086029;}
			else if(qp==6) {coords = {0.727272727272727, 0.090909090909091, 0.090909090909091}; w = 0.011645249086029;}
			else if(qp==7) {coords = {0.090909090909091, 0.727272727272727, 0.090909090909091}; w = 0.011645249086029;}
			else if(qp==8) {coords = {0.090909090909091, 0.090909090909091, 0.727272727272727}; w = 0.011645249086029;}
			else if(qp==9) {coords = {0.066550153573664, 0.433449846426336, 0.433449846426336}; w = 0.010949141561386;}
			else if(qp==10) {coords = {0.433449846426336, 0.066550153573664, 0.433449846426336}; w = 0.010949141561386;}
			else if(qp==11) {coords = {0.433449846426336, 0.433449846426336, 0.066550153573664}; w = 0.010949141561386;}
			else if(qp==12) {coords = {0.433449846426336, 0.066550153573664, 0.066550153573664}; w = 0.010949141561386;}
			else if(qp==13) {coords = {0.066550153573664, 0.433449846426336, 0.066550153573664}; w = 0.010949141561386;}
			else if(qp==14) {coords = {0.066550153573664, 0.066550153573664, 0.433449846426336}; w = 0.010949141561386;}
			else err_message("Quadrature point for 4th order quadrature of a Tet must be less than 15");
			break;
		}
	}
}

std::vector<double> Tet4::compute_shape(const std::vector<double>& coords) const
{
	double r = coords[0];
	double s = coords[1];
	double t = coords[2];
	std::vector<double> N = {1.0 - r - s - t,
							 r,
							 s,
							 t};
	return N;
}

std::vector<std::vector<double> > Tet4::compute_shape_grad(const std::vector<double>& coords) const
{
	std::vector<std::vector<double> > dN = {{-1.0, -1.0, -1.0},
											{1.0, 0.0, 0.0},
											{0.0, 1.0, 0.0},
											{0.0, 0.0, 1.0}};
	return dN;
}

std::vector<DenseMatrix<double> > Tet4::compute_shape_grad_grad(const std::vector<double>& coords) const
{
	DenseMatrix<double> ret(3,3);
	ret(0,0) = 0.0;	ret(0,1) = 0.0;	ret(0,2) = 0.0;
	ret(1,0) = 0.0;	ret(1,1) = 0.0;	ret(1,2) = 0.0;
	ret(2,0) = 0.0;	ret(2,1) = 0.0;	ret(2,2) = 0.0;
	std::vector<DenseMatrix<double> > d2N(n_nodes());
	for (id_type n=0; n<n_nodes(); ++n)
		d2N[n] = ret;
	return d2N;
}


bool Tet4::coords_inside(const std::vector<double>& rcoords)
{
	if (rcoords.size() != 3)
		err_message("Invalid coordinate size for coordinate inside check (TET4)!");
	
	if (rcoords[0] >= (0.0-_inside_tol) && rcoords[0] <= (1.0+_inside_tol))
	{
		if (rcoords[1] >= (0.0-_inside_tol) && rcoords[1] <= (1.0+_inside_tol))
		{
			if (rcoords[2] >= (0.0-_inside_tol) && rcoords[2] <= (1.0+_inside_tol))
			{
				if ((rcoords[0] + rcoords[1] + rcoords[2]) <= (1.0+_inside_tol))
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



// IGFEM MATERIAL BELOW HERE
//==============================================================================================================================

void Tet4::detection(std::vector<int>& node_detection, std::vector<Material*>& mats, bool isCohesive)
{
	if(node_detection.size() != n_nodes())
		err_message("The number of nodes used for element detection must be the same as the number of nodes in the element.");

	// NOTE: the mats structure is a pointer to the material for every inclusion in the mesh, with the 0th elment being a pointer to the base/matrix material
	// For this reason, when we index the mats array with the nodal detection values (incusion indices) we need to add 1
	//id_type n_inclusions = mats.size();
	
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
	
	
	
	
	// An interface cuts the element somehow
	else
	{
		_is_intersected = true;

		// Not cohesive
		if (!isCohesive)
		{
			if (n_unique_nodes ==2)
			{
				// If one of the materials is the matrix
				if (vec[0] < 0)
				{
					// Three nodes are located in one material and one in another
					if ((set_counts[0]==1) || (set_counts[0]==3))
					{
						if (n1e2 && n1e3)		// node zero is separate
							{_cutEdge = {0,3,2}; _int_elem_struct = {{4,5,6,0},{1,3,2,4},{2,4,3,5},{4,6,5,2}}; _int_elem_mat = {mats[n0+1],mats[n1+1],mats[n1+1],mats[n1+1]};}
							//{_cutEdge = {0, 3, 2}; _int_elem_struct = {{4, 5, 6, 0}, {1, 3, 2, 4, 5, 6}}; _int_elem_mat = {mats[n0+1], mats[n1+1]};}
						else if (n0e2 && n0e3)	// node one is separate
							{_cutEdge = {0,1,4}; _int_elem_struct = {{4,5,6,1},{0,2,3,4},{2,3,4,5},{4,6,5,3}}; _int_elem_mat = {mats[n1+1],mats[n0+1],mats[n0+1],mats[n0+1]};}
							//{_cutEdge = {0, 1, 4}; _int_elem_struct = {{4, 5, 6, 1}, {0, 2, 3, 4, 5, 6}}; _int_elem_mat = {mats[n1+1], mats[n0+1]};}
						else if (n0e1 && n0e3)	// node two is separate
							{_cutEdge = {1,2,5}; _int_elem_struct = {{4,5,6,2},{1,0,3,4},{0,3,4,5},{4,6,5,3}}; _int_elem_mat = {mats[n2+1],mats[n1+1],mats[n1+1],mats[n1+1]};}
							//{_cutEdge = {1, 2, 5}; _int_elem_struct = {{4, 5, 6, 2}, {1, 0, 3, 4, 5, 6}}; _int_elem_mat = {mats[n2+1], mats[n1+1]};}
						else if (n0e1 && n1e2)	// node three is separate
							{_cutEdge = {3,4,5}; _int_elem_struct = {{4,5,6,3},{0,1,2,4},{1,2,4,5},{4,6,5,2}}; _int_elem_mat = {mats[n3+1],mats[n1+1],mats[n1+1],mats[n1+1]};}
							//{_cutEdge = {3, 4, 5}; _int_elem_struct = {{4, 5, 6, 3}, {0, 1, 2, 4, 5, 6}}; _int_elem_mat = {mats[n3+1], mats[n0+1]};}
						else
							err_message("TETRAHEDRAL ELEMENT HAS 1 NODE SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
					}

					// Two nodes are in one material and two nodes are in another. Essentially amounts to an edge being cut off from the rest
					else if (set_counts[0]==2)
					{
						if (n0e1 && n2e3 && !n0e2)			// Edges 0 and 5 are split
							{_cutEdge = {1,2,3,4}; _int_elem_struct = {{0,1,4,7},{0,7,4,6},{0,6,4,5},{2,6,7,3},{2,6,4,7},{2,5,4,6}}; _int_elem_mat = {mats[n0+1],mats[n0+1],mats[n0+1],mats[n2+1],mats[n2+1],mats[n2+1]};}
							//{_cutEdge = {1, 2, 3, 4}; _int_elem_struct = {{0, 5, 6, 1, 4, 7}, {2, 5, 4, 3, 6, 7}}; _int_elem_mat = {mats[n0+1], mats[n2+1]};}
						else if (n1e2 && n0e3 && !n0e1)		// Edges 1 and 3 are split
							{_cutEdge = {0,4,5,2}; _int_elem_struct = {{1,2,4,5},{2,7,4,5},{2,7,5,6},{0,4,7,3},{4,5,7,3},{3,5,7,6}}; _int_elem_mat = {mats[n1+1],mats[n1+1],mats[n1+1],mats[n0+1],mats[n0+1],mats[n0+1]};}
							//{_cutEdge = {0, 4, 5, 2}; _int_elem_struct = {{0, 4, 7, 3, 5, 6}, {1, 4, 5, 2, 7, 6}}; _int_elem_mat = {mats[n0+1], mats[n1+1]};}
						else if (n0e2 && n1e3 && !n0e1)		// Edges 2 and 4 are split
							{_cutEdge = {0,1,5,3}; _int_elem_struct = {{0,4,2,7},{2,7,4,6},{2,4,5,6},{1,6,7,3},{1,6,4,7},{1,5,4,6}}; _int_elem_mat = {mats[n0+1],mats[n0+1],mats[n0+1],mats[n1+1],mats[n1+1],mats[n1+1]};}
							//{_cutEdge = {0, 1, 5, 3}; _int_elem_struct = {{1, 5, 4, 3, 6, 7}, {2, 5, 6, 0, 4, 7}}; _int_elem_mat = {mats[n1+1], mats[n2+1]};}
						else				// Interface is highly curved
							err_message("TETRAHEDRAL ELEMENT HAS 2 NODES SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
					}

					else
						err_message("Tet4 element detector: We'll never get here! (Hopefully)");

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
			else
			{
				// One of the materials is the matrix
				if (vec[0] < 0)
				{}

				// Otherwise 4 distinct inclusions cut this element
				else
				{}
			} // End 4 unique nodes



			// Actually do the work of assigning all of the variables
			_coh_elem_struct.clear();
			_enrichment_nodes.resize(_cutEdge.size());
		}




		// Cohesive
		else
		{
			// 2 unique nodal detection values
			if (n_unique_nodes == 2)
			{
				// If one of the materials is the matrix
				if (vec[0] < 0)
				{
					// Three nodes are located in one material and one in another
					if ((set_counts[0]==1) || (set_counts[0]==3))
					{
						if (n1e2 && n1e3)		// node zero is separate
						{
							_cutEdge = {0,3,2};
							_int_elem_mat = {mats[n0+1],mats[n1+1],mats[n1+1],mats[n1+1]};
							if (n0 < n1) {_int_elem_struct = {{4,5,6,0},{1,3,2,7},{2,7,3,8},{7,9,8,2}}; _coh_elem_struct = {{4,6,5,7,9,8}};}
							else {_int_elem_struct = {{7,8,9,0},{1,3,2,4},{2,4,3,5},{4,6,5,2}}; _coh_elem_struct = {{4,5,6,7,8,9}};}
							//_int_elem_mat = {mats[n0+1],mats[n1+1]};
							//if (n0 < n1) {_int_elem_struct = {{4,5,6,0},{1,3,2,7,8,9}}; _coh_elem_struct = {{4,6,5,7,9,8}};}
							//else {_int_elem_struct = {{7,8,9,0},{1,3,2,4,5,6}}; _coh_elem_struct = {{4,5,6,7,8,9}};}
						}
						else if (n0e2 && n0e3)	// node one is separate
							{
							_cutEdge = {0,1,4};
							_int_elem_mat = {mats[n1+1],mats[n0+1],mats[n0+1],mats[n0+1]};
							if (n1 < n0) {_int_elem_struct = {{4,5,6,1},{0,2,3,7},{2,3,7,8},{7,9,8,3}}; _coh_elem_struct = {{4,6,5,7,9,8}};}
							else {_int_elem_struct = {{7,8,9,1},{0,2,3,4},{2,3,4,5},{4,6,5,3}}; _coh_elem_struct = {{4,5,6,7,8,9}};}
							//_int_elem_mat = {mats[n1+1],mats[n0+1]};
							//if (n1 < n0) {_int_elem_struct = {{4,5,6,1},{0,2,3,7,8,9}}; _coh_elem_struct = {{4,6,5,7,9,8}};}
							//else {_int_elem_struct = {{4,5,6,1},{0,2,3,7,8,9}}; _coh_elem_struct = {{4,5,6,7,8,9}};}
						}
						else if (n0e1 && n0e3)	// node two is separate
							{
							_cutEdge = {1,2,5};
							_int_elem_mat = {mats[n2+1],mats[n1+1],mats[n1+1],mats[n1+1]};
							if (n2 < n0) {_int_elem_struct = {{4,5,6,2},{0,1,7,3},{0,3,7,8},{7,9,8,3}}; _coh_elem_struct = {{4,6,5,7,9,8}};}
							else {_int_elem_struct = {{7,8,9,2},{0,1,4,3},{0,3,4,5},{4,6,5,3}}; _coh_elem_struct = {{4,5,6,7,8,9}};}
							//_int_elem_mat = {mats[n2+1],mats[n1+1]};
							//if (n2 < n0) {_int_elem_struct = {{4,5,6,2},{1,0,3,4,5,6}}; _coh_elem_struct = {{4,6,5,7,9,8}};}
							//else {_int_elem_struct = {{7,8,9,2},{1,0,3,7,8,9}}; _coh_elem_struct = {{4,5,6,7,8,9}};}
						}
						else if (n0e1 && n1e2)	// node three is separate
							{
							_cutEdge = {3,4,5};
							_int_elem_mat = {mats[n3+1],mats[n1+1],mats[n1+1],mats[n1+1]};
							if (n3 < n1) {_int_elem_struct = {{4,5,6,3},{0,1,2,7},{1,2,7,8},{7,9,8,2}}; _coh_elem_struct = {{4,6,5,7,9,8}};}
							else {_int_elem_struct = {{7,8,9,3},{0,1,2,4},{1,2,4,5},{4,6,5,2}}; _coh_elem_struct = {{4,5,6,7,8,9}};}
							//_int_elem_mat = {mats[n3+1],mats[n1+1]};
							//if(n3 cc< n1) {_int_elem_struct = {{4,5,6,3},{0,1,2,4,5,6}}; _coh_elem_struct = {{4,6,5,7,9,8}};}
							//else {_int_elem_struct = {{7,8,9,3},{0,1,2,4,5,6}}; _coh_elem_struct = {{4,5,6,7,8,9}};}
						}
						else
							err_message("TETRAHEDRAL ELEMENT HAS 1 NODE SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
					}

					// Two nodes are in one material and two nodes are in another. Essentially amounts to an edge being cut off from the rest
					else if (set_counts[0]==2)
					{
						if (n0e1 && n2e3 && !n0e2)			// Edges 0 and 5 are split
						{
							_cutEdge = {1,2,3,4};
							_int_elem_mat = {mats[n0+1],mats[n0+1],mats[n0+1],mats[n2+1],mats[n2+1],mats[n2+1]};
							if (n0 < n2) {_int_elem_struct = {{0,1,4,7},{0,7,4,6},{0,6,4,5},{2,10,11,3},{2,10,8,11},{2,9,8,10}}; _coh_elem_struct = {{4,5,6,7,8,9,10,11}};}
							else {_int_elem_struct = {{0,1,8,11},{0,11,8,10},{0,10,8,9},{2,6,7,3},{2,6,4,7},{2,5,4,6}}; _coh_elem_struct = {{7,6,5,4,11,10,9,8}};}
							//_int_elem_mat = {mats[n0+1],mats[n2+1]};
							//if (n0 < n2) {_int_elem_struct = {{0,5,6,1,4,7},{2,9,8,3,10,11}}; _coh_elem_struct = {{4,5,6,7,8,9,10,11}};}
							//else {_int_elem_struct = {{0,9,10,1,8,11},{2,5,4,3,6,7}}; _coh_elem_struct = {{7,6,5,4,11,10,9,8}};}
						}
						else if (n1e2 && n0e3 && !n0e1)		// Edges 1 and 3 are split
						{
							_cutEdge = {0,4,5,2};
							_int_elem_mat = {mats[n1+1],mats[n1+1],mats[n1+1],mats[n3+1],mats[n3+1],mats[n3+1]};
							if (n1 < n3) {_int_elem_struct = {{1,2,4,5},{2,7,4,5},{2,7,5,6},{0,8,11,3},{8,9,11,3},{3,9,11,10}}; _coh_elem_struct = {{7,6,5,4,11,10,9,8}};}
							else {_int_elem_struct = {{1,2,8,9},{2,11,8,9},{2,11,9,10},{0,4,7,3},{4,5,7,3},{3,5,7,6}}; _coh_elem_struct = {{4,5,6,7,8,9,10,11}};}
							//_int_elem_mat = {mats[n1+1],mats[n3+1]};
							//if(n1 < n3) {_int_elem_struct = {{1,4,5,2,7,6},{0,8,11,3,9,10}}; _coh_elem_struct = {{7,6,5,4,11,10,9,8}};}
							//else {_int_elem_struct = {{1,8,9,2,11,10},{0,4,7,3,5,6}}; _coh_elem_struct = {{4,5,6,7,8,9,10,11}};}
						}
						else if (n0e2 && n1e3 && !n0e1)		// Edges 2 and 4 are split
						{
							_cutEdge = {0,1,5,3};
							_int_elem_mat = {mats[n0+1],mats[n0+1],mats[n0+1],mats[n1+1],mats[n1+1],mats[n1+1]};
							if (n0 < n1) {_int_elem_struct = {{0,4,2,7},{2,7,4,6},{2,4,5,6},{1,10,11,3},{1,10,8,11},{1,9,8,10}}; _coh_elem_struct = {{4,5,6,7,8,9,10,11}};}
							else {_int_elem_struct = {{0,8,2,11},{2,11,8,10},{2,8,9,10},{1,6,7,3},{1,6,4,7},{1,5,4,6}}; _coh_elem_struct = {{7,6,5,4,11,10,9,8}};}
							//_int_elem_mat = {mats[n0+1],mats[n1+1]};
							//if (n0 < n1) {_int_elem_struct = {{2,5,6,0,4,7},{1,9,8,3,10,11}; _coh_elem_struct = {{4,5,6,7,8,9,10,11}};}
							//else {_int_elem_struct = {{2,9,10,0,8,11},{1,5,4,3,6,7}; _coh_elem_struct = {{7,6,5,4,11,10,9,8}};}
						}
						else				// Interface is highly curved
							err_message("TETRAHEDRAL ELEMENT HAS 2 NODES SEPARATED IN AN UNKNOWN MANNER. MAY NEED REFINEMENT.");
					}

					else
						err_message("Tet4 element detector: We'll never get here! (Hopefully)");

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
			else
			{
				// One of the materials is the matrix
				if (vec[0] < 0)
				{}

				// Otherwise 4 distinct inclusions cut this element
				else
				{}
			} // End 4 unique nodes

			_enrichment_nodes.resize(_cutEdge.size() * 2);
			_coh_mat.resize(_coh_elem_struct.size());
		} // End is cohesive	
	} // End is intersected

	_detected = true;
} // detection




void Tet4::refinement_nodes(std::vector<std::vector<id_type> >& refine_nodes)
{
	std::vector<id_type> ids(n_nodes());
	for (id_type n=0; n<n_nodes(); ++n)
		ids[n] = _nodes[n]->get_id();
	refine_nodes = {{ids[0], ids[1]},
					{ids[1], ids[2]},
					{ids[2], ids[0]},
					{ids[0], ids[3]},
					{ids[1], ids[3]},
					{ids[2], ids[3]}};
}
void Tet4::refinement_structure(std::vector<std::vector<id_type> >& structure, std::vector<elem_type>& types)
{
	structure = {{0, 4, 6, 7},
				 {1, 5, 4, 8},
				 {2, 6, 5, 9},
				 {7, 8, 9, 3},
				 {4, 9, 6, 7},
				 {4, 5, 6, 9},
				 {4, 8, 5, 9},
				 {4, 8, 9, 7}};
	types = {TET4, TET4, TET4, TET4, TET4, TET4, TET4, TET4};
}
