/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#include "point1.h"

const id_type Point1::s_n_map[0][0] =
{};

const id_type Point1::e_n_map[0][0] =
{};

id_type Point1::side_nodes_map(id_type side, id_type node)
{
	err_message("Attempting to access the side of 0-D element.");
}

id_type Point1::edge_nodes_map(id_type edge, id_type node)
{
	err_message("Attempting to access the edge of 0-D element.");
}

Elem* Point1::build_side(id_type side)
{
	err_message("Attempting to build the side of a 0-D element.");
}



id_type Point1::n_q_points(int order) const
{
	return 1;
}
void Point1::q_point(std::vector<double>& coords, double& w, id_type qp, int order) const
{
	coords.clear();
	coords.push_back(0.0);
	w = 1.0;
}

std::vector<double> Point1::compute_shape(const std::vector<double>& coords) const
{
	return std::vector<double>(1, 1.0);
}

std::vector<std::vector<double> > Point1::compute_shape_grad(const std::vector<double>& coords) const
{
	return std::vector<std::vector<double> >(1, std::vector<double>(1, 0.0));
}

std::vector<DenseMatrix<double> > Point1::compute_shape_grad_grad(const std::vector<double>& coords) const
{
	return std::vector<DenseMatrix<double> >(1, DenseMatrix<double>(1,1));
}







// IGFEM MATERIAL BELOW HERE
//==============================================================================================================================

void Point1::detection(std::vector<int>& node_detection, std::vector<Material*>& inc_mats, bool isCohesive)
{
	err_message("Attempting to perform element detection on a 0-D element.");
} // detection



void Point1::refinement_nodes(std::vector<std::vector<id_type> >& refine_nodes)
{
	err_message("No refinement for Point1 elements.");
}
void Point1::refinement_structure(std::vector<std::vector<id_type> >& structure, std::vector<elem_type>& types)
{
	err_message("No refinement for Point1 elements.");
}
