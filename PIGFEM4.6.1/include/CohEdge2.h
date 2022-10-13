/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#ifndef _COH_EDGE2_H_
#define _COH_EDGE2_H_

#include "CohesiveElem.h"

// Predeclarations
//class CohesiveElem;

class CohEdge2 : public CohesiveElem
{
	public:
		CohEdge2();
		CohEdge2(std::vector<Node*> nodes, id_type id);
		
		// Element defintion
		id_type n_nodes() const {return 4;};
		id_type dim() const {return 1;};
		coh_elem_type get_type() const {return COHEDGE2;};

		// Quadrature definition
		id_type n_q_points(int order=3) const;
		void q_point(std::vector<double>& coords, double& w, id_type qp, int order=3) const;

		// Shape function definition
		double compute_shape(const std::vector<double>& coords, int n) const;
		std::vector<double> compute_shape_grad(const std::vector<double>& coords, int n) const;

		// Output helper function
		std::vector<id_type> plot_elem_ids();

	protected:

		// Defines the rotation of a general surace for this dimension element
		// void compute_rotation_general(const std::vector<double>& coords, const std::vector<std::vector<double> >& node_coords,
		// 									  DenseMatrix<double>& mat, double& dA) const;
};


inline
CohEdge2::CohEdge2()
{
	_nodes.resize(4);
}

inline
CohEdge2::CohEdge2(std::vector<Node*> nodes, id_type id)
	:CohesiveElem(nodes, id)
{
}

#endif
