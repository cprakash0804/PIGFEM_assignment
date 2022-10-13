/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#ifndef _COH_POINT1_H_
#define _COH_POINT1_H_

#include "CohesiveElem.h"

// Predeclarations
//class CohesiveElem;

class CohPoint1 : public CohesiveElem
{
	public:
		CohPoint1();
		CohPoint1(std::vector<Node*> nodes, id_type id);
		
		// Element defintion
		id_type n_nodes() const {return 2;};
		id_type dim() const {return 0;};
		coh_elem_type get_type() const {return COHPOINT1;};

		// Quadrature definition
		id_type n_q_points(int order=3) const;
		void q_point(std::vector<double>& coords, double& w, id_type qp, int order=3) const;

		// Shape function definition
		double compute_shape(const std::vector<double>& coords, int n) const;
		std::vector<double> compute_shape_grad(const std::vector<double>& coords, int n) const;

	protected:

		// Defines the rotation of a general surace for this dimension element
		// virtual void compute_rotation_general(const std::vector<double>& coords, const std::vector<std::vector<double> >& node_coords,
		// 								DenseMatrix<double>& mat, double& dA) const;
};


inline
CohPoint1::CohPoint1()
{
	_nodes.resize(4);
}

inline
CohPoint1::CohPoint1(std::vector<Node*> nodes, id_type id)
	:CohesiveElem(nodes, id)
{
}

#endif
