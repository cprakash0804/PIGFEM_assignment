/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#ifndef _QUAD4_H_
#define _QUAD4_H_

#include "elem.h"

class Quad4 : public Elem
{
	public:
		Quad4();
		Quad4(std::vector<Node*> nodes, id_type id);
		
		id_type n_nodes() const {return 4;};
		id_type n_edges() const {return 4;};
		id_type n_sides() const {return 4;};

		/**
		 * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
		 * element node numbers.
		*/
		static const id_type s_n_map[4][2];

		/**
		 * This maps the \f$ j^{th} \f$ node of the \f$ i^{th} \f$ side to
		 * element node numbers.
		*/
		static const id_type e_n_map[4][2];
		virtual id_type edge_nodes_map(id_type edge, id_type node);
		virtual id_type side_nodes_map(id_type side, id_type node);
		virtual id_type n_nodes_on_edge(id_type edge) {return 2;};
		virtual id_type n_nodes_on_side(id_type side) {return 2;};
		virtual Elem* build_side(id_type side);

		virtual id_type h_level() const {return 1;}; // 1=LINEAR, 2=QUADRATIC, etc
		virtual id_type dim() const {return 2;};
		virtual elem_type get_type() const {return QUAD4;};
		virtual id_type n_q_points(int order=3) const;
		virtual void q_point(std::vector<double>& coords, double& w, id_type qp, int order=3) const;
		virtual std::vector<double> compute_shape(const std::vector<double>& coords) const;
		virtual std::vector<std::vector<double> > compute_shape_grad(const std::vector<double>& coords) const;
		virtual std::vector<DenseMatrix<double> > compute_shape_grad_grad(const std::vector<double>& coords) const;

		virtual bool coords_inside(const std::vector<double>& rcoords);
		
		// FUNCTIONS HAVING TO DO WITH REFINEMENT
		//------------------------------------------------------------------

		virtual int n_refinement_elem() {return 4;};
		virtual void refinement_nodes(std::vector<std::vector<id_type> >& refine_nodes);
		virtual void refinement_structure(std::vector<std::vector<id_type> >& structure, std::vector<elem_type>& types);
		
		// IGFEM MATERIAL BELOW HERE
		virtual void detection(std::vector<int>& node_detection, std::vector<Material*>& mats, bool isCohesive);
};
inline
Quad4::Quad4()
{
	_nodes.resize(4);
}

inline
Quad4::Quad4(std::vector<Node*> nodes, id_type id)
	:Elem(nodes, id)
{
}

#endif
