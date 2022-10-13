/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#ifndef _Polygonal_H_
#define _Polygonal_H_

#include "elem.h"

class Polygon : public Elem
{
	private:
		id_type _nn; // number of nodes (and sides)

		// Function used to break the polygon into subtriangles in the parent coordinates
		void PolyTriangulate(std::vector<std::vector<double> >& p, std::vector<std::vector<id_type> >& Tri, std::vector<double> rcoords) const;

	public:
		Polygon(id_type nn);
		Polygon(std::vector<Node*> nodes, id_type id);
		
		id_type n_nodes() const {return _nn;};
		id_type n_edges() const {return _nn;};
		id_type n_sides() const {return _nn;};

		virtual id_type edge_nodes_map(id_type edge, id_type node);
		virtual id_type side_nodes_map(id_type side, id_type node);
		virtual id_type n_nodes_on_edge(id_type edge) {return 2;};
		virtual id_type n_nodes_on_side(id_type side) {return 2;};
		virtual Elem* build_side(id_type side);

		virtual id_type h_level() const {return 1;}; // 1=LINEAR, 2=QUADRATIC, etc
		virtual id_type dim() const {return 2;};
		virtual elem_type get_type() const {return POLYGON;};

		// Quadrature functions
		virtual id_type n_q_points(int order=3) const;
		virtual void q_point(std::vector<double>& coords, double& w, id_type qp, int order=5) const;

		// Shape function functions
		virtual std::vector<double> compute_shape(const std::vector<double>& coords) const;
		virtual std::vector<std::vector<double> > compute_shape_grad(const std::vector<double>& coords) const;
		virtual std::vector<DenseMatrix<double> > compute_shape_grad_grad(const std::vector<double>& coords) const;

		virtual bool coords_inside(const std::vector<double>& rcoords);

		// FUNCTIONS HAVING TO DO WITH REFINEMENT
		//------------------------------------------------------------------

		virtual int n_refinement_elem() {return 1;};
		virtual void refinement_nodes(std::vector<std::vector<id_type> >& refine_nodes);
		virtual void refinement_structure(std::vector<std::vector<id_type> >& structure, std::vector<elem_type>& types);
		
		// IGFEM MATERIAL BELOW HERE
		virtual void detection(std::vector<int>& node_detection, std::vector<Material*>& mats, bool isCohesive);
};
inline
Polygon::Polygon(id_type nn)
{
	_nodes.resize(nn);
	_nn = nn;
}

inline
Polygon::Polygon(std::vector<Node*> nodes, id_type id)
	:Elem(nodes, id)
{
}

#endif
