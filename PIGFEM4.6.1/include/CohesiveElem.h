/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _COH_ELEM_H_
#define _COH_ELEM_H_
#include <algorithm>
#include <vector>
#include <cmath>
#include "common.h"
#include "node.h"
#include "material.h"

// Predeclarations
//class CohPoint1;
//class CohEdge2;
//class CohTri3;
//class CohQuad4;

// Defines the cohesive elmeent base class
// This class basically is the same as a normal element but
//  all of its nodes are double nodes with a + and a - surface.
// There are also difference functions necessary for a cohesive element
//  than there are in an volumetric element. For instance, finding the 
//  normal vector and rotation matrix to the ttn "tangent-tanget-normal"
//  curvilinear coordinate system to the global xyz coordinate system.
// Also note, due to the cohesive nature of these elements, elements of 
//  dimension greater then 2 don't really make any sense, whereas elements
//  of dimension 0 actually make sense in this context
class CohesiveElem
{
	protected:
		CohesiveElem();
		CohesiveElem(std::vector<Node*> nodes, id_type id);
		std::vector<Node*> _nodes;
		id_type _global_id; // Don't know if I actualy need an id since each Cohesive element is associated with a volumetric element

	public:
		CohesiveElem(CohesiveElem & other_elem);
		virtual ~CohesiveElem();

		void clear();
		void copy(CohesiveElem & other_elem);

		void set_id(id_type id) {_global_id = id;};
		id_type get_id() {return _global_id;};
		virtual coh_elem_type get_type() const = 0;

		void set_nodes(std::vector<Node*>& nodes);
		void set_node(id_type i, Node* node);
		Node*& set_node(id_type i);
		std::vector<Node*> get_nodes() {return _nodes;};
		Node* get_node(id_type i);
		Node& operator()(id_type i);

		friend std::ostream& operator<<(std::ostream& os, const CohesiveElem& elem)
		{
			elem.print(os);
			return os;
		}

		void print(std::ostream & os) const;

		virtual id_type dim() const = 0;
		virtual id_type n_nodes() const = 0; // As yet undecided if this returns n or 2*n

		// Quadrature rules
		virtual id_type n_q_points(int order=3) const = 0;
		virtual void q_point(std::vector<double>& coords, double& w, id_type qp, int order=3) const = 0;

		// Shape function functions
		virtual double compute_shape(const std::vector<double>& coords, int n) const = 0;
		virtual std::vector<double> compute_shape_grad(const std::vector<double>& coords, int n) const = 0;
		virtual void ShapeFunctions(std::vector<double> rcoords, std::vector<double>& N, DenseMatrix<double>& Rot, double& J);
		void transform_gradients(std::vector<std::vector<double> >& parent_grad, std::vector<std::vector<double> >& dNdx, double& J);

		// Rotation/Orientation functions
		// Defines the rotation matrix from the xyz to the ttn coord system at the given quadrature point
		void compute_rotation(const std::vector<double>& coords, DenseMatrix<double>& mat, double& dA) const;
		void compute_deformed_rotation(const std::vector<double>& coords, const std::vector<double>& UEL,
									   DenseMatrix<double>& mat, double& dA) const;
		void RotationSensitivity(const std::vector<double>& rcoords, const std::vector<double>& coh_velocity, 
								 DenseMatrix<double>& dRot, std::vector<double>& v, double& div_v);

		// Output helper function
		virtual std::vector<id_type> plot_elem_ids();

		static CohesiveElem* build(coh_elem_type etype);

	protected:
		virtual void compute_rotation_general(const std::vector<double>& coords, const std::vector<std::vector<double> >& node_coords,
											  DenseMatrix<double>& mat, double& dA) const;
};









inline
void CohesiveElem::set_nodes(std::vector<Node*>& nodes)
{
	if(nodes.size() == n_nodes())
	{
		_nodes.clear();
		_nodes = nodes;
	}
	else
		err_message("Number of nodes in list must be equal to number of nodes in the element.");
}

inline
Node* CohesiveElem::get_node(id_type i)
{
	if(i<n_nodes())
		return _nodes[i];
	else
		err_message("Index of node must be less than the number of nodes.");
}
inline
void CohesiveElem::set_node(id_type i, Node* node)
{
	if(i<n_nodes())
		_nodes[i] = node;
	else
		err_message("Index of node must be less than the number of nodes.");
}
inline
Node*& CohesiveElem::set_node(id_type i)
{
	if(i<n_nodes())
		return _nodes[i];
	else
		err_message("Index of node must be less than the number of nodes.");
}



#endif

