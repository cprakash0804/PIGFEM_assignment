/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _ELEM_H_
#define _ELEM_H_

#include <algorithm>
#include <vector>
#include <iterator>
#include <set>
#include "common.h"
#include "node.h"
#include "CohesiveElem.h"
#include "material.h"
#include "DenseMatrix.h"
#include <gsl_vector.h> // These are for the inverse map
#include <gsl_multiroots.h>



class Elem
{
	protected:
		Elem();
		Elem(std::vector<Node*> nodes, id_type id);
		std::vector<Node*> _nodes;
		id_type _global_id;

		static double _inside_tol;
	public:
		Elem(Elem & other_elem);
		virtual ~Elem();

		void clear();
		void copy(Elem & other_elem);

		void set_id(id_type id) {_global_id = id;};
		id_type get_id() const {return _global_id;};
		virtual elem_type get_type() const = 0;

		void set_nodes(std::vector<Node*>& nodes);
		Node*& set_node(id_type i);
		void set_node(id_type i, Node* node);
		std::vector<Node*> get_nodes() {return _nodes;};
		Node* get_node(id_type i);
		Node& operator()(id_type i);
		
		virtual id_type edge_nodes_map(id_type edge, id_type node) = 0;
		virtual id_type side_nodes_map(id_type side, id_type node) = 0;
		virtual id_type n_nodes_on_edge(id_type edge) = 0;
		virtual id_type n_nodes_on_side(id_type side) = 0;

		virtual id_type h_level() const = 0; // 1=LINEAR, 2=QUADRATIC, etc
		virtual id_type dim() const = 0;
		virtual id_type n_nodes() const = 0;
		virtual id_type n_edges() const = 0;
		virtual id_type n_sides() const = 0;

		// Quadrature functions
		virtual id_type n_q_points(int order=3) const = 0;
		virtual void q_point(std::vector<double>& coords, double& w, id_type qp, int order=3) const = 0;

		// Shape function functions
		virtual std::vector<double> compute_shape(const std::vector<double>& coords) const = 0;
		virtual std::vector<std::vector<double> > compute_shape_grad(const std::vector<double>& coords) const = 0;
		virtual std::vector<DenseMatrix<double> > compute_shape_grad_grad(const std::vector<double>& coords) const = 0;
		virtual void ShapeFunctions(std::vector<double> rcoords, std::vector<double>& N, std::vector<std::vector<double> >& dNdx, double& J);
		void ShapeFunctionPartials(std::vector<double> rcoords, int child, const std::vector<double>& nodal_velocity,
								   std::vector<double>& dN_dd, std::vector<std::vector<double> >& d2N_dxdd,
								   std::vector<double>& v, double& div_v);
		void ChildShapeFunctionPartials(const std::vector<double>& rcoords, const std::vector<double>& child_nodal_velocity,
										std::vector<double>& dN_dd_child, std::vector<std::vector<double> >& d2N_dxdd_child,
										std::vector<double>& x, std::vector<double>& v, double& div_v);
		void transform_gradients(std::vector<std::vector<double> >& parent_grad, std::vector<std::vector<double> >& dNdx, double& J);
		
		friend std::ostream& operator<<(std::ostream& os, const Elem& elem)
		{
			elem.print(os);
			return os;
		}
		void print(std::ostream & os) const;
		
		static Elem* build(elem_type etype);
		virtual Elem* build_side(id_type side) = 0;

		// Function to compute the inverse map of the globalcoordinates
		// This will be used for elements with nonlinear shape function.
		// For simpler elements (triangle) explicit forms can be found so this will be overridden
		virtual std::vector<double> inverse_map(std::vector<double> gcoords, bool& inside);
		virtual bool coords_inside(const std::vector<double>& rcoords) = 0;


		// function that takes a vector of global node ids and returns a side that they have in common
		bool belongs_to_edge(id_type edge, id_type id);
		bool belongs_to_side(id_type side, id_type id);
		int get_edge_from_nodes(std::vector<id_type>& nodes);
		int get_side_from_nodes(std::vector<id_type>& nodes);

		// FUNCTIONS HAVING TO DO WITH REFINEMENT
		//------------------------------------------------------------------

		virtual int n_refinement_elem() = 0;
		virtual void refinement_nodes(std::vector<std::vector<id_type> >& refine_nodes) = 0;
		virtual void refinement_structure(std::vector<std::vector<id_type> >& structure, std::vector<elem_type>& types) = 0;
		
		
		// IGFEM MATERIAL BELOW HERE
		//------------------------------------------------------------------
		// Return element to pre-cut state
		void clearEnrichments();
		// Used to determine if/how an element is cut
		virtual void detection(std::vector<int>& node_detection, std::vector<Material*>& mats, bool isCohesive) = 0;

		bool is_intersected() {return _is_intersected;};
		
		// Generate the element's integration elmeents defined in the elemental detection routine
		void generate_integration_elem();
		
		// Function to return how many integration elements this elemnt has
		id_type n_integration_elem() {return _integrationElem.size();};

		// Function to return a pointer to an element's integration element
		Elem* get_integration_elem(id_type idx);

		id_type get_cut_edge(id_type e);
		id_type get_inclusion_from_enrich(id_type e);		
		
		// Returns the number of enrichment nodes stored on this element
		id_type n_enrich_nodes() const {return _enrichment_nodes.size();};
		
		// Returns a raw pointer to the requested enrichment node
		Node* get_enrich_node(id_type i);
		Node*& set_enrich_node(id_type i);
		void set_enrich_node(id_type i, Node* node);

		Material* get_int_elem_mat(id_type idx);

		//Functions used for copying an element
		bool& get_detected() {return _detected;};
		std::vector<short_id_type>& getIntersectedEdges() {return _cutEdge;};
		std::vector<id_type>& getInclusionNumbers() {return _enrichment_on_inclusion;};
		std::vector<std::vector<short_id_type> >& getIntegrationElemStructure() {return _int_elem_struct;};
		std::vector<Material*>& getIntegrationElemMaterials() {return _int_elem_mat;};
		std::vector<Elem*>& getIntegrationElem() {return _integrationElem;};
		std::vector<Node*>& getEnrichmentNodes() {return _enrichment_nodes;};
		bool& get_intersected() {return _is_intersected;};

		// Cohesive element functions
		id_type n_cohesive_elem() {return _cohesiveElem.size();};
		CohesiveElem* get_cohesive_elem(id_type idx);
		std::vector<CohesiveElem*>& getCohesiveElem() {return _cohesiveElem;};
		void set_cohesive_material(Material* mat, id_type n);
		Material* get_cohesive_material(id_type n);
		std::vector<std::vector<short_id_type> >& getCohesiveElemStructure() {return _coh_elem_struct;};


		
	protected:
		/*
		*  Integer to store the state of the element, whether or not it has been cut or not
		* 0 : Not cut by an interface. May be entirely in base material or entirely inside an incusion
		* >0 : Stores the type of cut that this element is experiencing. Numbering varies between the types of elements
		*	Ex: For a Quad4:
				1: One corner is dirrerent than the other 3
				2: Two adjacent nodes are different from the other two adjacent nodes
					Note: if the two nodes are not adjacent then the inclusion is highly curved and the element should ebe refined
		*/
		//int _elemDetector;
		bool _is_intersected;
		
		// Boolean to store whether or not the element has been detected yet
		bool _detected;

		/*
		*  Vector to store the list of edges that are cut in this element
		*  If no edges are cut, this vector will be empty
		*/
		std::vector<short_id_type> _cutEdge;

		/*
		* Integers storing the inclusions that each enrichment node lies on
		*/
		std::vector<id_type> _enrichment_on_inclusion;

		// This is for storing the structure of the integration elements prior to actually forming them with the enrchment node pointers
		std::vector<std::vector<short_id_type> > _int_elem_struct;

		// This is for storing the structure of the chesive elements prior to actually forming them with the enrchment node pointers
		std::vector<std::vector<short_id_type> > _coh_elem_struct;

		// This is for storing the pointers to the material that each integration element is made of
		std::vector<Material*> _int_elem_mat;

		// Material for the cohesive surface
		std::vector<Material*> _coh_mat;

		/*
		*  Once the elemnt has been cut, it is difficult to integrate the enrichment functions over the element
		*  For this reason the element is partitioned into several different easy to integrate elements
		*/
		std::vector<Elem*> _integrationElem;

		/*
		 * A vector of the cohesive elements associated with an intersected element
		 * These exist on the interface between the two materials and contain entirely enrichment nodes
		*/
		std::vector<CohesiveElem*> _cohesiveElem;
		
		// This is a list of the enriched nodes associated with the element. Determined via IGFEM
		std::vector<Node*> _enrichment_nodes;

};


// Structure used by GSL in the inverse map function
struct inv_map_params
{
	std::vector<double> goal; // The goal x, y, z position
	Elem* elem;				  // Pointer to the current element (Because I can't make inv_map_f, etc member functions)
};











inline
void Elem::set_nodes(std::vector<Node*>& nodes)
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
Node* Elem::get_node(id_type i)
{
	if(i<n_nodes())
		return _nodes[i];
	else
		err_message("Index of node must be less than the number of nodes.");
}
inline
Node*& Elem::set_node(id_type i)
{
	if(i<n_nodes())
		return _nodes[i];
	else
		err_message("Index of node must be less than the number of nodes.");
}
inline
void Elem::set_node(id_type i, Node* node)
{
	if(i<n_nodes())
		_nodes[i] = node;
	else
		err_message("Index of node must be less than the number of nodes.");
}

inline
Node*& Elem::set_enrich_node(id_type i)
{
	if(i<n_enrich_nodes())
		return _enrichment_nodes[i];
	else
		err_message("Index of enrichment node must be less than the number of enrichment nodes.");
}
inline
void Elem::set_enrich_node(id_type i, Node* node)
{
	if(i<n_enrich_nodes())
		_enrichment_nodes[i] = node;
	else
		err_message("Index of enrichment node must be less than the number of enrichment nodes.");
}


inline
Node* Elem::get_enrich_node(id_type i)
{
	if(i<n_enrich_nodes())
		return _enrichment_nodes[i];
	else
		err_message("Index of enrichment node must be less than the number of enrichment nodes.");
}

inline
Elem* Elem::get_integration_elem(id_type idx)
{
	if(idx<n_integration_elem())
		return _integrationElem[idx];
	else
		err_message("Index of the integration element must be less than the total number of integration elements.");
}

inline
Material* Elem::get_int_elem_mat(id_type idx)
{
	if(idx<n_integration_elem())
		return _int_elem_mat[idx];
	else
		err_message("Index of the integration element must be less than the total number of integration elements.");
}

inline
CohesiveElem* Elem::get_cohesive_elem(id_type idx)
{
	if(idx<n_cohesive_elem())
		return _cohesiveElem[idx];
	else
		err_message("Index of the cohesive element must be less than the total number of cohesive elements.");
}

inline
void Elem::set_cohesive_material(Material* mat, id_type n)
{
	if (n < _coh_mat.size())
		_coh_mat[n] = mat;
	else
		err_message("Attempted to add a cohesive material that is greater than the number of cohesive materials");
}

inline
Material* Elem::get_cohesive_material(id_type n)
{
	if (n < _coh_mat.size())
		return _coh_mat[n];
	else
		err_message("Attempted to access a cohesive material that is greater than the number of cohesive materials");
}




#endif
