/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated October 2016

##################################################################################
*/
#ifndef _POLYHEDRAL_INCLUSION_H_
#define _POLYHEDRAL_INCLUSION_H_


#include "Inclusion.h"


/*
 * This class defines an inclusion described by a general polyhedral made of triangles
 * The vertices of the polygon are the parameter for the shape
 * There is no ned to set the last index as being the same as the first
 * There is no restriction on the orientation or convexivity of the polygon
 */
class Polyhedron_Inclusion : public Inclusion
{
	private:
		std::vector<std::vector<double> > _facets;		// Sets of 9 doubles that each specify a facet
		std::vector<std::vector<double> > _normals;
		std::vector<double> _bounds;

		double coord(id_type facet, id_type vertex, id_type dir);
		id_type n_facets() {return _facets.size();};
		bool in_y_z_space(const std::vector<double>& P, id_type face);
		std::vector<double>& get_normal(id_type facet) {return _normals[facet];};

	public:
		Polyhedron_Inclusion();
		Polyhedron_Inclusion(Polyhedron_Inclusion & other_inc);
		virtual ~Polyhedron_Inclusion();

		// Here are several function that every type of inclusion must define
		// They involve setting and getting various parameters as well
		//	as detecting if a node is within an interface or not
		virtual void set_parameter(std::string name, double val);
		virtual double get_parameter(std::string name);
		virtual void set_vec_parameter(std::string name, std::vector<double>& val);
		virtual std::vector<double> get_vec_parameter(std::string name);
		virtual inclusion_type get_type() {return PLANE;};
		virtual std::vector<double> get_domain_lims();
		virtual id_type detect_node(Node* node);
		virtual std::vector<double> get_surface_normal(Node* node);
		void set_facets(const std::vector<std::vector<double> >& vert);
		void set_normals(const std::vector<std::vector<double> >& norm) {_normals = norm;};

		// Here are functions and strucures necessary for finding the intersection point of this inclusion
		// 	with a line defined by a set of nodes.
		// find_intersection will be a general purpose function to find the solution of a general nonlinear problem
		//	with int_f, int_df, and int_fdf defieing the problem.
		// find_intersection can be overridden if the intersection finding process is easier for a given inclusion type
		virtual std::vector<double> find_intersection(std::vector<Node>& nodes);
		virtual double int_f(double alpha, void *params);
		virtual double int_df(double alpha, void *params);
		virtual void int_fdf(double alpha, void *params, double *f, double *df);
		virtual std::vector<double> getIntersectionVelocity(std::vector<Node*>& nodes, std::vector<double> intersect, std::string param_name) {err_message("Shape velocity not implmented yet!");};

		virtual Inclusion* allocate_and_copy();
		virtual Inclusion* allocate();			// virtual function
		void copy(Polyhedron_Inclusion* other_inc);	// Overridden function

};

#endif
