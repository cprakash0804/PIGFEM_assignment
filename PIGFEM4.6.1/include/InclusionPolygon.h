/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated October 2016

##################################################################################
*/
#ifndef _POLYGON_INCLUSION_H_
#define _POLYGON_INCLUSION_H_


#include "Inclusion.h"


/*
 * This class defines an inclusion described by a general polygon
 * The vertices of the polygon are the parameter for the shape
 * There is no ned to set the last index as being the same as the first
 * There is no restriction on the orientation or convexivity of the polygon
 */
class Polygon_Inclusion : public Inclusion
{
	private:
		std::vector<double> _vertices;
		std::vector<double> _bounds;

		double& vx(const id_type& vertex);
		double& vy(const id_type& vertex);
		void set_vertices(const std::vector<double>& vert);
		id_type n_vertices() {return _vertices.size()/2;};

	public:
		Polygon_Inclusion();
		Polygon_Inclusion(Polygon_Inclusion & other_inc);
		virtual ~Polygon_Inclusion();

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
		void copy(Inclusion* other_inc);	// virtual function
};



#endif
