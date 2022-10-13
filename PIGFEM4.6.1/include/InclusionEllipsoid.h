/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#ifndef _ELLIPSOID_INCLUSION_H_
#define _ELLIPSOID_INCLUSION_H_


#include "Inclusion.h"

class Ellipsoid_Inclusion : public Inclusion
{
	private:
		std::vector<double> _center;
		double _alpha, _ca, _sa;	// first euler angle (x-axis)
		double _beta, _cb, _sb;	// second euler angle (y-axis)
		double _gamma, _cg, _sg;	// third euler angle (z-axis)
		double _a;		// major axis radius
		double _b;		// semi-major axis radius
		double _c;		// minor axis radius

	public:
		Ellipsoid_Inclusion();
		Ellipsoid_Inclusion(Ellipsoid_Inclusion & other_inc);
		~Ellipsoid_Inclusion();

		// Here are several function that every type of inclusion must define
		// They involve setting and getting various parameters as well
		//	as detecting if a node is within an interface or not
		virtual void set_parameter(std::string name, double val);
		virtual double get_parameter(std::string name);
		virtual void set_vec_parameter(std::string name, std::vector<double>& val);
		virtual std::vector<double> get_vec_parameter(std::string name);
		virtual inclusion_type get_type() {return ELLIPSOID;};
		virtual std::vector<double> get_domain_lims();
		virtual id_type detect_node(Node* node);
		virtual std::vector<double> get_surface_normal(Node* node);

		// Here are functions and strucures necessary for finding the intersection point of this inclusion
		// 	with a line defined by a set of nodes.
		// find_intersection will be a general purpose function to find the solution of a general nonlinear problem
		//	with int_f, int_df, and int_fdf defieing the problem.
		// find_intersection can be overridden if the intersection finding process is easier for a given inclusion type
		// virtual std::vector<double> find_intersection(std::vector<Node>& nodes);
		virtual double int_f(double alpha, void *params);
		virtual double int_df(double alpha, void *params);
		virtual void int_fdf(double alpha, void *params, double *f, double *df);
		virtual std::vector<double> getIntersectionVelocity(std::vector<Node*>& nodes, std::vector<double> intersect, std::string param_name);
		virtual SensitivityShapeParameter* getSensitivityParameter(std::string param_name, classification type);
		virtual std::vector<std::string> getSensitivityParameterNames();

		virtual Inclusion* allocate_and_copy();
		Inclusion* allocate();			// virtual function
		void copy(Inclusion* other_mat);	// virtual function
	private:
		void convert_to_local(double * global, double * local);

};



#endif
