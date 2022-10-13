/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _INCLUSION_H_
#define _INCLUSION_H_

#include <vector>
#include <algorithm>
#include <string>
#include "common.h"
#include "node.h"
#include "material.h"
#include <gsl_vector.h> // These are for the intersection finding
#include <gsl_roots.h>


// Predeclarations
class SensitivityShapeParameter;


// Template class to wrap the inherited objective member functions
// I don't really understand what this is doing at all...
template< typename F >
  class gsl_function_pp : public gsl_function {
  public:
  gsl_function_pp(const F& func) : _func(func) {
    function = &gsl_function_pp::invoke;
    params=this;
  }
  private:
  const F& _func;
  static double invoke(double x, void *params) {
    return static_cast<gsl_function_pp*>(params)->_func(x);
  }
};






class Inclusion
{
	protected:
		id_type _id;
		Material* _mat;

	public:
		Inclusion();
		Inclusion(Inclusion & other_inc);
		virtual ~Inclusion();
		
		Material* get_material() {return _mat;};
		void set_material(Material* mat) {_mat = mat;};
		id_type get_id() {return _id;};
		void set_id(id_type id) {_id = id;};

		// Here are several function that every type of inclusion must define
		// They involve setting and getting various parameters as well
		//	as detecting if a node is within an interface or not
		virtual void set_parameter(std::string name, double val) = 0;
		virtual double get_parameter(std::string name) = 0;
		virtual void set_vec_parameter(std::string name, std::vector<double>& val) = 0;
		virtual std::vector<double> get_vec_parameter(std::string name) = 0;
		virtual inclusion_type get_type() = 0;
		virtual std::vector<double> get_domain_lims() = 0;;
		virtual id_type detect_node(Node* node) = 0;
		virtual std::vector<double> get_surface_normal(Node* node) = 0;

		// Here are functions and strucures necessary for finding the intersection point of this inclusion
		// 	with a line defined by a set of nodes.
		// find_intersection will be a general purpose function to find the solution of a general nonlinear problem
		//	with int_f, int_df, and int_fdf defieing the problem.
		// find_intersection can be overridden if the intersection finding process is easier for a given inclusion type
		virtual std::vector<double> find_intersection(std::vector<Node>& nodes);
		virtual double int_f(double alpha, void *params) = 0;
		virtual double int_df(double alpha, void *params) = 0;
		virtual void int_fdf(double alpha, void *params, double *f, double *df) = 0;
		virtual std::vector<double> getIntersectionVelocity(std::vector<Node*>& nodes, std::vector<double> intersect, std::string param_name) = 0;
		virtual SensitivityShapeParameter* getSensitivityParameter(std::string param_name, classification type)  {err_message("Attempting to get the sensitivity parameter of an inclusion that does not support sensitivity!");};
		virtual std::vector<std::string> getSensitivityParameterNames() {err_message("Attempting to get the sensitivity parameter names of an inclusion that does not support sensitivity!");};

		// Functions necessary to copy an inclusion
		virtual Inclusion* allocate_and_copy() = 0;
		virtual Inclusion* allocate() = 0;              // Allocates a new derived type


		
};


/*
inline
Inclusion::Inclusion()
	: _id(0), _mat(NULL)
{
}

inline
Inclusion::Inclusion(Inclusion & other_inc)
{
	other_inc.set_id(_id);
	other_inc.set_material(_mat);
}

inline
Inclusion::~Inclusion()
{
}
*/



struct int_params
{
	std::vector<Node> nodes;
	Inclusion * pt_MyClass;
};

// Wrapper solution to passing class member functions.
// Involves storing a pointer to the current object in the params struct
// Very simple but I really don't like it very much...
// extern makes it so that when a new inclusion type references this header file the linker doesn't think this is a new definition of these functions
inline
extern double int_f_wrapper(double x, void *pp)
{
	int_params *p = (int_params *)pp;
	return p->pt_MyClass->int_f(x, p);
}
inline
extern double int_df_wrapper(double x, void *pp)
{
	int_params *p = (int_params *)pp;
	return p->pt_MyClass->int_df(x, p);
}
inline
extern void int_fdf_wrapper(double x, void *pp, double *f, double *df)
{
	int_params *p = (int_params *)pp;
	p->pt_MyClass->int_fdf(x, p, f, df);
}


#endif
