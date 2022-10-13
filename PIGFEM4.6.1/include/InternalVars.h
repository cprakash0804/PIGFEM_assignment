/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _INTERNAL_VARS_H_
#define _INTERNAL_VARS_H_

#include "common.h"
#include <vector>
#include <map>
#include <string>


// Predeclarations
class Problem;


class InternalVars
{
	private:

		/*
		 * Pointer to the problem object associated with this InternalVars object
		 */
		Problem* _prob;


		/*
		 * Boolean for whether or not this contains the derivative of internal variables or just regular internal vars
		 */
		bool _sensitivity;
		bool _has_int_vars;
		std::vector<bool> _elem_isvs;


		/*
		 * Storage for the internal variables
		 * First index is the local element number
		 * Second index is the element's quadrature point
		 * Third index is the internal variable index
		 */
		std::vector<std::vector<std::vector<double> > > _internal_vars;


		/*
		 * Internal variable maps
		 * Maps from the materials names to various values
		 *	First map contains initial values of each material's internal variables
		 *	Second map contains the maximum values of each material's internal variables
		 *	Third map contains how many internal variables each material has
		 */
		std::map<std::string, std::vector<double> > _internal_vars_init;
		std::map<std::string, std::vector<double> > _internal_vars_max;
		std::map<std::string, id_type> _internal_var_number;

		/*
		 *  Wrapper function to add all material's internal variables to the problem
		 */
		void add_material_internal_variables();


	public:


		/*
		 * The Big 3
		 * Constructor, needs the problem object associated with this object
		 * Destructor doesn't do anything
		 * Copy Constructor copies 
		 */
		InternalVars(Problem* prob, bool sensitivity);
		~InternalVars();
		InternalVars(const InternalVars& other);


		/*
		 * Used to copy over current internal variable data
		 */
		 void copy_data(std::vector<std::vector<std::vector<double> > >& storage);
		 void copy_data_in(std::vector<std::vector<std::vector<double> > >& storage);


		/*
		 * Return the problem object associated with this
		 */
		Problem* get_prob() const {return _prob;};


		/*
		 * Preallocates the private data members based on the Problem pointer
		 */
		void preallocate_internal_vars();
		void preallocate_internal_vars_object(std::vector<std::vector<std::vector<double> > >& storage);


		/*
		 * Get the initial values for all internal variables associated with a material name
		 */
		std::vector<double> get_internal_var_init(std::string name);


		/*
		 * Get the maximum values for all internal variables associated with a material name
		 */
		std::vector<double> get_internal_var_max(std::string name);


		/*
		 * Iterator definition
		 * Used to find all materials and their internal variables
		 */
		typedef std::map<std::string, id_type>::iterator internal_var_iterator;
		internal_var_iterator internal_var_begin() {return _internal_var_number.begin();};
		internal_var_iterator internal_var_end() {return _internal_var_number.end();};

		
		/*
		 * Return a writeable reference to the current internal variable value
		 */
		double& get_internal_var_local(id_type idx, id_type qp, id_type var);
		double& get_internal_var_global(id_type id, id_type qp, id_type var);
		double& get_internal_var(id_type id, id_type qp, id_type var) {return get_internal_var_global(id, qp, var);};


		/*
		 * Return a writable reference to the list of internal variables at a quadrature point
		 */
		std::vector<double>& get_internal_vars_local(id_type idx, id_type qp);
		std::vector<double>& get_internal_vars_global(id_type id, id_type qp);


		/*
		 * Return a writeable reference to all of the internal variables for an element
		 */
		std::vector<std::vector<double> >& get_elem_internal_vars_local(id_type idx);
		std::vector<std::vector<double> >& get_elem_internal_vars_global(id_type id);


		/*
		 * Return a writeable reference to the entire internal variable object (This is actually terrible)
		 */
		std::vector<std::vector<std::vector<double> > >& get_all_internal_vars() {return _internal_vars;};


		/*
		 * Returns whether or no this ISV object actually is associated with any ISVs
		 */
		bool haveISVs() {return _has_int_vars;};
		bool elemHasISVs(id_type l_elem) {return _elem_isvs[l_elem];};

};





#endif
