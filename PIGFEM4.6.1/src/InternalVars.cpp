/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "InternalVars.h"
#include "Problem.h"
#include "Mesh.h"
#include "elem.h"
#include "CohesiveElem.h"




/*
 * The Big 3
 * Constructor, needs the problem object associated with this object
 * Destructor doesn't do anything
 * Copy Constructor copies 
 */
InternalVars::InternalVars(Problem* prob, bool sensitivity)
	: _prob(prob), _sensitivity(sensitivity), _has_int_vars(false)
{
	preallocate_internal_vars();
}
InternalVars::~InternalVars()
{}
InternalVars::InternalVars(const InternalVars& other)
{
	_prob = other.get_prob();
}




/*
 * Used to copy current internal variable data intor storage struct
 */
 void InternalVars::copy_data(std::vector<std::vector<std::vector<double> > >& storage)
 {
 	for (id_type e=0; e<_internal_vars.size(); ++e)
 	{
 		if (e > storage.size())
 			err_message("Invalid number of elements in internal variable storage");

 		for (id_type qp=0; qp<_internal_vars[e].size(); ++qp)
 		{
 			if (qp > storage[e].size())
 				err_message("Invalid number of quadrature points in internal variable storage");

 			for (id_type iv=0; iv<_internal_vars[e][qp].size(); ++iv)
 			{
 				if (iv > storage[e][qp].size())
 					err_message("Invalid number of internal variables in internal variable storage");

 				storage[e][qp][iv] = _internal_vars[e][qp][iv];
 			}
 		}
 	}
 }
  void InternalVars::copy_data_in(std::vector<std::vector<std::vector<double> > >& storage)
 {
 	for (id_type e=0; e<_internal_vars.size(); ++e)
 	{
 		if (e > storage.size())
 			err_message("Invalid number of elements in internal variable storage");

 		for (id_type qp=0; qp<_internal_vars[e].size(); ++qp)
 		{
 			if (qp > storage[e].size())
 				err_message("Invalid number of quadrature points in internal variable storage");

 			for (id_type iv=0; iv<_internal_vars[e][qp].size(); ++iv)
 			{
 				if (iv > storage[e][qp].size())
 					err_message("Invalid number of internal variables in internal variable storage");

 				_internal_vars[e][qp][iv] = storage[e][qp][iv];
 			}
 		}
 	}
 }




/*
 *  Wrapper function to add all material's internal variables to the problem
 */
 
 
void InternalVars::add_material_internal_variables()
{
	// Get all of the materials of the mesh
	// Every mesh partition has a copy of every material for consiistency's sake
	std::vector<Material*> mats = _prob->get_mesh()->get_materials();

	// Loop over every material and add its internal variable structures to the problem
	for(id_type m=0; m<mats.size(); ++m)
	{
		Material* mat = mats[m];
		_internal_var_number.insert(std::pair<std::string, id_type>(mat->get_name(), mat->n_internal_vars()));
		_internal_vars_init.insert(std::pair<std::string, std::vector<double> >(mat->get_name(), mat->init_internal_vars()));
		_internal_vars_max.insert(std::pair<std::string, std::vector<double> >(mat->get_name(), mat->max_internal_vars()));
	}
}




/*
 * Preallocates the private data members based on the Problem pointer
 */
void InternalVars::preallocate_internal_vars()
{
	preallocate_internal_vars_object(_internal_vars);
	if (!_sensitivity)
		add_material_internal_variables();
}
void InternalVars::preallocate_internal_vars_object(std::vector<std::vector<std::vector<double> > >& storage)
{
	Mesh* mesh = _prob->get_mesh();

	storage.resize(mesh->n_local_elem());
	_elem_isvs.resize(mesh->n_local_elem());
	std::fill(_elem_isvs.begin(), _elem_isvs.end(), false);
	for(Mesh::element_iterator it=mesh->elements_begin(), end=mesh->elements_end(); it!=end; ++it) // Not this is element instead of active_element because of refienement and coarsening processes
	{
		id_type local_e = mesh->global_to_local_elem((*it)->get_id());
		Material* mat;
		if (!(*it)->is_intersected()) // non-intersected
		{
			// Get the material associated witht this element
			mat = mesh->get_element_material_global((*it)->get_id());
			if (mat->n_internal_vars() != 0)
				_has_int_vars = true;

			id_type nqp = (*it)->n_q_points();
			storage[local_e].resize(nqp);
			for(id_type qp=0; qp<nqp; ++qp)
			{
				if (!_sensitivity)
					storage[local_e][qp] = mat->init_internal_vars();
				else
					storage[local_e][qp].resize(mat->n_internal_vars() * _prob->n_SensitivityParameters(), 0.0);
				if (storage[local_e][qp].size() != 0)
					_elem_isvs[local_e] = true;
			}
		}



		else // intersected
		{
			// integration elements
			for(id_type ie=0; ie<(*it)->n_integration_elem(); ++ie)
			{
				Elem* int_el = (*it)->get_integration_elem(ie);
				mat = mesh->get_element_material_global((*it)->get_id(), ie);
				if (mat->n_internal_vars() != 0)
					_has_int_vars = true;

				id_type nqp = int_el->n_q_points();
				id_type curr_qp = storage[local_e].size();
				storage[local_e].resize(curr_qp+nqp);
				for(id_type qp=curr_qp; qp<(curr_qp+nqp); ++qp)
				{
					if (!_sensitivity)
						storage[local_e][qp] = mat->init_internal_vars();
					else
						storage[local_e][qp].resize(mat->n_internal_vars() * _prob->n_SensitivityParameters(), 0.0);
					if (storage[local_e][qp].size() != 0)
						_elem_isvs[local_e] = true;
				}
			}

			// cohesive elements
			for(id_type ce=0; ce<(*it)->n_cohesive_elem(); ++ce)
			{
				CohesiveElem* coh_el = (*it)->get_cohesive_elem(ce);
				mat = (*it)->get_cohesive_material(ce);
				if (mat->n_internal_vars() != 0)
					_has_int_vars = true;

				id_type nqp = coh_el->n_q_points();
				id_type curr_qp = storage[local_e].size();
				storage[local_e].resize(curr_qp+nqp);
				for(id_type qp=curr_qp; qp<(curr_qp+nqp); ++qp)
				{
					if (!_sensitivity)
						storage[local_e][qp] = mat->init_internal_vars();
					else
						storage[local_e][qp].resize(mat->n_internal_vars() * _prob->n_SensitivityParameters(), 0.0);
					if (storage[local_e][qp].size() != 0)
						_elem_isvs[local_e] = true;
				}
			}
		}
	}
}




/*
 * Get the initial values for all internal variables associated with a material name
 */
std::vector<double> InternalVars::get_internal_var_init(std::string name)
{
	std::map<std::string, std::vector<double> >::iterator it = _internal_vars_init.find(name);
	if(it != _internal_vars_init.end())
		return it->second;
	else
		err_message("Invlaid material name.");
}




/*
 * Get the maximum values for all internal variables associated with a material name
 */
std::vector<double> InternalVars::get_internal_var_max(std::string name)
{
	std::map<std::string, std::vector<double> >::iterator it = _internal_vars_max.find(name);
	if(it != _internal_vars_max.end())
		return it->second;
	else
		err_message("Invlaid material name.");
}



		
/*
 * Return a writeable reference to the current internal variable value
 */
double& InternalVars::get_internal_var_local(id_type idx, id_type qp, id_type var)
{
	if(idx >= _internal_vars.size())
		err_message("Internal variables have not been allocated yet. Please allocate before attempting to set an internal variable.");
	if(qp >= _internal_vars[idx].size())
		err_message("The given quadrature point exceeds the number of quadrature point for the given element.");
	if(var >= _internal_vars[idx][qp].size())
	{
		char  buf[100];
		sprintf(buf, "The given internal variable: (%d, %d, %d) does not exist!", idx, qp, var);
		err_message( buf );
	}
	
	return _internal_vars[idx][qp][var];
}
double& InternalVars::get_internal_var_global(id_type id, id_type qp, id_type var)
{
	id_type idx = _prob->get_mesh()->global_to_local_elem(id);
	return get_internal_var_local(idx, qp, var);
}




/*
 * Return a writable reference to the list of internal variables at a quadrature point
 */
std::vector<double>& InternalVars::get_internal_vars_local(id_type idx, id_type qp)
{
	if(idx >= _internal_vars.size())
		err_message("Internal variables have not been allocated yet. Please allocate before attempting to set an internal variable.");
	if(qp >= _internal_vars[idx].size())
		err_message("The given quadrature point exceeds the number of quadrature point for the given element.");

	return _internal_vars[idx][qp];
}
std::vector<double>& InternalVars::get_internal_vars_global(id_type id, id_type qp)
{
	id_type idx = _prob->get_mesh()->global_to_local_elem(id);
	return get_internal_vars_local(idx, qp);
}




/*
 * Return a writeable reference to all of the internal variables for an element
 */
std::vector<std::vector<double> >& InternalVars::get_elem_internal_vars_local(id_type idx)
{
	if(idx >= _internal_vars.size())
		err_message("Internal variables have not been allocated yet. Please allocate before attempting to set an internal variable.");
	else
		return _internal_vars[idx];
}
std::vector<std::vector<double> >& InternalVars::get_elem_internal_vars_global(id_type id)
{
	id_type idx = _prob->get_mesh()->global_to_local_elem(id);
	return get_elem_internal_vars_local(idx);
}
