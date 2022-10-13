/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "NodalData.h"
#include "Mesh.h"




/*
 * The Big 3.
 * Constructor requires a Problem pointer to be built
 * Destructor and Copy Constructor don't do anything
 */
NodalData::NodalData(Mesh* mesh)
	: _nndof(0), _preallocated(false), _mesh(mesh)
{
}




/*
 * Destructor doesn't do anything (vectors deallocate themselves)
 */
NodalData::~NodalData()
{}




/*
 * Copy constructor copies over problem and data
 */
NodalData::NodalData(const NodalData& other)
{
	_mesh = other.get_mesh();
	if (other.preallocated())
		preallocate_storage(other.nndof());

	_data.copy(other._data);
	_enrich_data.copy(other._enrich_data);
}




/*
 * Function preallocates storage based on the current state
 * of the Problem. May be called after mesh refinement to
 * add storage
 */
void NodalData::preallocate_storage(id_type nndof)
{
	_nndof = nndof;
	_data.resize(_mesh->n_local_nodes(), nndof);
	_enrich_data.resize(_mesh->n_local_enrich_nodes(), nndof);
	_preallocated = true;
}




/*
 * Get reference to the data at a local node id
 * Can be used to get or set data
 */
double& NodalData::get_value_local(id_type idx, id_type dof)
{
	if(idx < ENRICH_START)
	{
		if(idx < _mesh->n_local_nodes())
		{
			if(dof < _nndof)
				return _data(idx, dof);
			else
				err_message("Attempting to access the data for a degree of freedom not existing on the given node.");
		}
		else
			err_message("Attempting to access the data for a node not existing on the local mesh.");
	}
	else if((idx-ENRICH_START) < _mesh->n_local_enrich_nodes())
	{
		if(dof < _nndof)
			return _enrich_data(idx-ENRICH_START, dof);
		else
			err_message("Attempting to access the data for a degree of freedom not existing on the given node.");
	}
	else
		err_message("Attempting to access the data for a node not existing on the global mesh.");
}




/*
 * Get reference to the data at a global node id
 * Can be used to get or set data
 */
double& NodalData::get_value_global(id_type id, id_type dof)
{
	id_type idx = _mesh->global_to_local_node(id);
	return get_value_local(idx, dof);
}


/*
 * Get reference to the data at a local node id
 * Can be used to get or set data
 */
double NodalData::get_value_local(id_type idx, id_type dof) const
{
	if(idx < ENRICH_START)
	{
		if(idx < _mesh->n_local_nodes())
		{
			if(dof < _nndof)
				return _data(idx, dof);
			else
				err_message("Attempting to access the data for a degree of freedom not existing on the given node.");
		}
		else
			err_message("Attempting to access the data for a node not existing on the local mesh.");
	}
	else if((idx-ENRICH_START) < _mesh->n_local_enrich_nodes())
	{
		if(dof < _nndof)
			return _enrich_data(idx-ENRICH_START, dof);
		else
			err_message("Attempting to access the data for a degree of freedom not existing on the given node.");
	}
	else
		err_message("Attempting to access the data for a node not existing on the global mesh.");
}




/*
 * Get reference to the data at a global node id
 * Can be used to get or set data
 */
double NodalData::get_value_global(id_type id, id_type dof) const
{
	id_type idx = _mesh->global_to_local_node(id);
	return get_value_local(idx, dof);
}