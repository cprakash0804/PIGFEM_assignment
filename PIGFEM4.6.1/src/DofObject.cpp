/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "DofObject.h"
#include "Mesh.h"
#include "Problem.h"
#include "BoundaryObject.h"
#include "elem.h"
#include "CohesiveElem.h"
#include "Utilities.h"

DofObject::DofObject()
	: _prob(NULL), _mesh(NULL)
{
}

DofObject::DofObject(const DofObject& other)
{
	clear();
	copy(other);
}

void DofObject::copy(const DofObject& other)
{
	_prob = other._prob;
	_mesh = other._mesh;
	_nndof = other._nndof;
	_node_global_dofs = other._node_global_dofs;
	_node_local_dofs = other._node_local_dofs;
	_global_local_dof = other._global_local_dof;
	_local_global_dof = other._local_global_dof;
	_global_dof_to_node = other._global_dof_to_node;
	_n_local_owned_free_dofs = other._n_local_owned_free_dofs;
	_n_local_owned_const_dofs = other._n_local_owned_const_dofs;
	_n_local_free_dofs = other._n_local_free_dofs;
	_n_local_const_dofs = other._n_local_const_dofs;
	_n_global_free_dofs = other._n_global_free_dofs;
	_n_global_const_dofs = other._n_global_const_dofs;
	 _dof_distributed = other._dof_distributed;
	_enrich_node_global_dofs = other._enrich_node_global_dofs;
	_enrich_node_local_dofs = other._enrich_node_local_dofs;
}

void DofObject::clear()
{
	_nndof = 0;
	_node_global_dofs.clear();
	_node_local_dofs.clear();
	_global_local_dof.clear();
	_local_global_dof.clear();
	_global_dof_to_node.clear();
	_n_local_owned_free_dofs = 0;
	_n_local_owned_const_dofs = 0;
	_n_local_free_dofs = 0;
	_n_local_const_dofs = 0;
	_n_global_free_dofs = 0;
	_n_global_const_dofs = 0;
	 _dof_distributed = false;
	_enrich_node_global_dofs.clear();
	_enrich_node_local_dofs.clear();
}


void DofObject::attachProblem(Problem* prob)
{
	_prob = prob;
	_mesh = _prob->get_mesh(); // for convenience
}

// Returns a pair containing the globa node number of local dof# corresponding to the global dof passed in
std::pair<id_type, id_type> DofObject::get_node_from_dof(id_type dof)
{
	auto it = _global_dof_to_node.find(dof);
	if(it != _global_dof_to_node.end())
		return it->second;
	else
		err_message("Attempting to access a node-dof pair not existing on the local partition.");
}

// Attempts to find the global dof number associated with the node with the given global node id and the dof'th dof on the node
id_type DofObject::get_global_dof(id_type id, id_type dof)
{
	if(!_dof_distributed) err_message("Must distribute degrees of freedom prior to calling get_global_dof.");
	
	if(dof<_nndof)
	{
		if(id < ENRICH_START) // A normal node
		{
			if(id < _mesh->n_global_nodes())
				return _node_global_dofs[_mesh->global_to_local_node(id)][dof];
			else
				err_message("Attempted to access global dof of a node not in the global mesh.");
		}
		else if((id-ENRICH_START) < _mesh->n_global_enrich_nodes()) // An enriched node
			return _enrich_node_global_dofs[_mesh->global_to_local_node(id)-ENRICH_START][dof];
		else
			err_message("Attempted to access global dof of a node not in the global mesh.");
	}
	else
		err_message("The global dof requested for the given global node id does not exist on the given node.");
}


// Attempts to find the local dof number associated with the node with the given global node id and the dof'th dof on the node
id_type DofObject::get_local_dof(id_type id, id_type dof)
{
	if(!_dof_distributed) err_message("Must distribute degrees of freedom prior to calling get_global_dof.");
	
	if(dof<_nndof)
	{
		if(id < ENRICH_START) // A normal node
		{
			if(id < _mesh->n_global_nodes())
				return _node_local_dofs[_mesh->global_to_local_node(id)][dof];
			else
				err_message("Attempted to access local dof of a node not in the global mesh.");
		}
		else if((id-ENRICH_START) < _mesh->n_global_enrich_nodes()) // An enriched node
			return _enrich_node_local_dofs[_mesh->global_to_local_node(id)-ENRICH_START][dof];
		else
			err_message("Attempted to access local dof of a node not in the global mesh.");
	}
	else
		err_message("The local dof requested for the given global node id does not exist on the given node.");
}


// Returns the global dofs associated with the local node idx
std::vector<id_type>& DofObject::get_nodal_global_dofs_local(id_type idx)
{
	if(!_dof_distributed) err_message("Must distribute degrees of freedom prior to calling get_nodal_global_dofs_local.");
	
	if(idx < ENRICH_START)
	{
		if(idx < _mesh->n_local_nodes())
			return _node_global_dofs[idx];
		else
			err_message("Attempted to access global dofs of a node not in the local mesh.");
	}
	else if((idx-ENRICH_START) < _mesh->n_local_enrich_nodes())
		return _enrich_node_global_dofs[idx-ENRICH_START];
	else
		err_message("Attempted to access global dofs of a node not in the local mesh.");
}


// Returns the global dofs associated with the global node id
std::vector<id_type>& DofObject::get_nodal_global_dofs_global(id_type id)
{
	id_type idx = _mesh->global_to_local_node(id);
	return get_nodal_global_dofs_local(idx);
}


// Returns the local dofs associated with the local node idx
std::vector<id_type>& DofObject::get_nodal_local_dofs_local(id_type idx)
{
	if(!_dof_distributed) err_message("Must distribute degrees of freedom prior to calling get_nodal_local_dofs_local.");
	
	if(idx < ENRICH_START)
	{
		if(idx < _mesh->n_local_nodes())
			return _node_local_dofs[idx];
		else
			err_message("Attempted to access local dofs of a node not in the local mesh.");
	}
	else if((idx-ENRICH_START) < _mesh->n_local_enrich_nodes())
		return _enrich_node_local_dofs[idx-ENRICH_START];
	else
		err_message("Attempted to access local dofs of a node not in the local mesh.");
}


// Returns the local dofs associated with the global node id
std::vector<id_type>& DofObject::get_nodal_local_dofs_global(id_type id)
{
	id_type idx = _mesh->global_to_local_node(id);
	return get_nodal_local_dofs_local(idx);
}


// Fills the passed in vector with all of the global dofs associated with all of the nodes of the local element idx
void DofObject::get_elem_global_dofs(Elem* el, std::vector<id_type> & vec)
{
	if(!_dof_distributed) err_message("Must distribute degrees of freedom prior to calling get_elem_global_dofs.");

	// Get rid of all previous entries in the dof vector
	vec.clear();

	// Grab a reference to the elem_node table for this element as well as a pointer to the element itself
	std::vector<id_type>& elem_node = _mesh->get_elem_node_global(el->get_id());

	// Add regualar dofs to the vector first
	for(id_type i=0; i<elem_node.size(); ++i)
	{
		id_type l_node = elem_node[i];
		vec.insert(vec.end(), _node_global_dofs[l_node].begin(), _node_global_dofs[l_node].end());
	}
	// Now add enriched dofs to the vector
	for(id_type i=0; i<el->n_enrich_nodes(); ++i)
	{
		id_type l_node = _mesh->global_to_local_node(el->get_enrich_node(i)->get_id()) - ENRICH_START;
		vec.insert(vec.end(), _enrich_node_global_dofs[l_node].begin(), _enrich_node_global_dofs[l_node].end());
	}
}
// Fills the passed in vector with all of the global dofs associated with all of the nodes of the local element idx
void DofObject::get_elem_global_dofs(CohesiveElem* el, std::vector<id_type> & vec)
{
	if(!_dof_distributed) err_message("Must distribute degrees of freedom prior to calling get_elem_global_dofs.");

	// Get rid of all previous entries in the dof vector
	vec.clear();

	// Add regualar dofs to the vector first
	for(id_type i=0; i<el->n_nodes(); ++i)
	{
		id_type l_node = _mesh->global_to_local_node(el->get_node(i)->get_id()) - ENRICH_START; // Assume cohesive elements will always be constructed from enrichment nodes
		vec.insert(vec.end(), _enrich_node_global_dofs[l_node].begin(), _enrich_node_global_dofs[l_node].end());
	}
}














void DofObject::distribute_dofs(id_type nndof)
{
	// CLear out any old distribution
	clear();

	_nndof = nndof;

	id_type partial_sum_free, partial_sum_pres;
	determine_local_info(partial_sum_free, partial_sum_pres);

	id_type local_dof_free, local_dof_pres, local_dof;
	std::set<id_type> primary_pnodes;
	set_local_owned_normal_dofs(partial_sum_free, partial_sum_pres,
								local_dof_free, local_dof_pres, local_dof,
								primary_pnodes);

	if (!_mesh->serial())
		communicate_normal_dofs(local_dof, primary_pnodes);

	if (_mesh->periodic())
		set_periodic_nodes(primary_pnodes);

	if (_mesh->IGFEM())
	{
		set_local_owned_enrich_dofs(partial_sum_free, partial_sum_pres,
									local_dof_free, local_dof_pres,local_dof,
									primary_pnodes);

		if (!_mesh->serial())
			communicate_enrich_dofs(local_dof);
	}

	Finalize();
	nndof++;
	nndof--;
}



void DofObject::determine_local_info(id_type& partial_sum_free, id_type& partial_sum_pres)
{
	BoundaryObject* boundary = _prob->get_boundary();

	// Normal nodes
	_n_local_owned_free_dofs = 0;
	_n_local_owned_const_dofs = 0;
	for(Mesh::node_iterator it=_mesh->nodes_begin(), end=_mesh->nodes_end(); it!=end; ++it)
	{
		if (_mesh->check_node_responsibility( (*it)->get_id() ))
		{
			id_type n_const = boundary->n_dirichlet_on_node((*it)->get_id());

			if(n_const > _nndof)  // Something went wrong
				err_message("More contrained dofs than available degrees of freedom on a node.");
				
			_n_local_owned_free_dofs += (_nndof-n_const);
			_n_local_owned_const_dofs += n_const;
		}
	}

	// Enriched nodes
	for(Mesh::enrich_node_iterator it=_mesh->enrich_nodes_begin(), end=_mesh->enrich_nodes_end(); it!=end; ++it)
	{
		if (_mesh->check_node_responsibility( (*it)->get_id() ))
		{
			id_type n_const = boundary->n_dirichlet_on_node((*it)->get_id());

			if(n_const > _nndof)  // Something went wrong
				err_message("More contrained dofs than available degrees of freedom on a node.");
				
			_n_local_owned_free_dofs += (_nndof-n_const);
			_n_local_owned_const_dofs += n_const;
		}
	}

	if (_mesh->serial())
	{
		_n_global_free_dofs = _n_local_owned_free_dofs;
		_n_global_const_dofs = _n_local_owned_const_dofs;
		partial_sum_free = 0;
		partial_sum_pres = 0;
	}
	else
	{
		MPI_Scan(&_n_local_owned_free_dofs, &partial_sum_free, 1, MPI_ID, MPI_SUM, _mesh->get_comm());
		MPI_Allreduce(&_n_local_owned_free_dofs, &_n_global_free_dofs, 1, MPI_ID, MPI_SUM, _mesh->get_comm());
		partial_sum_free -= _n_local_owned_free_dofs;   // Get total number BEFORE this rank
		MPI_Scan(&_n_local_owned_const_dofs, &partial_sum_pres, 1, MPI_ID, MPI_SUM, _mesh->get_comm());
		MPI_Allreduce(&_n_local_owned_const_dofs, &_n_global_const_dofs, 1, MPI_ID, MPI_SUM, _mesh->get_comm());
		partial_sum_pres -= _n_local_owned_const_dofs; // Get total number BEFORE this rank
	}
}



void DofObject::set_local_owned_normal_dofs(const id_type& partial_sum_free, const id_type& partial_sum_pres,
											id_type& local_dof_free, id_type& local_dof_pres, id_type& local_dof,
											std::set<id_type>& primary_pnodes)
{
	local_dof_free=0;
	local_dof_pres=0;
	local_dof = 0; // Doesn't matter if the dof is free or constrained in the local numbering
	if (_mesh->periodic()) // Preallocate space for local dofs
		_local_global_dof.resize((_mesh->n_local_non_periodic_nodes()+_mesh->n_local_enrich_nodes()) * _nndof);
	else
		_local_global_dof.resize((_mesh->n_local_nodes()+_mesh->n_local_enrich_nodes()) * _nndof);

	// Assign dofs to the regular nodes first
	for (Mesh::node_iterator it=_mesh->nodes_begin(), end=_mesh->nodes_end(); it!=end; ++it)
	{
		id_type node_id = (*it)->get_id();
		id_type l_node = _mesh->global_to_local_node(node_id);
		_node_global_dofs.push_back(std::vector<id_type>(_nndof));
		_node_local_dofs.push_back(std::vector<id_type>(_nndof));
		if (_mesh->check_node_responsibility( node_id )) // If I own this node
		{
			for (id_type j=0; j<_nndof; ++j)
			{
				// Add a global degree of freedom
				if ( _prob->get_boundary()->has_dirichlet_bc(node_id, j) )  // Have a constrained dof
				{
					_node_global_dofs[l_node][j] = _n_global_free_dofs + partial_sum_pres + local_dof_pres;   // All constrained dofs are numbered after free dofs
					local_dof_pres++;
				}
				else
				{
					_node_global_dofs[l_node][j] = partial_sum_free + local_dof_free;
					local_dof_free++;
				}
				// Add the corresponding local degree of freedom
				_node_local_dofs[l_node][j] = local_dof;
				
				// Add entries to various maps from global to local and such
				_local_global_dof[local_dof] = _node_global_dofs[l_node][j]; // Create the portion of the _local_global_dof map that corresponds to owned dofs
				_global_local_dof[_node_global_dofs[l_node][j]] = local_dof; // Insert the localy owned global dofs into the _global_local_dof map
				//_global_dof_to_node[_node_global_dofs[l_node][j]] = std::pair<id_type, id_type>(node_id, j);  // Add an entry to the global dof-> global node_id/dof# map
				local_dof++;
			}

			// If this was a periodic node then I need to set the dofs for the other nodes later
			if (_mesh->node_periodic(node_id))
				primary_pnodes.insert(node_id);
		}
	}
}



void DofObject::communicate_normal_dofs(id_type& local_dof, std::set<id_type>& primary_pnodes)
{
	id_type msg_size = 1 + _nndof;   // number of id_types that need to be sent per node
	std::map<int, id_type> n_expected;
	std::map<int, std::vector<id_type> > dof_send_lists;
	for (Mesh::partition_iterator it=_mesh->partition_interface_begin(), end=_mesh->partition_interface_end(); it!=end; ++it)
	{
		id_type id = it->first;
		int owner = _mesh->get_node_owner_global(id);
		if (_mesh->check_node_responsibility( id )) // If I set the dofs for this node
		{
			for (id_type i=0; i<(*it).second.size(); ++i)
			{
				int part = (*it).second[i];

				if (part!=_mesh->get_rank())
				{
					// Make sure I have an empty vector to insert into
					if (dof_send_lists.find(part)==dof_send_lists.end())
						dof_send_lists.insert(std::pair<int, std::vector<id_type> >(part, std::vector<id_type>()));

					dof_send_lists[part].push_back(id);
					for (id_type k=0; k<_nndof; ++k)
						dof_send_lists[part].push_back(_node_global_dofs[_mesh->global_to_local_node(id)][k]);
				}
			}
		}
		else if (owner != _mesh->get_rank())
		{
			if (_mesh->node_periodic(id))
			{
				if (_mesh->primary_periodic_node(id))
					n_expected[owner]++;
			}
			else
				n_expected[owner]++;
		}
	}
	
	// Send all of the dofs
	int n_sends = dof_send_lists.size();
	MPI_Request reqs[n_sends];
	int n_sent = 0;
	for(auto it=dof_send_lists.begin(), end=dof_send_lists.end(); it!=end; ++it)
		MPI_Isend(&(*it).second[0], (*it).second.size(), MPI_ID, (*it).first, 0, _mesh->get_comm(), &reqs[n_sent++]);

	// Recieve all the dofs from normal nodes
	for(auto it=n_expected.begin(), end=n_expected.end(); it!=end; ++it)
	{
		int part = (*it).first;
		std::vector<id_type> dof_recv_vec;		Utilities::RecieveUnknown(dof_recv_vec, part, 0, MPI_ID, _mesh->get_comm());
		id_type n_recvd = dof_recv_vec.size()/msg_size;
		if((*it).second != n_recvd)  // Safety check which hopefully shouldn't happen. Mostly for debugging
			err_message("During dof distribution, number of dofs recieved did not match number of dofs expected.");
		
		// Put these dofs in my dof data structure
		for (id_type n=0; n<n_recvd; ++n)
		{
			id_type n_id = dof_recv_vec[n*msg_size];

			// If this is a periodic node and isn't the primary one then I want to skip it
			if (_mesh->node_periodic(n_id))
			{
				if (!_mesh->primary_periodic_node(n_id))
					continue;
				else
					primary_pnodes.insert(n_id);
			}

			id_type l_node = _mesh->global_to_local_node(n_id);
			for (id_type d=0; d<_nndof; ++d)
			{
				id_type global_dof = dof_recv_vec[n*msg_size + d + 1];
				_node_global_dofs[l_node][d] = global_dof;
				_node_local_dofs[l_node][d] = local_dof;
				_local_global_dof[local_dof] = global_dof; // Create the portion of the _local_global_dof map that corresponds to remote dofs
				_global_local_dof[global_dof] = local_dof; // Insert the remotely owned global dofs into the _global_local_dof map
				//_global_dof_to_node[ _node_global_dofs[l_node][idx-1] ] = std::pair<id_type, id_type>(node_id, idx-1);  // Add an entry to the global dof-> global node_id/dof# map
				local_dof++;
			}
		}
	}

	// MPI Cleanup
	MPI_Waitall(n_sends, reqs, MPI_STATUSES_IGNORE);
}



void DofObject::set_periodic_nodes(const std::set<id_type>& primary_pnodes)
{
	// Loop over all of the primary periodic nodes and set the rest of the set's dofs
	for (auto it=primary_pnodes.begin(), end=primary_pnodes.end(); it!=end; ++it)
	{
		id_type n_id = *it;
		id_type l_node = _mesh->global_to_local_node(n_id);
		std::vector<id_type> set = _mesh->get_node_periodic_nodeset(n_id);
		for (id_type n=1; n<set.size(); ++n)
		{
			id_type l_node2 = _mesh->global_to_local_node(set[n]);
			_node_global_dofs[l_node2] = _node_global_dofs[l_node];
			_node_local_dofs[l_node2] = _node_local_dofs[l_node];
		}
	}
}




void DofObject::set_local_owned_enrich_dofs(const id_type& partial_sum_free, const id_type& partial_sum_pres,
											id_type& local_dof_free, id_type& local_dof_pres, id_type& local_dof,
											std::set<id_type>& primary_pnodes)
{
	for(Mesh::enrich_node_iterator it=_mesh->enrich_nodes_begin(), end=_mesh->enrich_nodes_end(); it!=end; ++it)
	{
		id_type node_id = (*it)->get_id();
		id_type l_node = _mesh->global_to_local_node(node_id);
		id_type l_e_node = l_node - ENRICH_START;
		_enrich_node_global_dofs.push_back(std::vector<id_type>(_nndof));
		_enrich_node_local_dofs.push_back(std::vector<id_type>(_nndof));
		if(_mesh->check_node_responsibility(node_id)) // If I own this node
		{
			for(id_type j=0; j<_nndof; ++j)
			{
				// Add a global degree of freedom
				if ( _prob->get_boundary()->has_dirichlet_bc(node_id, j) )  // Have a constrained dof
				{
					_enrich_node_global_dofs[l_e_node][j] = _n_global_free_dofs + partial_sum_pres + local_dof_pres;   // All constrained dofs are numbered after free dofs
					local_dof_pres++;
				}
				else
				{
					_enrich_node_global_dofs[l_e_node][j] = partial_sum_free + local_dof_free;
					local_dof_free++;
				}
				// Add the corresponding local degree of freedom
				_enrich_node_local_dofs[l_e_node][j] = local_dof;
				
				// Add entries to various maps from global to local and such
				_local_global_dof[local_dof] = _enrich_node_global_dofs[l_e_node][j]; // Create the portion of the _local_global_dof map that corresponds to owned dofs
				_global_local_dof[_enrich_node_global_dofs[l_e_node][j]] = local_dof; // Insert the localy owned global dofs into the _global_local_dof map
				//_global_dof_to_node[_enrich_node_global_dofs[l_e_node][j]] = std::pair<id_type, id_type>(node_id, j);  // Add an entry to the global dof-> global node_id/dof# map
				local_dof++;
			}
		}
	}
}




void DofObject::communicate_enrich_dofs(id_type& local_dof)
{
	// Loop over all of the enrichmnt nodes on the interface between 2 partitions
	id_type msg_size = 1 + _nndof;
	std::map<int, id_type> n_e_expected;
	std::map<int, std::vector<id_type> > dof_e_send_lists;
	for(Mesh::enrich_partition_iterator it=_mesh->enrich_partition_interface_begin(), end=_mesh->enrich_partition_interface_end(); it!=end; ++it)
	{
		id_type id = it->first;
		id_type l_node = _mesh->global_to_local_node(id);
		id_type l_e_node = l_node - ENRICH_START;
		int owner = _mesh->get_node_owner_local(l_node);
		if(owner==_mesh->get_rank())
		{
			for(id_type i=0; i<it->second.size(); ++i)
			{
				int part = it->second[i];
				if(part!=_mesh->get_rank())
				{
					// Make sure I have an empty vector to insert into
					if(dof_e_send_lists.find(part)==dof_e_send_lists.end())
						dof_e_send_lists.insert(std::pair<int, std::vector<id_type> >(part, std::vector<id_type>()));

					dof_e_send_lists[part].push_back( id );
					for(id_type j=0; j<_nndof; ++j)
						dof_e_send_lists[part].push_back(_enrich_node_global_dofs[l_e_node][j]);
				}
			}
		}
		else
			n_e_expected[owner]++;
	}

	// Send enriched dofs
	int n_sends_e = dof_e_send_lists.size();
	MPI_Request reqs_e[n_sends_e];
	int n_sent_e = 0;
	for(auto it=dof_e_send_lists.begin(), end=dof_e_send_lists.end(); it!=end; ++it)
		MPI_Isend(&(*it).second[0], (*it).second.size(), MPI_ID, (*it).first, 1, _mesh->get_comm(), &reqs_e[n_sent_e++]);

	// Recieve all the dofs from enrichment nodes
	for(auto it=n_e_expected.begin(), end=n_e_expected.end(); it!=end; ++it)
	{
		int part = (*it).first;
		std::vector<id_type> dof_recv_vec;		Utilities::RecieveUnknown(dof_recv_vec, part, 1, MPI_ID, _mesh->get_comm());
		id_type n_recvd = dof_recv_vec.size()/msg_size;
		if((*it).second != n_recvd)  // Safety check which hopefully shouldn't happen. Mostly for debugging
			err_message("During dof distribution, number of enriched dofs recieved did not match number of dofs expected.");
		
		// Put these dofs in my dof data structure
		for (id_type n=0; n<n_recvd; ++n)
		{
			id_type node_id = dof_recv_vec[n*msg_size];
			id_type l_e_node = _mesh->global_to_local_node(node_id) - ENRICH_START;
			for (id_type d=0; d<_nndof; ++d)
			{
				id_type global_dof = dof_recv_vec[n*msg_size + d + 1];
				_enrich_node_global_dofs[l_e_node][d] = global_dof;
				_enrich_node_local_dofs[l_e_node][d] = local_dof;
				_local_global_dof[local_dof] = global_dof; // Create the portion of the _local_global_dof map that corresponds to remote dofs
				_global_local_dof[global_dof] = local_dof; // Insert the remotely owned global dofs into the _global_local_dof map
				//_global_dof_to_node[ _enrich_node_global_dofs[l_e_node][idx-1] ] = std::pair<id_type, id_type>(node_id, idx-1);  // Add an entry to the global dof-> global node_id/dof# map
				local_dof++;
			}
		}
	}

	// MPI Cleanup
	MPI_Waitall(n_sends_e, reqs_e, MPI_STATUSES_IGNORE);
}



void DofObject::Finalize()
{
	// Figure out the number of free and constrained degrees of freedom that I have the nodes for, not necessarily just the ones that I own
	// So the sum of _n_local_free_dofs over all processors won't equal _n_global_free_dofs, it will be greater
	_n_local_const_dofs = 0;
	_n_local_free_dofs = 0;
	for(id_type i=0; i<_mesh->n_local_nodes(); ++i)
	{
		for(id_type j=0; j<_nndof; ++j)
		{
			if(_node_global_dofs[i][j] < _n_global_free_dofs)
				_n_local_free_dofs++;
			else
				_n_local_const_dofs++;
		}
	}
	_n_local_free_dofs += _mesh->n_local_enrich_nodes() * _nndof; // Because eriched dofs are always free
	
	_dof_distributed = true;
}
