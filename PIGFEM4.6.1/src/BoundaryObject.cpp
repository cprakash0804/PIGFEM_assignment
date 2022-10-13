/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated September 2016

##################################################################################
*/
#include "BoundaryObject.h"
#include "Mesh.h"
#include "Problem.h"
#include "node.h"
#include <algorithm>

BoundaryObject::BoundaryObject()
	: _prob(NULL)
{
}


void BoundaryObject::clear()
{
	_prob = NULL;
	_dirichlet_bcs.clear();
	_individuals.clear();
	_dirichlet_nodeset_bcs.clear();
}


// Function to apply an unapplied BCs
void BoundaryObject::apply_BCs()
{
	// Clear any old boundary conditions
	_dirichlet_bcs.clear();
	Mesh* mesh = _prob->get_mesh();

	// Apply all BCs applied to a nodeset
	for (auto it=_dirichlet_nodeset_bcs.begin(), end=_dirichlet_nodeset_bcs.end(); it!=end; ++it)
	{
		for (id_type i=0; i<it->second.size(); ++i)
		{
			id_type dof = it->second[i].first;
			double val = it->second[i].second;
			// Grab a refernece to the set
			std::set<id_type> set = mesh->get_nodeset(it->first);
			for (auto it2=set.begin(), end2=set.end(); it2!=end2; ++it2)
			{
				id_type node = *it2;
				_dirichlet_bcs[node][dof] = val;
			}
		}
	}

	// All BCs applied to individual nodes
	for (auto it=_individuals.begin(), end=_individuals.end(); it!=end; ++it)
		for (auto it2=it->second.begin(), end2=it->second.end(); it2!=end2; ++it2)
			_dirichlet_bcs[it->first][it2->first] = it2->second;

	// If I want to lock down any direction completely
	for (auto it=_lock_down_directions.begin(), end=_lock_down_directions.end(); it!=end; ++it)
	{
		// Apply bcs to all normal nodes
		for (Mesh::node_iterator it2=mesh->nodes_begin(), end2=mesh->nodes_end(); it2!=end2; ++it2)
		{
			id_type id = (*it2)->get_id();
			if (mesh->check_node_responsibility( id ))
				_dirichlet_bcs[id][*it] = 0.0;
		}

		// Apply bcs to all enriched nodes
		for (Mesh::enrich_node_iterator it2=mesh->enrich_nodes_begin(), end2=mesh->enrich_nodes_end(); it2!=end2; ++it2)
		{
			id_type id = (*it2)->get_id();
			if (mesh->check_node_responsibility( id ))
				_dirichlet_bcs[id][*it] = 0.0;
		}
	}
}







/*
 * Functions to add a prescribed dof to the mesh
 * 1. Add a single dirichlet bc so a single node and single dof
 * 2. Add a vector of dirichlet bcs to a vector of nodes and vector of dofs
 * 3. Add a single dirichlet bc to a vector of nodes on a single dof
 * 4. Add a vector of dirichlet bcs to a single nodeset (number) and a vector of dofs (FIXME: Don't know how to make agnostic of partitioning)
 * 5. Add a vector of dirichlet bcs to a single nodeset (name) and a vector of dofs (FIXME: Don't know how to make agnostic of partitioning)
 * 6. Add a single dirichlet bc to a single nodeset (number) and a single dof
 * 7. Add a single dirichlet bc to a single nodeset (name) and a single dof
*/

// 1. Add a single dirichlet bc so a single node and single dof
void BoundaryObject::set_dirichlet_bc(id_type node, id_type dof, double val) // FIXME: should check if dof is less than the number of dofs existsing for it
{
	if (dof > _prob->nndof())
		err_message("Degree-of-freedom " << dof << " does nt exist in this problem");

	Mesh* mesh = _prob->get_mesh();
	if ( mesh->node_in_local_mesh(node) ) // This mesh has a copy of this node
	{
		if ( mesh->own_node(node) ) // I own this node
		{
			if (_individuals.find(node)==_individuals.end()) // Make sure there at least an empty map to insert into
				_individuals[node] = std::map<id_type, double>();

			_individuals[node][dof] = val;
		}
	}
}
// 2. Add a vector of dirichlet bcs to a vector of nodes and vector of dofs
void BoundaryObject::set_dirichlet_bcs(std::vector<id_type>& nodes, std::vector<id_type>& dofs, std::vector<double>& vals)
{
	if(nodes.size()!=dofs.size() || dofs.size()!=vals.size())
		err_message("Vectors for creation of a boundary condition must be the same size.");
	
	for(id_type i=0; i<nodes.size(); ++i)
		set_dirichlet_bc(nodes[i], dofs[i], vals[i]);
}
// 3. Add a single dirichlet bc to a vector of nodes on a single dof
void BoundaryObject::set_dirichlet_bcs(std::vector<id_type>& nodes, id_type dof, double val)
{
	for(id_type i=0; i<nodes.size(); ++i)
		set_dirichlet_bc(nodes[i], dof, val);
}
/*
// 4. Add a vector of dirichlet bcs to a single nodeset (number) and a vector of dofs
void Mesh::set_dirichlet_bcs_from_nodeset(id_type set, std::vector<id_type>& dofs, std::vector<double>& vals)
{
	if(!_init) err_message("Must initiailize the mesh prior to calling set_dirichlet_bcs_from_nodeset.");
	
	if(set>=_nodesets.size())
		err_message("Please select a valid nodeset.");
	
	if(_nodesets[set].size()!=dofs.size() || dofs.size()!=vals.size())
		err_message("Vectors for creation of a boundary condition must be the same size.");
	
	NamedSet::iterator it = _nodesets[set].begin();
	NamedSet::iterator end = _nodesets[set].end();
	int i = 0;
	for(; it!=end; ++it)
	{
		set_dirichlet_bc((*it), dofs[i], vals[i]);
		i++;
	}
}
*/
/*
// 5. Add a vector of dirichlet bcs to a single nodeset (name) and a vector of dofs
void BoundaryObject::set_dirichlet_bcs_from_nodeset(std::string set_name, std::vector<id_type>& dofs, std::vector<double>& vals)
{	
	if(_nodesets[set_name].size()!=dofs.size() || dofs.size()!=vals.size())
		err_message("Vectors for creation of a boundary condition must be the same size.");

	if(_nodesets.find(set_name) == _nodesets.end())
		err_message("Please select a valid nodeset.");
	else
	{
		if(_nodesets[set_name].size() != dofs.size())
			err_message("There must be only one BC per node in the local mesh for the given nodeset.");

		int i=0;
		for(auto it=_nodesets[set_name].begin(), end=_nodesets[set_name].end(); it!=end; ++it, ++i)
			set_dirichlet_bc((*it), dofs[i], vals[i]);
	}
}
*/
// 7. Add a single dirichlet bc to a single nodeset (name) and a single dof
void BoundaryObject::set_dirichlet_bcs_from_nodeset(std::string set_name, id_type dof, double val)
{
	std::transform(set_name.begin(), set_name.end(), set_name.begin(), ::toupper); // Capitilize the name
	if (dof > _prob->nndof())
		err_message("Degree-of-freedom " << dof << " does nt exist in this problem");

	if ( _prob->get_mesh()->nodesetExists(set_name) )
		_dirichlet_nodeset_bcs[set_name].push_back(std::pair<id_type, double>(dof, val));
	else
		err_message("Nodeset " << set_name << " does not exist");
}


// Function to determine if a given node and dof has a boundary condition
bool BoundaryObject::has_dirichlet_bc(id_type node, id_type dof)
{
	if(_dirichlet_bcs.find(node)==_dirichlet_bcs.end()) // Node either isn't in the local mesh or it doesn't have any bcs on it
		return false;
	else
	{
		if(_dirichlet_bcs[node].find(dof)==_dirichlet_bcs[node].end()) // Node doesn't have a bc on that dof
			return false;
		else
			return true;
	}
}

// Function to get the prescribed value. Maybe should just define an iterator?
double BoundaryObject::get_dirichlet_bc(id_type node, id_type dof, bool& exists)
{
	exists = has_dirichlet_bc(node, dof);
	if(exists)
		return _dirichlet_bcs[node][dof];
	else
		return 0.0; // Depend on checking exists to see if this is actually a boundary condition
}

id_type BoundaryObject::n_dirichlet_on_node(id_type node)
{
	if(_dirichlet_bcs.find(node)!=_dirichlet_bcs.end())
		return _dirichlet_bcs[node].size();
	else
		return 0;
}
