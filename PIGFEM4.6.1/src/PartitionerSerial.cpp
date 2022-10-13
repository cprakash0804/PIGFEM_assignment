/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "PartitionerSerial.h"
#include "Mesh.h"
#include "elem.h"
#include "node.h"
#include "Utilities.h"



/*
 * The Big 3
 * The constructor takes the mesh that will be partitioned
 * Destructor doesn't do anything
 * Copy constructor copies over the mesh pointer
 */
PartitionerSerial::PartitionerSerial(Mesh* mesh)
	: Partitioner(mesh)
{}
PartitionerSerial::~PartitionerSerial()
{}
PartitionerSerial::PartitionerSerial(const PartitionerSerial& other)
{
	_mesh = other.get_mesh();
}




/*
 * The main function that is called to partition the mesh
 * Takes mesh that exists on one processor and partitions
 * it to n processors based on the METIS decomposition
 */
void PartitionerSerial::partition()
{
	// Running in serial, not much to do
	if (_mesh->serial())
		handle_serial_case();

	// Parallel, mesh actually has to be partitioned
	else
	{
		// Get what elements and nodes go to what processors
		std::vector<idx_t> epart, npart;
		get_decomposition(epart, npart);

		// Something to determine the node send lists
		std::vector<id_type> n_elem_per_proc, n_node_per_proc, eind_size;
		std::vector<std::vector<id_type> > node_to_parts;
		std::vector<std::set<id_type> > part_node_id_lists, periodic_sets;
		generate_node_lists(epart, n_elem_per_proc, n_node_per_proc, eind_size, node_to_parts, part_node_id_lists, periodic_sets);

		// Communicate all nodal information
		communicate_nodes(npart, n_node_per_proc, part_node_id_lists);

		// Communicate all element information
		communicate_elements(epart, n_elem_per_proc, eind_size);

		// Communicate the periodic nodesets associated with a periodic mesh
		if (_mesh->periodic())
			commmunicate_periodic_nodesets(part_node_id_lists, periodic_sets);

		// Communicate the information about the mesh partition niterfaces
		communicate_partition_interface(part_node_id_lists);

		// // Communicate all boundary conditions
		// communicate_boundary_conditions(npart);

		// Communicate the nodesets
		communicate_nodesets(node_to_parts);

		// Communicate the element sets
		communicate_elemsets(epart);

		// Communicate the sidesetsd
		communicate_sidesets(epart);

		// Communicate all material information (CURRENTLY BROKEN)
		//communicate_materials(epart, n_elem_per_proc);

		// Determine all global mesh information
		determine_global_info();
	}
}




/*
 * Fills in any mesh member variables that need to be
 * managed in the case of a serial mesh
 */
void PartitionerSerial::handle_serial_case() const
{
	_mesh->_node_owners.resize(_mesh->n_local_nodes());
	std::fill(_mesh->_node_owners.begin(), _mesh->_node_owners.end(), _mesh->_rank);
	_mesh->_n_global_nodes = _mesh->_nodes.size();
	_mesh->_n_global_elem = _mesh->_elem.size();
	_mesh->_n_local_owned_nodes = _mesh->_nodes.size();
	_mesh->_partitioned = true;
}




/*
 * Uses METIS to get the decomposition of the mesh
 */
void PartitionerSerial::get_decomposition(std::vector<idx_t>& epart, std::vector<idx_t>& npart) const
{
	// Only actually do this if this is proc 0
	if (_mesh->get_rank() == 0)
	{
		if (!_mesh->periodic())
		{
			// Set up some variables
			idx_t ne, nn, ncommon, nparts, objval;
			ne = _mesh->n_local_elem();
			nn = _mesh->n_local_nodes();
			ncommon = 1;
			nparts = _mesh->n_ranks();

			// Set up the mesh connectivity vectors
			std::vector<idx_t> eptr, eind;
			eptr.resize(ne + 1);
			int eind_size = 0;
			for(idx_t i=0; i<ne; ++i)
				eind_size += _mesh->get_elem_local(i)->n_nodes();
			eind.resize(eind_size);
			id_type curr_ind = 0;
			for(int i=0; i<ne; ++i)
			{
				for(id_type j=0; j<_mesh->_elem_node[i].size(); ++j)
					eind[curr_ind+j] = _mesh->_elem_node[i][j];

				curr_ind += _mesh->_elem_node[i].size();
				eptr[i+1] = curr_ind;
			}

			// Actually call the METIS partitioning routine	
			epart.resize(ne);
			npart.resize(nn);
			int ierr = METIS_PartMeshDual(&ne, &nn, eptr.data(), eind.data(), NULL, NULL, &ncommon, &nparts, NULL, NULL, &objval, epart.data(), npart.data());
			if (ierr)
			{
				// Gives me options to output error messages here
			}
		}
		else // We are working with a periodic mesh
		{
			// Set up some variables
			idx_t ne, nn, nnnp, ncommon, nparts, objval;
			ne = _mesh->n_local_elem();
			nn = _mesh->n_local_nodes();
			nnnp = _mesh->n_local_non_periodic_nodes();
			ncommon = 1;
			nparts = _mesh->n_ranks();

			// Create a map from nodal ids to the mesh local non-periodic node numbering
			std::unordered_map<id_type, id_type> id_to_non_periodic_local(nnnp);
			id_type count = 0;
			for (id_type n=0; n<(id_type)nn; ++n)
			{
				if (_mesh->node_periodic(n)) // This node does lie along a periodic boundary
				{
					if (_mesh->primary_periodic_node(n)) // If this is the first element in the periodic nodeset then I will use this node
					{
						id_to_non_periodic_local[n] = count;
						count++;
					}
					else
					{
						std::vector<id_type> set = _mesh->get_node_periodic_nodeset(n);
						id_to_non_periodic_local[n] = id_to_non_periodic_local[set[0]];
					}
				}
				else
				{
					id_to_non_periodic_local[n] = count;
					count++;
				}
			}

			// Set up the mesh connectivity vectors
			std::vector<idx_t> eptr, eind;
			eptr.resize(ne + 1);
			int eind_size = 0;
			for(idx_t i=0; i<ne; ++i)
				eind_size += _mesh->get_elem_local(i)->n_nodes();
			eind.resize(eind_size);
			id_type curr_ind = 0;
			for(int i=0; i<ne; ++i)
			{
				Elem* el = _mesh->get_elem_local(i);
				for(id_type n=0; n<el->n_nodes(); ++n)
				{
					id_type n_id = el->get_node(n)->get_id();
					if (_mesh->node_periodic(n_id)) // This node does lie along a periodic boundary
					{
						std::vector<id_type> periodic_set = _mesh->get_node_periodic_nodeset(n_id);
						eind[curr_ind+n] = id_to_non_periodic_local[periodic_set[0]]; // Just use the one that was stored first (should be the lowest id too)
					}
					else
						eind[curr_ind+n] = id_to_non_periodic_local[n_id];
				}

				curr_ind += _mesh->_elem_node[i].size();
				eptr[i+1] = curr_ind;
			}

			// Actually call the METIS partitioning routine	
			epart.resize(ne);
			std::vector<idx_t> npart_temp(nnnp);
			int ierr = METIS_PartMeshDual(&ne, &nnnp, eptr.data(), eind.data(), NULL, NULL, &ncommon, &nparts, NULL, NULL, &objval, epart.data(), npart_temp.data());
			if (ierr)
			{
				// Gives me options to output error messages here
			}

			// Fill in the npart vector wih the right owners now
			npart.resize(nn);
			for (id_type n=0; n<(id_type)nn; ++n)
			{
				if (_mesh->node_periodic(n)) // This node does lie along a periodic boundary
				{
					std::vector<id_type> periodic_set = _mesh->get_node_periodic_nodeset(n);
					npart[n] = npart_temp[id_to_non_periodic_local[periodic_set[0]]]; // Just use the one that was stored first (should be the lowest id too)
				}
				else
					npart[n] = npart_temp[id_to_non_periodic_local[n]];
			}
		}
	}
}




/*
 * Generate several variables that will be used in several contexts
 *	n_elem_per_proc: self explanatory... Number of elements per processor
 *	n_nodes_per_proc: The number of nodes that each processor will recieve
 *		(Not necessarily the same as the number of nodes assigned by METIS)
 *	eind_size: The total number of node points for all elements in each mesh partition
 *	node_to_parts: A vector for every node in the mesh containing which partitions it is being sent to
 *	part_node_id_lists: A set of unique node ids for each mesh partition. Each set contains all of the
 		nodes necessary to build all of the elements that the mesh partition owns.
 */
void PartitionerSerial::generate_node_lists(const std::vector<idx_t>& epart,
											std::vector<id_type>& n_elem_per_proc, std::vector<id_type>& n_nodes_per_proc,
											std::vector<id_type>& eind_size, std::vector<std::vector<id_type> >& node_to_parts,
											std::vector<std::set<id_type> >& part_node_id_lists,
											std::vector<std::set<id_type> >& periodic_sets) const
{
	if (_mesh->get_rank() == 0)
	{
		// Resize things appropriately
		id_type nparts = _mesh->n_ranks();
		n_elem_per_proc.resize(nparts);
		n_nodes_per_proc.resize(nparts);
		eind_size.resize(nparts);
		part_node_id_lists.resize(nparts);
		node_to_parts.resize(_mesh->n_local_nodes());

		// Loop over all elements, counting them and building the node id sets
		for(id_type e=0; e<_mesh->n_local_elem(); ++e)
		{
			n_elem_per_proc[epart[e]]++;
			Elem* el = _mesh->get_elem_local(e);
			eind_size[epart[e]] += el->n_nodes();
			for(id_type n=0; n<el->n_nodes(); ++n)
				part_node_id_lists[epart[e]].insert( el->get_node(n)->get_id() );  // Will add unique nodes to the set
		}

		if(_mesh->periodic())
		{
			// Determine which periodic sets each mesh partition has a copy of
			periodic_sets.resize(nparts);
			for (id_type p=0; p< part_node_id_lists.size(); ++p)
				for (auto it=part_node_id_lists[p].begin(), end=part_node_id_lists[p].end(); it!=end; ++it)
					if (_mesh->node_periodic(*it)) // this is a periodic node
						periodic_sets[p].insert(_mesh->_periodic_id_to_set[*it]);
			// Add all of the nodes in each periodic set to every partition that has a copy of it
			for (id_type p=0; p<nparts; ++p)
			{
				for (auto it=periodic_sets[p].begin(), end=periodic_sets[p].end(); it!=end; ++it)
				{
					std::vector<id_type> set = _mesh->_periodic_nodesets[*it];
					part_node_id_lists[p].insert(set.begin(), set.end());
				}
			}
		}

		for(id_type i=0; i<(nparts); ++i)
			n_nodes_per_proc[i] = part_node_id_lists[i].size();

		for(id_type rank=0; rank<nparts; ++rank)
			for(auto it=part_node_id_lists[rank].begin(), end=part_node_id_lists[rank].end(); it!=end; ++it)
				node_to_parts[*it].push_back(rank);
	}
}
	



/*
 * Scatter all of the nodal information from proc 0 to the other procs
 */
void PartitionerSerial::communicate_nodes(const std::vector<idx_t>& npart, const std::vector<id_type>& n_nodes_per_proc,
										  const std::vector<std::set<id_type> >& part_node_id_lists) const
{
	// Send the info
	if (_mesh->get_rank() == 0)
	{
		id_type nparts = _mesh->n_ranks();
		std::vector<MPI_Request> reqs(3*(nparts-1));
		std::vector<std::vector<id_type> > node_id_lists(nparts);
		std::vector<std::vector<double> > node_coord_lists(nparts);
		std::vector<std::vector<int> > node_owner_lists(nparts);
		for(id_type i=0; i<(nparts); ++i)
		{
			// Preallocate space
			node_id_lists[i].resize( n_nodes_per_proc[i] );
			node_coord_lists[i].resize( n_nodes_per_proc[i]*3 );
			node_owner_lists[i].resize( n_nodes_per_proc[i] );

			int count = 0;
			for(auto j=part_node_id_lists[i].begin(), end=part_node_id_lists[i].end();j!=end; ++j)
			{
				node_id_lists[i][count] = _mesh->get_node_local(*j)->get_id();
				node_owner_lists[i][count] = npart[*j];
				node_coord_lists[i][3*count] = (*(_mesh->get_node_local(*j)))(0);
				node_coord_lists[i][3*count+1] = (*(_mesh->get_node_local(*j)))(1);
				node_coord_lists[i][3*count+2] = (*(_mesh->get_node_local(*j)))(2);
				count++;
			}
			// Now that the vector has been generated, might as well send it (non-blocking send)
			if (i!=0)
			{
				MPI_Isend(node_id_lists[i].data(), node_id_lists[i].size(), MPI_ID, i, 0, _mesh->get_comm(), &reqs[(i-1)+(0*(nparts-1))]);
				MPI_Isend(node_coord_lists[i].data(), node_coord_lists[i].size(), MPI_DOUBLE, i, 1, _mesh->get_comm(), &reqs[(i-1)+(1*(nparts-1))]);
				MPI_Isend(node_owner_lists[i].data(), node_owner_lists[i].size(), MPI_INT, i, 2, _mesh->get_comm(), &reqs[(i-1)+(2*(nparts-1))]);
			}
		}

		// Set up my partition
		build_nodes(node_id_lists[0], node_coord_lists[0], node_owner_lists[0]);

		// MPI Cleanup
		MPI_Waitall(3*(nparts-1), reqs.data(), MPI_STATUSES_IGNORE);
	}

	// Recieve the info
	else
	{
		std::vector<id_type> ids;	Utilities::RecieveUnknown(ids, 0, 0, MPI_ID, _mesh->get_comm());
		std::vector<double> coords;	Utilities::RecieveUnknown(coords, 0, 1, MPI_DOUBLE, _mesh->get_comm());
		std::vector<int> owners;	Utilities::RecieveUnknown(owners, 0, 2, MPI_INT, _mesh->get_comm());

		// Set up nodes
		build_nodes(ids, coords, owners);
	}
}
void PartitionerSerial::build_nodes(const std::vector<id_type>& ids, const std::vector<double>& coords, const std::vector<int>& owners) const
{
	for(id_type i=0; i<_mesh->n_local_nodes(); ++i)
		delete _mesh->_nodes[i];
	_mesh->_nodes.clear();
	_mesh->_global_to_local_node.clear();
	_mesh->_node_owners.clear();
	_mesh->_global_to_local_node.rehash(std::ceil(ids.size() / _mesh->_global_to_local_node.max_load_factor()));
	for(id_type i=0; i<ids.size(); ++i)
		_mesh->add_node(coords[i*3], coords[i*3+1], coords[i*3+2],
						ids[i], owners[i]);
	_mesh->_n_local_owned_nodes = count(_mesh->_node_owners.begin(), _mesh->_node_owners.end(), _mesh->get_rank());
}




/*
 * Scatter all of the elemental structure information from proc 0 to the other procs
 */
void PartitionerSerial::communicate_elements(const std::vector<idx_t>& epart, const std::vector<id_type>& n_elem_per_proc,
											 const std::vector<id_type>& eind_size) const
{
	// Send the informaion
	if (_mesh->get_rank() == 0)
	{
		id_type nparts = _mesh->n_ranks();
		std::vector<MPI_Request> reqs(4*(nparts-1));
		std::vector<std::vector<id_type> > eptr_send_lists(nparts, std::vector<id_type>());
		std::vector<std::vector<id_type> > eind_send_lists(nparts, std::vector<id_type>());
		std::vector<std::vector<id_type> > elem_id_lists(nparts);
		std::vector<std::vector<elem_type> > elem_type_lists(nparts);
		std::vector<id_type> curr_inds(nparts);
		std::vector<id_type> n_elem_done(nparts);
		for(id_type i=0; i<(nparts); ++i)
		{
			eptr_send_lists[i].resize(n_elem_per_proc[i] + 1);
			eind_send_lists[i].resize(eind_size[i]);
			elem_id_lists[i].resize(n_elem_per_proc[i]);
			elem_type_lists[i].resize(n_elem_per_proc[i]);
		}
		for(id_type i=0; i<_mesh->n_local_elem(); ++i)
		{
			Elem* el = _mesh->get_elem_local( i );
			id_type part = epart[i];

			// Elem-node connectivity
			for(id_type n=0; n<el->n_nodes(); ++n)
				eind_send_lists[part][curr_inds[part] + n] = _mesh->_elem_node[i][n];

			curr_inds[part] = curr_inds[part] + el->n_nodes();
			eptr_send_lists[part][n_elem_done[part]+1] = curr_inds[part];

			// IDs and types
			elem_id_lists[part][n_elem_done[part]] = el->get_id();
			elem_type_lists[part][n_elem_done[part]] = el->get_type();

			n_elem_done[part]++;
		}
		for(id_type i=1; i<(nparts); ++i)
		{
			MPI_Isend(eptr_send_lists[i].data(), eptr_send_lists[i].size(), MPI_ID, i, 0, _mesh->get_comm(), &reqs[(i-1)+(0*(nparts-1))]);
			MPI_Isend(eind_send_lists[i].data(), eind_send_lists[i].size(), MPI_ID, i, 1, _mesh->get_comm(), &reqs[(i-1)+(1*(nparts-1))]);
			MPI_Isend(elem_id_lists[i].data(), elem_id_lists[i].size(), MPI_ID, i, 2, _mesh->get_comm(), &reqs[(i-1)+(2*(nparts-1))]);
			MPI_Isend(elem_type_lists[i].data(), elem_type_lists[i].size(), MPI_INT, i, 3, _mesh->get_comm(), &reqs[(i-1)+(3*(nparts-1))]);
		}

		// Set up my partition
		build_elements(elem_id_lists[0], elem_type_lists[0], eptr_send_lists[0], eind_send_lists[0]);

		// MPI Cleanup
		MPI_Waitall(4*(nparts-1), reqs.data(), MPI_STATUSES_IGNORE);
	}

	// Recieve the information
	else
	{
		std::vector<id_type> eptr;		Utilities::RecieveUnknown(eptr, 0, 0, MPI_ID, _mesh->get_comm());
		std::vector<id_type> eind;		Utilities::RecieveUnknown(eind, 0, 1, MPI_ID, _mesh->get_comm());
		std::vector<id_type> ids;		Utilities::RecieveUnknown(ids, 0, 2, MPI_ID, _mesh->get_comm());
		std::vector<elem_type> types;	Utilities::RecieveUnknown(types, 0, 3, MPI_INT, _mesh->get_comm());
		if (ids.size()!=types.size() || ids.size()!=(eptr.size()-1))
			err_message("Inconsistent number of elements recieved.");

		// Create my elements!
		build_elements(ids, types, eptr, eind);
	}
}
void PartitionerSerial::build_elements(const std::vector<id_type>& ids, const std::vector<elem_type>& types,
									   const std::vector<id_type>& eptr, const std::vector<id_type>& eind) const
{
	for(id_type i=0; i<_mesh->n_local_elem(); ++i)
		delete _mesh->_elem[i];
	_mesh->_elem.clear();
	_mesh->_global_to_local_elem.clear();
	_mesh->_elem_node.clear();
	_mesh->_global_to_local_elem.rehash( std::ceil(ids.size() / _mesh->_global_to_local_elem.max_load_factor()) );
	for(id_type i=0; i<ids.size(); ++i)
	{
		id_type start = eptr[i];
		id_type stop = eptr[i+1];
		std::vector<id_type> node_ids(eind.begin()+start, eind.begin()+stop);
		_mesh->add_elem(types[i], node_ids, ids[i]);
	}
}




/*
 * Communicates the periodic nodesets of a periodic mesh
 */
void PartitionerSerial::commmunicate_periodic_nodesets(const std::vector<std::set<id_type> >& part_node_id_lists,
													   const std::vector<std::set<id_type> >& periodic_sets)
{
	// Send the informaion
	if (_mesh->get_rank() == 0)
	{
		id_type nparts = _mesh->n_ranks();
		std::vector<MPI_Request> reqs((nparts-1)*2);
		std::vector<std::vector<id_type> > periodic_ptr(nparts, std::vector<id_type>(1, 0));
		std::vector<std::vector<id_type> > periodic_ind(nparts);
		std::vector<id_type> curr_inds(nparts);
		for (id_type p=0; p< part_node_id_lists.size(); ++p)
		{
			for (auto it=periodic_sets[p].begin(), end=periodic_sets[p].end(); it!=end; ++it)
			{
				std::vector<id_type> set = _mesh->get_periodic_nodeset(*it);
				periodic_ind[p].insert(periodic_ind[p].end(), set.begin(), set.end());
				curr_inds[p] += set.size();
				periodic_ptr[p].push_back(curr_inds[p]);
			}
			if (p!=0)
			{
				MPI_Isend(periodic_ptr[p].data(), periodic_ptr[p].size(), MPI_ID, p, 0, _mesh->get_comm(), &reqs[(p-1)*2]);
				MPI_Isend(periodic_ind[p].data(), periodic_ind[p].size(), MPI_ID, p, 1, _mesh->get_comm(), &reqs[(p-1)*2+1]);
			}
		}

		// Set up my partition
		build_periodicity(periodic_ptr[0], periodic_ind[0]);

		// MPI Cleanup
		MPI_Waitall(2*(nparts-1), reqs.data(), MPI_STATUSES_IGNORE);
	}

	// Recieve the information
	else
	{
		std::vector<id_type> pptr;		Utilities::RecieveUnknown(pptr, 0, 0, MPI_ID, _mesh->get_comm());
		std::vector<id_type> pind;		Utilities::RecieveUnknown(pind, 0, 1, MPI_ID, _mesh->get_comm());

		// Create my periodicity structures!
		build_periodicity(pptr, pind);
	}
}
void PartitionerSerial::build_periodicity(const std::vector<id_type> periodic_ptr, const std::vector<id_type> periodic_ind)
{
	_mesh->_periodic_nodesets.clear();
	_mesh->_periodic_id_to_set.clear();

	id_type count = 0;
	for (id_type pset=0; pset<(periodic_ptr.size()-1); ++pset)
	{
		id_type start = periodic_ptr[pset];
		id_type stop = periodic_ptr[pset+1];
		_mesh->_periodic_nodesets.push_back(std::vector<id_type>(periodic_ind.begin()+start, periodic_ind.begin()+stop));

		for (id_type n=start; n<stop; ++n)
			_mesh->_periodic_id_to_set[periodic_ind[n]] = count;

		count++;
	}
}



/*
 * Scatter all of the information about the interface between
 * mesh partitions from proc 0 to the other procs
 */
void PartitionerSerial::communicate_partition_interface(const std::vector<std::set<id_type> >& part_node_id_lists) const
{
	if (_mesh->get_rank() == 0)
	{
		id_type nparts = _mesh->n_ranks();
		std::vector<MPI_Request> reqs(2*(nparts-1));
		std::vector<std::vector<id_type> > interface_nodes;
		std::vector<std::vector<int> > part_idx_map;
		id_type idx = 0;
		for(id_type i=0; i<(nparts-1); ++i)
		{
			part_idx_map.push_back(std::vector<int>(nparts-1-i));
			for(id_type j=(i+1); j<(nparts); ++j)
			{
				part_idx_map[i][j-(i+1)] = idx;
				interface_nodes.push_back( std::vector<id_type>(std::max(part_node_id_lists[i].size(), part_node_id_lists[j].size())) );
				std::vector<id_type>::iterator iter;
				iter = std::set_intersection(part_node_id_lists[i].begin(), part_node_id_lists[i].end(), part_node_id_lists[j].begin(), part_node_id_lists[j].end(), interface_nodes[idx].begin());
				interface_nodes[idx].resize(iter-interface_nodes[idx].begin());
				idx++;
			}
		}
		std::vector<std::vector<id_type> > pptr_send_lists;
		std::vector<std::vector<id_type> > pind_send_lists;
		for(id_type i=0; i<(nparts); ++i)
		{
			pptr_send_lists.push_back(std::vector<id_type>(nparts+1));
			pind_send_lists.push_back(std::vector<id_type>());
			int curr_ind = 0;
			for(id_type j=0; j<(nparts); ++j)
			{
				pptr_send_lists[i][j] = curr_ind;
				if (i<j) idx = part_idx_map[i][j-(i+1)];
				else if (i>j) idx = part_idx_map[j][i-(j+1)];
				else idx = nparts*(nparts-1)/2 + 1;
				if (idx < nparts*(nparts-1)/2) // Actual set intersection
				{
					for(id_type k=0; k<interface_nodes[idx].size(); ++k)
					{
						pind_send_lists[i].push_back(interface_nodes[idx][k]);
						curr_ind++;
					}
				}
				// else do nothing
			}
			pptr_send_lists[i][nparts] = curr_ind;

			// Now that the vectors are made might as well send them
			if (i!=0)
			{
				MPI_Isend(pptr_send_lists[i].data(), pptr_send_lists[i].size(), MPI_ID, i, 0, _mesh->get_comm(), &reqs[(i-1)+(0*(nparts-1))]);
				MPI_Isend(pind_send_lists[i].data(), pind_send_lists[i].size(), MPI_ID, i, 1, _mesh->get_comm(), &reqs[(i-1)+(1*(nparts-1))]);
			}
		}

		// Build my partition interface
		build_partition_interface(pptr_send_lists[0], pind_send_lists[0]);

		// MPI cleanup
		MPI_Waitall(2*(nparts-1), reqs.data(), MPI_STATUSES_IGNORE);
	}

	else
	{
		std::vector<id_type> pptr_recv_vec;		Utilities::RecieveUnknown(pptr_recv_vec, 0, 0, MPI_ID, _mesh->get_comm());
		std::vector<id_type> pind_recv_vec;		Utilities::RecieveUnknown(pind_recv_vec, 0, 1, MPI_ID, _mesh->get_comm());

		// Build my partition interface
		build_partition_interface(pptr_recv_vec, pind_recv_vec);
	}
}
void PartitionerSerial::build_partition_interface(const std::vector<id_type>& pptr, const std::vector<id_type>& pind) const
{
	// Set up the _node_partition_interface map
	_mesh->_node_partition_interface.clear();
	for(id_type interface=0; interface<(pptr.size()-1); ++interface)
	{
		int start = pptr[interface];
		int end = pptr[interface+1];
		for(int j=start; j<end; ++j)
		{
			int node = pind[j];
			if (_mesh->_node_partition_interface.find(node)==_mesh->_node_partition_interface.end()) // node isn't already in the map (add a new vector)
			{
				std::vector<int> v(1);
				v[0] = interface;
				_mesh->_node_partition_interface.insert(std::pair<id_type, std::vector<int> >(node, v));
			}
			else // Node lies on yet another partition (add another entry to the existing vector)
			{
				_mesh->_node_partition_interface[node].push_back(interface);
			}
		}
	}
	for(auto it=_mesh->_node_partition_interface.begin(), end=_mesh->_node_partition_interface.end(); it!=end; ++it)
		(*it).second.push_back(_mesh->get_rank()); // Stored as a pair in the map data structure
}




// /*
//  * Scatter all of the boundary condition information from proc 0 to the other procs
//  */
// void PartitionerSerial::communicate_boundary_conditions(const std::vector<idx_t>& npart) const
// {
// 	if (_mesh->get_rank() == 0)
// 	{
// 		id_type nparts = _mesh->n_ranks();
// 		std::vector<MPI_Request> reqs(2*(nparts-1));
// 		std::vector<std::vector<id_type> > bcs_Node_and_dof_send_lists(nparts);
// 		std::vector<std::vector<double> > bcs_vals_send_lists(nparts);
// 		for(auto it=_mesh->_boundary->dirichlet_begin(), end=_mesh->_boundary->dirichlet_end(); it!=end; ++it)
// 		{
// 			id_type node_id = (*it).first;
// 			int send_part = npart[node_id];
// 			for(auto it2=(*it).second.begin(), end2=(*it).second.end(); it2!=end2; ++it2)
// 			{
// 				id_type dof = (*it2).first;
// 				bcs_Node_and_dof_send_lists[send_part].push_back(node_id);
// 				bcs_Node_and_dof_send_lists[send_part].push_back(dof);
// 				bcs_vals_send_lists[send_part].push_back((*it2).second);
// 			}	
// 		}
// 		// Send the information about boundary conditions
// 		for(id_type i=1; i<(nparts); ++i) // Send to all processes other than 0
// 		{
// 			MPI_Isend(bcs_Node_and_dof_send_lists[i].data(), bcs_Node_and_dof_send_lists[i].size(), MPI_ID, i, 0, _mesh->get_comm(), &reqs[(i-1)+(0*(nparts-1))]);
// 			MPI_Isend(bcs_vals_send_lists[i].data(), bcs_vals_send_lists[i].size(), MPI_DOUBLE, i, 1, _mesh->get_comm(), &reqs[(i-1)+(1*(nparts-1))]);
// 		}

// 		// Build my partition
// 		build_boundary_conditions(bcs_Node_and_dof_send_lists[0], bcs_vals_send_lists[0]);

// 		// MPI cleanup
// 		MPI_Waitall(2*(nparts-1), reqs.data(), MPI_STATUSES_IGNORE);
// 	}

// 	// Recieve the info
// 	else
// 	{
// 		std::vector<id_type> bc_node_and_dof_recv_vec;	Utilities::RecieveUnknown(bc_node_and_dof_recv_vec, 0, 0, MPI_ID, _mesh->get_comm());
// 		std::vector<double> bc_vals_recv_vec;			Utilities::RecieveUnknown(bc_vals_recv_vec, 0, 1, MPI_DOUBLE, _mesh->get_comm());

// 		// Build my partition
// 		build_boundary_conditions(bc_node_and_dof_recv_vec, bc_vals_recv_vec);
// 	}
// }
// void PartitionerSerial::build_boundary_conditions(const std::vector<id_type>& bc_node_dofs, const std::vector<double>& bc_vals) const
// {
// 	_mesh->_boundary->_dirichlet_bcs.clear();
// 	_mesh->_boundary->attach_mesh(_mesh); // Make sure we're working with the right mesh
// 	// Set up the boundary condition information
// 	for(id_type i=0; i<bc_vals.size(); ++i)
// 	{
// 		double val = bc_vals[i];
// 		id_type node = bc_node_dofs[i*2];
// 		id_type dof = bc_node_dofs[i*2+1];
// 		_mesh->_boundary->set_dirichlet_bc(node, dof, val);
// 	}
// }




/*
 * Scatter all of the nodeset information from proc 0 to the other procs
 */
void PartitionerSerial::communicate_nodesets(const std::vector<std::vector<id_type> >& node_to_parts) const
{
	if (_mesh->get_rank() == 0)
	{
		id_type nparts = _mesh->n_ranks();
		std::vector<MPI_Request> reqs(4*(nparts-1));
		std::vector<id_type> nodeset_name_ptr_send_list(1, 0);
		std::vector<char> nodeset_name_ind_send_list;
		std::vector<std::vector<id_type> > nodeset_id_ptr_send_lists(nparts, std::vector<id_type>(1, 0));
		std::vector<std::vector<id_type> > nodeset_id_ind_send_lists(nparts);
		id_type curr_ind_name = 0;
		std::vector<id_type> curr_inds_ids(nparts, 0);
		// Send every nodeset name to every partition. This is the only way for the user to use the code agnostic of partitioning
		for(auto it1=_mesh->nodesets_begin(), end1=_mesh->nodesets_end(); it1!=end1; ++it1)
		{
			// Send every nodeset name to every partition. This is the only way for the user to use the code agnostic of partitioning
			std::string name = (*it1).first;
			for(id_type j=0; j<name.length(); ++j)
			{
				nodeset_name_ind_send_list.push_back( name[j] );     // Store the individual characters of the name
				curr_ind_name++;
			}
			nodeset_name_ptr_send_list.push_back( curr_ind_name ); // Store one after the ending index of the name

			// Actually put all the node ids that each partition will be recieving into the lists
			for(auto it2=(*it1).second.begin(), end2=(*it1).second.end(); it2!=end2; ++it2) // Iterate through all the nodes of this set
			{
				id_type id = (*it2);
				for(id_type r=0; r<node_to_parts[id].size(); ++r) // Iterate through all of the partitions that have a copy of this node
				{
					int part = node_to_parts[id][r];
					nodeset_id_ind_send_lists[part].push_back( id ); // Actually add the id to the send list
					curr_inds_ids[part]++; // Update the current index for the current partition
				}
			}
			for(id_type part=0; part<nparts; ++part)
				nodeset_id_ptr_send_lists[part].push_back( curr_inds_ids[part] );
		}
		for(id_type i=1; i<nparts; ++i) // Send to all processes other than 0
		{
			MPI_Isend(nodeset_name_ptr_send_list.data(), nodeset_name_ptr_send_list.size(), MPI_ID, i, 0, _mesh->get_comm(), &reqs[(i-1)+(0*(nparts-1))]);
			MPI_Isend(nodeset_name_ind_send_list.data(), nodeset_name_ind_send_list.size(), MPI_CHAR, i, 1, _mesh->get_comm(), &reqs[(i-1)+(1*(nparts-1))]);
			MPI_Isend(nodeset_id_ptr_send_lists[i].data(), nodeset_id_ptr_send_lists[i].size(), MPI_ID, i, 2, _mesh->get_comm(), &reqs[(i-1)+(2*(nparts-1))]);
			MPI_Isend(nodeset_id_ind_send_lists[i].data(), nodeset_id_ind_send_lists[i].size(), MPI_ID, i, 3, _mesh->get_comm(), &reqs[(i-1)+(3*(nparts-1))]);
		}

		// Build my nodesets
		build_nodesets(nodeset_name_ptr_send_list, nodeset_name_ind_send_list,
					   nodeset_id_ptr_send_lists[0], nodeset_id_ind_send_lists[0]);

		// MPI cleanup
		MPI_Waitall(4*(nparts-1), reqs.data(), MPI_STATUSES_IGNORE);
	}

	// Recieve the info
	else
	{
		std::vector<id_type> nodeset_name_ptr_recv_vec;	Utilities::RecieveUnknown(nodeset_name_ptr_recv_vec, 0, 0, MPI_ID, _mesh->get_comm());
		std::vector<char> nodeset_name_ind_recv_vec;	Utilities::RecieveUnknown(nodeset_name_ind_recv_vec, 0, 1, MPI_CHAR, _mesh->get_comm());
		std::vector<id_type> nodeset_ids_ptr_recv_vec;	Utilities::RecieveUnknown(nodeset_ids_ptr_recv_vec, 0, 2, MPI_ID, _mesh->get_comm());
		std::vector<id_type> nodeset_ids_ind_recv_vec;	Utilities::RecieveUnknown(nodeset_ids_ind_recv_vec, 0, 3, MPI_ID, _mesh->get_comm());
		if (nodeset_name_ptr_recv_vec.size() != nodeset_ids_ptr_recv_vec.size())
			err_message("Number of nodeset names recieved did not match number of nodeset sets of ids.");

		// Build my nodesets
		build_nodesets(nodeset_name_ptr_recv_vec, nodeset_name_ind_recv_vec,
					   nodeset_ids_ptr_recv_vec, nodeset_ids_ind_recv_vec);
	}
}
void PartitionerSerial::build_nodesets(const std::vector<id_type>& nodeset_name_ptr, const std::vector<char>& nodeset_name_ind,
									   const std::vector<id_type>& nodeset_ids_ptr, const std::vector<id_type>& nodeset_ids_ind) const
{
	_mesh->_nodesets.clear();
	std::map<std::string, std::set<id_type> > nodesets;
	for(id_type i=0; i<(nodeset_name_ptr.size()-1); ++i)
	{
		id_type name_start = nodeset_name_ptr[i];
		id_type name_end = nodeset_name_ptr[i+1];
		id_type ids_start = nodeset_ids_ptr[i];
		id_type ids_end = nodeset_ids_ptr[i+1];
		std::string name;
		name.assign(&nodeset_name_ind[name_start], name_end-name_start); // Assigns the string starting from the index name_start of length (name_end-name_start) to name
		nodesets.insert(std::pair<std::string, std::set<id_type> >
						(name, std::set<id_type>(nodeset_ids_ind.begin()+ids_start, nodeset_ids_ind.begin()+ids_end)));
		_mesh->add_nodeset_total(name, nodesets[name]); // Note, this call essentiall clears the local nodesets variable every time because it uses a swap
	}
}




/*
 * Scatter all of the element set information from proc 0 to the other procs
 */
void PartitionerSerial::communicate_elemsets(const std::vector<idx_t>& epart) const
{
	if (_mesh->get_rank() == 0)
	{
		id_type nparts = _mesh->n_ranks();
		std::vector<MPI_Request> reqs(4*(nparts-1));
		std::vector<id_type> elemset_name_ptr_send_list(1, 0);
		std::vector<char> elemset_name_ind_send_list;
		std::vector<std::vector<id_type> > elemset_id_ptr_send_lists(nparts, std::vector<id_type>(1, 0));
		std::vector<std::vector<id_type> > elemset_id_ind_send_lists(nparts);
		id_type curr_ind_name = 0;
		std::vector<id_type> curr_inds_ids(nparts, 0);
		// Send every nodeset name to every partition. This is the only way for the user to use the code agnostic of partitioning
		for(auto it1=_mesh->_elemsets.begin(), end1=_mesh->_elemsets.end(); it1!=end1; ++it1)
		{
			// Send every nodeset name to every partition. This is the only way for the user to use the code agnostic of partitioning
			std::string name = (*it1).first;
			for(id_type j=0; j<name.length(); ++j)
			{
				elemset_name_ind_send_list.push_back( name[j] );     // Store the individual characters of the name
				curr_ind_name++;
			}
			elemset_name_ptr_send_list.push_back( curr_ind_name ); // Store one after the ending index of the name

			// Actually put all the node ids that each partition will be recieving into the lists
			for(auto it2=(*it1).second.begin(), end2=(*it1).second.end(); it2!=end2; ++it2) // Iterate through all the nodes of this set
			{
				id_type id = (*it2);
				int part = epart[id];
				elemset_id_ind_send_lists[part].push_back( id ); // Actually add the id to the send list
				curr_inds_ids[part]++; // Update the current index for the current partition
			}
			for(id_type part=0; part<nparts; ++part)
				elemset_id_ptr_send_lists[part].push_back( curr_inds_ids[part] );
		}
		for(id_type i=1; i<nparts; ++i) // Send to all processes other than 0
		{
			MPI_Isend(elemset_name_ptr_send_list.data(), elemset_name_ptr_send_list.size(), MPI_ID, i, 0, _mesh->get_comm(), &reqs[(i-1)+(0*(nparts-1))]);
			MPI_Isend(elemset_name_ind_send_list.data(), elemset_name_ind_send_list.size(), MPI_CHAR, i, 1, _mesh->get_comm(), &reqs[(i-1)+(1*(nparts-1))]);
			MPI_Isend(elemset_id_ptr_send_lists[i].data(), elemset_id_ptr_send_lists[i].size(), MPI_ID, i, 2, _mesh->get_comm(), &reqs[(i-1)+(2*(nparts-1))]);
			MPI_Isend(elemset_id_ind_send_lists[i].data(), elemset_id_ind_send_lists[i].size(), MPI_ID, i, 3, _mesh->get_comm(), &reqs[(i-1)+(3*(nparts-1))]);
		}

		// Build my partition
		build_elemsets(elemset_name_ptr_send_list, elemset_name_ind_send_list,
					   elemset_id_ptr_send_lists[0], elemset_id_ind_send_lists[0]);

		MPI_Waitall(4*(nparts-1), reqs.data(), MPI_STATUSES_IGNORE);
	}

	// Recieve the info
	else
	{
		std::vector<id_type> elemset_name_ptr_recv_vec;	Utilities::RecieveUnknown(elemset_name_ptr_recv_vec, 0, 0, MPI_ID, _mesh->get_comm());
		std::vector<char> elemset_name_ind_recv_vec;	Utilities::RecieveUnknown(elemset_name_ind_recv_vec, 0, 1, MPI_CHAR, _mesh->get_comm());
		std::vector<id_type> elemset_ids_ptr_recv_vec;	Utilities::RecieveUnknown(elemset_ids_ptr_recv_vec, 0, 2, MPI_ID, _mesh->get_comm());
		std::vector<id_type> elemset_ids_ind_recv_vec;	Utilities::RecieveUnknown(elemset_ids_ind_recv_vec, 0, 3, MPI_ID, _mesh->get_comm());
		if (elemset_name_ptr_recv_vec.size() != elemset_ids_ptr_recv_vec.size())
			err_message("Number of elemset names recieved did not match number of elemset sets of ids.");

		// Build my partition
		build_elemsets(elemset_name_ptr_recv_vec, elemset_name_ind_recv_vec,
					   elemset_ids_ptr_recv_vec, elemset_ids_ind_recv_vec);
	}
}
void PartitionerSerial::build_elemsets(const std::vector<id_type>& elemset_name_ptr, const std::vector<char>& elemset_name_ind,
									   const std::vector<id_type>& elemset_ids_ptr, const std::vector<id_type>& elemset_ids_ind) const
{
	_mesh->_elemsets.clear();
	// Set up my element sets
	for(id_type i=0; i<(elemset_name_ptr.size()-1); ++i)
	{
		id_type name_start = elemset_name_ptr[i];
		id_type name_end = elemset_name_ptr[i+1];
		id_type ids_start = elemset_ids_ptr[i];
		id_type ids_end = elemset_ids_ptr[i+1];
		std::string name;
		name.assign(&elemset_name_ind[name_start], name_end-name_start); // Assigns the string starting from the index name_start of length (name_end-name_start) to name
		_mesh->_elemsets.insert(std::pair<std::string, std::set<id_type> >(name, std::set<id_type>()));
		for(id_type id=ids_start; id<ids_end; ++id)
			_mesh->_elemsets[name].insert(elemset_ids_ind[id]);
	}
}




/*
 * Scatter all of the side set information from proc 0 to the other procs
 */
void PartitionerSerial::communicate_sidesets(const std::vector<idx_t>& epart) const
{
	if (_mesh->get_rank() == 0)
	{
		id_type nparts = _mesh->n_ranks();
		std::vector<MPI_Request> reqs(4*(nparts-1));
		std::vector<id_type> sideset_name_ptr_send_list(1, 0);
		std::vector<char> sideset_name_ind_send_list;
		std::vector<std::vector<id_type> > sideset_id_ptr_send_lists(nparts, std::vector<id_type>(1, 0));
		std::vector<std::vector<id_type> > sideset_id_ind_send_lists(nparts);
		id_type curr_ind_name = 0;
		std::vector<id_type> curr_inds_ids(nparts, 0);
		// Send every nodeset name to every partition. This is the only way for the user to use the code agnostic of partitioning
		for(auto it=_mesh->sidesets_begin(), end=_mesh->sidesets_end(); it!=end; ++it)
		{
			// Send every nodeset name to every partition. This is the only way for the user to use the code agnostic of partitioning
			std::string name = (*it).first;
			for(id_type j=0; j<name.length(); ++j)
			{
				sideset_name_ind_send_list.push_back( name[j] );     // Store the individual characters of the name
				curr_ind_name++;
			}
			sideset_name_ptr_send_list.push_back( curr_ind_name ); // Store one after the ending index of the name

			// Actually put all the node ids that each partition will be recieving into the lists
			for(auto it2=(*it).second.begin(), end2=(*it).second.end(); it2!=end2; ++it2) // Iterate through all the nodes of this set
			{
				id_type id = (*it2).first;
				id_type side = (*it2).second;
				int part = epart[id];
				sideset_id_ind_send_lists[part].push_back( id ); // Actually add the id to the send list
				sideset_id_ind_send_lists[part].push_back( side ); // Actually add the side to the send list
				curr_inds_ids[part] += 2; // Update the current index for the current partition
			}
			for(id_type part=0; part<nparts; ++part)
				sideset_id_ptr_send_lists[part].push_back( curr_inds_ids[part] );
		}
		for(id_type i=1; i<nparts; ++i) // Send to all processes other than 0
		{
			MPI_Isend(sideset_name_ptr_send_list.data(), sideset_name_ptr_send_list.size(), MPI_ID, i, 0, _mesh->get_comm(), &reqs[(i-1)+(0*(nparts-1))]);
			MPI_Isend(sideset_name_ind_send_list.data(), sideset_name_ind_send_list.size(), MPI_CHAR, i, 1, _mesh->get_comm(), &reqs[(i-1)+(1*(nparts-1))]);
			MPI_Isend(sideset_id_ptr_send_lists[i].data(), sideset_id_ptr_send_lists[i].size(), MPI_ID, i, 2, _mesh->get_comm(), &reqs[(i-1)+(2*(nparts-1))]);
			MPI_Isend(sideset_id_ind_send_lists[i].data(), sideset_id_ind_send_lists[i].size(), MPI_ID, i, 3, _mesh->get_comm(), &reqs[(i-1)+(3*(nparts-1))]);
		}

		// Build my sidesets
		build_sidesets(sideset_name_ptr_send_list, sideset_name_ind_send_list,
					   sideset_id_ptr_send_lists[0], sideset_id_ind_send_lists[0]);

		MPI_Waitall(4*(nparts-1), reqs.data(), MPI_STATUSES_IGNORE);
	}

	// Recieve the info
	else
	{
		std::vector<id_type> sideset_name_ptr_recv_vec;	Utilities::RecieveUnknown(sideset_name_ptr_recv_vec, 0, 0, MPI_ID, _mesh->get_comm());
		std::vector<char> sideset_name_ind_recv_vec;	Utilities::RecieveUnknown(sideset_name_ind_recv_vec, 0, 1, MPI_CHAR, _mesh->get_comm());
		std::vector<id_type> sideset_ids_ptr_recv_vec;	Utilities::RecieveUnknown(sideset_ids_ptr_recv_vec, 0, 2, MPI_ID, _mesh->get_comm());
		std::vector<id_type> sideset_ids_ind_recv_vec;	Utilities::RecieveUnknown(sideset_ids_ind_recv_vec, 0, 3, MPI_ID, _mesh->get_comm());
		if (sideset_name_ptr_recv_vec.size() != sideset_ids_ptr_recv_vec.size())
			err_message("Number of sideset names recieved did not match number of sideset sets of ids.");

		// Build my sidesets
		build_sidesets(sideset_name_ptr_recv_vec, sideset_name_ind_recv_vec,
					   sideset_ids_ptr_recv_vec, sideset_ids_ind_recv_vec);

	}
}
void PartitionerSerial::build_sidesets(const std::vector<id_type>& sideset_name_ptr, const std::vector<char>& sideset_name_ind,
									   const std::vector<id_type>& sideset_ids_ptr, const std::vector<id_type>& sideset_ids_ind) const
{
	_mesh->_sidesets.clear();


	std::map<std::string, std::set<std::pair<id_type, id_type> > > sidesets;
	for(id_type i=0; i<(sideset_name_ptr.size()-1); ++i)
	{
		id_type name_start = sideset_name_ptr[i];
		id_type name_end = sideset_name_ptr[i+1];
		id_type ids_start = sideset_ids_ptr[i];
		id_type ids_end = sideset_ids_ptr[i+1];
		std::string name;
		name.assign(&sideset_name_ind[name_start], name_end-name_start); // Assigns the string starting from the index name_start of length (name_end-name_start) to name
		sidesets.insert(std::pair<std::string, std::set<std::pair<id_type, id_type> > >(name, std::set<std::pair<id_type, id_type> >()));
		for(id_type id=ids_start; id<ids_end; id+=2)
			sidesets[name].insert(std::pair<id_type, id_type>(sideset_ids_ind[id], sideset_ids_ind[id+1])); // Not as easy to just use a different constructor here :(
		_mesh->add_sideset_total(name, sidesets[name]);
	}
}




/*
 * Scatter all of the material information from proc 0 to the other procs
 *	NOTE: THIS IS BROKEN THIS SHOULD BE COMMENTED OUT IN THE PARTITION FUNCTION
 */
void PartitionerSerial::communicate_materials(const std::vector<idx_t> epart, const std::vector<id_type> n_elem_per_proc) const
{
	if (_mesh->get_rank() == 0)
	{
		id_type nparts = _mesh->n_ranks();
		std::vector<MPI_Request> reqs(2*(nparts-1));
		std::vector<char> material_send_list;
		std::vector<std::vector<int> > elem_material_send_lists(nparts);
		std::vector<id_type> n_elem_done(nparts);
		// Send every material to every partition. This is the only way for the user to use the code agnostic of partitioning
		for(id_type m=0; m<_mesh->_materials.size(); ++m)
		{
			std::vector<char> vec(MAX_MATERIAL_SIZE);
			char* buf = &vec[0];
			_mesh->_materials[m]->pack(buf);   // Packs the material into a contiguous char buffer
			material_send_list.insert(material_send_list.end(), vec.begin(), vec.end());  // Append the current material vector to the end of the current partition's material send lists
		}
		// Preallocate memory for the elemental material ids for each partition
		for(id_type part=0; part<(nparts); ++part)
			elem_material_send_lists[part].resize(n_elem_per_proc[part]);
		// Loop over all elements and add the elemental material id to the correct partition send list
		for(id_type i=0; i<epart.size(); ++i)
		{
			int part = epart[i];
			elem_material_send_lists[part][n_elem_done[part]] = _mesh->_elem_to_material[i];
			n_elem_done[part]++;
		}
		for(id_type i=1; i<(nparts); ++i) // Send to all processes other than 0
		{
			MPI_Isend(material_send_list.data(), material_send_list.size(), MPI_CHAR, i, 0, _mesh->get_comm(), &reqs[(i-1)+(0*(nparts-1))]);
			MPI_Isend(elem_material_send_lists[i].data(), elem_material_send_lists[i].size(), MPI_INT, i, 1, _mesh->get_comm(), &reqs[(i-1)+(1*(nparts-1))]);
		}

		// Build my materials
		build_materials(material_send_list, elem_material_send_lists[0]);

		MPI_Waitall(2*(nparts-1), reqs.data(), MPI_STATUSES_IGNORE);
	}

	// Recieve the info
	else
	{
		std::vector<char> material_recv_vec;		Utilities::RecieveUnknown(material_recv_vec, 0, 0, MPI_CHAR, _mesh->get_comm());
		std::vector<int> elem_material_recv_vec;	Utilities::RecieveUnknown(elem_material_recv_vec, 0, 1, MPI_INT, _mesh->get_comm());

		// Build my materials
		build_materials(material_recv_vec, elem_material_recv_vec);
	}
}
void PartitionerSerial::build_materials(std::vector<char>& mat_vec, std::vector<int>& elem_mat_vec) const
{
	char* buf = &mat_vec[0];
	size_t offset;

	// Add any actual materials
	int nmat = mat_vec.size() / MAX_MATERIAL_SIZE;
	for(int i=0; i<nmat; ++i)
	{
		offset = MAX_MATERIAL_SIZE*i;
		Material* mat = Material::unpack(buf+offset);
		_mesh->add_material( mat ); // Copies the material into a new object
		delete mat;
	}

	// Set up the elem-to-material vector
	_mesh->_elem_to_material.swap(elem_mat_vec);

	// Assigns the materials assigned boolean variable
	bool all_set = true;
	for(id_type i=0; i<_mesh->_elem_to_material.size(); ++i)
	{
		if (_mesh->_elem_to_material[i]<0)
		{
			all_set = false;
			break;
		}
	}
	if (all_set)
		_mesh->_materials_assigned = true;
}




/*
 * Perform reductions on all of the mesh partitions to determine global information
 */
void PartitionerSerial::determine_global_info() const
{
	// Might as well do the global reduction on number of nodes and elements here
	MPI_Allreduce(&_mesh->_n_local_owned_nodes, &_mesh->_n_global_nodes, 1, MPI_ID, MPI_SUM, _mesh->get_comm());
	id_type n_elem = _mesh->n_local_elem();
	MPI_Allreduce(&n_elem, &_mesh->_n_global_elem, 1, MPI_ID, MPI_SUM, _mesh->get_comm());

	// Create a set of all of the processors that I share part of a partition interface withc
	for (Mesh::partition_iterator it=_mesh->partition_interface_begin(), end=_mesh->partition_interface_end(); it!=end; ++it)
		_mesh->_proc_neighbors.insert(it->second.begin(), it->second.end());

	// Set boolean
	_mesh->_partitioned = true;

	// If the mesh is periodic, determine some extra info about it
	if (_mesh->periodic())
	{
		for (id_type n=0; n<_mesh->n_local_nodes(); ++n)
		{
			id_type node_id = _mesh->get_node_local(n)->get_id();
			if (!_mesh->node_periodic(node_id))
			{
				_mesh->_n_local_non_periodic_nodes++;
				if (_mesh->own_node_local(n))
					_mesh->_n_local_owned_non_periodic_nodes++;
			}
			else
			{
				if (_mesh->primary_periodic_node(node_id)) // Only count primary periodic nodes
				{
					_mesh->_n_local_non_periodic_nodes++;
					if (_mesh->own_node_local(n))
						_mesh->_n_local_owned_non_periodic_nodes++;
				}
			}
		}
		MPI_Allreduce(&_mesh->_n_local_owned_non_periodic_nodes, &_mesh->_n_global_non_periodic_nodes, 1, MPI_ID, MPI_SUM, _mesh->get_comm());
	}
}