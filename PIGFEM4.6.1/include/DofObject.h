/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _DOFS_H_
#define _DOFS_H_
#include "common.h"
#include <map>
#include <vector>
#include <set>
#include <unordered_map>

// Forward Declaraions
class Mesh;
class Problem;
class Elem;
class CohesiveElem;

class DofObject
{
	protected:
		// The mesh object associated with this DOfObject
		Problem* _prob;
		Mesh* _mesh;

		// INFO ABOUT DEGREES OF FREEDOM, BOTH GLOBAL AND LOCAL
		//----------------------------------------------------------------------------------------------------------------------------------------------

		// Stores number of nodal dofs
		id_type _nndof;

		/*
		 * Data structure to store the global dofs associated with each local node
		 * _node_global_dofs[i][j] correspnds to the jth global dof # of the ith local node
		*/
		std::vector<std::vector<id_type> > _node_global_dofs;
		
		/*
		 * Data structure to store the local dofs associated with each local node
		 * _node_local_dofs[i][j] correspnds to the jth local dof # of the ith local node
		*/
		std::vector<std::vector<id_type> > _node_local_dofs;
		
		/*
		 * Structure to map from the global dofs contained on this mesh to the local dof numbers
		*/
		std::map<id_type, id_type> _global_local_dof;
		
		/*
		 * Structure to map from the local dof numbers to the global dof numbers
		 * Only the first n_local_owned_dofs() actually belong to the process. HIgher local dofs are actually owned by other processes
		*/
		std::vector<id_type> _local_global_dof;
		
		/*
		 * Map that goes from global dof to a pair representing (global node#, dof# on that node)
		 * At the moment this only contains global dofs that I own. We'l see if I need to change this
		 * FIXME: change this. and also make his an unordered_map
		*/
		std::unordered_map<id_type, std::pair<id_type, id_type> > _global_dof_to_node;

		id_type _n_local_owned_free_dofs;
		id_type _n_local_owned_const_dofs;
		id_type _n_local_free_dofs;
		id_type _n_local_const_dofs;
		id_type _n_global_free_dofs;
		id_type _n_global_const_dofs;

		// Boolean that actually stores whether or not the dofs have been distributed yet
		bool _dof_distributed;

		// IGFEM information
		//------------------------------------------------------------------------------------------------------------------------
		/*
		 * Data structure to store the global dofs associated with each local node
		 * _node_global_dofs[i][j] correspnds to the jth global dof # of the ith local node
		*/
		std::vector<std::vector<id_type> > _enrich_node_global_dofs;
		
		/*
		 * Data structure to store the local dofs associated with each local node
		 * _node_local_dofs[i][j] correspnds to the jth local dof # of the ith local node
		*/
		std::vector<std::vector<id_type> > _enrich_node_local_dofs;

	public:

		DofObject();
		DofObject(const DofObject& other);
		virtual ~DofObject() {};
		void clear();
		void copy(const DofObject& other);

		void attachProblem(Problem* prob);

		id_type nndof() {return _nndof;};

		id_type n_local_owned_dofs() {return _n_local_owned_free_dofs+_n_local_owned_const_dofs;};
		id_type n_local_owned_free_dofs() {return _n_local_owned_free_dofs;};
		id_type n_local_owned_const_dofs() {return _n_local_owned_const_dofs;};
		id_type n_local_dofs() {return _n_local_free_dofs+_n_local_const_dofs;};
		id_type n_local_free_dofs() {return _n_local_free_dofs;};
		id_type n_local_const_dofs() {return _n_local_const_dofs;};
		id_type n_global_dofs() {return _n_global_free_dofs+_n_global_const_dofs;};
		id_type n_global_free_dofs() {return _n_global_free_dofs;};
		id_type n_global_const_dofs() {return _n_global_const_dofs;};

		// Returns a pair containing the globa node number of local dof# corresponding to the global dof passed in
		std::pair<id_type, id_type> get_node_from_dof(id_type dof);

		// Attempts to find the global dof number associated with the node with the given global node id and the dof'th dof on the node
		virtual id_type get_global_dof(id_type node, id_type dof);
		
		// Attempts to find the local dof number associated with the node with the given global node id and the dof'th dof on the node
		virtual id_type get_local_dof(id_type node, id_type dof);
		
		// Returns the global dofs associated with the local node idx
		virtual std::vector<id_type>& get_nodal_global_dofs_local(id_type idx);
		
		// Returns the global dofs associated with the global node id
		virtual std::vector<id_type>& get_nodal_global_dofs_global(id_type id);
		
		// Returns the local dofs associated with the local node idx
		virtual std::vector<id_type>& get_nodal_local_dofs_local(id_type idx);
		
		// Returns the local dofs associated with the global node id
		virtual std::vector<id_type>& get_nodal_local_dofs_global(id_type id);
		
		// Fills the passed in vector with all of the global dofs associated with all of the nodes of the local element idx
		virtual void get_elem_global_dofs(Elem* el, std::vector<id_type> & vec);
		virtual void get_elem_global_dofs(CohesiveElem* el, std::vector<id_type> & vec);

		// Distribute the dofs associated with each node accross all partitions
		virtual void distribute_dofs(id_type nndof);

	private:

		void determine_local_info(id_type& partial_sum_free, id_type& partial_sum_const);

		void set_local_owned_normal_dofs(const id_type& partial_sum_free, const id_type& partial_sum_const,
										 id_type& local_dof_free, id_type& local_dof_const, id_type& local_dof,
										 std::set<id_type>& primary_pnodes);

		void communicate_normal_dofs(id_type& local_dof, std::set<id_type>& primary_pnodes);

		void set_periodic_nodes(const std::set<id_type>& primary_pnodes);

		void set_local_owned_enrich_dofs(const id_type& partial_sum_free, const id_type& partial_sum_const,
										 id_type& local_dof_free, id_type& local_dof_const, id_type& local_dof,
										 std::set<id_type>& primary_pnodes);

		void communicate_enrich_dofs(id_type& local_dof);

		void Finalize();
};




#endif
