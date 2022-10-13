/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated September 2016

##################################################################################
*/
#ifndef _BCS_H_
#define _BCS_H_
#include "common.h"
#include <set>
#include <map>
#include <string>


// Forward declaration of mesh class
class Problem;

class BoundaryObject
{
	protected:

		// The mesh that this BoundaryObject is associated with
		Problem* _prob;

		typedef std::map<id_type, std::map<id_type, double> > dirichlet_t;
		
		
		
		
		// INFORMATION ABOUT BOUNDARY CONDITIONS
		//----------------------------------------------------------------------------------------------------------------------------------------------
		/*
		 * Structure to hold the dirichlet bcs associated with the local mesh
		 * _dirichlet_bcs[i][j] corresponds to to the node with global id "i" and the j dof on that node (0-2 for 3D elasticity)
		*/
		dirichlet_t _dirichlet_bcs;
		dirichlet_t _individuals;

		// Structures containing all of the BCs added before applcation
		std::map<std::string, std::vector<std::pair<id_type, double> > > _dirichlet_nodeset_bcs;

		std::set<id_type> _lock_down_directions;

		// Define friend classes
		friend class Mesh;
		friend class Partitioner;
		friend class PartitionerSerial;


	public:

		// Constructor
		BoundaryObject();

		// Used to add a mesh object to a BoundaryObject that has been cleared
		void attachProblem(Problem* prob) {_prob = prob;};

		// Clears the boundary object
		void clear();

		






		// Function to apply BCs that were assciated with a nodeset or a sideset
		void apply_BCs();

		/*
		* Functions to add a prescribed dof to the mesh
		* 1. Add a single dirichlet bc so a single node and single dof
		* 2. Add a vector of dirichlet bcs to a vector of nodes and vector of dofs
		* 3. Add a single dirichlet bc to a vector of nodes on a single dof
		* 4. Add a vector of dirichlet bcs to a single nodeset (number) and a vector of dofs
		* 5. Add a vector of dirichlet bcs to a single nodeset (name) and a vector of dofs
		* 6. Add a single dirichlet bc to a single nodeset (number) and a single dof
		* 7. Add a single dirichlet bc to a single nodeset (name) and a single dof
		*/
		void set_dirichlet_bc(id_type node, id_type dof, double val);
		void set_dirichlet_bcs(std::vector<id_type>& nodes, std::vector<id_type>& dofs, std::vector<double>& vals);
		void set_dirichlet_bcs(std::vector<id_type>& nodes, id_type dof, double val);
		//void set_dirichlet_bcs_from_nodeset(id_type set, std::vector<id_type>& dofs, std::vector<double>& vals);
		//void set_dirichlet_bcs_from_nodeset(id_type set, id_type dof, double val);
		//void set_dirichlet_bcs_from_nodeset(std::string set_name, std::vector<id_type>& dofs, std::vector<double>& vals);
		void set_dirichlet_bcs_from_nodeset(std::string set_name, id_type dof, double val);

		void lock_down_direction(id_type dir) {_lock_down_directions.insert(dir);};
		
		// Function to determine if a given node and dof has a boundary condition
		bool has_dirichlet_bc(id_type node, id_type dof);

		id_type n_dirichlet_on_node(id_type node);

		// Function to get the prescribed value 
		double get_dirichlet_bc(id_type node, id_type dof, bool& exists);

		// Boundary condition iterator
		typedef dirichlet_t::iterator dirichlet_iterator;
		dirichlet_iterator dirichlet_begin() {return _dirichlet_bcs.begin();};
		dirichlet_iterator dirichlet_end() {return _dirichlet_bcs.end();};
};


#endif
