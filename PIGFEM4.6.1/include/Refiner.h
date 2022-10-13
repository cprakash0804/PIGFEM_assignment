/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _REFINER_H_
#define _REFINER_H_

#include "common.h"
#include <set>
#include <vector>
#include <unordered_map>
#include <map>

// Predeclarations
class Problem;
class Mesh;


/*
 * This is an abstract base class to refine a mesh based upon
 * certain mesh refinement criteria (geometrical, damage,
 * anything problem specific). This class implements the ability
 * to refine a general set of active elements and interpolate the
 * the current solution. To instantiate a refiner a criteria function
 * must be defined.
 */
class Refiner
{
	public:


		/*
		 * The Big 3
		 * The constructor takes the problem that will be partitioned
		 * Destructor doesn't do anything
		 * Copy constructor copies over the problem pointer
		 */
		Refiner(Problem* prob);
		Refiner();
		~Refiner();
		Refiner(const Refiner& other);


		/*
		 * The main interface with this class. Will be called to refine the mesh.
		 * Can be overridden for specific implementation in derived classes
		 */
		virtual void refine();


		/*
		 * Returns the Problem that is being refined
		 */
		Problem* get_prob() const {return _prob;};


		/*
		 * Returns the Mesh that is being refined
		 */
		Mesh* get_mesh() const;


	protected:


		/*
		 * The problem that contains the mesh that is to be refined
		 */
		Problem* _prob;


		/*
		 * This is the function that must be implemented to generate a set of elements to refine
		 */
		virtual void generate_refine_set(std::set<id_type>& elements) = 0;


		/*
		 * If the mesh hasn't beeen refined before, set some things
		 */
		void initialize_mesh_refinement();


		/*
		 * Given a set of new nodes, returns the owner of the refinemenet node generated from them
		 */
		int get_new_refine_node_owner(const std::vector<id_type>& node_group);


		/*
		 * Refine a general set of active elements (I think these need to be global ids)
		 */
		void refine_elements(const std::set<id_type>& refined_elements);


		/*
		 * Will be called after the mesh is refined to interpolate the current solution field
		 * Handles solution field as well as quadrature point internal variables
		 */
		void interpolate_solution();


		/*
		 * Enforces the restriction that neighboring elements cannot be more than 1 p-level different
		 */
		void enforce_p_level_constraint(std::set<id_type>& refined_elements);


		/*
		* Creates a set of node groups. Each node group represents a new node that will be added to the mesh.
		*/
		void generate_new_node_groups(const std::set<id_type>& refined_elements,
									  std::set<std::vector<id_type> >& new_node_groups);


		/*
		 * If this is being run in parallel then communicate the node groups
		 * that exist entirely along partition interfaces. Assign owners to these node groups
		 */
		void communicate_new_node_groups(std::set<std::vector<id_type> >& new_node_groups);


		/*
		 * Create all of the new nodes from the node groups
		 */
		void create_new_nodes(const std::set<std::vector<id_type> >& new_node_groups);


		/*
		 * Actually find the coordinates of the new nodes and create the node objects
		 */
		void allocate_new_nodes(const std::set<std::vector<id_type> >& new_node_groups);


		/*
		 * Set the node ids for the nodes that I own
		 */
		void set_new_node_ids(const std::set<std::vector<id_type> >& new_node_groups);


		/*
		 * In parallel, communicate the global node ids of new nodes appearing on partition interfaces
		 * Maintains consistent node numbering across mesh partitions
		 */
		void communicate_new_node_ids(const std::set<std::vector<id_type> >& new_node_groups);


		/*
		 * Add the newly created nodes to the appropriate node sets
		 */
		void add_new_nodes_to_nodesets(const std::set<std::vector<id_type> >& new_node_groups);


		/*
		 * Detect all of the new nodes (in an IGFEM mesh)
		 */
		void detect_new_nodes(const std::set<std::vector<id_type> >& new_node_groups);


		/*
		 * Actually create the new elements from the given set of elements to be refined
		 */
		void create_new_elements(const std::set<id_type>& refined_elements);


		/*
		 * Gather any global information regarding the global mesh
		 */
		void determine_global_information();


		/*
		 * Adds to the set refined_elements all those elements that break the p-level restriction
		 * starting from the elements in the new_elements set. Also adds lists of new remote elements
		 * that should be refined in parallel
		 */
		void find_local_p_level_additions(std::set<id_type>& refined_elements, std::set<id_type>& new_elements,
										  std::map<int, std::vector<id_type> >& remote_additions);


		/*
		 * Communicates to all neighboring mesh partitions those elements which should be refined
		 * based on criteria present on remote partitions. Stores all of these new elements in the
		 * new_local_elements set
		 */
		void communicate_new_remote_elements(const std::set<id_type>& refined_elements, std::set<id_type>& new_local_elements,
											 std::map<int, std::vector<id_type> >& remote_additions);
};


#endif