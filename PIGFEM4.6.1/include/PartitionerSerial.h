/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "Partitioner.h"
#include "metis.h"
#include <vector>
#include <set>



/*
 * Derived class that takes a mesh that exists on
 * one processor and partitiones it to n processers
 * based on the METIS decomposition
 */
class PartitionerSerial : public Partitioner
{
	public:


		/*
		 * The main function that is called to partition the mesh
		 */
		virtual void partition();


		/*
		 * The Big 3
		 * The constructor takes the mesh that will be partitioned
		 * Destructor doesn't do anything
		 * Copy constructor copies over the mesh pointer
		 */
		PartitionerSerial(Mesh* mesh);
		~PartitionerSerial();
		PartitionerSerial(const PartitionerSerial& other);


	private:


		/*
		 * Fills in any mesh member variables that need to be
		 * managed in the case of a serial mesh
		 */
		void handle_serial_case() const;


		/*
		 * Uses METIS to get the decomposition of the mesh
		 */
		void get_decomposition(std::vector<idx_t>& epart, std::vector<idx_t>& npart) const;


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
		void generate_node_lists(const std::vector<idx_t>& epart,
								 std::vector<id_type>& n_elem_per_proc, std::vector<id_type>& n_nodes_per_proc,
								 std::vector<id_type>& eind_size, std::vector<std::vector<id_type> >& node_to_parts,
								 std::vector<std::set<id_type> >& part_node_id_lists,
								 std::vector<std::set<id_type> >& periodic_sets) const;


		/*
		 * Scatter all of the nodal information from proc 0 to the other procs
		 */
		void communicate_nodes(const std::vector<idx_t>& npart, const std::vector<id_type>& n_nodes_per_proc,
							   const std::vector<std::set<id_type> >& part_node_id_lists) const;
		void build_nodes(const std::vector<id_type>& ids, const std::vector<double>& coords, const std::vector<int>& owners) const;


		/*
		 * Scatter all of the elemental structure information from proc 0 to the other procs
		 */
		void communicate_elements(const std::vector<idx_t>& epart, const std::vector<id_type>& n_elem_per_proc,
								  const std::vector<id_type>& eind_size) const;
		void build_elements(const std::vector<id_type>& ids, const std::vector<elem_type>& types,
							const std::vector<id_type>& eptr, const std::vector<id_type>& eind) const;


		/*
		 * Communicates the periodic nodesets of a periodic mesh
		 */
		void commmunicate_periodic_nodesets(const std::vector<std::set<id_type> >& part_node_id_lists,
											const std::vector<std::set<id_type> >& periodic_sets);
		void build_periodicity(const std::vector<id_type> periodic_ptr, const std::vector<id_type> periodic_ind);


		/*
		 * Scatter all of the information about the interface between
		 * mesh partitions from proc 0 to the other procs
		 */
		void communicate_partition_interface(const std::vector<std::set<id_type> >& part_node_id_lists) const;
		void build_partition_interface(const std::vector<id_type>& pptr, const std::vector<id_type>& pind) const;


		// /*
		//  * Scatter all of the boundary condition information from proc 0 to the other procs
		//  */
		// void communicate_boundary_conditions(const std::vector<idx_t>& npart) const;
		// void build_boundary_conditions(const std::vector<id_type>& bc_node_dofs, const std::vector<double>& bc_vals) const;


		/*
		 * Scatter all of the nodeset information from proc 0 to the other procs
		 */
		void communicate_nodesets(const std::vector<std::vector<id_type> >& node_to_parts) const;
		void build_nodesets(const std::vector<id_type>& nodeset_name_ptr, const std::vector<char>& nodeset_name_ind,
							const std::vector<id_type>& nodeset_ids_ptr, const std::vector<id_type>& nodeset_ids_ind) const;


		/*
		 * Scatter all of the element set information from proc 0 to the other procs
		 */
		void communicate_elemsets(const std::vector<idx_t>& epart) const;
		void build_elemsets(const std::vector<id_type>& elemset_name_ptr, const std::vector<char>& elemset_name_ind,
							const std::vector<id_type>& elemset_ids_ptr, const std::vector<id_type>& elemset_ids_ind) const;


		/*
		 * Scatter all of the side set information from proc 0 to the other procs
		 */
		void communicate_sidesets(const std::vector<idx_t>& epart) const;
		void build_sidesets(const std::vector<id_type>& sideset_name_ptr, const std::vector<char>& sideset_name_ind,
							const std::vector<id_type>& sideset_ids_ptr, const std::vector<id_type>& sideset_ids_ind) const;


		/*
		 * Scatter all of the material information from proc 0 to the other procs
		 *	NOTE: THIS IS BROKEN THIS SHOULD BE COMMENTED OUT IN THE PARTITION FUNCTION
		 */
		void communicate_materials(const std::vector<idx_t> epart, const std::vector<id_type> n_elem_per_proc) const;
		void build_materials(std::vector<char>& mat_vec, std::vector<int>& elem_mat_vec) const;


		/*
		 * Perform reductions on all of the mesh partitions to determine global information
		 */
		void determine_global_info() const;
};