/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _MESH_H_
#define _MESH_H_

//#include <cstdlib>   // Do I need this?
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <functional> // for std::hash
#include <string>
#include "common.h"
#include "mpi.h" // kinda an important one here...
#include "petsc.h"
// #include "node.h"
// #include "elem.h"
// #include "material.h"
// #include "Inclusion.h"
// #include "Utilities.h"
// #include "mpi.h" // kinda an important one here...

// Forward Declarations
template<class T>
class DenseMatrix;
class BoundaryObject;
class DofObject;
class Elem;
class Node;
class Inclusion;
class Material;


namespace std {

	template <class T> struct hash<std::vector<T> >
 	{
		std::size_t operator()(const std::vector<T>& in) const
		{
			// Compute individual hash values for first,
			// second and third and combine them using XOR
			// and bit shifting:

			std::hash<id_type> hasher;
			size_t seed = hasher(in[0]);
			for (size_t i=1; i<in.size(); ++i)
				seed ^= hasher(in[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);

			return seed;
		}
	};

	template <class T> struct hash<std::pair<T, T> >
 	{
		std::size_t operator()(const std::pair<T, T>& in) const
		{
			// Compute individual hash values for first,
			// second and third and combine them using XOR
			// and bit shifting:

			std::hash<id_type> hasher;
			size_t seed = hasher(in.first);
			seed ^= hasher(in.second);

			return seed;
		}
	};
}



class Mesh /* inherit from charm object here??? */
{
	protected:
	
		// INFORMATION ABOUT THE MESH'S NODES, ELEMENTS, AND THEIR CONNECTIVITY
		//----------------------------------------------------------------------------------------------------------------------------------------------
		/*
		 * Vector of pointers to the local elements
		*/
		std::vector<Elem*> _elem;

		/*
		 * Vector of pointers to the local nodes
		 *	NOTE: some of these nodes may be shared with other chares/processors but a local pointer is stored here
		*/
		std::vector<Node*> _nodes;

		/*
		 * Vector of partition owners of the local nodes
		 *	NOTE: nodes along boundaries might not be owned by this proc but a copy is stored here
		*/
		std::vector<int> _node_owners;
		id_type _n_local_owned_nodes;

		/*
		 * Map from the global to the local node numbers
		*/
		std::unordered_map<id_type, id_type> _global_to_local_node;
		
		/*
		 * Map from the global to the local element numbers
		*/
		std::unordered_map<id_type, id_type> _global_to_local_elem;

		/*
		 * Data structure containing map from local element numbers to its local node numbers
		*/
		std::vector< std::vector<id_type> > _elem_node;

		/*
		 * Data structure containing map from local node numbers to its GLOBAL element numbers
		*/
		std::vector< std::vector<id_type> > _node_elem;	
		
		/*
		 * Data structure containing map from global node number to a map from partition number to its GLOBAL element numbers
		 * EXAMPLE: TO access the vector of element number for global node 912 on partition 3, you would access _remote_node_elem[912][3].
		 *			To access the global element number of the 4th element in that list you would access _remote_node_elem[912][3][3] (because of zero-index in vector)
		*/
		std::unordered_map< id_type, std::map<int, std::vector<id_type> > > _remote_node_elem;

		id_type _n_global_nodes;
		id_type _n_global_elem;
		
		
		
		
		// INFORMATION ABOUT INTERFACES BETWEEN PARTITIONS OF THE MESH AND THE NODES THAT LIE ON THAT PARTITION
		//----------------------------------------------------------------------------------------------------------------------------------------------
		/*
		 * Map from the global node numbers to vector of partition interfaces they lie on
		 * If a node isn't in the map, that means it doesn't lie on a partition interface
		*/
		std::unordered_map<id_type, std::vector<int> > _node_partition_interface;

		// Set of all processes that neighbor this one
		std::set<int> _proc_neighbors;
		
		
		
		
	
		// INFORMATION ABOUT GROUPS OF NODES AND ELEMENTS
		//----------------------------------------------------------------------------------------------------------------------------------------------		
		/*
		 * Data structure containing global element ids.
		 * Each elemset will contain a list of element numbers associated with it
		 * Useful for applying materials
		*/
		//std::vector<NamedSet> _elemsets;
		std::map<std::string, std::set<id_type> > _elemsets;
		
		
		
		
		
		// INFORMATION ABOUT MATERIALS
		//----------------------------------------------------------------------------------------------------------------------------------------------
		/*
		 * Data structure to contain all the materials in the mesh
		*/
		std::vector<Material*> _materials;

		// Data structure that maps from the material's name to the material's index
		std::map<std::string, int> _mat_name_to_idx;
		
		/*
		 * Data structure to map from the local element numbers ot the indicies of the material in the material vector
		*/
		std::vector<int> _elem_to_material;




		// INFORMATION ABOUT PERIODICITY
		//----------------------------------------------------------------------------------------------------------------------------------------------
		/*
		 * Data structure that contains the sets of nodes that represent the same dofs
		 */
		std::vector<std::vector<id_type> > _periodic_nodesets;

		/*
		 * Map from the periodic node id to the local vector of he periodic nodeset
		 */
		std::unordered_map<id_type, id_type> _periodic_id_to_set;

		// Number of nodes that will actually represent dofs
		id_type _n_global_non_periodic_nodes;
		id_type _n_local_non_periodic_nodes;
		id_type _n_local_owned_non_periodic_nodes;

		/*
		 * Function to build the periodic nodesets associated with a serial generated mesh
		 */
		void build_periodic_nodesets(id_type Nx, id_type Ny, id_type Nz);




		// INFORMATION ABOUT THE QUADRATURE POINTS
		//----------------------------------------------------------------------------------------------------------------------------------------------
		// Storage of shape function and Their gradients
		// storage order goes:
		//	-first index: local active element index number
		//	-second index: elemental quadrature point
		//	-third index: element local node number
		//	-fourth index: 
		std::vector<std::vector<std::vector<double> > > _shape;
		std::vector<std::vector<std::vector<std::vector<double> > > > _shape_grad;
		std::vector<std::vector<double> > _J;
		std::vector<std::vector<double> > _W;
		std::vector<std::vector<DenseMatrix<double> > > _Rot;
		std::vector<double> _volumetric_qps;
		
		
		
		
		// OUTSIDE OBJETCS THAT CONTAIN INFORMATION RELEVANT TO THE PROBLEM
		//----------------------------------------------------------------------------------------------------------------------------------------------
		typedef std::map<std::string, std::set<id_type> > nodeset_t;
		typedef std::map<std::string, std::set<std::pair<id_type, id_type> > > sideset_t;

		/*
		 * Data structure containing global node ids.
		 * Each nodeset will contain a list of node numbers associated with it
		 * Useful for applying BCs
		*/
		//std::vector<NamedSet> _nodesets;
		nodeset_t _nodesets;

		 /*
		  *	Data structure to hold sideset information
		  * Contains pairs of gloabl element numbers and the element local sides that are part of the set
		 */
		sideset_t _sidesets;

		// Bounding box on th mesh
		std::vector<double> _minimum_bounds;
		std::vector<double> _maximum_bounds;




		// INFORMATION ABOUT GENERAL MESH PARAMETERS AND BOOLEANS TO TRACK WHAT HAS BEEN DONE TO THE MESH
		//----------------------------------------------------------------------------------------------------------------------------------------------
		// Number of dimensions contained in this mesh (1, 2, or 3)
		id_type _dim;
		
		// MPI communicator for this mesh (should usually be MPI_COMM_WORLD)
		MPI_Comm _communicator;
		
		// Stores the MPI rank associated with this mesh partition
		int _rank;
		
		// Stores the how many ranks are associated with this mesh
		int _nranks;

		// Stores whether or not the mesh has been initialized
		bool _init;

		// Stores whether or not the mesh has been partitioned
		bool _partitioned;
		
		// Stores whether or not the _node_elem table has been generated
		bool _node_elem_generated;
		
		// Stores whether or not every element in the mesh has a material assigned to it
		bool _materials_assigned;
		
		// Stores if the mesh is serial
		bool _serial;

		// Stores whether or not this mesh has IGFEM components
		bool _IGFEM;

		// Stores whether or no the mesh was generated here or if it was read from Abaqus
		bool _was_generated;
		id_type _Nx, _Ny, _Nz;

		// Stores whether or not this mesh is periodic
		bool _periodic;
		
		// Command Line parameters
		int *_argc;
		char*** _argv;


		// Friend class declaration so partitioners can access private members
		friend class Partitioner;
		friend class PartitionerSerial;
		friend class DofObject;
		friend class BoundaryObject;
		friend class Refiner;
		friend class RefinerGeometricalInclusion;
		friend class RefinerGeometricalInclusionTraversal;



		
	public:
		/*
		 * Default constructor
		*/
		Mesh(int *argc, char***argv);

		/*
		 * Destructor
		*/
		~Mesh();
		
		// Clears all data structures and reinitializs the mesh
		void clear();

		/*
		 * Copy constructor
		*/
		Mesh(Mesh & other_mesh);
		
		// Routine to initialize anything that needs to be initialized before running stuff
		PetscErrorCode init(int *argc, char***argv);

		// Routine to finalize all data structures
		void Finalize();


		// Functions to manage the boundary
		// Add a new nodeset. To do this you need a nam and a vector of the global node ids that will belong to the nodeset
		void add_nodeset(std::string name, std::set<id_type>& nodes);

		// Add a new nodeset without any checks so its faster
		void add_nodeset_total(std::string name, std::set<id_type>& nodes);

		// Adds an individual value to the given nodeset
		void add_to_nodeset(std::string name, id_type val);

		// Add a new sideset. To do this you need a name and a vector of the global element ids and side ids that will belong to the elemset
		void add_sideset(std::string name, std::set<std::pair<id_type, id_type> >& sides);

		// Add a new sideset without checks so its faster
		void add_sideset_total(std::string name, std::set<std::pair<id_type, id_type> >& sides);

		// Adds an individual value to the given sideset
		void add_to_sideset(std::string name, std::pair<id_type, id_type> val);

		// Function to return the given nodeset
		std::set<id_type> get_nodeset(id_type idx);
		
		// Function to return the given nodeset
		std::set<id_type> get_nodeset(std::string name);

		// Function to return the given sideset
		std::set<std::pair<id_type, id_type> > get_sideset(id_type idx);
		
		// Function to return the given sideset
		std::set<std::pair<id_type, id_type> > get_sideset(std::string name);

		// Function to return the number of nodesets
		id_type n_nodesets() {return _nodesets.size();};

		// Function to return the number of sidesets
		id_type n_sidesets() {return _sidesets.size();};

		bool nodesetExists(std::string name);
		bool sidesetExists(std::string name);


		// Set iterator definitions
		typedef nodeset_t::iterator nodeset_iterator;
		typedef sideset_t::iterator sideset_iterator;
		nodeset_iterator nodesets_begin() {return _nodesets.begin();};
		nodeset_iterator nodesets_end() {return _nodesets.end();};
		sideset_iterator sidesets_begin() {return _sidesets.begin();};
		sideset_iterator sidesets_end() {return _sidesets.end();};



		// Function to return the MPI rank on which I belong
		bool own_node(id_type id) {return own_node_global(id);};
		bool own_node_global(id_type id);
		bool own_node_local(id_type idx);

		// Functions to return simple values that control mesh behavoir
		MPI_Comm get_comm() const {return _communicator;};
		int get_rank() const {return _rank;};
		id_type n_ranks() const {return _nranks;};
		bool serial() const {return _serial;};
		bool IGFEM() const {return _IGFEM;};
		bool generated() const {return _was_generated;};
		id_type nElemInDim(id_type dim) const;
		bool periodic() const {return _periodic;};

		// Some useful functions
		id_type dim() const {return _dim;};
		id_type n_nodes() const {return n_global_nodes();};
		id_type n_local_nodes() const {return _nodes.size();};
		id_type n_local_owned_nodes() const {return _n_local_owned_nodes;};
		id_type n_global_nodes() const {return _n_global_nodes;};
		id_type n_elem() const {return n_global_elem();};
		id_type n_local_elem() const {return _elem.size();};
		id_type n_global_elem() const {return _n_global_elem;};
		id_type n_local_non_periodic_nodes() const {return _n_local_non_periodic_nodes;};
		id_type n_local_owned_non_periodic_nodes() const {return _n_local_owned_non_periodic_nodes;};
		id_type n_global_non_periodic_nodes() const {return _n_global_non_periodic_nodes;};
		
		id_type global_to_local_node(id_type id);
		id_type global_to_local_elem(id_type id);
		std::vector<id_type>& get_elem_node_local(id_type idx);
		std::vector<id_type>& get_elem_node_global(id_type id);

		// Return true if the given node id is on a  periodic boundary
		bool node_periodic(id_type id);

		void set_periodic(bool periodic);

		// Returns the set of node ids that is associated with this periodic node id
		std::vector<id_type> get_node_periodic_nodeset(id_type id);

		// Returns the idx'th periodic nodeset
		std::vector<id_type> get_periodic_nodeset(id_type idx) {return _periodic_nodesets[idx];};

		// Returns how many periodic nodesets I have
		id_type n_periodic_sets() const {return _periodic_nodesets.size();};

		// Returns whether or not the given node id is the primary node of its periodic nodeset (the lowest id)
		bool primary_periodic_node(id_type id);

		// Returns true if I own this node and if its in a periodic nodeset, if this is the primary node
		bool check_node_responsibility(id_type id);

		// Iterates through the elements of the mesh and stores al the shape functions, gradients, gauss weihgts, and mapping jacobians at all o fthe geuss points of the local mesh
		void store_shape_functions();
		double& get_W(id_type local_e, id_type qp);
		double& get_J(id_type local_e, id_type qp);
		std::vector<double>& get_shape(id_type local_e, id_type qp);
		std::vector<std::vector<double> >& get_shape_grad(id_type local_e, id_type qp);
		DenseMatrix<double>& get_rot_matrix(id_type local_e, id_type qp);
		id_type n_volumetric_qps(id_type local_e);
		std::vector<double> get_global_coords_undeformed(id_type l_elem, id_type local_qp);
		
	
		/*
		* A function to add a new material to the mesh
		* To use this you should create a local instance of the material and then pass in the address to this function.
		* You must define all parameters and the name of the material. This function will handle the assignment of material ids
		* Ex:
		*	LinearElasticIsotropicMaterial material;
		*	material.set_parameter("E", 69000000000);
		*	material.set_parameter("nu", 0.334);
		*	material.set_name("Aluminum");
		*	mesh.add_material(&material);
		*/
		void add_material(Material* mat);

		// Function to Determine how many materials have been added to the mesh
		id_type n_materials() {return _materials.size();};
		
		// Function that returns the raw pointer to the idx'th material object. This will allow for modification of current materials in the mesh
		Material* get_material(id_type idx);
		
		// Function that returns the raw pointer to the material with the matching name. This will allow for modification of current materials in the mesh
		Material* get_material(std::string name);

		// Function that returns the index in the materials vector that corresponds to the material name (return -1 for an invalid material)
		int find_mat_number(std::string name);

		// Function that just returns whatever is in the _elem_to_material entry
		id_type get_mat_num_local(id_type idx) {return _elem_to_material[idx];};

		// Function to return a reference to the full material vector
		std::vector<Material*>& get_materials() {return _materials;};

		// Function that a user should call if they wish to remain ignorant of mesh partitioning. Simply calls get_element_material_global
		Material* get_elem_material(id_type id, id_type int_el = 0);

		// Function that returns a pointer to the material contained in the element with global id id
		Material* get_element_material_global(id_type id, id_type int_el = 0);

		// Function that returns a pointer to the material contained in the local element idx
		Material* get_element_material_local(id_type idx, id_type int_el = 0);
		
		// Add a new elemset. To do this you need a name and a vector of the global element ids that will belong to the elemset
		void add_elemset(std::string name, std::set<id_type>& elems);

		void add_to_elemset(std::string name, id_type val);
		
		// Function to return the given elemset
		std::set<id_type> get_elemset(id_type idx);
		
		// Function to return the given elemset
		std::set<id_type> get_elemset(std::string name);

		// Function to return the number of elemsets
		id_type n_elemsets() {return _elemsets.size();};

		// Set iterator definitions
		typedef std::map<std::string, std::set<id_type> >::iterator elemset_iterator;
		elemset_iterator elemsets_begin() {return _elemsets.begin();};
		elemset_iterator elemsets_end() {return _elemsets.end();};
		
		/*
		 * Functions to assign materials to elements
		 *	1. Assigns the material with the given material_id to every element of the mesh
		 *	2. Assigns the material with the given material_name to every element of the mesh
		 *	3. Assigns the material with the given material_id to the elements with the given list of global element ids
		 *	4. Assigns the material with the given material_name to the elements with the given list of global element ids
		 *	6. Assigns the material with the given material_name to the elements in the elemset with the given name
		*/
		void set_material(id_type id);
		void set_material(std::string name);
		void set_material(id_type id, std::vector<id_type>& elems);
		void set_material(std::string name, std::vector<id_type>& elems);
		void set_material_from_elemset(std::string mat_name, std::string elem_name);
		
		// Function that makes it easy to determine if a vector of global node ids are all in the local mesh.
		// Returns true if all ids exist in the local mesh. Returns false if any node id does not exist in the local mesh.
		bool nodes_in_local_mesh(std::vector<id_type>& nodes);
		
		// Function that makes it easy to determine if a vector of global node ids are all in the local mesh.
		// Returns true if all ids exist in the local mesh. Returns false if any node id does not exist in the local mesh.
		bool nodes_in_local_mesh(std::set<id_type>& nodes);
		
		// Nice wrapper function to determine if a global node id exists in the local mesh
		bool node_in_local_mesh(id_type id);

		// Function that makes it easy to determine if a vector of global element ids are all in the local mesh.
		// Returns true if all ids exist in the local mesh. Returns false if any element id does not exist in the local mesh.
		bool elems_in_local_mesh(std::vector<id_type>& elems);
		
		// Function that makes it easy to determine if a vector of global element ids are all in the local mesh.
		// Returns true if all ids exist in the local mesh. Returns false if any element id does not exist in the local mesh.
		bool elems_in_local_mesh(std::set<id_type>& elems);

		// Nice wrapper function to determine if a global elem id exists in the local mesh
		bool elem_in_local_mesh(id_type id);

		// Finds all of the partitions that the list of global node ids has in common
		std::vector<int> partitions_in_common(const std::vector<id_type>& nodes);
		
		// Function to determine if the given global node id is located on a partition interface
		bool node_on_part_interface(id_type id);

		// Function to determine if a group of nodes are all on a partition interface
		bool nodes_on_part_interface(const std::vector<id_type>& node_group);

		// Get the vector of partitions that have a copy of this node
		std::vector<int> get_node_parts(id_type id);
		
		// Function to determine how many patitions have a copy of this node. Return 1 if the node is not n a partition interface.
		id_type n_parts_with_node(id_type id);

		// Used to iterate through which partitions have a copy of this node. Usually would be called after calling n_parts_with_node(id).
		id_type get_part_int_from_node(id_type id, id_type interface);
		
		// Function returns the raw pointer to the local element idx. Good for element modification. Use carefully.
		Elem* get_elem_local(id_type idx);

		// Function returns the raw pointer to the global element id. Good for element modification. Use carefully.
		Elem* get_elem_global(id_type id);
		
		// Function returns the raw pointer to the local node idx. Good for node modification. Use carefully.
		Node* get_node_local(id_type idx);
		
		// Function returns the raw pointer to the global node id. Good for node modification. Use carefully.
		Node* get_node_global(id_type id);

		// Function that return the mesh partition number that corresponds to the owner of the given global node id
		int get_node_owner_global(id_type id);
		
		// Function that return the mesh partition number that corresponds to the owner of the given local node number
		int get_node_owner_local(id_type idx);
		
		// Returns the full node-elem result, which is both the node_elem_local and remote_node_elem. Based on global node number
		std::vector<id_type> get_node_elem(id_type id);

		// Returns the global element numbers of all local elements neighboring the node with the given local node id
		// NOTE: local elements means that this will only return element ids tha exist on the local processor. To get element ids on another processor use get_remote_node_elem.
		//		 Should really abstract that away so that there's a function to get all neighboring global element ids regardless of which processor they're on.
		std::vector<id_type>& get_node_elem_local(id_type idx);
		
		// Returns the global element numbers of all elements neighboring the node with the given global node id
		// NOTE: local elements means that this will only return element ids tha exist on the local processor. To get element ids on another processor use get_remote_node_elem.
		//		 Should really abstract that away so that there's a function to get all neighboring global element ids regardless of which processor they're on.
		std::vector<id_type>& get_node_elem_global(id_type id);
		
		// Returns the global element numbers of all remote elements neighboring the node with the given global node id
		// NOTE: if this node is not on a boundary then this function returns an empty vector
		std::vector<id_type> get_remote_node_elem(id_type id, int rank);

		// Find the global elements (local or remote) that the list of global nodes has in common
		std::vector<id_type> elem_in_common(const std::vector<id_type>& nodes);

		std::vector<id_type> local_elem_in_common(const std::vector<id_type>& nodes);

		std::map<int, std::vector<id_type> > remote_elem_in_common(const std::vector<id_type>& nodes);
		
		void delete_elem(id_type id);
		Elem* add_elem(elem_type etype, std::vector<id_type>& node_ids, id_type global_id);
		Node* add_node(double x, double y, double z, id_type global_id, int owner);
		void set_mesh_dimension(id_type dim) {_dim = dim;};
		
		
		// Element and node iterator definitions
		typedef std::vector<Elem*>::iterator element_iterator;
		typedef std::vector<Node*>::iterator node_iterator;
		
		element_iterator elements_begin() {return _elem.begin();};
		element_iterator elements_end() {return _elem.end();};
		node_iterator nodes_begin() {return _nodes.begin();};
		node_iterator nodes_end() {return _nodes.end();};

		// Partition interface iterators
		typedef std::unordered_map<id_type, std::vector<int> >::iterator partition_iterator;
		partition_iterator partition_interface_begin() {return _node_partition_interface.begin();};
		partition_iterator partition_interface_end() {return _node_partition_interface.end();};
		
		
		
		
		
		
		
		
		// MOST USEFUL FUNCTIONS GO HERE. THESE FUNCTIONS WILL NORMALLY BE CALLED BY THE USER
		//----------------------------------------------------------------------------------------------------------------------------------------------
		
		// Generate a mesh from some parameters
		void generate_mesh(elem_type e_type, double min_x, double max_x, id_type Nx, double min_y=0.0, double max_y=0.0, id_type Ny=0, double min_z=0.0, double max_z=0.0, id_type Nz=0);
		void getDomainLimits(std::vector<double>& mins, std::vector<double>& maxes);

		// Rotate the mesh about the origin. Angles are about the fixed z-y-x coordinate axes
		void rotate_mesh(double z_angle, double y_angle=0.0, double x_angle=0.0);

		// Create the node-elem table
		void generate_node_elem();
		
		// Read a mesh from an input file (Abaqus)
		void read_mesh(const std::string& name);









	// EVERYTHING FROM HERE ON IS RELATED TO IGFEM
	// ============================================================================================================================================

	private:
		// INFORMATION ABOUT REGULAR AND ENRICHMENT NODES
		//-----------------------------------------------------------------------------------------------------------------------
		
		// Structure to contain the list of all local enrichment node
		std::vector<Node*> _enrich_nodes;
		
		// Map from the global to local enrichment node ids
		std::unordered_map<id_type, id_type> _global_to_local_enrich_node;
		
		// Structure that contains the owner of all of the enricment nodes
		std::vector<int> _enrich_owners;
		
		// Vector containing the nodal detection of each local node 
		// negative if it belongs to a subdomain material, index of the inclusion otherwise
		std::vector<int> _node_detect;
		
		// Stores number of global enrichment nodes
		id_type _n_global_enrich_nodes;

		// Stores the number of enrichment nodes that this processor owns
		id_type _n_local_owned_enrich_nodes;

		// Map that contains lists of interface with copies of the given global enrichment node id
		std::unordered_map<id_type, std::vector<int> > _enrich_nodes_interface_lists;

		// Data structure containing map from local node numbers to its GLOBAL element numbers
		std::vector< std::vector<id_type> > _enrich_node_elem;

		// Data structure to contain the local node numbers of all elements that have been intersected somehow
		std::vector<id_type> _intersected_elem;

		// maps from the local enrichment nodes to the edge node pair that it is on
		// FIXME: Is it better to store this the other way around?
		std::vector<std::pair<id_type, id_type> > _enrich_node_to_edge;
		std::vector<id_type> _enrich_node_to_inclusion;

		// INFORMATION ABOUT INCLUSIONS AND THEIR MATERIAL PROPERTIES
		//-----------------------------------------------------------------------------------------------------------------------
		
		// Vector of pointers to all of the inclusion objects
		std::vector<Inclusion*> _inclusions;

		// INFORMATION ABOUT THE COHESIVE NATURE OF THE MESH
		//-----------------------------------------------------------------------------------------------------------------------

		// Points to the materials in the _materials vector that are cohesive
		std::vector<id_type> _cohesive_mats;
		std::vector<std::pair<int, int> > _cohesive_to_mat_pair;

		// MISCELLANEOUS STUFF
		//-----------------------------------------------------------------------------------------------------------------------
		
		// Various Booleans
		bool _nodes_detected;
		bool _elements_detected;
		bool _is_cohesive; // Stores whether or not we are using a cohesive 




	public:

		// FUNCTIONS THAT ARE UNIQUE TO IGFEM (MAINLY RELATING TO INCLUSIONS AND ENRICHMENT NODES)
		//--------------------------------------------------------------------------------------------------

		// Function which resets the mesh to the Pre-IGFEM state
		void clearEnrichments();

		// Function that returns whether or not the mesh has cohesive materials
		bool is_cohesive() {return _is_cohesive;};

		void add_cohesive_material_pair(std::string coh_mat, std::string mat1, std::string mat2);

		// A function to adda new inclusion to the mesh
		void add_inclusion(Inclusion* inc);

		// gets the idx'th inclusion in the list
		Inclusion* get_inclusion(id_type idx);
		
		id_type n_inclusions() {return _inclusions.size();};
		
		// Analyze all nodes and determine which inclusion they lie in
		void analyze_nodes();
		
		// Analyze all elements and determine if/how they are intersected by the inclusions 
		void analyze_elements();
		
		// Function to add all enrichment nodes to the mesh, keeping consisten across processors
		void add_enrichments();

		id_type n_local_enrich_nodes() {return _enrich_nodes.size();};

		id_type n_global_enrich_nodes() {return _n_global_enrich_nodes;};

		id_type n_local_owned_enrich_nodes() {return _n_local_owned_enrich_nodes;};

		std::pair<id_type, id_type> get_enriched_edge_local(id_type idx);
		id_type get_enriched_inclusion_local(id_type idx);

		// Function to compute how far along an edge the given enrichment node is
		double compute_enrich_node_interpolation(id_type id);

		// Iterators
		typedef std::vector<Node*>::iterator enrich_node_iterator;
		enrich_node_iterator enrich_nodes_begin() {return _enrich_nodes.begin();};
		enrich_node_iterator enrich_nodes_end() {return _enrich_nodes.end();};
		typedef std::unordered_map<id_type, std::vector<int> >::iterator enrich_partition_iterator;
		enrich_partition_iterator enrich_partition_interface_begin() {return _enrich_nodes_interface_lists.begin();};
		enrich_partition_iterator enrich_partition_interface_end() {return _enrich_nodes_interface_lists.end();};

		// DEBUGGING FUNCTIONS
		std::vector<int> get_nodal_detection() {return _node_detect;};
		Node* get_enrich_node(id_type idx) {return _enrich_nodes[idx];}; // This can be done with get_node but oh well

		std::pair<id_type, id_type> get_edge_from_enrich_node(id_type enrich_id);

		int get_nodal_detection_local(id_type idx);


	private:

		// Helper function for nodal detection
		int detect_node(Node* node);

		// Helper function for nodal detection
		void reprocess_nodes(std::unordered_map<id_type, Inclusion*>& nodes_to_reprocess);

		// Helper function for actually moving the Nodal position
		void modify_node_location(Node* node, Inclusion* inc);

		// Helper function for adding enrichments
		void add_enrich_node_to_elem(std::pair<id_type, id_type> edge, id_type inclusion, id_type n, id_type n_initial);




	// EVERYTHING FROM HERE ON IS RELATED TO MESH REFINEMENT
	// ============================================================================================================================================

	private:

		// Defines how many times the mesh has been refined
		unsigned char _mesh_topology; 

		// Vector that is n_local_active_elem() long and contains pointers ot all of the local active elements (no children)
		std::vector<Elem*> _active_elem;

		// Vector that us n_local_elem() long and contains booleans for every element in the mesh
		// true if the local element is active and false otherwise. Used for telling if any arbitrary element is active quickly
		std::vector<bool> _elem_activity;

		// Vector of the p-level of all of the elements in the mesh
		std::vector<unsigned char> _p_level;

		// Vector that contains the parentage for the elements of the mesh (-1 for no parents)
		// (Local or global? Not sure yet. Probably local)
		std::vector<int> _elem_parents;

		// Maps from a vector of global node ids to the associated refienement node
		std::unordered_map<std::vector<id_type>, id_type> _refine_node_neighbors;

		// Store the global ids of the hangining nodes
		std::set<id_type> _edge_hanging_nodes;
		std::set<id_type> _face_hanging_nodes;

		// Stores the constraints used for each hanging node in a list of pairs
		// The first entry is the global node id of a node that this hanging node should be constrained with
		// The second entry is the coefficient in the constraint equation corresponding to that node
		std::unordered_map<id_type, std::vector<std::pair<id_type, double> > > _hanging_node_constraints;
		std::unordered_map<id_type, std::vector<id_type> > _hanging_node_neighbors;

		// Number storage
		id_type _n_global_active_elem;
		// Store the number of local non-refinement nodes
		id_type _n_original_nodes;
		// Store number of base parenet elements
		id_type _n_original_elements;

		// In the case of an enrichment occuring along an edge that has a hanging node on it, maps from the parent (neighbor nodes) to the refined nodes
		//std::unordered_map<std::pair<id_type, id_type>, std::pair<id_type, id_type> > parent_to_refined_edge;e

	public:

		// Get how many times the mesh has been refined
		unsigned char get_mesh_topology() {return _mesh_topology;};

		id_type n_local_active_elem() {return _active_elem.size();};
		id_type n_global_active_elem() {return _n_global_active_elem;};

		// iterators
		element_iterator active_elements_begin() {return _active_elem.begin();};
		element_iterator active_elements_end() {return _active_elem.end();};

		// Searches the _refine_node_neighbors structure for maatching global node id for the given neighbors
		id_type get_refine_node_from_neighbors(const std::vector<id_type>& neighbors);

		typedef std::unordered_map<id_type, std::vector<std::pair<id_type, double> > >::iterator hanging_node_iterator;
		hanging_node_iterator hanging_nodes_begin() {return _hanging_node_constraints.begin();};
		hanging_node_iterator hanging_nodes_end() {return _hanging_node_constraints.end();};

		// Gets the parent node neighbors (nodes that it was built from) for each of the hanging nodes
		std::vector<id_type>& get_hanging_node_neighbors(id_type id);

		// Returns whether or not the global node id is a hanging node
		bool is_hanging_node(id_type id);

		// Function to analyze the hanging node constraints and determine the constraint equations
		void set_constraint_values();

		id_type n_hanging_nodes() {return _edge_hanging_nodes.size() + _face_hanging_nodes.size();};

		std::vector<std::pair<id_type, double> >& get_hanging_node_constraint(id_type id);

		id_type n_pre_refinement_nodes() const {return _n_original_nodes;};

		unsigned char get_p_level_local(id_type idx);

	private:

		// Function to detect which of my local nodes are hanging
		void detect_hanging_nodes();

		void set_active_elements();

		// If any of the enrichment nodes lie on an enriched edge or face, this detects that and lists whichglobal enrichment nodes are associated with each hanging node
		void find_hanging_enriched_neighbors(std::unordered_map<id_type, id_type>& edges_map, std::unordered_map<id_type, std::vector<id_type> >& faces_map);
};

#endif
