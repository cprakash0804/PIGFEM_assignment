/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#ifndef _Writer_VTK_NODAL_H_
#define _Writer_VTK_NODAL_H_

#include "Writer.h"
#include "Utilities.h"
#include <map>
#include <unordered_map>




/*
 * This class is used to write the Problem data
 * to an ASCII text file following the legacy
 * VTK file format described at:
 * http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 */
class Writer_VTK_Nodal : public Writer
{
	public:


		/*
		 * The Big 3
		 * Constructor need the Problem
		 * Destructor doesn't do anything
		 * Copy Constructor copies the problem pointer and stores the mesh
		 */
		Writer_VTK_Nodal(Problem* prob);
		Writer_VTK_Nodal();
		virtual ~Writer_VTK_Nodal();
		Writer_VTK_Nodal(const Writer_VTK_Nodal& other);


		/*
		 * Main function used to write problem to any inherited file type
		 */
		virtual void write(std::string filename, double curr_t);


		/*
		 * Stores the mesh in this Writer's internal storage
		 * Can be recalled if the mesh is adapted
		 */
		virtual void store_mesh();

		virtual void setParameter(std::string param, bool val);
		virtual bool getBoolParameter(std::string param);


	protected:


		/*
		 * Defining the map from my elmeent types to the VTK numbering system
		 */
		static std::map<elem_type, unsigned char> create_map()
		{
			std::map<elem_type, unsigned char> m;
        	m[POINT1] = 1;
			m[EDGE2] = 3;
			m[EDGE3] = 21;
			m[TRI3] = 5;
			m[TRI6] = 22;
			m[QUAD4] = 9;
			m[QUAD8] = 23;
			m[TET4] = 10;
			m[TET10] = 24;
			m[HEX8] = 12;
			m[HEX20] = 25; // HEX27 doesn't have a type...
			m[PRISM6] = 13;// None of the other prism types have VTK types...
			m[POLYGON] = 7; // General polygon
			return m;
        }
        static std::map<elem_type, unsigned char> elem_type_to_VTK;



        static std::map<coh_elem_type, elem_type> create_cohesive_type_map()
        {
        	std::map<coh_elem_type, elem_type> m;
        	m[COHPOINT1] = EDGE2;
        	m[COHEDGE2] = QUAD4;
        	m[COHTRI3] = PRISM6;
        	m[COHQUAD4] = HEX8;
        	return m;
        }
        static std::map<coh_elem_type, elem_type> cohesive_elem_to_plot_elem;



        /*
		 *  Whether or not to output the cohesive element contributions
		 */
		bool _output_cohesive;


		/*
		 * Function to compute the associated cohesive damage for a given cohesive element
		 */
		float compute_cohesive_damage(CohesiveElem* coh_el, Material* mat, id_type l_elem, id_type initial_qp);


		/*
		 * Private data members that store the mesh in a useful format for this writer
		 */
		std::vector<id_type> _node_id_list;
		std::vector<float> _node_coord_list;
		std::vector<id_type> _eptr, _eind, _elem_mats;
		std::vector<unsigned char> _elem_types;
		std::vector<float> _elem_vols;
		std::unordered_map<id_type, std::vector<float> > _int_elem_vols; // Store the voluem fo the child elements for the intersected elements
		id_type _n_elem, _n_points;


		/*
         * The actual function that will do all of the writing (maybe by calling other functions)
         */
        virtual void writeFromStream(std::ofstream& myfile, double curr_t);


		/*
		 * Stores the node information in the private data members used to store nodes
		 */
		void store_nodes(std::vector<id_type>& id_list, std::vector<float>& coords_list);


		/*
		 * Stores the element information in the private data members
		 */
		void store_elements(std::vector<id_type>& eptr, std::vector<id_type>& eind,
							std::vector<unsigned char>& types, std::vector<id_type>& mats);


		/*
		 * Writes the nodes and element information
		 */
		void write_mesh(std::ofstream& myfile);


		/*
		 * Gather and write data associated with the nodes (solution, external load)
		 */
		void write_node_info(std::ofstream& myfile);


		/*
		 * Gather and write data associated with the elements (stress, strain, internal variables, material)
		 */
		void write_element_info(std::ofstream& myfile);


		/*
		 * Function that communicates all of the existing mesh data to proc 0
		 */
		void communicate_mesh(std::vector<std::vector<id_type> >& remote_node_ids, std::vector<std::vector<float> >& remote_coords,
							  std::vector<std::vector<id_type> >& remote_eptr, std::vector<std::vector<id_type> >& remote_eind, std::vector<std::vector<unsigned char> >& remote_etype);


		/*
		 * Communicated solution and rhs data to proc 0
		 */
		void communicate_node_data(const std::vector<float>& sol_list, const std::vector<float>& rhs_list,
								   const std::vector<std::vector<float> >& avg_strain, const std::vector<std::vector<float> >& avg_stress,
								   std::vector<std::vector<float> >& remote_sol, std::vector<std::vector<float> >& remote_rhs,
								   std::vector<std::vector<std::vector<float> > >& remote_strain, std::vector<std::vector<std::vector<float> > >& remote_stress);


		/*
		 * Communicate, materials, internal variables, stresses and strains to proc 0
		 */
		void communicate_element_data(const std::vector<id_type>& _elem_mats, const std::map<std::string, std::vector<float> >& internal_vars_to_avgs,
									  std::vector<std::vector<id_type> >& remote_mats, std::vector<std::vector<std::vector<float> > >& remote_int_vars);


		/*
		 * Gathers info about the solution and rhs in a serialized form
		 */
		void gather_nodal_data(std::vector<float>& sol_list, std::vector<float>& rhs_list,
							   std::vector<std::vector<float> >& strain_list, std::vector<std::vector<float> >& stress_list);
		void project_stress_strain(std::vector<std::vector<float> >& strain_list, std::vector<std::vector<float> >& stress_list);
		void fill_local_stress_strain_contributions(std::vector<std::vector<std::vector<float> > >& node_strain_contributions,
													std::vector<std::vector<std::vector<float> > >& node_stress_contributions,
													std::vector<std::vector<float> >& node_volume_contributions);
		void communicate_stress_strain_contributions(std::vector<std::vector<std::vector<float> > >& node_strain_contributions,
													 std::vector<std::vector<std::vector<float> > >& node_stress_contributions,
													 std::vector<std::vector<float> >& node_volume_contributions);


		/*
		 * Gathers info about the internal variables, stress, and strain in a serialized form
		 */
		void gather_element_data(std::map<std::string, std::vector<float> >& internal_vars_to_avgs);


		/*
		 * Virtual function to open the stream in ascii or binary mode
		 */
		virtual void open(std::ofstream& myfile, std::string& filename) = 0;


		/*
		 * Writes the standard (ASCII/binary) header for a legacy VTK file
		 */
		virtual void write_header(std::ofstream& myfile, double curr_t) = 0;


		/*
		 * Writes all of the nodal points to the VTK format. For serial ouput simply pass in empty remote lists
		 */
		virtual void write_nodes(std::ofstream& myfile,
							     std::vector<id_type>& local_to_VTK,
							     const std::vector<id_type>& id_list, const std::vector<float>& coords_list,
							     const std::vector<std::vector<id_type> >& remote_id_lists, const std::vector<std::vector<float> >& remote_coords_lists) = 0;


		/*
		 * Writes all of the elements to the VTK format. For serial ouput simply pass in empty remote list
		 */
		virtual void write_elements(std::ofstream& myfile,
									const std::vector<id_type>& local_to_VTK,
									const std::vector<id_type>& eptr, const std::vector<id_type>& eind, const std::vector<unsigned char>& elem_to_type,
									const std::vector<std::vector<id_type> >& remote_eptr, const std::vector<std::vector<id_type> >& remote_eind,
									const std::vector<std::vector<unsigned char> >& remote_etype) = 0;


		/*
		 * Writes all of the nodal data to the VTK format. For serial ouput simply pass in empty remote lists
		 */
		virtual void write_nodal_data(std::ofstream& myfile,
									  const std::vector<float>& sol_list, const std::vector<float>& rhs_list,
									  const std::vector<std::vector<float> >& avg_strain, const std::vector<std::vector<float> >& avg_stress,
									  const std::vector<std::vector<float> >& remote_sol_lists, const std::vector<std::vector<float> >& remote_rhs_lists,
									  const std::vector<std::vector<std::vector<float> > >& remote_strain, const std::vector<std::vector<std::vector<float> > >& remote_stress) = 0;


		/*
		 * Writes all of the element data to the VTK format. For serial ouput simply pass in empty remote lists
		 */
		virtual void write_element_data(std::ofstream& myfile,
										const std::vector<id_type>& elem_mats, const std::map<std::string, std::vector<float> >& internal_vars_to_avgs,
										const std::vector<std::vector<id_type> >& remote_mats, const std::vector<std::vector<std::vector<float> > >& remote_int_vars) = 0;
};





#endif