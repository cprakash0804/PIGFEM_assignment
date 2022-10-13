/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#ifndef _Writer_ASCII_VTK_NODAL_H_
#define _Writer_ASCII_VTK_NODAL_H_

#include "Writer_VTK_Nodal.h"
#include "Utilities.h"
#include <map>




/*
 * This class is used to write the Problem data
 * to an ASCII text file following the legacy
 * VTK file format described at:
 * http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 */
class Writer_ASCII_VTK_Nodal : public Writer_VTK_Nodal
{
	public:


		/*
		 * The Big 3
		 * Constructor need the Problem
		 * Destructor doesn't do anything
		 * Copy Constructor copies the problem pointer and stores the mesh
		 */
		Writer_ASCII_VTK_Nodal(Problem* prob);
		virtual ~Writer_ASCII_VTK_Nodal();
		Writer_ASCII_VTK_Nodal(const Writer_ASCII_VTK_Nodal& other);


		/*
		 * Virtual function to open the stream in ascii or binary mode
		 */
		virtual void open(std::ofstream& myfile, std::string& filename);


	private:
		

		/*
		 * Writes the standard ASCII header for a legacy VTK file
		 */
		virtual void write_header(std::ofstream& myfile, double curr_t);


		/*
		 * Writes all of the nodal points to the VTK format. For serial ouput simply pass in empty remote lists
		 */
		virtual void write_nodes(std::ofstream& myfile,
							     std::vector<id_type>& local_to_VTK,
							     const std::vector<id_type>& id_list, const std::vector<float>& coords_list,
							     const std::vector<std::vector<id_type> >& remote_id_lists, const std::vector<std::vector<float> >& remote_coords_lists);


		/*
		 * Writes all of the elements to the VTK format. For serial ouput simply pass in empty remote list
		 */
		virtual void write_elements(std::ofstream& myfile,
									const std::vector<id_type>& local_to_VTK,
									const std::vector<id_type>& eptr, const std::vector<id_type>& eind, const std::vector<unsigned char>& elem_to_type,
									const std::vector<std::vector<id_type> >& remote_eptr, const std::vector<std::vector<id_type> >& remote_eind,
									const std::vector<std::vector<unsigned char> >& remote_etype);


		/*
		 * Writes all of the nodal data to the VTK format. For serial ouput simply pass in empty remote lists
		 */
		virtual void write_nodal_data(std::ofstream& myfile,
									  const std::vector<float>& sol_list, const std::vector<float>& rhs_list,
									  const std::vector<std::vector<float> >& avg_strain, const std::vector<std::vector<float> >& avg_stress,
									  const std::vector<std::vector<float> >& remote_sol_lists, const std::vector<std::vector<float> >& remote_rhs_lists,
									  const std::vector<std::vector<std::vector<float> > >& remote_strain, const std::vector<std::vector<std::vector<float> > >& remote_stress);


		/*
		 * Writes all of the element data to the VTK format. For serial ouput simply pass in empty remote lists
		 */
		virtual void write_element_data(std::ofstream& myfile,
										const std::vector<id_type>& elem_mats, const std::map<std::string, std::vector<float> >& internal_vars_to_avgs,
										const std::vector<std::vector<id_type> >& remote_mats, const std::vector<std::vector<std::vector<float> > >& remote_int_vars);
};





#endif