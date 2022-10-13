/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#include "Writer_ASCII_VTK_Nodal.h"
#include "Problem.h"
#include "Mesh.h"
#include "InternalVars.h"
#include "NodalData.h"
#include "Utilities.h"
#include "mpi.h"



/*
 * The Big 3
 * Constructor need the Problem
 * Destructor doesn't do anything
 * Copy Constructor copies the problem pointer
 */
Writer_ASCII_VTK_Nodal::Writer_ASCII_VTK_Nodal(Problem* prob)
	: Writer_VTK_Nodal(prob)
{}
Writer_ASCII_VTK_Nodal::~Writer_ASCII_VTK_Nodal()
{}
Writer_ASCII_VTK_Nodal::Writer_ASCII_VTK_Nodal(const Writer_ASCII_VTK_Nodal& other)
{
	_prob = other.get_prob();
	store_mesh();
}




/*
 * Virtual function to open the stream in ascii or binary mode
 */
void Writer_ASCII_VTK_Nodal::open(std::ofstream& myfile, std::string& filename)
{
	myfile.open(filename.c_str(), std::ofstream::out);
	if (!myfile.good())
	{
		std::string out = "Error opening the output file " + filename;
		err_message( out.data() );
	}
}




/*
 * Writes the standard ASCII header for a legacy VTK file
 */
void Writer_ASCII_VTK_Nodal::write_header(std::ofstream& myfile, double curr_t)
{
	Mesh* mesh = _prob->get_mesh();
	myfile << "# vtk DataFile Version 3.0\n";
	if (mesh->IGFEM())
		myfile << "IGFEM " << mesh->dim() << "D simulation: t=" << curr_t << "\n";
	else
		myfile << "FEM " << mesh->dim() << "D simulation: t=" << curr_t << "\n";
	myfile << "ASCII\n\n";

	// Print the Dataset
	myfile << "DATASET UNSTRUCTURED_GRID\n";
}




/*
 * Writes all of the nodal points to the VTK format. For serial ouput simply pass in empty remote lists
 */
void Writer_ASCII_VTK_Nodal::write_nodes(std::ofstream& myfile,
								   std::vector<id_type>& local_to_VTK,
								   const std::vector<id_type>& id_list, const std::vector<float>& coords_list,
								   const std::vector<std::vector<id_type> >& remote_id_lists, const std::vector<std::vector<float> >& remote_coords_lists)
{
	Mesh* mesh = _prob->get_mesh();

	// Only print if this is rank 0
	if (mesh->get_rank()==0)
	{
		id_type n_ranks = remote_id_lists.size() + 1;

		// Set up the mapping vector from global ids to VTK numbering
		local_to_VTK.resize(mesh->n_global_nodes() + mesh->n_global_enrich_nodes());

		// Print out the nodes themselves
		id_type curr_node = 0;
		myfile << "POINTS " << (mesh->n_global_nodes()+mesh->n_global_enrich_nodes()) << " float\n";
		for(id_type n=0; n<id_list.size(); ++n)
		{
			id_type id = id_list[n];
			if(id >= ENRICH_START)
				id = id - ENRICH_START + mesh->n_global_nodes();
			local_to_VTK[id] = curr_node;
			myfile << coords_list[n*3] << " " << coords_list[n*3+1] << " " << coords_list[n*3+2] << "\n";
			curr_node++;
		}
		for(id_type part=1; part<n_ranks; ++part)
		{
			for(id_type n=0; n<remote_id_lists[part-1].size(); ++n)
			{
				id_type id = remote_id_lists[part-1][n];
				if(id >= ENRICH_START)
					id = id - ENRICH_START + mesh->n_global_nodes();
				local_to_VTK[id] = curr_node;
				myfile << remote_coords_lists[part-1][n*3] << " " << remote_coords_lists[part-1][n*3+1] << " " << remote_coords_lists[part-1][n*3+2] << "\n";
				curr_node++;
			}
		}
	}
}




/*
 * Writes all of the elements to the VTK format. For serial ouput simply pass in empty remote list
 */
void Writer_ASCII_VTK_Nodal::write_elements(std::ofstream& myfile,
									  const std::vector<id_type>& local_to_VTK,
									  const std::vector<id_type>& eptr, const std::vector<id_type>& eind, const std::vector<unsigned char>& elem_to_type,
									  const std::vector<std::vector<id_type> >& remote_eptr, const std::vector<std::vector<id_type> >& remote_eind,
									  const std::vector<std::vector<unsigned char> >& remote_etype)
{
	Mesh* mesh = _prob->get_mesh();

	// Only print if this is rank 0
	if (mesh->get_rank()==0)
	{
		id_type n_ranks = mesh->n_ranks();

		id_type n_elem = eptr.size() - 1;
		id_type n_points = eind.size();
		for(id_type part=1; part<mesh->n_ranks(); ++part)
		{
			n_elem += (remote_eptr[part-1].size() - 1);
			n_points += remote_eind[part-1].size();
		}
		// Print all the element-node table information
		myfile << "\nCELLS " << n_elem << " " << (n_elem+n_points) << "\n";
		for(id_type e=0; e<(eptr.size()-1); ++e)
		{
			id_type start = eptr[e];
			id_type end = eptr[e+1];
			myfile << (end-start) << " ";
			for(id_type idx=start; idx<end; ++idx)
			{
				id_type id = eind[idx];
				if(id >= ENRICH_START)
					id = id - ENRICH_START + mesh->n_global_nodes();
				myfile << local_to_VTK[id];
				if(idx==(end-1)) myfile << "\n";
				else myfile << " ";
			}
		}
		for(id_type part=1; part<n_ranks; ++part)
		{
			for(id_type e=0; e<(remote_eptr[part-1].size()-1); ++e)
			{
				id_type start = remote_eptr[part-1][e];
				id_type end = remote_eptr[part-1][e+1];
				myfile << (end-start) << " ";
				for(id_type idx=start; idx<end; ++idx)
				{
					id_type id = remote_eind[part-1][idx];
					if(id >= ENRICH_START)
						id = id - ENRICH_START + mesh->n_global_nodes();
					myfile << local_to_VTK[id];
					if(idx==(end-1)) myfile << "\n";
					else myfile << " ";
				}
			}
		}

		// Print out the cell types data
		myfile << "\nCELL_TYPES " << n_elem << "\n";
		for(id_type e=0; e<elem_to_type.size(); ++e)
			myfile << static_cast<int>(elem_to_type[e]) << "\n";
		for(id_type part=1; part<n_ranks; ++part)
		{
			for(id_type e=0; e<remote_etype[part-1].size(); ++e)
				myfile << static_cast<int>(remote_etype[part-1][e]) << "\n";
		}
	}
}




/*
 * Writes all of the nodal data to the VTK format. For serial ouput simply pass in empty remote lists
 */
void Writer_ASCII_VTK_Nodal::write_nodal_data(std::ofstream& myfile,
									  const std::vector<float>& sol_list, const std::vector<float>& rhs_list,
									  const std::vector<std::vector<float> >& avg_strain, const std::vector<std::vector<float> >& avg_stress,
									  const std::vector<std::vector<float> >& remote_sol_lists, const std::vector<std::vector<float> >& remote_rhs_lists,
									  const std::vector<std::vector<std::vector<float> > >& remote_strain, const std::vector<std::vector<std::vector<float> > >& remote_stress)
{
	// Get the mesh
	Mesh* mesh = _prob->get_mesh();

	if (mesh->get_rank() == 0)
	{
		id_type n_ranks = mesh->n_ranks();
		id_type nndof = (_prob->get_classification()==STRUCTURAL) ? 3 : _prob->nndof(); // Set the number of entries for each solution component

		// Print out solution information
		myfile << "\nPOINT_DATA " << (mesh->n_global_nodes() + mesh->n_global_enrich_nodes()) << "\n";
		if(_prob->get_classification()==STRUCTURAL)
			myfile << "VECTORS displacement float " << "\n";
		else
		{
			myfile << "SCALARS solution float " << nndof << "\n"; // I don't actually know if this is right
			myfile << "LOOKUP_TABLE default\n";
		}
		for(id_type n=0; n<(sol_list.size()/nndof); ++n)
		{
			for(id_type d=0; d<nndof; ++d)
			{
				myfile << sol_list[n*nndof + d];
				if(d!=(nndof-1)) myfile << " ";
			}
			myfile << "\n";
		}
		for(id_type part=1; part<n_ranks; ++part)
		{
			for(id_type n=0; n<(remote_sol_lists[part-1].size()/nndof); ++n)
			{
				for(id_type d=0; d<nndof; ++d)
				{
					myfile << remote_sol_lists[part-1][n*nndof + d];
					if(d!=(nndof-1)) myfile << " ";
				}
				myfile << "\n";
			}
		}

		// Print out external load information
		if(_prob->get_classification()==STRUCTURAL)
			myfile << "VECTORS external_force float " << "\n";
		else
		{
			myfile << "SCALARS external_load float " << nndof << "\n"; // I don't actually know if this is right
			myfile << "LOOKUP_TABLE default\n";
		}
		for(id_type n=0; n<(rhs_list.size()/nndof); ++n)
		{
			for(id_type d=0; d<nndof; ++d)
			{
				myfile << rhs_list[n*nndof + d];
				if(d!=(nndof-1)) myfile << " ";
			}
			myfile << "\n";
		}
		for(id_type part=1; part<n_ranks; ++part)
		{
			for(id_type n=0; n<(remote_rhs_lists[part-1].size()/nndof); ++n)
			{
				for(id_type d=0; d<nndof; ++d)
				{
					myfile << remote_rhs_lists[part-1][n*nndof + d];
					if(d!=(nndof-1)) myfile << " ";
				}
				myfile << "\n";
			}
		}

		// Print out stresses and strains
		if (_prob->get_classification() == STRUCTURAL)
		{
			id_type voigt = avg_strain.size();
			std::vector<std::string> voigt_map;
			if (voigt == 1)
				voigt_map = {"11"};
			else if (voigt == 3)
				voigt_map = {"11", "22", "12"};
			else if (voigt == 6)
				voigt_map = {"11", "22", "33", "23", "13", "12"};
			else
				err_message("Invalid number of average strains/stresses");

			// Output strains
			for (id_type v=0; v<voigt; ++v)
			{
				myfile << "SCALARS e" << voigt_map[v] << " float\n";
				myfile << "LOOKUP_TABLE default\n";
				for (id_type n=0; n<avg_strain[v].size(); ++n)
					myfile << avg_strain[v][n] << "\n";

				for (id_type part=1; part<n_ranks; ++part)
					for (id_type n=0; n<remote_strain[part-1][v].size(); ++n)
						myfile << remote_strain[part-1][v][n] << "\n";
			}

			// Output stresses
			for (id_type v=0; v<voigt; ++v)
			{
				myfile << "SCALARS s" << voigt_map[v] << " float\n";
				myfile << "LOOKUP_TABLE default\n";
				for (id_type n=0; n<avg_stress[v].size(); ++n)
					myfile << avg_stress[v][n] << "\n";

				for (id_type part=1; part<n_ranks; ++part)
					for (id_type n=0; n<remote_stress[part-1][v].size(); ++n)
						myfile << remote_stress[part-1][v][n] << "\n";
			}

			// Output von Mises stresses
			myfile << "SCALARS von_mises float\n";
			myfile << "LOOKUP_TABLE default\n";
			for (id_type n=0; n<avg_stress[0].size(); ++n)
			{
				std::vector<float> node_stress(voigt);
				for (id_type v=0; v<voigt; ++v)
					node_stress[v] = avg_stress[v][n];
				myfile << Utilities::von_mises(node_stress) << "\n";
			}

			for (id_type part=1; part<n_ranks; ++part)
				for (id_type n=0; n<remote_stress[part-1][0].size(); ++n)
				{
					std::vector<float> node_stress(voigt);
					for (id_type v=0; v<voigt; ++v)
						node_stress[v] = remote_stress[part-1][v][n];
					myfile << Utilities::von_mises(node_stress) << "\n";
				}
		}
	}
}




/*
 * Writes all of the element data to the VTK format. For serial ouput simply pass in empty remote lists
 */
void Writer_ASCII_VTK_Nodal::write_element_data(std::ofstream& myfile,
										  const std::vector<id_type>& elem_mats, const std::map<std::string, std::vector<float> >& internal_vars_to_avgs,
										  const std::vector<std::vector<id_type> >& remote_mats, const std::vector<std::vector<std::vector<float> > >& remote_int_vars)
{
	// Get the mesh
	Mesh* mesh = _prob->get_mesh();

	if (mesh->get_rank() == 0)
	{
		id_type n_ranks = mesh->n_ranks();
		id_type n_elem = elem_mats.size();
		for (id_type p=0; p<remote_mats.size(); ++p)
			n_elem += remote_mats[p].size();

		// Print out the element materials
		// Print out the cell types data
		myfile << "\nCELL_DATA " << n_elem << "\n";
		myfile << "SCALARS materials int\n";
		myfile << "LOOKUP_TABLE default\n";
		for(id_type e=0; e<elem_mats.size(); ++e)
			myfile << elem_mats[e] << "\n";
		for(id_type part=1; part<n_ranks; ++part)
		{
			for(id_type e=0; e<remote_mats[part-1].size(); ++e)
				myfile << remote_mats[part-1][e] << "\n";
		}

		// Print out partitioning information
		myfile << "\nSCALARS partitions int\n";
		myfile << "LOOKUP_TABLE default\n";
		for(id_type e=0; e<elem_mats.size(); ++e)
			myfile << 0 << "\n";
		for(id_type part=1; part<n_ranks; ++part)
			for(id_type e=0; e<remote_mats[part-1].size(); ++e)
				myfile << part << "\n";

		// Print out the internal variable information
		id_type iv = 0;
		for (auto it=internal_vars_to_avgs.begin(), end=internal_vars_to_avgs.end(); it!=end; ++it)
		{
			myfile << "\nSCALARS " << it->first << " float\n";
			myfile << "LOOKUP_TABLE default\n";
			const std::vector<float>& elem_avgs = it->second; // Get a reference
			for (id_type e=0; e<elem_avgs.size(); ++e)
				myfile << elem_avgs[e] << "\n";

			for(id_type part=1; part<n_ranks; ++part)
				for(id_type e=0; e<remote_int_vars[part-1][iv].size(); ++e)
					myfile << remote_int_vars[part-1][iv][e] << "\n";
			iv++;
		}
	}
}