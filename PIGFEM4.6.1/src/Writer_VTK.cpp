/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#include "Writer_VTK.h"
#include "Problem.h"
#include "Mesh.h"
#include "InternalVars.h"
#include "NodalData.h"
#include "Assembler.h"
#include "SubscaleModel.h"
#include "BodyLoad.h"
#include "Utilities.h"
#include "CohesiveElem.h"
#include "elem.h"
#include "mpi.h"


std::map<elem_type, unsigned char> Writer_VTK::elem_type_to_VTK = Writer_VTK::create_map();
std::map<coh_elem_type, elem_type> Writer_VTK::cohesive_elem_to_plot_elem = Writer_VTK::create_cohesive_type_map();



/*
 * The Big 3
 * Constructor need the Problem
 * Destructor doesn't do anything
 * Copy Constructor copies the problem pointer
 */
Writer_VTK::Writer_VTK(Problem* prob)
	: Writer(prob), _output_cohesive(false)
{}
Writer_VTK::Writer_VTK()
{}
Writer_VTK::~Writer_VTK()
{}
Writer_VTK::Writer_VTK(const Writer_VTK& other)
{
	_prob = other.get_prob();
	store_mesh();
}







void Writer_VTK::setParameter(std::string param, bool val)
{
	std::transform(param.begin(), param.end(), param.begin(), ::toupper); // Capitilize the name
	if (param=="COH" || param=="COHESIVE")
		_output_cohesive = val;
	else
		err_message("Invalid parameter name!");
}
bool Writer_VTK::getBoolParameter(std::string param)
{
	std::transform(param.begin(), param.end(), param.begin(), ::toupper); // Capitilize the name
	if (param=="COH" || param=="COHESIVE")
		return _output_cohesive;
	else
		err_message("Invalid parameter name!");
}












/*
 * Function to compute the associated cohesive damage for a given cohesive element
 */
float Writer_VTK::compute_cohesive_damage(CohesiveElem* coh_el, Material* mat,
									  id_type l_elem, id_type initial_qp)
{
	Mesh* mesh = _prob->get_mesh();

	// Assemble the cohesive element solution vector
	id_type nndof = _prob->nndof();
	std::vector<double> elem_sol(coh_el->n_nodes() * nndof);
	for (id_type n=0; n<coh_el->n_nodes(); ++n)
	{
		id_type id = coh_el->get_node(n)->get_id();
		for (id_type d=0; d<nndof; ++d)
			elem_sol[d + n*nndof] = _prob->get_solution()->get_value_global(id, d);
	}

	// Loop over the quadrature points
	id_type nqp = coh_el->n_q_points();
	std::vector<double> avg_delta(_prob->get_mesh()->dim());
	DenseMatrix<double> Nmat, Bcr;
	for(id_type qp=0; qp<nqp; ++qp)
	{
		// Get the quadrature point in local coordinates ( I hate having to do this but storing this means storing the quadrature points for all of the other elements as weel)
		// Maybe I could jus make a map for the cohesive elements?
		// I really just need to rethink how to do the cohesive stuff
		std::vector<double> rcoords;
		double W;
		coh_el->q_point(rcoords, W, qp);

		//Assemble the N matrix to compute the cohesive opening
		ProblemUtilities::Assemble_Cohesive_N_Mat(Nmat, mesh->get_shape(l_elem, initial_qp+qp), avg_delta.size());

		// Get the rotation matrix and differential area
		double dA;
		DenseMatrix<double> Rot_xyz_to_ntt;
		coh_el->compute_rotation(rcoords, Rot_xyz_to_ntt, dA);
		//coh_el->compute_deformed_rotation(rcoords, elem_sol, Rot_xyz_to_ntt, dA);

		// Compute opening in the rotated coordinate frame
		Bcr = Rot_xyz_to_ntt * Nmat;
		std::vector<double> delta = Bcr * elem_sol;
		for (id_type i=0; i<delta.size(); ++i)
			avg_delta[i] += delta[i];
	}
	for (id_type i=0; i<avg_delta.size(); ++i)
		avg_delta[i] /= nqp;


	id_type zone = mat->get_cohesive_damage_zone(avg_delta);
	if (zone == 0) // Compressive
		return 0.2;
	else if (zone == 1) // Increasing Curve
		return 0.4;
	else if (zone == 2)	// Decreasing Curve
		return 0.6;
	else if (zone == 3)	// Totally Failed
		return 0.8;
	else
		err_message("Unknown cohesive damage zone");
}
















/*
 * Stores the mesh in this Writer's internal storage
 * Can be recalled if the mesh is adapted
 */
void Writer_VTK::store_mesh()
{
	// Clear any old mesh
	_node_id_list.clear();
	_node_coord_list.clear();
	_eptr.clear();
	_eind.clear();
	_elem_types.clear();
	_elem_mats.clear();

	// Store the node information
	store_nodes(_node_id_list, _node_coord_list);

	// Store the element information
	store_elements(_eptr, _eind, _elem_types, _elem_mats);
}




/*
 * Stores the node information in the private data members used to store nodes
 */
void Writer_VTK::store_nodes(std::vector<id_type>& id_list, std::vector<float>& coords_list)
{
	Mesh* mesh = _prob->get_mesh();
	id_list.resize(mesh->n_local_owned_nodes()+mesh->n_local_owned_enrich_nodes());
	coords_list.resize(id_list.size() * 3);
	id_type curr_node = 0;
	for(Mesh::node_iterator it=mesh->nodes_begin(), end=mesh->nodes_end(); it!=end; ++it)
	{
		id_type l_node = mesh->global_to_local_node((*it)->get_id());
		if(mesh->own_node_local(l_node))
		{
			id_list[curr_node] = (*it)->get_id();
			coords_list[3*curr_node] = (*(*it))(0);   coords_list[3*curr_node+1] = (*(*it))(1);   coords_list[3*curr_node+2] = (*(*it))(2);
			curr_node++;
		}
	}
	for(Mesh::enrich_node_iterator it=mesh->enrich_nodes_begin(), end=mesh->enrich_nodes_end(); it!=end; ++it)
	{
		id_type l_node = mesh->global_to_local_node((*it)->get_id());
		if(mesh->own_node_local(l_node))
		{
			id_list[curr_node] = (*it)->get_id();
			coords_list[3*curr_node] = (*(*it))(0);   coords_list[3*curr_node+1] = (*(*it))(1);   coords_list[3*curr_node+2] = (*(*it))(2);
			curr_node++;
		}
	}
}




/*
 * Stores the element information in the private data members
 */
void Writer_VTK::store_elements(std::vector<id_type>& eptr, std::vector<id_type>& eind,
								std::vector<unsigned char>& types, std::vector<id_type>& mats)
{
	Mesh* mesh = _prob->get_mesh();
	// Get the total number of elements (each integration element counts as 1 element)
	_n_elem = 0;
	_n_points = 0;
	for(Mesh::element_iterator it=mesh->active_elements_begin(), end=mesh->active_elements_end(); it!=end; ++it)
	{
		if((*it)->is_intersected())
		{
			for(id_type ie=0; ie<(*it)->n_integration_elem(); ++ie)
			{
				_n_elem++;
				_n_points += (*it)->get_integration_elem(ie)->n_nodes();
			}
			// If I also want to output the cohesive elements, store them here
			if (_output_cohesive && mesh->is_cohesive())
			{
				for (id_type ce=0; ce<(*it)->n_cohesive_elem(); ++ce)
				{
					_n_elem++;
					_n_points += (*it)->get_cohesive_elem(ce)->n_nodes();
				}
			}
		}
		else
		{
			_n_elem++;
			_n_points += (*it)->n_nodes();
		}
	}	
	eptr.resize(_n_elem+1); eptr[0] = 0;
	eind.resize(_n_points);
	types.resize(_n_elem); // Vector that contains the element to VTK type mapping
	mats.resize(_n_elem);

	id_type curr_ind = 0;
	id_type curr_elem = 0;
	std::vector<Material*>& materials = mesh->get_materials();
	for(Mesh::element_iterator it=mesh->active_elements_begin(), end=mesh->active_elements_end(); it!=end; ++it)
	{
		if((*it)->is_intersected())
		{
			for(id_type ie=0; ie<(*it)->n_integration_elem(); ++ie)
			{
				Elem* int_el = (*it)->get_integration_elem(ie);

				// Get the material and type information for this integration element
				Material* mat = (*it)->get_int_elem_mat(ie);
				mats[curr_elem] = std::find(materials.begin(), materials.end(), mat) - materials.begin();
				types[curr_elem] = elem_type_to_VTK[int_el->get_type()];

				// Print out elem_node information for this integration element
				for(id_type n=0; n<int_el->n_nodes(); ++n)
				{
					eind[curr_ind] = int_el->get_node(n)->get_id();
					curr_ind++;
				}

				// Update the current element
				curr_elem++;
				// Set the node ptr value
				eptr[curr_elem] = curr_ind;
			}

			// If I also want to output the cohesive elements, store them here
			if (_output_cohesive && mesh->is_cohesive())
			{
				for (id_type ce=0; ce<(*it)->n_cohesive_elem(); ++ce)
				{
					CohesiveElem* coh_el = (*it)->get_cohesive_elem(ce);

					// Get the material and type information for this cohesive element
					Material* mat = (*it)->get_cohesive_material(ce);
					mats[curr_elem] = std::find(materials.begin(), materials.end(), mat) - materials.begin();
					types[curr_elem] = elem_type_to_VTK[ cohesive_elem_to_plot_elem[coh_el->get_type()] ];

					// Print out elem_node information for this integration element
					std::vector<id_type> ids = coh_el->plot_elem_ids(); // Special function because cohesive elements are weird
					for(id_type n=0; n<ids.size(); ++n)
					{
						eind[curr_ind] = ids[n];
						curr_ind++;
					}

					// Update the current element
					curr_elem++;
					// Set the node ptr value
					eptr[curr_elem] = curr_ind;
				}
			}
		}
		else
		{
			id_type l_id = mesh->global_to_local_elem((*it)->get_id());

			// Get the material and type information for this element
			id_type mat_num = mesh->get_mat_num_local(l_id);
			mats[curr_elem] = mat_num;
			types[curr_elem] = elem_type_to_VTK[(*it)->get_type()];

			// Print out elem_node information
			for(id_type n=0; n<(*it)->n_nodes(); ++n)
			{
				eind[curr_ind] = (*it)->get_node(n)->get_id();
				curr_ind++;
			}

			// Update the current element
			curr_elem++;
			// Set the node ptr value
			eptr[curr_elem] = curr_ind;
		}
	}
}




/*
 * Main function used to write problem to any inherited file type
 */
void Writer_VTK::write(std::string filename, double curr_t)
{
	// Open the file
	std::ofstream myfile;
	open(myfile, filename); // Open in binary or ascii depending on the inherited object

	// Write everything from the general function
	writeFromStream(myfile, curr_t);

	myfile.close();
}



/*
 * The actual function that will do all of the writing (maybe by calling other functions)
 */
void Writer_VTK::writeFromStream(std::ofstream& myfile, double curr_t)
{
	// Write the VTK header
	write_header(myfile, curr_t);

	// Write the mesh (nodes and elements)
	write_mesh(myfile);

	// Gather and write data associated with the nodes (solution, external load)
	write_node_info(myfile);

	// Gather and write data associated with the elements (stress, strain, internal variables, material)
	write_element_info(myfile);
}




/*
 * Writes the nodes and element information
 */
void Writer_VTK::write_mesh(std::ofstream& myfile)
{
	// This vector maps from the local id numbering to the order the nodes are output to file
	std::vector<id_type> local_to_VTK;

	// Get the information from other mesh partitions if this in parallel
	std::vector<std::vector<id_type> > remote_node_ids, remote_eptr, remote_eind;
	std::vector<std::vector<unsigned char> > remote_etype;
	std::vector<std::vector<float> > remote_coords; 
	communicate_mesh(remote_node_ids, remote_coords, remote_eptr, remote_eind, remote_etype);

	// Write all of the nodes
	write_nodes(myfile, local_to_VTK,
				_node_id_list, _node_coord_list,
				remote_node_ids, remote_coords);

	// Write all of the elements
	write_elements(myfile, local_to_VTK,
				   _eptr, _eind, _elem_types,
				   remote_eptr, remote_eind, remote_etype);
}




/*
 * Gather and write data associated with the nodes (solution, external load)
 */
void Writer_VTK::write_node_info(std::ofstream& myfile)
{
	// Gather all of the data from the nodes
	std::vector<float> sol_list, rhs_list;
	gather_nodal_data(sol_list, rhs_list);

	// Gather info from other processes to process 0
	std::vector<std::vector<float> > remote_sol, remote_rhs;
	communicate_node_data(sol_list, rhs_list,
						  remote_sol, remote_rhs);

	// Write all of the nodal data to the VTK format
	write_nodal_data(myfile,
					 sol_list, rhs_list,
					 remote_sol, remote_rhs);
}




/*
 * Gather and write data associated with the elements (stress, strain, internal variables, material)
 */
void Writer_VTK::write_element_info(std::ofstream& myfile)
{
	// Gather all of the data from the nodes
	std::vector<std::vector<float> > avg_strain, avg_stress;
	std::map<std::string, std::vector<float> > internal_vars_to_avgs;
	gather_element_data(internal_vars_to_avgs, avg_strain, avg_stress);

	// Gather info from other processes to process 0
	std::vector<std::vector<std::vector<float> > > remote_int_vars, remote_strain, remote_stress;
	std::vector<std::vector<id_type> > remote_mats;
	communicate_element_data(_elem_mats, internal_vars_to_avgs, avg_strain, avg_stress,
							 remote_mats, remote_int_vars, remote_strain, remote_stress);

	// Writes all of the elemental data to the VTK format. For serial ouput simply pass in empty remote list
	write_element_data(myfile,
					   _elem_mats, internal_vars_to_avgs, avg_strain, avg_stress,
					   remote_mats, remote_int_vars, remote_strain, remote_stress);
}




/*
 * Function that communicates all of the existing mesh data to proc 0
 */
void Writer_VTK::communicate_mesh(std::vector<std::vector<id_type> >& remote_node_ids, std::vector<std::vector<float> >& remote_coords,
										std::vector<std::vector<id_type> >& remote_eptr, std::vector<std::vector<id_type> >& remote_eind, std::vector<std::vector<unsigned char> >& remote_etype)
{
	Mesh* mesh = _prob->get_mesh();

	// No need to communicate anything if this is serial
	if (mesh->serial())
		return;

	else
	{
		// Not proc 0, send all of my info
		if (mesh->get_rank() != 0)
		{
			MPI_Send(_node_id_list.data(), _node_id_list.size(), MPI_ID, 0, 0, mesh->get_comm());
			MPI_Send(_node_coord_list.data(), _node_coord_list.size(), MPI_FLOAT, 0, 1, mesh->get_comm());
			MPI_Send(_eptr.data(), _eptr.size(), MPI_ID, 0, 2, mesh->get_comm());
			MPI_Send(_eind.data(), _eind.size(), MPI_ID, 0, 3, mesh->get_comm());
			MPI_Send(_elem_types.data(), _elem_types.size(), MPI_UNSIGNED_CHAR, 0, 4, mesh->get_comm());
		}

		// Proc 0, recieve all info
		else
		{
			int nproc = mesh->n_ranks();
			remote_node_ids.resize(nproc-1);
			remote_coords.resize(nproc-1);
			remote_eptr.resize(nproc-1);
			remote_eind.resize(nproc-1);
			remote_etype.resize(nproc-1);
			for (int part=1; part<nproc; ++part)
			{
				Utilities::RecieveUnknown(remote_node_ids[part-1], part, 0, MPI_ID, mesh->get_comm());
				Utilities::RecieveUnknown(remote_coords[part-1], part, 1, MPI_FLOAT, mesh->get_comm());
				Utilities::RecieveUnknown(remote_eptr[part-1], part, 2, MPI_ID, mesh->get_comm());
				Utilities::RecieveUnknown(remote_eind[part-1], part, 3, MPI_ID, mesh->get_comm());
				Utilities::RecieveUnknown(remote_etype[part-1], part, 4, MPI_UNSIGNED_CHAR, mesh->get_comm());
			}
		}
	}
}




/*
 * Communicated solution and rhs data to proc 0
 */
void Writer_VTK::communicate_node_data(const std::vector<float>& sol_list, const std::vector<float>& rhs_list,
											std::vector<std::vector<float> >& remote_sol, std::vector<std::vector<float> >& remote_rhs)
{
	Mesh* mesh = _prob->get_mesh();

	// No need to communicate anything if this is serial
	if (mesh->serial())
		return;

	else
	{
		// Not proc 0, send all of my info
		if (mesh->get_rank() != 0)
		{
			MPI_Send(sol_list.data(), sol_list.size(), MPI_FLOAT, 0, 0, mesh->get_comm());
			MPI_Send(rhs_list.data(), rhs_list.size(), MPI_FLOAT, 0, 1, mesh->get_comm());
		}

		// Proc 0, recieve all info
		else
		{
			int nproc = mesh->n_ranks();
			remote_sol.resize(nproc-1);
			remote_rhs.resize(nproc-1);
			for (int part=1; part<nproc; ++part)
			{
				Utilities::RecieveUnknown(remote_sol[part-1], part, 0, MPI_FLOAT, mesh->get_comm());
				Utilities::RecieveUnknown(remote_rhs[part-1], part, 1, MPI_FLOAT, mesh->get_comm());
			}
		}
	}
}




/*
 * Communicate, materials, internal variables, stresses and strains to proc 0
 */
void Writer_VTK::communicate_element_data(const std::vector<id_type>& _elem_mats, const std::map<std::string, std::vector<float> >& internal_vars_to_avgs,
												const std::vector<std::vector<float> >& avg_strain, const std::vector<std::vector<float> >& avg_stress,
												std::vector<std::vector<id_type> >& remote_mats, std::vector<std::vector<std::vector<float> > >& remote_int_vars,
												std::vector<std::vector<std::vector<float> > >& remote_strain, std::vector<std::vector<std::vector<float> > >& remote_stress)
{
	Mesh* mesh = _prob->get_mesh();

	// No need to communicate anything if this is serial
	if (mesh->serial())
		return;

	else
	{
		// Not proc 0, send all of my info
		if (mesh->get_rank() != 0)
		{
			MPI_Send(_elem_mats.data(), _elem_mats.size(), MPI_ID, 0, 0, mesh->get_comm());

			// Send internal variable info
			int tag = 1;
			for (auto it=internal_vars_to_avgs.begin(), end=internal_vars_to_avgs.end(); it!=end; ++it)
			{
				MPI_Send(it->second.data(), it->second.size(), MPI_FLOAT, 0, tag, mesh->get_comm());
				tag++;
			}

			// If its structural, send stress and strain info
			if (_prob->get_classification() == STRUCTURAL)
			{
				std::vector<float> strains, stresses;
				for (id_type e=0; e<avg_strain.size(); ++e)
				{
					strains.insert(strains.end(), avg_strain[e].begin(), avg_strain[e].end());
					stresses.insert(stresses.end(), avg_stress[e].begin(), avg_stress[e].end());
				}
				MPI_Send(strains.data(), strains.size(), MPI_FLOAT, 0, tag, mesh->get_comm());
				MPI_Send(stresses.data(), stresses.size(), MPI_FLOAT, 0, tag+1, mesh->get_comm());
			}
		}

		// Proc 0, recieve all info
		else
		{
			int nproc = mesh->n_ranks();
			remote_mats.resize(nproc-1);
			remote_int_vars.resize(nproc-1);
			remote_strain.resize(nproc-1);
			remote_stress.resize(nproc-1);
			for (int part=1; part<nproc; ++part)
			{
				Utilities::RecieveUnknown(remote_mats[part-1], part, 0, MPI_ID, mesh->get_comm());

				// Recieve internal variable info
				int tag = 1;
				remote_int_vars[part-1].resize(internal_vars_to_avgs.size());
				for (id_type iv=0; iv<internal_vars_to_avgs.size(); ++iv)
				{
					Utilities::RecieveUnknown(remote_int_vars[part-1][iv], part, tag, MPI_FLOAT, mesh->get_comm());
					tag++;
				}

				// Recieve stresses and strains
				if(_prob->get_classification() == STRUCTURAL)
				{
					std::vector<float> strains, stresses;
					Utilities::RecieveUnknown(strains, part, tag, MPI_FLOAT, mesh->get_comm());
					Utilities::RecieveUnknown(stresses, part, tag+1, MPI_FLOAT, mesh->get_comm());
					if (strains.size() != stresses.size())
						err_message("Invalid number of stress and strain components recieved.");

					id_type nelem = remote_mats[part-1].size();
					int stride = strains.size() / nelem;
					remote_strain[part-1].resize(nelem);
					remote_stress[part-1].resize(nelem);
					for (id_type e=0; e<nelem; ++e)
					{
						id_type start = e*stride;
						id_type end = (e+1)*stride;
						remote_strain[part-1][e] = std::vector<float>(strains.begin()+start, strains.begin()+end);
						remote_stress[part-1][e] = std::vector<float>(stresses.begin()+start, stresses.begin()+end);
					}
				}
			}
		}
	}
}





/*
 * Gathers info about the solution and rhs in a serialized form
 */
void Writer_VTK::gather_nodal_data(std::vector<float>& sol_list, std::vector<float>& rhs_list)
{
	// Get the mesh
	Mesh* mesh = _prob->get_mesh();

	// Convert the enrichment solution to the displacement solution at the enrichment nodes
	NodalData enrich_sol_plot(_prob->get_mesh());
	enrich_sol_plot.preallocate_storage(_prob->nndof());
	NodalData enrich_load_plot(_prob->get_mesh());
	enrich_load_plot.preallocate_storage(_prob->nndof());
	if(mesh->IGFEM())
		_prob->interpolate_nodal_enrichment(enrich_sol_plot, enrich_load_plot);

	// Get the macroscopic strain if this is a subscale simulations
	std::vector<double> curr_macro_strain;
	if (_prob->subscale())
	{
		if (_prob->linear())
			curr_macro_strain = (*(_prob->body_loads_begin()))->get_vec_parameter("macro strain"); // NOTE: Assume that the macro strain is the first and only body load here which it should be
		else
			curr_macro_strain = _prob->get_assembler()->get_vec_parameter("current strain");
	}

	id_type nndof = _prob->nndof();
	int nndof_write = (_prob->get_classification()==STRUCTURAL) ? 3 : nndof; // Set the number of entries for each solution component
	sol_list.resize(_node_id_list.size() * nndof_write);
	rhs_list.resize(_node_id_list.size() * nndof_write);
	id_type curr_node = 0;
	for(Mesh::node_iterator it=mesh->nodes_begin(), end=mesh->nodes_end(); it!=end; ++it)
	{
		id_type l_node = mesh->global_to_local_node((*it)->get_id());
		if(mesh->own_node_local(l_node))
		{
			for(id_type d=0; d<nndof; ++d)
			{
				sol_list[curr_node*nndof_write+d] = _prob->get_solution()->get_value_local(l_node, d);
				rhs_list[curr_node*nndof_write+d] = _prob->get_external_load()->get_value_local(l_node, d);
			}
			if (_prob->get_classification()==STRUCTURAL) // Doing this here saves a lot of copying in the actual output functions (for binary)
				for (id_type d=nndof; d<3; ++d)
				{
					sol_list[curr_node*nndof_write+d] = 0.0;
					rhs_list[curr_node*nndof_write+d] = 0.0;
				}
			if (_prob->subscale())
			{
				std::vector<double> coords = (*it)->get_coords();
				std::vector<double> macro_disp = _prob->get_subscale_model()->compute_macro_displacement(curr_macro_strain, coords);
				for(id_type d=0; d<nndof; ++d)
					sol_list[curr_node*nndof_write+d] += macro_disp[d];
			}
			curr_node++;
		}
	}
	if(mesh->IGFEM())
	{
		for(Mesh::enrich_node_iterator it=mesh->enrich_nodes_begin(), end=mesh->enrich_nodes_end(); it!=end; ++it)
		{
			id_type l_node = mesh->global_to_local_node((*it)->get_id());
			if(mesh->own_node_local(l_node))
			{
				for(id_type d=0; d<nndof; ++d)
				{
					sol_list[curr_node*nndof_write+d] = enrich_sol_plot.get_value_local(l_node, d);
					rhs_list[curr_node*nndof_write+d] = enrich_load_plot.get_value_local(l_node, d);
				}
				if (_prob->get_classification()==STRUCTURAL)
					for (id_type d=nndof; d<3; ++d)
					{
						sol_list[curr_node*nndof_write+d] = 0.0;
						rhs_list[curr_node*nndof_write+d] = 0.0;
					}
				if (_prob->subscale())
				{
					std::vector<double> coords = (*it)->get_coords();
					std::vector<double> macro_disp = _prob->get_subscale_model()->compute_macro_displacement(curr_macro_strain, coords);
					for(id_type d=0; d<nndof; ++d)
						sol_list[curr_node*nndof_write+d] += macro_disp[d];
				}
				curr_node++;
			}
		}
	}
}




/*
 * Gathers info about the internal variables, stress, and strain in a serialized form
 */
void Writer_VTK::gather_element_data(std::map<std::string, std::vector<float> >& internal_vars_to_avgs,
										   std::vector<std::vector<float> >& avg_strain, std::vector<std::vector<float> >& avg_stress)
{
	// Get the mesh
	Mesh * mesh = _prob->get_mesh();

	// Create map data structures of all the internal variables in the mesh that we're actually going to print
	std::vector<Material*>& materials = mesh->get_materials();
	std::map<std::string, float> print_internal_vars_to_init;
	for (id_type m=0; m<materials.size(); ++m)
	{
		std::vector<std::string> names = materials[m]->internal_vars_name();
		std::vector<bool> print = materials[m]->internal_vars_print();
		std::vector<double> init = materials[m]->init_internal_vars();

		for (id_type iv=0; iv<names.size(); ++iv)
			if(print[iv])
				print_internal_vars_to_init[names[iv]] = init[iv];
	}	

	// If I'm ouputting cohesive elements and I don't already have a damage internal variable then I need to add it
	if (_output_cohesive)
		if (print_internal_vars_to_init.find("damage") == print_internal_vars_to_init.end())
			print_internal_vars_to_init.insert(std::pair<std::string, float>("damage", 0.0));

	// Set variable sizes appropriately
	for (auto it=print_internal_vars_to_init.begin(), end=print_internal_vars_to_init.end(); it!=end; ++it) // A vector of length n_elem for each internal variable that we are printing
		internal_vars_to_avgs.insert( std::pair<std::string, std::vector<float> >(it->first, std::vector<float>(_n_elem)) );
	if (_prob->get_classification() == STRUCTURAL)
	{
		avg_strain.resize(_n_elem);
		avg_stress.resize(_n_elem);
	}

	// Loop over element and store the internal variables
	// If the problem is structural also stroe stress and strain
	id_type curr_elem = 0;
	for(Mesh::element_iterator it=mesh->active_elements_begin(), end=mesh->active_elements_end(); it!=end; ++it)
	{
		// Get the local element number
		id_type l_elem = mesh->global_to_local_elem((*it)->get_id());

		// If the problem is structural then I will fill in the stresses and strains
		std::vector<std::vector<double> > strains;
		std::vector<std::vector<double> > stresses;
		if (_prob->get_classification() == STRUCTURAL)
			_prob->compute_stress_strain(strains, stresses, (*it));
			

		if((*it)->is_intersected())
		{
			id_type initial_qp = 0;
			for(id_type ie=0; ie<(*it)->n_integration_elem(); ++ie)
			{
				// Store the averages of the internal variables
				// First store the initial values of all of them. Then we will overwrite the ones that actually exist in this element
				Material* mat = mesh->get_element_material_global((*it)->get_id(), ie);
				for (auto it2=print_internal_vars_to_init.begin(), end2=print_internal_vars_to_init.end(); it2!=end2; ++it2)
				{
					internal_vars_to_avgs[it2->first][curr_elem] = it2->second;
					if (it2->first == "damage")
						if (mesh->get_element_material_local(l_elem, ie)->get_name() != "Epoxy")	// FIXME: I really don't like this hardcoding
							internal_vars_to_avgs[it2->first][curr_elem] -= 0.2;
				}
				std::vector<std::string> names = mat->internal_vars_name();
				std::vector<bool> print = mat->internal_vars_print();
				id_type nqp = (*it)->get_integration_elem(ie)->n_q_points();
				for (id_type iv=0; iv<names.size(); ++iv)
				{
					if (print[iv])
					{
						internal_vars_to_avgs[names[iv]][curr_elem] = 0; // Set initial value to 0 so we can sum and then divide by the number of quadrature points
						for (id_type qp=0; qp<nqp; ++qp)
							internal_vars_to_avgs[names[iv]][curr_elem] += _prob->get_internal_vars()->get_internal_var_global((*it)->get_id(), initial_qp+qp, iv);
						internal_vars_to_avgs[names[iv]][curr_elem] /= nqp;
					}
				}

				// Compute the averages of the strains and stresses
				if (_prob->get_classification() == STRUCTURAL)
				{
					id_type voigt = strains[0].size();
					avg_strain[curr_elem].resize(voigt);
					avg_stress[curr_elem].resize(voigt);
					for (id_type v=0; v<voigt; ++v)
						for (id_type qp=0; qp<nqp; ++qp)
						{
							avg_strain[curr_elem][v] += strains[initial_qp+qp][v];
							avg_stress[curr_elem][v] += stresses[initial_qp+qp][v];
						}
					for (id_type v=0; v<voigt; ++v)
					{
						avg_strain[curr_elem][v] /= nqp;
						avg_stress[curr_elem][v] /= nqp;
					}
				}

				// Update the current initial quadrature point
				initial_qp += nqp;
				// Update the current element
				curr_elem++;
			}

			// If I also want to output the cohesive elements, store them here
			if (_output_cohesive && mesh->is_cohesive())
			{
				for (id_type ce=0; ce<(*it)->n_cohesive_elem(); ++ce)
				{
					CohesiveElem* coh_el = (*it)->get_cohesive_elem(ce);

					// Store the averages of the internal variables
					// First store the initial values of all of them. Then we will overwrite the ones that actually exist in this element (on;y damage)
					for (auto it2=print_internal_vars_to_init.begin(), end2=print_internal_vars_to_init.end(); it2!=end2; ++it2)
						internal_vars_to_avgs[it2->first][curr_elem] = it2->second;

					internal_vars_to_avgs["damage"][curr_elem] = compute_cohesive_damage(coh_el, (*it)->get_cohesive_material(ce), l_elem, initial_qp);

					// Compute the averages of the strains and stresses (Just set to 0 here...)
					if (_prob->get_classification() == STRUCTURAL)
					{
						id_type voigt = strains[0].size();
						avg_strain[curr_elem].resize(voigt);
						avg_stress[curr_elem].resize(voigt);
					}

					// Update the current initial quadrature point
					initial_qp += coh_el->n_q_points();
					// Update the current element
					curr_elem++;
				}
			}
		}
		else
		{
			// Store the averages of the internal variables
			// First store the initial values of all of them. Then we will overwrite the ones that actually exist in this element
			Material* mat = mesh->get_element_material_global((*it)->get_id());
			for (auto it2=print_internal_vars_to_init.begin(), end2=print_internal_vars_to_init.end(); it2!=end2; ++it2)
			{
				internal_vars_to_avgs[it2->first][curr_elem] = it2->second;
				if (it2->first == "damage")
					if (mesh->get_element_material_local(l_elem)->get_name() != "Epoxy")	// FIXME: I really don't like this hardcoding
						internal_vars_to_avgs[it2->first][curr_elem] -= 0.2;
			}
			std::vector<std::string> names = mat->internal_vars_name();
			std::vector<bool> print = mat->internal_vars_print();
			id_type nqp = (*it)->n_q_points();
			for (id_type iv=0; iv<names.size(); ++iv)
			{
				if (print[iv])
				{
					internal_vars_to_avgs[names[iv]][curr_elem] = 0; // Set initial value to 0 so we can sum and then divide by the number of quadrature points
					for (id_type qp=0; qp<nqp; ++qp)
						internal_vars_to_avgs[names[iv]][curr_elem] += _prob->get_internal_vars()->get_internal_var_global((*it)->get_id(), qp, iv);
					internal_vars_to_avgs[names[iv]][curr_elem] /= nqp;
				}
			}

			// Compute the averages of the strains and stresses
			if (_prob->get_classification() == STRUCTURAL)
			{
				id_type voigt = strains[0].size();
				avg_strain[curr_elem].resize(voigt);
				avg_stress[curr_elem].resize(voigt);
				for (id_type v=0; v<voigt; ++v)
					for (id_type qp=0; qp<nqp; ++qp)
					{
						avg_strain[curr_elem][v] += strains[qp][v];
						avg_stress[curr_elem][v] += stresses[qp][v];
					}
				for (id_type v=0; v<voigt; ++v)
				{
					avg_strain[curr_elem][v] /= nqp;
					avg_stress[curr_elem][v] /= nqp;
				}
			}

			// Update the current element
			curr_elem++;
		}
	}
}
