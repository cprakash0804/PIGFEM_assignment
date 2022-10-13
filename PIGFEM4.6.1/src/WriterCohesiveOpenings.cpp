/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated April 2017

##################################################################################
*/
#include "WriterCohesiveOpenings.h"
#include "Problem.h"
#include "Mesh.h"
#include "elem.h"
#include "CohesiveElem.h"
#include "Inclusion.h"
#include "NodalData.h"
#include "Utilities.h"

/*
 * The Big 3
 * Constructor need the Problem
 * Destructor doesn't do anything
 * Copy Constructor copies the problem pointer
 */
WriterCohesiveOpenings::WriterCohesiveOpenings(Problem* prob)
	: Writer(prob)
{}
WriterCohesiveOpenings::WriterCohesiveOpenings()
{}
WriterCohesiveOpenings::~WriterCohesiveOpenings()
{}
WriterCohesiveOpenings::WriterCohesiveOpenings(const WriterCohesiveOpenings& other)
{
	_prob = other.get_prob();
	store_mesh();
}






/*
 * Main function used to write problem to any inherited file type
 */
void WriterCohesiveOpenings::write(std::string filename, double curr_t)
{
	std::ofstream myfile;
	if (_prob->get_mesh()->get_rank() == 0)
	{
		myfile.open(filename.c_str(), std::ofstream::out);
		if (!myfile.good())
		{
			std::string out = "Error opening the output file " + filename;
			err_message( out.data() );
		}
		myfile.precision(16);
	}
	
	writeFromStream(myfile, curr_t);

	if (_prob->get_mesh()->get_rank() == 0)
		myfile.close();
}




/*
 * The actual function that will do all of the writing (maybe by calling other functions)
 */
void WriterCohesiveOpenings::writeFromStream(std::ofstream& myfile, double curr_t)
{
	Mesh* mesh = _prob->get_mesh();
	if (mesh->dim() != 2)
		err_message("Cohesive opening output is only defined for 2D!");
	
	id_type n_inc = mesh->n_inclusions();
	std::vector<std::vector<double> > centers(n_inc);
	for (id_type i=0; i<n_inc; ++i)
		centers[i] = mesh->get_inclusion(i)->get_vec_parameter("center");

	std::vector<std::pair<id_type, std::pair<double, std::vector<double> > > > data; // structure: (inc, (theta, <delta>))
	for (Mesh::element_iterator it=mesh->active_elements_begin(), end=mesh->active_elements_end(); it!=end; ++it)
	{
		Elem* el = (*it);
		if (el->is_intersected())
		{
			id_type l_elem = mesh->global_to_local_elem(el->get_id());
			id_type curr_qp = mesh->n_volumetric_qps(l_elem);

			// Build the elemental solution vector
			std::vector<double> elem_sol;
			ProblemUtilities::formElementalVector(elem_sol, el, _prob->get_solution());

			// Get the cohesive element structure (Needed for mapping the cohesive matrix back to the elemental matrix)
			std::vector<std::vector<short_id_type> >& coh_elem_struct = el->getCohesiveElemStructure();

			// Loop over the cohesive elements
			for(id_type ce=0; ce<el->n_cohesive_elem(); ++ce)
			{
				// Get a pointer to the current cohesive element
				CohesiveElem* coh_el = el->get_cohesive_elem(ce);

				// Get the normal and tangential deltac
				double dcn, dct;
				Material* mat = el->get_cohesive_material(ce);
				if (mat->get_type()==OP_COHESIVE || mat->get_type()==OP_COHESIVE_NO_UNLOADING)
				{
					dcn = mat->get_parameter("DELTA_C");
					dct = dcn;
				}
				else if (mat->get_type()==XN_COHESIVE || mat->get_type()==XN_COHESIVE_NO_UNLOADING)
				{
					dcn = mat->get_parameter("DELTA_CN");
					dct = mat->get_parameter("DELTA_CT");
				}
				else
					err_message("Unknown cohesive material!");

				// Get the nodal solution object for the current cohesive element
				std::vector<double> coh_U_curr(coh_el->n_nodes() * 2); // The vectorized form of the previous structure
				for(id_type n=0; n<coh_el->n_nodes(); ++n)
				{
					id_type node = coh_elem_struct[ce][n];
					for(id_type d=0; d<2; ++d)
						coh_U_curr[n*2 + d] = elem_sol[node*2 + d];
				}

				std::vector<double> N = {0.5, 0.5}; // Assuming here that the cohesive elem is a CohEdge2 and we're just looking at the middle

				// Compute the theta value
				id_type curr_inc = el->get_inclusion_from_enrich(coh_elem_struct[ce][0] - el->n_nodes());
				std::vector<double> coords(2, 0.0);
				for (id_type d=0; d<2; ++d)
					for (id_type n=0; n<2; ++n)
						coords[d] += N[n] * (*(mesh->get_node_global(coh_el->get_node(n)->get_id())))(d);
				double theta = atan2(coords[1]-centers[curr_inc][1], coords[0]-centers[curr_inc][0]);

				// Compute the delta vector
				std::vector<double> delta(2, 0.0);
				ProblemUtilities::cohesiveDeltaFast(delta, N, coh_U_curr, mesh->get_rot_matrix(l_elem,curr_qp));
				delta[0] /= dct;
				delta[1] /= dcn;
				std::pair<double, std::vector<double> > temp = std::make_pair(theta, delta);
				std::pair<id_type, std::pair<double, std::vector<double> > > temp2 = std::make_pair(curr_inc, temp);
				data.push_back(temp2);

				id_type nqp = coh_el->n_q_points();
				curr_qp += nqp;
			}
		}
	}

	// if this is in parallel send the infor to proc 0
	if (!mesh->serial())
	{
		if (mesh->get_rank() != 0)
		{
			std::vector<id_type> incs(data.size());
			std::vector<double> delta_data(3*data.size());
			for (id_type i=0; i<data.size(); ++i)
			{
				incs[i] = data[i].first;
				delta_data[3*i] = data[i].second.first;
				delta_data[3*i+1] = data[i].second.second[0];
				delta_data[3*i+2] = data[i].second.second[1];
			}
			MPI_Send(incs.data(), incs.size(), MPI_ID, 0, 0, mesh->get_comm());
			MPI_Send(delta_data.data(), delta_data.size(), MPI_DOUBLE, 0, 1, mesh->get_comm());
		}
		else
		{
			for (id_type p=1; p<mesh->n_ranks(); ++p)
			{
				std::vector<id_type> incs;		Utilities::RecieveUnknown(incs, p, 0, MPI_ID, mesh->get_comm());
				std::vector<double> delta_data;	Utilities::RecieveUnknown(delta_data, p, 1, MPI_DOUBLE, mesh->get_comm());
				if ((incs.size() * 3) != delta_data.size())
					err_message("Mismatch between data found in the cohesive openings output!");
				for (id_type i=0; i<incs.size(); ++i)
				{
					std::vector<double> delta = {delta_data[i*3+1], delta_data[i*3+2]};
					std::pair<double, std::vector<double> > temp = std::make_pair(delta_data[3*i], delta);
					std::pair<id_type, std::pair<double, std::vector<double> > > temp2 = std::make_pair(incs[i], temp);
					data.push_back(temp2);
				}
			}
		}
	}

	// Sort the data into a good order
	std::sort(data.begin(), data.end());

	// Print the data
	if (mesh->serial() || (!mesh->serial() && mesh->get_rank()==0))
		for (id_type i=0; i<data.size(); ++i)
			myfile << data[i].first << "," << data[i].second.first << "," << data[i].second.second[0] << "," << data[i].second.second[1] << "\n";
}