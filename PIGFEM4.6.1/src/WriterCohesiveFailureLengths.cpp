/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated April 2017

##################################################################################
*/
#include "WriterCohesiveFailureLengths.h"
#include "Problem.h"
#include "Mesh.h"
#include "elem.h"
#include "NodalData.h"
#include <numeric>

/*
 * The Big 3
 * Constructor need the Problem
 * Destructor doesn't do anything
 * Copy Constructor copies the problem pointer
 */
WriterCohesiveFailureLengths::WriterCohesiveFailureLengths(Problem* prob)
	: Writer(prob), _individual(false)
{}
WriterCohesiveFailureLengths::WriterCohesiveFailureLengths()
{}
WriterCohesiveFailureLengths::~WriterCohesiveFailureLengths()
{}
WriterCohesiveFailureLengths::WriterCohesiveFailureLengths(const WriterCohesiveFailureLengths& other)
{
	_prob = other.get_prob();
	store_mesh();
}





void WriterCohesiveFailureLengths::setParameter(std::string param, bool val)
{
	std::transform(param.begin(), param.end(), param.begin(), ::toupper); // Capitilize the name
	if (param=="INDIVIDUAL" || param=="EACH")
		_individual = val;
	else
		err_message("Invalid parameter name!");
}
bool WriterCohesiveFailureLengths::getBoolParameter(std::string param)
{
	std::transform(param.begin(), param.end(), param.begin(), ::toupper); // Capitilize the name
	if (param=="INDIVIDUAL" || param=="EACH")
		return _individual;
	else
		err_message("Invalid parameter name!");
}




/*
 * Stores the mesh in this Writer's internal storage
 * Can be recalled if the mesh is adapted
 */
void WriterCohesiveFailureLengths::store_mesh()
{
	_local_lengths.clear();
	Mesh* mesh = _prob->get_mesh();
	for (Mesh::element_iterator it=mesh->active_elements_begin(), end=mesh->active_elements_end(); it!=end; ++it)
	{
		Elem* el = *it;
		if (el->is_intersected())
		{
			for (id_type ce=0; ce<el->n_cohesive_elem(); ++ce)
			{
				std::vector<double> coords1 = el->get_cohesive_elem(ce)->get_node(0)->get_coords();
				std::vector<double> coords2 = el->get_cohesive_elem(ce)->get_node(1)->get_coords();
				std::vector<double> vec = {coords2[0]-coords1[0], coords2[1]-coords1[1]};
				double dist = sqrt(vec[0]*vec[0] + vec[1]*vec[1]);
				_local_lengths.push_back( dist );
			}
		}
	}
}








/*
 * Main function used to write problem to any inherited file type
 */
void WriterCohesiveFailureLengths::writeConsecutive(std::string filename, double curr_t, id_type step)
{
	std::ofstream myfile;
	if (_prob->get_mesh()->get_rank() == 0)
	{
		if (step==0)
			myfile.open(filename.c_str(), std::ofstream::out);
		else
			myfile.open(filename.c_str(), std::ofstream::app);
		if (!myfile.good())
		{
			std::string out = "Error opening the output file " + filename;
			err_message( out.data() );
		}
		myfile.precision(16);
	}
	
	// If this is the first step, output all of the lengths of the cohesive 
	if (step == 0)
	{
		// Store the local lengths
		store_mesh();

		// Write the lengths of each interface or just the total length
		writeLengths(myfile);
	}

	// Write all of the failed cohesive lengths for this step
	writeFromStream(myfile, curr_t);

	if (_prob->get_mesh()->get_rank() == 0)
		myfile.close();
}



/*
 * Used to write all of the cohesive langths in the very first step
 */
void WriterCohesiveFailureLengths::writeLengths(std::ofstream& myfile)
{
	Mesh* mesh = _prob->get_mesh();
	if (mesh->get_rank() == 0)
	{
		if (_individual)
		{
			std::vector<std::vector<double> > other_lengths(mesh->n_ranks() - 1);
			for (id_type p=1; p<mesh->n_ranks(); ++p)
				Utilities::RecieveUnknown(other_lengths[p-1], p, 0, MPI_DOUBLE, mesh->get_comm());

			for (id_type i=0; i<_local_lengths.size(); ++i)
				myfile << _local_lengths[i] << ",";
			for (id_type p=1; p<mesh->n_ranks(); ++p)
				for (id_type i=0; i<other_lengths[p-1].size(); ++i)
					myfile << other_lengths[p-1][i] << ",";
		}

		else
		{
			double local_total = std::accumulate(_local_lengths.begin(), _local_lengths.end(), 0.0);
			MPI_Status stat;
			double remote_data;
			for (id_type p=1; p<mesh->n_ranks(); ++p)
			{
				MPI_Recv(&remote_data, 1, MPI_DOUBLE, p, 0, _prob->get_mesh()->get_comm(), &stat);
				local_total += remote_data;
			}
			myfile << local_total;
		}

		myfile << "\n";
	}
	else
	{
		if (_individual)
			MPI_Send(_local_lengths.data(), _local_lengths.size(), MPI_DOUBLE, 0, 0, mesh->get_comm());
		else
		{
			double local_total = std::accumulate(_local_lengths.begin(), _local_lengths.end(), 0.0);
			MPI_Send(&local_total, 1, MPI_DOUBLE, 0, 0, mesh->get_comm());
		}
	}
}



/*
 * The actual function that will do all of the writing (maybe by calling other functions)
 */
void WriterCohesiveFailureLengths::writeFromStream(std::ofstream& myfile, double curr_t)
{
	Mesh* mesh = _prob->get_mesh();
	std::vector<double> local_deltas(_local_lengths.size(), 0.0);
	id_type nndof  = _prob->nndof();
	id_type count = 0;
	for (Mesh::element_iterator it=mesh->active_elements_begin(), end=mesh->active_elements_end(); it!=end; ++it)
	{
		Elem* el = *it;
		if (el->is_intersected())
		{
			id_type l_elem = mesh->global_to_local_elem(el->get_id());
			id_type curr_qp = mesh->n_volumetric_qps(l_elem);
			for (id_type ce=0; ce<el->n_cohesive_elem(); ++ce)
			{
				CohesiveElem* coh_el = el->get_cohesive_elem(ce);
				std::vector<double> coh_U_curr(coh_el->n_nodes() * nndof);
				for (id_type n=0; n<coh_el->n_nodes(); ++n)
				{
					for (id_type d=0; d<nndof; ++d)
						coh_U_curr[n*nndof + d] = _prob->get_solution()->get_value_global(coh_el->get_node(n)->get_id(), d);
				}

				id_type nqp = coh_el->n_q_points();
				double delta_eff_avg = 0;
				for (id_type qp=0; qp<nqp; ++qp)
				{
					// Compute the current opening vector in the ntt coordinate system
					std::vector<double> delta_ntt;
					ProblemUtilities::cohesiveDeltaFast(delta_ntt, mesh->get_shape(l_elem,curr_qp), coh_U_curr, mesh->get_rot_matrix(l_elem, curr_qp));

					// Depending on the cohesive law, we might want a different measure of faiure
					delta_eff_avg += el->get_cohesive_material(ce)->getNormalizedCohesiveFailure( delta_ntt );

					curr_qp++;
				}
				delta_eff_avg /= nqp;
				local_deltas[count] = delta_eff_avg;
				count++;
			}
		}
	}
	if (mesh->get_rank() == 0)
	{
		if (_individual)
		{
			std::vector<std::vector<double> > other_deltas(mesh->n_ranks() - 1);
			for (id_type p=1; p<mesh->n_ranks(); ++p)
				Utilities::RecieveUnknown(other_deltas[p-1], p, 1, MPI_DOUBLE, mesh->get_comm());

			for (id_type i=0; i<local_deltas.size(); ++i)
				myfile << local_deltas[i] << ",";
			for (id_type p=1; p<mesh->n_ranks(); ++p)
				for (id_type i=0; i<other_deltas[p-1].size(); ++i)
					myfile << other_deltas[p-1][i] << ",";
		}
		else
		{
			std::vector<double> failed_lengths(5);
			for (id_type i=0; i<_local_lengths.size(); ++i)
				for (id_type j=0; j<failed_lengths.size(); ++j)
					if (local_deltas[i] > double(j+1))
						failed_lengths[j] += _local_lengths[i];

			MPI_Status stat;
			std::vector<double> remote_data(failed_lengths.size());
			for (id_type p=1; p<mesh->n_ranks(); ++p)
			{
				MPI_Recv(remote_data.data(), remote_data.size(), MPI_DOUBLE, p, 1, _prob->get_mesh()->get_comm(), &stat);
				for (id_type j=0; j<failed_lengths.size(); ++j)
					failed_lengths[j] += remote_data[j];
			}
			for (id_type j=0; j<failed_lengths.size(); ++j)
			{
				myfile << failed_lengths[j];
				if ((j+1) != failed_lengths.size())
					myfile << ",";
			}
		}

		myfile << "\n";
	}
	else
	{
		if (_individual)
			MPI_Send(local_deltas.data(), local_deltas.size(), MPI_DOUBLE, 0, 1, mesh->get_comm());
		else
		{
			std::vector<double> failed_lengths(5);
			for (id_type i=0; i<_local_lengths.size(); ++i)
				for (id_type j=0; j<failed_lengths.size(); ++j)
					if (local_deltas[i] > double(j+1))
						failed_lengths[j] += _local_lengths[i];
			MPI_Send(failed_lengths.data(), failed_lengths.size(), MPI_DOUBLE, 0, 1, mesh->get_comm());
		}
	}
}