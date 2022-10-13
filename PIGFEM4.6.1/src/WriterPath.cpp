/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated April 2017

##################################################################################
*/
#include "WriterPath.h"
#include "Problem.h"
#include "Mesh.h"
#include "elem.h"
#include "NodalData.h"
#include "material.h"
#include "KDTree.h"

/*
 * The Big 3
 * Constructor need the Problem
 * Destructor doesn't do anything
 * Copy Constructor copies the problem pointer
 */
WriterPath::WriterPath(Problem* prob)
	: Writer(prob), _output_stress(true)
{}
WriterPath::WriterPath()
{}
WriterPath::~WriterPath()
{}
WriterPath::WriterPath(const WriterPath& other)
{
	_prob = other.get_prob();
	store_mesh();
}




void WriterPath::setParameter(std::string param, bool val)
{
	std::transform(param.begin(), param.end(), param.begin(), ::toupper); // Capitilize the name
	if (param=="STRESS" || param=="OUTPUT STRESS" || param=="OUTPUT_STRESS")
		_output_stress = val;
	else
		err_message("Invalid parameter name!");
}
bool WriterPath::getBoolParameter(std::string param)
{
	std::transform(param.begin(), param.end(), param.begin(), ::toupper); // Capitilize the name
	if (param=="STRESS" || param=="OUTPUT STRESS" || param=="OUTPUT_STRESS")
		return _output_stress;
	else
		err_message("Invalid parameter name!");
}





void WriterPath::setParameter(std::string param, std::vector<double> val)
{
	std::transform(param.begin(), param.end(), param.begin(), ::toupper); // Capitilize the name
	if (param=="POINT" || param=="PATH")
	{
		if (val.size() == _prob->get_mesh()->dim())
			_points.push_back(val);
		else
			err_message("Included points for path output must be the same dimension as the problem!");
	}
	else
		err_message("Invalid parameter name!");
}
std::vector<double> WriterPath::getVecParameter(std::string param)
{
	id_type idx = std::stoi(param);
	if (idx < 0 || idx >= _points.size())
		err_message("Invalid point selection");
	else
		return _points[idx];
}




/*
 * Stores the mesh in this Writer's internal storage
 * Can be recalled if the mesh is adapted
 */
void WriterPath::store_mesh()
{
	_N.clear();
	_N.resize(_points.size());
	_dN.clear();
	_dN.resize(_points.size());
	_elements.clear();
	_elements.resize(_points.size());
	_mats.clear();
	_mats.resize(_points.size());
	_on_partition.clear();
	_on_partition.resize(_points.size(), 0);

	// Create a KDTree of the mesh's nodes
	Mesh* mesh = _prob->get_mesh();
	std::vector<std::vector<double> > nodes(mesh->n_local_nodes(), std::vector<double>());
	for (id_type n=0; n<mesh->n_local_nodes(); ++n)
		nodes[n] = mesh->get_node_local(n)->get_coords();
	KDTree* tree = new KDTree(nodes);

	// Do an inverse lookup to find the element that each point lies in along with the associated shape functions
	for (id_type p=0; p<_points.size(); ++p)
	{
		// Perform a nearest-neighbor lookup to find the node that is closest to this query point
		id_type nearest_node_local = tree->closest_point(_points[p]);

		// Loop over all of this node's neighboring elements and compute the inverse map to see if this point lies within that element
		std::vector<id_type>& node_elem  = mesh->get_node_elem_local(nearest_node_local);
		bool found_elem = false;
		for (id_type e=0; e<node_elem.size(); ++e)
		{
			bool inside = false;
			Elem* el = mesh->get_elem_global(node_elem[e]);
			std::vector<double> r_coords = el->inverse_map(_points[p], inside);
			if (inside)
			{
				_elements[p] = el->get_id();
				_on_partition[p] = 1;
				found_elem = true;

				// Get the shape functions of the normal nodes
				double J; // Don't actually care about this but I need it for the function call
				el->ShapeFunctions(r_coords, _N[p], _dN[p], J);

				if (el->is_intersected())
				{
					// Find the correct child element as well
					bool found_child = false;
					for (id_type ie=0; ie<el->n_integration_elem(); ++ie)
					{
						bool inside_child = false;
						Elem* int_el = el->get_integration_elem(ie);
						if (p==39)
							std::cout << "Placeholder\n";
						std::vector<double> rc_coords = int_el->inverse_map(_points[p], inside_child);
						if (inside_child)
						{
							found_child = true;
							// Compute the enrichment
							_mats[p] = mesh->get_element_material_global(el->get_id(), ie);
							std::vector<double> child_N(el->n_nodes());
							std::vector<std::vector<double> > child_dN(el->n_nodes());
							int_el->ShapeFunctions(rc_coords, child_N, child_dN, J);
							// Add the enrichment functions to the end of the vectors
							_N[p].resize(el->n_nodes() + el->n_enrich_nodes());
							for(id_type en=0; en<el->n_enrich_nodes(); ++en)
							{
								_dN[p].push_back(std::vector<double>(_points[0].size(), 0.0));
								id_type enrich_id = el->get_enrich_node(en)->get_id();
								for(id_type n=0; n<int_el->n_nodes(); ++n)
								{
									id_type int_node_id = (*int_el)(n).get_id();
									if (int_node_id == enrich_id)
									{
										_N[p][el->n_nodes()+en] = child_N[n];
										_dN[p][el->n_nodes()+en] = child_dN[n];
										break;
									}
								}
							}
							break;
						}
						else
							continue;
					}
					if (!found_child)
						err_message("Unable to find the correct child element for a point in the path output!");
				} // end is intersected
				else
					_mats[p] = mesh->get_element_material_global(el->get_id());

				break;
			} // end found parent element
			else
				continue;
		} // End loop over neighboring elements
	} // end points loop

	// Do a reduction across all processors to make sure that at least one partition found every point
	std::vector<int> res(_points.size(), 0);
	MPI_Reduce(_on_partition.data(), res.data(), _on_partition.size(), MPI_INT, MPI_SUM, 0, mesh->get_comm());
	if (mesh->get_rank()==0)
		for (id_type i=0; i<res.size(); ++i)
			if (res[i] == 0)
			{
				std::stringstream ss;
				ss << i;
				std::string out = "Unable to find point " + ss.str();
				err_message( out.data() );
			}
}








/*
 * Main function used to write problem to any inherited file type
 */
void WriterPath::write(std::string filename, double curr_t)
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
void WriterPath::writeFromStream(std::ofstream& myfile, double curr_t)
{
	Mesh* mesh = _prob->get_mesh();
	id_type nndof = _prob->nndof();

	// Loop over all of the query points and compute the solution
	std::vector<id_type> voigt = {1,3,6};
	id_type dim = mesh->dim();
	std::vector<std::vector<double> > sol(_points.size(), std::vector<double>(nndof, 0.0));
	std::vector<std::vector<double> > stress(_points.size(), std::vector<double>(voigt[dim-1], 0.0));
	for (id_type p=0; p<_points.size(); ++p)
	{
		if (_on_partition[p])
		{
			std::vector<double> elem_sol;
			ProblemUtilities::formElementalVector(elem_sol, mesh->get_elem_global(_elements[p]), _prob->get_solution());

			// Compute the actual solution
			for (id_type d=0; d<nndof; ++d)
				for (id_type n=0; n<_N[p].size(); ++n)
					sol[p][d] += _N[p][n] * elem_sol[n*nndof + d];

			// If I'm supposed to compute the stress, compute it here
			if (_output_stress)
			{
				// I'm gonna ignore internal variables for now because it shouldn't matter for the problems we're interested in right now
				Material::input_params input;
				input.dim = mesh->dim();
				input.delta_t = 0.0; // Not actually doing anything
				std::vector<double> strain;
				_prob->stress_strain_kernel(strain, stress[p],
											_dN[p], _mats[p], input, elem_sol);
			}
		} // End loop over neighboring elements
	} // end points loop


	// If this is parallel, send all the points to partition 0
	if (!mesh->serial())
	{
		id_type n_each = (_output_stress)? (nndof+voigt[dim-1]) : nndof;
		if (mesh->get_rank() != 0)
		{
			id_type n_local = std::count(_on_partition.begin(), _on_partition.end(), 1);
			std::vector<id_type> local_points(n_local);
			std::vector<double> point_data(n_local * n_each);
			id_type count = 0;
			for (id_type p=0; p<_points.size(); ++p)
			{
				if (_on_partition[p])
				{
					local_points[count] = p;
					for (id_type d=0; d<nndof; ++d)
						point_data[n_each * count + d] = sol[p][d];
					if (_output_stress)
						for (id_type s=0; s<voigt[dim-1]; ++s)
							point_data[n_each * count + nndof + s] = stress[p][s];
					count++;
				}
			}
			MPI_Send(local_points.data(), local_points.size(), MPI_ID, 0, 0, mesh->get_comm());
			MPI_Send(point_data.data(), point_data.size(), MPI_DOUBLE, 0, 1, mesh->get_comm());
		}
		else
		{
			for (id_type p=1; p<mesh->n_ranks(); ++p)
			{
				std::vector<id_type> remote_points;		Utilities::RecieveUnknown(remote_points, p, 0, MPI_ID, mesh->get_comm());
				std::vector<double> remote_data;		Utilities::RecieveUnknown(remote_data, p, 1, MPI_DOUBLE, mesh->get_comm());
				if (remote_data.size() != (remote_points.size() * n_each))
					err_message("Mismatch between number of path points and path data sizes!");

				for (id_type pt=0; pt<remote_points.size(); ++pt)
				{
					_on_partition[remote_points[pt]] = 1;
					for (id_type d=0; d<nndof; ++d)
						sol[remote_points[pt]][d] = remote_data[pt*n_each + d];
					if (_output_stress)
						for (id_type s=0; s<voigt[dim-1]; ++s)
							stress[remote_points[pt]][s] = remote_data[n_each * pt + nndof + s];
				}
			}
		}
	}


	// Write all points to file
	if (mesh->serial() || (!mesh->serial() && mesh->get_rank()==0))
	{
		// Actually output the info
		for (id_type p=0; p<_points.size(); ++p)
		{
			for (id_type d=0; d<_points[p].size(); ++d)
				myfile << _points[p][d] << ",";
			if (!_output_stress)
			{
				for (id_type d=0; d<nndof; ++d)
				{
					myfile << sol[p][d];
					if ((d+1) == nndof)
						myfile << "\n";
					else
						myfile << ",";
				}
			}
			else
			{
				for (id_type d=0; d<nndof; ++d)
					myfile << sol[p][d] << ",";
				for (id_type s=0; s<stress[p].size(); ++s)
				{
					myfile << stress[p][s];
					if ((s+1) == stress[p].size())
						myfile << "\n";
					else
						myfile << ",";
				}
			}
		}
	}
}