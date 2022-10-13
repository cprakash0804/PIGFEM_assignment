/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated November 2016

##################################################################################
*/
#include "SubscaleModel.h"
#include "Problem.h"
#include "Mesh.h"
#include "NodalData.h"
#include "InternalVars.h"
#include "BoundaryObject.h"
#include "elem.h"
#include "node.h"
#include "Utilities.h"





SubscaleModel::SubscaleModel()
	: _strain_set(false)
{
}
SubscaleModel::SubscaleModel(int *argc, char***argv)
	: _strain_set(false)
{
}
SubscaleModel::~SubscaleModel()
{
	Mesh* mesh = _prob->get_mesh();
	delete _prob;
	delete mesh;
}




void SubscaleModel::homogenize_stress_strain(const std::vector<double>& curr_macro_strain, std::vector<double>& avg_stress, std::vector<double>& avg_perturbation_strain)
{
	Mesh* mesh = get_mesh();
	avg_stress.clear();
	avg_stress.resize(curr_macro_strain.size());
	avg_perturbation_strain.clear();
	avg_perturbation_strain.resize(curr_macro_strain.size());
	double volume_local = 0;
	std::vector<double> avg_stress_local(avg_stress.size()), perturbation_strain_local(avg_perturbation_strain.size());
	id_type nndof = _prob->nndof();

	for (Mesh::element_iterator it=mesh->active_elements_begin(), end=mesh->active_elements_end(); it!=end; ++it)
	{
		Elem* el = (*it);
		id_type l_elem = mesh->global_to_local_elem(el->get_id());

		// Assemble the current displacement solution
		id_type nn = el->n_nodes() + el->n_enrich_nodes();
		id_type ndof = nn*nndof;
		std::vector<double> elem_sol(ndof);
		for(id_type n=0; n<nn; ++n)
		{
			id_type l_node;
			if(n<el->n_nodes())
				l_node = mesh->global_to_local_node(el->get_node(n)->get_id());
			else
				l_node = mesh->global_to_local_node(el->get_enrich_node(n-el->n_nodes())->get_id());

			for(id_type d=0; d<nndof; ++d)
				elem_sol[n*nndof + d] = _prob->get_solution()->get_value_local(l_node, d);
		}

		// Set up the material input object
		Material::input_params input;
		input.dim = mesh->dim();
		input.delta_t = 0.0; // Not actually doing anything
		if (_prob->get_classification() == STRUCTURAL)
		{
			std::vector<int> switch_dim = {1, 3, 6};
			input.strain.resize(switch_dim[input.dim - 1]); // I don't know if I need this here. Just a placeholder in case I don't wanna compute the actual strain
			input.delta_t = 0.0;
			if (_prob->get_parameter("plane_strain") == 0.0)
				input.plane_strain = false;
			else
				input.plane_strain = true;
		}

		// Set up quadrature point vectors
		std::vector<double> stress(curr_macro_strain.size()), strain(curr_macro_strain.size());

		if(!el->is_intersected())
		{
			// Get a pointer to the material object associated with this non-intersected element
			Material* curr_mat = mesh->get_element_material_global(el->get_id());

			// Loop over all of the quadrature points for this element
			id_type nqp = el->n_q_points();
			for(id_type qp=0; qp<nqp; ++qp)
			{
				std::vector<double> int_vars_copy = _prob->get_internal_vars()->get_internal_vars_local(l_elem, qp);
				input.internal_vars = &int_vars_copy;
				stress_strain_kernel(strain, stress,
									 mesh->get_shape_grad(l_elem, qp),
									 curr_mat, input, elem_sol,
									 curr_macro_strain);

				// Update current stress and volume
				double& J = mesh->get_J(l_elem, qp);
				double& W = mesh->get_W(l_elem, qp);
				for (id_type i=0; i<avg_stress_local.size(); ++i)
					avg_stress_local[i] += W * J* stress[i];
				for (id_type i=0; i<perturbation_strain_local.size(); ++i)
					perturbation_strain_local[i] += W * J * strain[i];
				volume_local += W * J;
			}
		}
		else
		{
			// Loop over all of the integration elements for this intersected element
			id_type curr_qp = 0; // Can't just use the qp f the integration element to get the internal vars because I store them sequentially in the Problem object
			for(id_type ie=0; ie<el->n_integration_elem(); ++ie)
			{
				// Get the material of the current integration element directly from the element
				Material* curr_mat = mesh->get_element_material_global(el->get_id(), ie);

				// Loop over the quadrature points
				id_type nqp = el->get_integration_elem(ie)->n_q_points();
				for(id_type qp=0; qp<nqp; ++qp)
				{
					std::vector<double> int_vars_copy = _prob->get_internal_vars()->get_internal_vars_local(l_elem, curr_qp);
					input.internal_vars = &int_vars_copy;
					stress_strain_kernel(strain, stress,
										 mesh->get_shape_grad(l_elem, curr_qp),
										 curr_mat, input, elem_sol,
										 curr_macro_strain);

					// Update current stress and volume
					double& J = mesh->get_J(l_elem, curr_qp);
					double& W = mesh->get_W(l_elem, curr_qp);
					for (id_type i=0; i<avg_stress_local.size(); ++i)
						avg_stress_local[i] += W * J * stress[i];
					for (id_type i=0; i<perturbation_strain_local.size(); ++i)
						perturbation_strain_local[i] += W * J * strain[i];
					volume_local += W * J;

					// Update the current qp
					curr_qp++;
				}
			}
		}
	}

	// Global reductions
	double volume;
	if (_prob->get_mesh()->serial())
	{
		avg_stress = avg_stress_local;
		volume = volume_local;
	}
	else
	{
		MPI_Reduce(&volume_local, &volume, 1, MPI_DOUBLE, MPI_SUM, 0, _prob->get_mesh()->get_comm());
		MPI_Reduce(avg_stress_local.data(), avg_stress.data(), avg_stress.size(), MPI_DOUBLE, MPI_SUM, 0, _prob->get_mesh()->get_comm());
		MPI_Reduce(perturbation_strain_local.data(), avg_perturbation_strain.data(), avg_perturbation_strain.size(), MPI_DOUBLE, MPI_SUM, 0, _prob->get_mesh()->get_comm());
	}

	// Do the homogenization
	if (_prob->get_mesh()->get_rank() == 0)
	{
		for (id_type i=0; i<avg_stress.size(); ++i)
		{
			avg_stress[i] /= volume;
			avg_perturbation_strain[i] /= volume;
		}
	}
}

void SubscaleModel::stress_strain_kernel(std::vector<double>& strain, std::vector<double>& stress,
										 const std::vector<std::vector<double> >& shape_grad,
										 Material* mat, Material::input_params& input,
										 const std::vector<double>& elem_U_curr,
										 const std::vector<double>& curr_macro_strain)
{
	ProblemUtilities::SmallStrainFast(strain, shape_grad, elem_U_curr);

	// Set up input object (add in the macroscopic strain)
	input.strain = strain;
	for (id_type i=0; i<strain.size(); ++i)
		input.strain[i] += curr_macro_strain[i];

	// Compute the constitutive relations for this quadrature point
	Material::output_params* output = mat->Constitutive(input);   // Computes the current stress (as well as the current stiffness matrix but oh well)
	stress = output->stress;
}


void SubscaleModel::write_stress_strain(std::string filename, id_type step_iter, const std::vector<double>& current_macro_strain,
										const std::vector<double>& stress, const std::vector<double>& perturbation_strain)
{
	// Add the file extension if it isn't already there
	if (filename.find(".csv", filename.length()-4) == std::string::npos)
		filename = filename + ".csv";

	// Open file for writing or appending base on the load step
	std::ofstream myfile;
	if (step_iter==0)
		myfile.open(filename.c_str(), std::ofstream::out);
	else
		myfile.open(filename.c_str(), std::ofstream::app);
	myfile.precision(16);

	if (_prob->get_mesh()->get_rank() == 0)
	{
		for (id_type i=0; i<current_macro_strain.size(); ++i)
			myfile << current_macro_strain[i] << ",";
		for (id_type i=0; i<stress.size(); ++i)
			myfile << stress[i] << ",";
		for (id_type i=0; i<perturbation_strain.size(); ++i)
		{
			myfile << perturbation_strain[i];
			if (i < (perturbation_strain.size()-1))
				myfile << ",";
			else
				myfile << "\n";
		}
	}
}



void SubscaleModel::init()
{
	// Assign BCs
	BoundaryObject* boundary = _prob->get_boundary();
	for(unsigned int d=0; d<get_mesh()->dim(); ++d)
		boundary->set_dirichlet_bc(0, d, 0);

	// Determine the lowest value of each coordinate for macro-displacement purposes
	determine_origin();

	// initialize the problem
	_prob->init();
}
void SubscaleModel::solve()
{
	if (!_strain_set)
		err_message("Please set a macro-scale strain before attempting to solve a subscale model");

 	_prob->solve_problem();
}





void SubscaleModel::determine_origin()
{
	// Determine the lowest vlues of coordinates on my processor
	std::vector<double> mins(_mesh->dim(), 99999999999999999999999999.9);
	for (Mesh::node_iterator it=_mesh->nodes_begin(), end=_mesh->nodes_end(); it!=end; ++it)
	{
		for (id_type d=0; d<_mesh->dim(); ++d)
		{
			double coord = (*(*it))(d);
			if (coord < mins[d])
				mins[d] = coord;
		}
	}

	// If this is parallel perform a global reduction
	if (_mesh->serial())
		_origin = mins;
	else
	{
		_origin.resize(_mesh->dim());
		MPI_Allreduce(mins.data(), _origin.data(), _origin.size(), MPI_DOUBLE, MPI_MIN, _mesh->get_comm());
	}
}


std::vector<double> SubscaleModel::compute_macro_displacement(const std::vector<double>& curr_macro_strain, const std::vector<double> coords)
{
	// 1D
	if (curr_macro_strain.size() == 1)
	{
		std::vector<double> ret(coords.size());
		ret[0] = (coords[0] - _origin[0]) * curr_macro_strain[0];
		return ret;
	}

	// 2D
	else if (curr_macro_strain.size() == 3)
	{
		std::vector<double> ret(coords.size());
		ret[0] = (coords[0] - _origin[0]) * curr_macro_strain[0] + (coords[1] - _origin[1]) * curr_macro_strain[2]; // Not 100% sure about the shear strains. Might need a 1/2
		ret[1] = (coords[0] - _origin[0]) * curr_macro_strain[2] + (coords[1] - _origin[1]) * curr_macro_strain[1];
		return ret;
	}

	// 3D
	else if (curr_macro_strain.size() == 6)
	{
		std::vector<double> ret(coords.size());
		ret[0] = (coords[0] - _origin[0]) * curr_macro_strain[0] + (coords[1] - _origin[1]) * curr_macro_strain[5] + (coords[2] - _origin[2]) * curr_macro_strain[4]; // Not 100% sure about the shear strains. Might need a 1/2
		ret[1] = (coords[0] - _origin[0]) * curr_macro_strain[5] + (coords[1] - _origin[1]) * curr_macro_strain[1] + (coords[2] - _origin[2]) * curr_macro_strain[3];
		ret[2] = (coords[0] - _origin[0]) * curr_macro_strain[4] + (coords[1] - _origin[1]) * curr_macro_strain[3] + (coords[2] - _origin[2]) * curr_macro_strain[2];
		return ret;
	}

	else
		err_message("Invalid macroscopic strain input!");
}