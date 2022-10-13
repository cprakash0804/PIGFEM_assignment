/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#include "SensitivityShapeParameter.h"
#include "Problem.h"
#include "Assembler.h"
#include "Mesh.h"
#include "Inclusion.h"
#include "NodalData.h"
#include "elem.h"
#include "CohesiveElem.h"
#include "DofObject.h"
#include "InternalVars.h"
#include "SensitivitySolver.h"
#include "Utilities.h"
#include "gsl_cblas.h"

SensitivityShapeParameter::SensitivityShapeParameter()
	: _velocity(NULL)
{
}


SensitivityShapeParameter::~SensitivityShapeParameter()
{
	if (_velocity != NULL)
		delete _velocity;
}


void SensitivityShapeParameter::init()
{
	Mesh* mesh = _prob->get_mesh();
	_assemble_elem.resize(mesh->n_local_elem());
	std::fill(_assemble_elem.begin(), _assemble_elem.end(), false);
	std::vector<id_type> ids = {_inclusion_id};
	set_shape_element_assemblies(_prob, ids);
	compute_velocity(_prob);
}



void SensitivityShapeParameter::compute_velocity(Problem* prob)
{
	Mesh* mesh = prob->get_mesh();

	// Preallocate the velocity nodal data strucuture
	_velocity = new NodalData(mesh);
	_velocity->preallocate_storage(mesh->dim());
	
	// Loop over all of the enrichment nodes and determine the velocity of their intersection
	for (Mesh::enrich_node_iterator it=mesh->enrich_nodes_begin(), end=mesh->enrich_nodes_end(); it!=end; ++it)
	{
		id_type l_node = mesh->global_to_local_node((*it)->get_id());
		id_type l_inclusion = mesh->get_enriched_inclusion_local(l_node);

		// Actually on the correct inclusion
		if (mesh->get_inclusion(l_inclusion)->get_id() == _inclusion_id)
		{
			std::pair<id_type, id_type> e = mesh->get_enriched_edge_local(l_node);
			std::vector<Node*> edge = {mesh->get_node_global(e.first), mesh->get_node_global(e.second)};	// Pointers to the nodes f the parent edge
			std::vector<double> intersection = mesh->get_node_local(l_node)->get_coords(); // Coordinates of the intersection

			// Compute the actual intersection velocity
			std::vector<double> v = mesh->get_inclusion(l_inclusion)->getIntersectionVelocity(edge, intersection, _param_name);
			for (id_type d=0; d<mesh->dim(); ++d)
				_velocity->get_value_local(l_node, d) = v[d];
		}
	}



	// Now that I've compute the nodal velocities, look over all of the intersected elements and store their shape function partials fields
	_dshape.resize(mesh->n_local_elem());
	_dshape_grad.resize(mesh->n_local_elem());
	_div_v.resize(mesh->n_local_elem());
	_dRot.resize(mesh->n_local_elem());
	double w; // temp
	for (Mesh::element_iterator it=mesh->active_elements_begin(), end=mesh->active_elements_end(); it!=end; ++it)
	{
		Elem* el =(*it);
		id_type local_e = mesh->global_to_local_elem(el->get_id());

		// Fill in placeholder entries here
		if (!el->is_intersected())
		{
			id_type nqp = el->n_q_points();
			_dshape[local_e].resize(nqp);
			_dshape_grad[local_e].resize(nqp);
			_div_v[local_e].resize(nqp);
		}

		// Otherwise have to actually do shape function computation
		else
		{
			id_type ngqp = 0;
			for (id_type ie=0; ie<el->n_integration_elem(); ++ie)
				ngqp += el->get_integration_elem(ie)->n_q_points();
			id_type nvqp = ngqp;
			for (id_type ce=0; ce<el->n_cohesive_elem(); ++ce)
				ngqp += el->get_cohesive_elem(ce)->n_q_points();

			_dshape[local_e].resize(ngqp);
			_dshape_grad[local_e].resize(ngqp);
			_div_v[local_e].resize(ngqp);
			_dRot[local_e].resize(ngqp - nvqp);

			// Assemble the elemental velocity solution
			std::vector<double> nodal_velocity;
			ProblemUtilities::formElementalVector(nodal_velocity, el, _velocity);

			id_type local_qp = 0;;
			for (id_type ie=0; ie<(*it)->n_integration_elem(); ++ie)
			{
				Elem* int_el = el->get_integration_elem(ie);
				id_type nqp = int_el->n_q_points();
				for (id_type qp=0; qp<nqp; ++qp)
				{
					// First get the integration point in child element referenece coordinates 
					std::vector<double> child_rcoords;
					int_el->q_point(child_rcoords, w, qp);

					std::vector<double> v;
					el->ShapeFunctionPartials(child_rcoords, ie, nodal_velocity, _dshape[local_e][local_qp],
											  _dshape_grad[local_e][local_qp], v, _div_v[local_e][local_qp]);
					local_qp++;
				}
			}

			// Now do the cohesive elements
			std::vector<std::vector<short_id_type> >& coh_elem_struct = el->getCohesiveElemStructure();
			for(id_type ce=0; ce<el->n_cohesive_elem(); ++ce)
			{
				CohesiveElem* coh_el = el->get_cohesive_elem(ce);

				// Get the nodal solution object for the current cohesive element
				std::vector<double> coh_velocity(coh_el->n_nodes() *  mesh->dim()); // The vectorized form of the previous structure
				for(id_type n=0; n<coh_el->n_nodes(); ++n)
				{
					id_type node = coh_elem_struct[ce][n];
					for(id_type d=0; d< mesh->dim(); ++d)
						coh_velocity[n* mesh->dim() + d] = nodal_velocity[node* mesh->dim() + d];
				}

				id_type nqp = coh_el->n_q_points();
				for (id_type qp=0; qp<nqp; ++qp)
				{
					std::vector<double> coh_rcoords;
					coh_el->q_point(coh_rcoords, w, qp);
					std::vector<double> v;
					coh_el->RotationSensitivity(coh_rcoords, coh_velocity, _dRot[local_e][local_qp-nvqp], v, _div_v[local_e][local_qp]);
					local_qp++;
				}
			}
		}
	}

	// If the assembler the problem used stored the B-matrices, then I should store the B-matrics as well
	Assembler* assem = prob->get_assembler();
	if (assem->storedBmats())
	{
		Mesh* mesh = prob->get_mesh();
		_dBmats.resize( mesh->n_local_active_elem() );
		for (Mesh::element_iterator it=mesh->active_elements_begin(), end=mesh->active_elements_end(); it!=end; ++it)
		{
			id_type l_elem = mesh->global_to_local_elem((*it)->get_id());
			if (assembleElem(l_elem))
			{
				id_type nqp = mesh->n_volumetric_qps(l_elem);
				_dBmats[l_elem].resize( nqp );
				for (id_type qp=0; qp<nqp; ++qp)
					assem->fillBmat(_dBmats[l_elem][qp], get_shape_grad(l_elem, qp));
			}
		}
	}
}

double& SensitivityShapeParameter::get_div_v(id_type local_e, id_type qp)
{
	if (local_e < _div_v.size())
	{
		if (qp < _div_v[local_e].size())
			return _div_v[local_e][qp];
		else
			err_message("Invalid Quadrature point.");
	}
	else
		err_message("Invalid element number.");
}
std::vector<double>& SensitivityShapeParameter::get_shape(id_type local_e, id_type qp)
{
	if (local_e < _dshape.size())
	{
		if (qp < _dshape[local_e].size())
			return _dshape[local_e][qp];
		else
			err_message("Invalid Quadrature point.");
	}
	else
		err_message("Invalid element number.");
}
std::vector<std::vector<double> >& SensitivityShapeParameter::get_shape_grad(id_type local_e, id_type qp)
{
	if (local_e < _dshape_grad.size())
	{
		if (qp < _dshape_grad[local_e].size())
			return _dshape_grad[local_e][qp];
		else
			err_message("Invalid Quadrature point.");
	}
	else
		err_message("Invalid element number.");
}
DenseMatrix<double>& SensitivityShapeParameter::get_rot_matrix(id_type local_e, id_type qp)
{
	if (local_e < _dRot.size())
	{
		id_type nvqp = _dshape[local_e].size() - _dRot[local_e].size();
		qp -= nvqp;
		if (qp >= 0 && qp < _dRot[local_e].size())
			return _dRot[local_e][qp];
		else
			err_message("Invalid Quadrature point.");
	}
	else
		err_message("Invalid element number.");
}
DenseMatrix<double>& SensitivityShapeParameter::getBmat(id_type local_e, id_type qp)
{
	if (local_e < _dBmats.size())
	{
		if (qp < _dBmats[local_e].size())
			return _dBmats[local_e][qp];
		else
			err_message("Invalid Quadrature point.");
	}
	else
		err_message("Invalid element number.");
}


















void SensitivityShapeParameter::assemble_elem_dPdd(Elem* el, std::vector<double>& dP_el, std::vector<double>& elem_U_curr, Problem* prob)
{
	// Get the mesh and dof objects
	Mesh* mesh = prob->get_mesh();
	Assembler* assem = prob->get_assembler();
	id_type nndof = prob->get_dofs()->nndof();

	// Get the local element number
	id_type l_elem = mesh->global_to_local_elem(el->get_id());

	// Define the initial input parameters
	Material::input_params input;
	Material::sensitivity_input_params sens_input;
	input.dim = mesh->dim();
	sens_input.dim = mesh->dim();
	sens_input.sensitivity_mat_name = "#SHAPE_SENSITIVITY#";
	sens_input.sensitivity_param_name = _param_name;
	sens_input.parameter_id = _id;

	if (prob->get_classification() == STRUCTURAL)
	{
		std::vector<int> switch_dim = {1, 3, 6};
		input.strain.resize(switch_dim[input.dim - 1]); // I don't know if I need this here. Just a placeholder in case I don't wanna compute the actual strain
		sens_input.strain.resize(switch_dim[input.dim - 1]);
		if (prob->get_parameter("plane_strain") == 0.0)
			input.plane_strain = false;
		else
			input.plane_strain = true;
		sens_input.plane_strain = input.plane_strain;
	}

	// Copy over the current internal state valiables into a new object so that when I call the Constitutive call it doesn't update the ISVs
	std::vector<std::vector<double> > temp_ISVs = prob->get_internal_vars()->get_elem_internal_vars_local(l_elem);


	// If this is a non-intersected element
	// (Note, this should never ever be actually called because the velocity in a non-intersected element is 0. Should be caught in the _assemble_elem vector)
	if (!el->is_intersected())
	{
		// Get a pointer to the material object associated with this non-intersected element
		Material* curr_mat = mesh->get_element_material_global(el->get_id());
		bool mat_has_isvs = (curr_mat->n_internal_vars() != 0);

		// Only need to compute anything for a non-intersected element if said element has ISVs
		if (mat_has_isvs)
		{
			std::vector<double> dP_el_temp = dP_el;
			// Loop over all of the quadrature points for this element
			id_type nqp = el->n_q_points();
			for(id_type qp=0; qp<nqp; ++qp)
			{
				input.internal_vars = &(temp_ISVs[qp]);
				sens_input.internal_vars = &(prob->get_internal_vars()->get_internal_vars_local(l_elem, qp));
				sens_input.internal_vars_sensitivity = &(prob->get_sensitivity_solver()->get_internal_vars()->get_internal_vars_local(l_elem, qp));
				if (assem->storedBmats()) // Call the B-matrix variant
				{
					KernelVolumetric (dP_el_temp,
									  mesh->get_shape(l_elem,qp), assem->getBmat(l_elem,qp),
									  get_shape(l_elem,qp), getBmat(l_elem,qp), get_div_v(l_elem, qp),
									  curr_mat, input, sens_input,
									  elem_U_curr, assembleElem(l_elem)); // assembleElem() should always be false here
				}
				else
				{
					KernelVolumetric (dP_el_temp,
									  mesh->get_shape(l_elem,qp), mesh->get_shape_grad(l_elem,qp),
									  get_shape(l_elem,qp), get_shape_grad(l_elem,qp), get_div_v(l_elem, qp),
									  curr_mat, input, sens_input,
									  elem_U_curr, assembleElem(l_elem));
				}
				double& J = mesh->get_J(l_elem,qp);
				double& W = mesh->get_W(l_elem,qp);
				double wgt = J * W;
				cblas_daxpy(dP_el.size(), wgt, dP_el_temp.data(), 1, dP_el.data(), 1);
			}
		}
	}
	
	// Otherwise, this element is intersected somehow
	else
	{
		std::vector<double> dP_el_temp = dP_el;

		// Loop over all of the integration elements for this intersected element
		id_type curr_qp = 0;
		for(id_type ie=0; ie<el->n_integration_elem(); ++ie)
		{
			// Get the material of the current integration element directly from the element
			Material* curr_mat = mesh->get_element_material_global(el->get_id(), ie);
			bool mat_has_isvs = (curr_mat->n_internal_vars() != 0);

			id_type nqp = el->get_integration_elem(ie)->n_q_points();
			if (mat_has_isvs || assembleElem(l_elem))
			{
				// Loop over the quadrature points
				for(id_type qp=0; qp<nqp; ++qp)
				{
					// Get the internal variable list associated with this quadrature point
					input.internal_vars = &(temp_ISVs[qp]);
					sens_input.internal_vars = &(prob->get_internal_vars()->get_internal_vars_local(l_elem, curr_qp));
					sens_input.internal_vars_sensitivity = &(prob->get_sensitivity_solver()->get_internal_vars()->get_internal_vars_local(l_elem, curr_qp));
					if (assem->storedBmats()) // Call the B-matrix variant
					{
						KernelVolumetric (dP_el_temp,
										  mesh->get_shape(l_elem,curr_qp), assem->getBmat(l_elem,curr_qp),
										  get_shape(l_elem,curr_qp), getBmat(l_elem,curr_qp), get_div_v(l_elem, curr_qp),
										  curr_mat, input, sens_input,
										  elem_U_curr, assembleElem(l_elem));
					}
					else
					{
						KernelVolumetric (dP_el_temp,
										  mesh->get_shape(l_elem,curr_qp), mesh->get_shape_grad(l_elem,curr_qp),
										  get_shape(l_elem,curr_qp), get_shape_grad(l_elem,curr_qp), get_div_v(l_elem, curr_qp),
										  curr_mat, input, sens_input,
										  elem_U_curr, assembleElem(l_elem));
					}
					double& J = mesh->get_J(l_elem,curr_qp);
					double& W = mesh->get_W(l_elem,curr_qp);
					double wgt = J * W;
					cblas_daxpy(dP_el.size(), wgt, dP_el_temp.data(), 1, dP_el.data(), 1);
					// Update the current qp
					curr_qp++;
				}
			}
			else
			{
				curr_qp += nqp;
				continue;
			}
		}


		// If the problem is cohesive then I need to add the contribution from the cohesive opening
		// This really doesn't fit int the general model of the code but I can't think of another way to include it right now...
		if (mesh->is_cohesive())
		{
			// Get the cohesive element structure (Needed for mapping the cohesive matrix back to the elemental matrix)
			std::vector<std::vector<short_id_type> >& coh_elem_struct = el->getCohesiveElemStructure();

			// Loop over the cohesive elements
			for(id_type ce=0; ce<el->n_cohesive_elem(); ++ce)
			{
				// Get a pointer to the current cohesive element
				CohesiveElem* coh_el = el->get_cohesive_elem(ce);

				// Get the material of the current cohesive element directly from the element
				Material* curr_mat = el->get_cohesive_material(ce);
				bool mat_has_isvs = (curr_mat->n_internal_vars() != 0);

				id_type nqp = coh_el->n_q_points();
				if (mat_has_isvs || assembleElem(l_elem))
				{
					// Initialize the cohesive contributions
					id_type n_dof_coh = coh_el->n_nodes() * nndof;
					std::vector<double> dP_el_coh(n_dof_coh);
					std::vector<double> dP_el_coh_temp = dP_el_coh;

					// Get the nodal solution object for the current cohesive element
					std::vector<double> coh_U_curr(coh_el->n_nodes() * nndof);
					for(id_type n=0; n<coh_el->n_nodes(); ++n)
					{
						id_type node = coh_elem_struct[ce][n];
						for(id_type d=0; d<nndof; ++d)
							coh_U_curr[n*nndof + d] = elem_U_curr[node*nndof + d];
					}

					// Loop over the quadrature points
					for(id_type qp=0; qp<nqp; ++qp)
					{
						input.internal_vars = &(temp_ISVs[curr_qp]);
						sens_input.internal_vars = &(prob->get_internal_vars()->get_internal_vars_local(l_elem, curr_qp));
						sens_input.internal_vars_sensitivity = &(prob->get_sensitivity_solver()->get_internal_vars()->get_internal_vars_local(l_elem, curr_qp));
						KernelCohesive (dP_el_coh_temp,
										mesh->get_shape(l_elem,curr_qp), mesh->get_rot_matrix(l_elem,curr_qp),
										get_shape(l_elem,curr_qp), get_rot_matrix(l_elem,curr_qp), get_div_v(l_elem, curr_qp),
										curr_mat, input, sens_input,
										coh_U_curr, assembleElem(l_elem));

						double& J = mesh->get_J(l_elem,curr_qp);
						double& W = mesh->get_W(l_elem,curr_qp);
						double wgt = J * W;
						cblas_daxpy(dP_el_coh.size(), wgt, dP_el_coh_temp.data(), 1, dP_el_coh.data(), 1);
						// Update the quadrature point
						curr_qp++;
					}




					// Add the cohesive contributions to the elemental stiffness matrix and internal force vector
					for(id_type n1=0; n1<coh_elem_struct[ce].size(); ++n1)
					{
						id_type row_node = coh_elem_struct[ce][n1];
						for(id_type d1=0; d1<nndof; ++d1)
						{
							id_type row_dof = row_node*nndof + d1;
							id_type row_dof_local = n1*nndof + d1;
							// Internal force contribution
							dP_el[row_dof] += dP_el_coh[row_dof_local];
						}
					} // End cohesive contribution update
				}
				else
				{
					curr_qp += nqp;
					continue;
				}
			} // End cohesive elem loop
		} // End is coheisve if-statement

	} // End if intersected
}


