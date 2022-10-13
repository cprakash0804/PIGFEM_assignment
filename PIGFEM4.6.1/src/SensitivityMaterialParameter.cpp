/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated January 2017

##################################################################################
*/
#include "SensitivityMaterialParameter.h"
#include "elem.h"
#include "CohesiveElem.h"
#include "DofObject.h"
#include "Problem.h"
#include "Mesh.h"
#include "Assembler.h"
#include "DenseMatrix.h"
#include "NodalData.h"
#include "InternalVars.h"
#include "SensitivitySolver.h"
#include "Utilities.h"
#include "Solver.h"
#include "gsl_cblas.h"
#include <iostream>










void SensitivityMaterialParameter::init()
{
	Mesh* mesh = _prob->get_mesh();
	_assemble_elem.resize(mesh->n_local_active_elem());
	std::fill(_assemble_elem.begin(), _assemble_elem.end(), false);

	std::vector<std::string> name = {_mat_name};
	set_material_element_assemblies(_prob, name);
}










void SensitivityMaterialParameter::assemble_elem_dPdd(Elem* el, std::vector<double>& dP_el, std::vector<double>& elem_U_curr, Problem* prob)
{
	// Get the mesh and dof objects
	Mesh* mesh = prob->get_mesh();
	Assembler* assem = prob->get_assembler();
	id_type nndof = prob->get_dofs()->nndof();

	// Get the local element number
	id_type l_elem = mesh->global_to_local_elem(el->get_id());

	// Define the initial input parameters
	Material::sensitivity_input_params input;
	input.dim = mesh->dim();
	input.sensitivity_mat_name = _mat_name;
	input.sensitivity_param_name = _param_name;
	input.parameter_id = _id;
	if (prob->get_classification() == STRUCTURAL)
	{
		std::vector<int> switch_dim = {1, 3, 6};
		input.strain.resize(switch_dim[input.dim - 1]); // I don't know if I need this here. Just a placeholder in case I don't wanna compute the actual strain
		if (prob->get_parameter("plane_strain") == 0.0)
			input.plane_strain = false;
		else
			input.plane_strain = true;
	}


	// If this is a non-intersected element
	if (!el->is_intersected())
	{	
		// Get a pointer to the material object associated with this non-intersected element
		Material* curr_mat = mesh->get_element_material_global(el->get_id());
		bool mat_has_isvs = (curr_mat->n_internal_vars() != 0);
		bool right_material = (curr_mat->get_name()==_mat_name);

		// Loop over all of the quadrature points for this element
		if (mat_has_isvs || right_material)
		{
			std::vector<double> dP_el_temp = dP_el;
			id_type nqp = el->n_q_points();
			for(id_type qp=0; qp<nqp; ++qp)
			{
				// Get the internal variable list associated with this quadrature point
				input.internal_vars = &(prob->get_internal_vars()->get_internal_vars_local(l_elem, qp));
				input.internal_vars_sensitivity = &(prob->get_sensitivity_solver()->get_internal_vars()->get_internal_vars_local(l_elem, qp));
				if (assem->storedBmats()) // Call the B-matrix variant
				{
					KernelVolumetric(dP_el_temp,
									 mesh->get_shape(l_elem,qp), assem->getBmat(l_elem,qp),
									 curr_mat, input,
									 elem_U_curr);
				}
				else
				{
					KernelVolumetric(dP_el_temp,
									 mesh->get_shape(l_elem,qp), mesh->get_shape_grad(l_elem,qp),
									 curr_mat, input,
									 elem_U_curr);
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
		id_type curr_qp = 0; // Can't just use the qp f the integration element to get the internal vars because I store them sequentially in the Problem object
		for(id_type ie=0; ie<el->n_integration_elem(); ++ie)
		{
			// Get the material of the current integration element directly from the element
			Material* curr_mat = mesh->get_element_material_global(el->get_id(), ie);
			bool mat_has_isvs = (curr_mat->n_internal_vars() != 0);
			bool right_material = (curr_mat->get_name()==_mat_name);

			// Loop over the quadrature points
			id_type nqp = el->get_integration_elem(ie)->n_q_points();
			if (mat_has_isvs || right_material)
			{
				for(id_type qp=0; qp<nqp; ++qp)
				{
					// Get the internal variable list associated with this quadrature point
					input.internal_vars = &(prob->get_internal_vars()->get_internal_vars_local(l_elem, curr_qp));
					input.internal_vars_sensitivity = &(prob->get_sensitivity_solver()->get_internal_vars()->get_internal_vars_local(l_elem, curr_qp));
					if (assem->storedBmats()) // Call the B-matrix variant
					{
						KernelVolumetric(dP_el_temp,
										 mesh->get_shape(l_elem,qp), assem->getBmat(l_elem,qp),
										 curr_mat, input,
										 elem_U_curr);
					}
					else
					{
						KernelVolumetric(dP_el_temp,
										 mesh->get_shape(l_elem,qp), mesh->get_shape_grad(l_elem,qp),
										 curr_mat, input,
										 elem_U_curr);
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
				bool right_material = (curr_mat->get_name()==_mat_name);

				// Loop over the quadrature points
				id_type nqp = coh_el->n_q_points();
				if (mat_has_isvs || right_material)
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

					for(id_type qp=0; qp<nqp; ++qp)
					{
						input.internal_vars = &(prob->get_internal_vars()->get_internal_vars_local(l_elem, curr_qp));
						input.internal_vars_sensitivity = &(prob->get_sensitivity_solver()->get_internal_vars()->get_internal_vars_local(l_elem, curr_qp));
						KernelCohesive (dP_el_coh_temp,
										mesh->get_shape(l_elem,curr_qp), mesh->get_rot_matrix(l_elem,curr_qp),
										curr_mat, input,
										coh_U_curr);

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











