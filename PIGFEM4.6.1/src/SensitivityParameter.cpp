/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#include "SensitivityParameter.h"
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
#include <algorithm>


SensitivityParameter::SensitivityParameter()
	: _id(0), _prob(NULL)
{
}

void SensitivityParameter::init()
{
	_assemble_elem.resize(_prob->get_mesh()->n_local_active_elem());
	std::fill(_assemble_elem.begin(), _assemble_elem.end(), false);
}


void SensitivityParameter::set_material_element_assemblies(Problem* prob, std::vector<std::string>& names)
{
	std::sort(names.begin(), names.end());
	Mesh* mesh = prob->get_mesh();
	for(Mesh::element_iterator it=mesh->active_elements_begin(), end=mesh->active_elements_end(); it!=end; ++it)
	{
		Elem* el = (*it);
		std::vector<std::string> elem_names;
		if (!el->is_intersected())
		{
			elem_names.push_back( mesh->get_element_material_global(el->get_id())->get_name() );
			std::vector<std::string> intersection;
			std::set_intersection(names.begin(), names.end(),
                        		  elem_names.begin(), elem_names.end(),
                        		  std::back_inserter(intersection));
			if (intersection.size() != 0)
			{
				id_type l_elem = mesh->global_to_local_elem((*it)->get_id());
				_assemble_elem[l_elem] = true;
			}
		}

		else
		{
			for (id_type ie=0; ie<el->n_integration_elem(); ++ie)
				elem_names.push_back( mesh->get_element_material_global(el->get_id(), ie)->get_name() );

			if (mesh->is_cohesive())
			{
				for (id_type ce=0; ce<el->n_cohesive_elem(); ++ce)
					elem_names.push_back( el->get_cohesive_material(ce)->get_name() );
			}

			std::sort(elem_names.begin(), elem_names.end());
			std::vector<std::string> intersection;
			std::set_intersection(names.begin(), names.end(),
                        		  elem_names.begin(), elem_names.end(),
                        		  std::back_inserter(intersection));
			if (intersection.size() != 0)
			{
				id_type l_elem = mesh->global_to_local_elem((*it)->get_id());
				_assemble_elem[l_elem] = true;
			}
		}
	}
}
void SensitivityParameter::set_shape_element_assemblies(Problem* prob, std::vector<id_type>& ids)
{
	std::sort(ids.begin(), ids.end());
	Mesh* mesh = prob->get_mesh();
	for(Mesh::element_iterator it=mesh->active_elements_begin(), end=mesh->active_elements_end(); it!=end; ++it)
	{
		Elem* el = (*it);
		if (el->is_intersected())
		{
			std::vector<id_type> elem_ids = el->getInclusionNumbers();
			std::sort(elem_ids.begin(), elem_ids.end());
			std::vector<id_type> intersection;
			std::set_intersection(ids.begin(), ids.end(),
                        		  elem_ids.begin(), elem_ids.end(),
                        		  std::back_inserter(intersection));
			if (intersection.size() != 0)
			{
				id_type l_elem = mesh->global_to_local_elem((*it)->get_id());
				_assemble_elem[l_elem] = true;
			}
		}
	}
}



PetscErrorCode SensitivityParameter::assemble_dPdd(pvector& dPdd, Problem* prob)
{
	// Some initializations
	PetscErrorCode ierr;
	Mesh* mesh = prob->get_mesh();
	InternalVars* int_vars = prob->get_internal_vars();
	DofObject* dofs = prob->get_dofs();
	id_type ngfd = prob->get_dofs()->n_global_free_dofs();

	// Assemble the partial derivative of internal load
	std::vector<double> dP_el;
	std::vector<id_type> elem_dofs;
	std::vector<double> elem_sol;
	for(Mesh::element_iterator it=mesh->active_elements_begin(), end=mesh->active_elements_end(); it!=end; ++it)
	{
		id_type l_elem = mesh->global_to_local_elem((*it)->get_id());
		if (_assemble_elem[l_elem] || int_vars->elemHasISVs(l_elem))
		{
			// Get the global dofs for this element
			dofs->get_elem_global_dofs((*it), elem_dofs); // elem_dofs will now contain a vector with all the global dofs associated with the element
			dP_el.resize(elem_dofs.size());
			std::fill(dP_el.begin(), dP_el.end(), 0.0);

			// Assemble the current solution vector
			ProblemUtilities::formElementalVector(elem_sol, (*it), prob->get_solution());

			// Assemble the elemental contribution
			assemble_elem_dPdd((*it), dP_el, elem_sol, prob);

			// Assemble elemental vector into the global vector
			set_dP(dPdd, dP_el, elem_dofs, ngfd);
		}
	}
	
	// PETSc communication routines
	ierr = VecAssemblyBegin(dPdd.f);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(dPdd.p);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(dPdd.f);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(dPdd.p);CHKERRQ(ierr);
	
	return ierr;
}




void SensitivityParameter::updateSensitivityISVs(Problem* prob, NodalData* dU_dd)
{
	// Some initializations
	Mesh* mesh = prob->get_mesh();
	InternalVars* int_vars = prob->get_internal_vars();

	// Assemble the partial derivative of internal load
	std::vector<double> elem_sol;
	std::vector<double> delem_sol;
	for(Mesh::element_iterator it=mesh->active_elements_begin(), end=mesh->active_elements_end(); it!=end; ++it)
	{
		if (int_vars->elemHasISVs(mesh->global_to_local_elem((*it)->get_id())))
		{
			// Assemble the current solution vector
			ProblemUtilities::formElementalVector(elem_sol, (*it), prob->get_solution());
			ProblemUtilities::formElementalVector(delem_sol, (*it), dU_dd);

			// Assemble the elemental contribution
			updateSensitivityElemISVs((*it), elem_sol, delem_sol, prob);
		}
	}
}


void SensitivityParameter::updateSensitivityElemISVs(Elem* el, std::vector<double>& elem_U_curr, std::vector<double>& delem_U_curr, Problem* prob)
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
	input.parameter_id = _id;
	input.delta_t = prob->get_solver()->getDeltaT();

	// If this is a non-intersected element
	if (!el->is_intersected())
	{	
		Material* curr_mat = mesh->get_element_material_global(el->get_id());
		bool mat_has_isvs = (curr_mat->n_internal_vars() != 0);
		// Only actually run the Kernel if there are ISVs to update
		if (mat_has_isvs)
		{
			// Get a pointer to the material object associated with this non-intersected element
			
			// Loop over all of the quadrature points for this element
			id_type nqp = el->n_q_points();
			for(id_type qp=0; qp<nqp; ++qp)
			{
				// Get the internal variable list associated with this quadrature point
				input.internal_vars = &(prob->get_internal_vars()->get_internal_vars_local(l_elem, qp));
				input.internal_vars_sensitivity = &(prob->get_sensitivity_solver()->get_internal_vars()->get_internal_vars_local(l_elem, qp));
				if (assem->storedBmats()) // Call the B-matrix variant
				{
					KernelVolumetricUpdate (mesh->get_shape(l_elem,qp), assem->getBmat(l_elem,qp),
											curr_mat, input,
											elem_U_curr, delem_U_curr);
				}
				else
				{
					KernelVolumetricUpdate (mesh->get_shape(l_elem,qp), mesh->get_shape_grad(l_elem,qp),
											curr_mat, input,
											elem_U_curr, delem_U_curr);
				}
			}
		}
	}
	
	// Otherwise, this element is intersected somehow
	else
	{
		// Loop over all of the integration elements for this intersected element
		id_type curr_qp = 0; // Can't just use the qp f the integration element to get the internal vars because I store them sequentially in the Problem object
		for(id_type ie=0; ie<el->n_integration_elem(); ++ie)
		{
			Material* curr_mat = mesh->get_element_material_global(el->get_id(), ie);
			bool mat_has_isvs = (curr_mat->n_internal_vars() != 0);

			// Only actually run the Kernel if there are ISVs to update
			id_type nqp = el->get_integration_elem(ie)->n_q_points();
			if (mat_has_isvs)
			{
				// Loop over the quadrature points
				for(id_type qp=0; qp<nqp; ++qp)
				{
					// Get the internal variable list associated with this quadrature point
					input.internal_vars = &(prob->get_internal_vars()->get_internal_vars_local(l_elem, curr_qp));
					input.internal_vars_sensitivity = &(prob->get_sensitivity_solver()->get_internal_vars()->get_internal_vars_local(l_elem, curr_qp));
					if (assem->storedBmats()) // Call the B-matrix variant
					{
						KernelVolumetricUpdate (mesh->get_shape(l_elem,curr_qp), assem->getBmat(l_elem,curr_qp),
												curr_mat, input,
												elem_U_curr, delem_U_curr);
					}
					else
					{
						KernelVolumetricUpdate (mesh->get_shape(l_elem,curr_qp), mesh->get_shape_grad(l_elem,curr_qp),
												curr_mat, input,
												elem_U_curr, delem_U_curr);
					}

					// Update the current qp
					curr_qp++;
				}
			}
			else
				curr_qp += nqp;
		}


		// If the problem is cohesive then I need to add the contribution from the cohesive opening
		if (mesh->is_cohesive())
		{
			// Get the cohesive element structure (Needed for mapping the cohesive matrix back to the elemental matrix)
			std::vector<std::vector<short_id_type> >& coh_elem_struct = el->getCohesiveElemStructure();

			// Loop over the cohesive elements
			for(id_type ce=0; ce<el->n_cohesive_elem(); ++ce)
			{
				// Get a pointer to the current cohesive element
				CohesiveElem* coh_el = el->get_cohesive_elem(ce);

				Material* curr_mat = el->get_cohesive_material(ce);
				bool mat_has_isvs = (curr_mat->n_internal_vars() != 0);

				id_type nqp = coh_el->n_q_points();
				if (mat_has_isvs)
				{
					// Get the nodal solution object for the current cohesive element
					std::vector<double> coh_U_curr(coh_el->n_nodes() * nndof);
					std::vector<double> dcoh_U_curr(coh_el->n_nodes() * nndof);
					for(id_type n=0; n<coh_el->n_nodes(); ++n)
					{
						id_type node = coh_elem_struct[ce][n];
						for(id_type d=0; d<nndof; ++d)
						{
							coh_U_curr[n*nndof + d] = elem_U_curr[node*nndof + d];
							dcoh_U_curr[n*nndof + d] = delem_U_curr[node*nndof + d];
						}
					}

					// Get the material of the current cohesive element directly from the element
					

					// Loop over the quadrature points
					for(id_type qp=0; qp<nqp; ++qp)
					{
						input.internal_vars = &(prob->get_internal_vars()->get_internal_vars_local(l_elem, curr_qp));
						input.internal_vars_sensitivity = &(prob->get_sensitivity_solver()->get_internal_vars()->get_internal_vars_local(l_elem, curr_qp));
						KernelCohesiveUpdate (mesh->get_shape(l_elem,curr_qp), mesh->get_rot_matrix(l_elem,curr_qp),
											  curr_mat, input,
											  coh_U_curr, dcoh_U_curr);

						// Update the quadrature point
						curr_qp++;
					}
				}
				else
					curr_qp += nqp;
			} // End cohesive elem loop
		} // End is coheisve if-statement
	} // End if intersected
}


























PetscErrorCode SensitivityParameter::assemble_dPdd_mat(Mat* dPfdd, Mat* dPpdd, std::vector<SensitivityParameter*> parameters, Problem* prob)
{
	return 0;
}
void SensitivityParameter::assemble_elem_dPdd_mat(Elem* el, DenseMatrix<double>& dP_el, std::vector<double>& elem_sol, Problem* prob)
{
}

void SensitivityParameter::updateSensitivityISVs_mat(Problem* prob, Mat* dUf_dd, Mat* dUp_dd)
{
	// // Some initializations
	// Mesh* mesh = prob->get_mesh();
	// InternalVars* int_vars = prob->get_internal_vars();
	// DofObject* dofs = prob->get_dofs();

	// std::vector<id_type> elem_dofs;
	// std::vector<double> elem_sol;
	// std::vector<std::vector<double> > delem_sol; // A n_parameters by ndof matrix (So each elemenal parameter derivative is contiguous in memory)
	// for(Mesh::element_iterator it=mesh->active_elements_begin(), end=mesh->active_elements_end(); it!=end; ++it)
	// {
	// 	if (int_vars->elemHasISVs(mesh->global_to_local_elem((*it)->get_id())))
	// 	{
	// 		// Assemble the current solution vector
	// 		dofs->get_elem_global_dofs((*it), elem_dofs);
	// 		ProblemUtilities::formElementalVector(elem_sol, (*it), prob->get_solution());
	// 		ProblemUtilities::formElementalMatrix(delem_sol, elem_dofs, dUf_dd, dUp_dd);

	// 		// Assemble the elemental contribution
	// 		updateSensitivityElemISVs_mat((*it), elem_sol, delem_sol, prob);
	// 	}
	// }
}
void SensitivityParameter::updateSensitivityElemISVs_mat(Elem* el, std::vector<double>& elem_sol, std::vector<std::vector<double> >& delem_sol, Problem* prob)
{
}










PetscErrorCode SensitivityParameter::set_dP(pvector& dPdd, std::vector<double>& vec, std::vector<id_type>& dofs, id_type ngfd)
{
	id_type nd = dofs.size();
	PetscErrorCode ierr(0);
	for (id_type i=0; i<nd; ++i)
	{
		if (vec[i]!=0.0)
		{
			if (dofs[i] < ngfd)
				{ierr = VecSetValue(dPdd.f, dofs[i], vec[i], ADD_VALUES);CHKERRQ(ierr);}
			else
				{ierr = VecSetValue(dPdd.p, dofs[i]-ngfd, vec[i], ADD_VALUES);CHKERRQ(ierr);}
		}
	}
	return ierr;
}