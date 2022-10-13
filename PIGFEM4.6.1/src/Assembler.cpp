/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#include "Assembler.h"
#include "Mesh.h"
#include "DofObject.h"
#include "BoundaryObject.h"
#include "Problem.h"
#include "NodalData.h"
#include "InternalVars.h"
#include "Solver.h"
#include "Utilities.h"
#include "BodyLoad.h"
#include "elem.h"
#include "gsl_cblas.h"





Assembler::~Assembler()
{
	_prob = NULL;
	_Bmats.clear();
}

void Assembler::clear()
{
	_prob = NULL;
	_Bmats.clear();
}






// ===========================================================================================================================
// LINEAR STIFFNESS MATRIX ASSEMBLY
// ===========================================================================================================================

PetscErrorCode Assembler::assemble_linear(pmatrix& K)
{
	// Some initializations
	PetscErrorCode ierr;
	Mesh* mesh = _prob->get_mesh();
	DofObject* dofs = _prob->get_dofs();
	
	// Assemble the stiffness matrix
	DenseMatrix<double> K_el;
	std::vector<id_type> elem_dofs;
	// double max = 0.0;
	for(Mesh::element_iterator it=mesh->active_elements_begin(), end=mesh->active_elements_end(); it!=end; ++it)
	{
		// Get the global dofs for this element
		dofs->get_elem_global_dofs((*it), elem_dofs); // elem_dofs will now contain a vector with all the global dofs associated with the element

		// Assemble the eleental stiffness matrices
		assemble_elem_linear((*it), K_el);

		set_matrix(K, K_el, elem_dofs);

		// for(id_type i=0; i<K_el.n_rows(); ++i)
		// 	for(id_type j=0; j<K_el.n_cols(); ++j)
		// 		max = std::max(max, std::abs(K_el(i,j)));
	}

	// // Apply the hanging node constraints using the penalty method
	// if (_prob->get_enforcement_method() == PENALTY)
	// 	apply_penalty_constraint(Kff, Kfp, Kpf, Kpp, max);
	
	
	ierr = MatAssemblyBegin(K.ff, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(K.fp, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(K.pf, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(K.pp, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K.ff, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K.fp, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K.pf, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(K.pp, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	return ierr;
}




void Assembler::assemble_elem_linear(Elem* el,
									 DenseMatrix<double>& K_el)
{
	Mesh* mesh = _prob->get_mesh();
	DofObject* dofs = _prob->get_dofs();

	// Define the initial input parameters
	Material::input_params input;
	input.dim = mesh->dim();
	input.delta_t = 0.0;
	if (_prob->get_classification() == STRUCTURAL)
	{
		std::vector<int> switch_dim = {1, 3, 6};
		input.strain.resize(switch_dim[input.dim - 1]); // I don't know if I need this here. Just a placeholder in case I don't wanna compute the actual strain
		if (_prob->get_parameter("plane_strain") == 0.0)
			input.plane_strain = false;
		else
			input.plane_strain = true;
	}

	id_type l_elem = mesh->global_to_local_elem(el->get_id());

	// Initialize the Stiffness matrix
	id_type nn = el->n_nodes() + el->n_enrich_nodes();
	id_type ndof = nn*dofs->nndof();
	K_el.clear();
	K_el.resize(ndof, ndof);
	DenseMatrix<double> K_temp = K_el;

	// If this is a nonintersected element
	if (!el->is_intersected())
	{		
		// Get a pointer to the material object associated with this non-intersected element
		Material* curr_mat = mesh->get_element_material_global(el->get_id());

		// Loop over all of the quadrature points for this element and update the stiffness matrix using the Kernel function for each point
		id_type nqp = el->n_q_points();
		for(id_type qp=0; qp<nqp; ++qp)
		{
			if (storedBmats())
			{
				KernelLinear (K_temp,
							  mesh->get_shape(l_elem, qp), getBmat(l_elem, qp),
							  curr_mat, input);
			}
			else
			{
				KernelLinear (K_temp,
							  mesh->get_shape(l_elem, qp), mesh->get_shape_grad(l_elem, qp),
							  curr_mat, input);
			}
			double& J = mesh->get_J(l_elem,qp);
			double& W = mesh->get_W(l_elem,qp);
			double wgt = J * W;
			K_el.AXPY(wgt, K_temp);
		}
	}
	
	// Otherwise, this element is intersected somehow
	else
	{
		id_type local_qp = 0;
		for(id_type ie=0; ie<el->n_integration_elem(); ++ie)
		{
			// Get the material of the current integration element using a mesh function
			Material* curr_mat = mesh->get_element_material_global(el->get_id(), ie); // Get the material pointer for the current integration element

			// Loop over the quadrature points
			id_type nqp = el->get_integration_elem(ie)->n_q_points();
			for(id_type qp=0; qp<nqp; ++qp)
			{
				if (storedBmats())
				{
					KernelLinear (K_temp,
								  mesh->get_shape(l_elem, qp), getBmat(l_elem, local_qp),
								  curr_mat, input);
				}
				else
				{
					KernelLinear (K_temp,
								  mesh->get_shape(l_elem, qp), mesh->get_shape_grad(l_elem, local_qp),
								  curr_mat, input);
				}
				double& J = mesh->get_J(l_elem,local_qp);
				double& W = mesh->get_W(l_elem,local_qp);
				double wgt = J * W;
				K_el.AXPY(wgt, K_temp);
				local_qp++;
			}
		}
	}
}


















// ===========================================================================================================================
// NONLINEAR STIFFNESS MATRIX AND INTERNAL LOAD ASSEMBLY
// ===========================================================================================================================

PetscErrorCode Assembler::assemble_nonlinear(pmatrix& K, pvector& P_int,
											 double delta_t, std::vector<std::vector<std::vector<double> > >& update_ISVs,
											 NodalData* solution, bool assembleFunc, bool assembleJac)
{
	if (!assembleFunc && !assembleJac)
		return 0;

	// Some initializations
	PetscErrorCode ierr;
	Mesh* mesh = _prob->get_mesh();
	DofObject* dofs = _prob->get_dofs();

	// Zero the entries in preparation for assembly
	if (assembleJac)
	{
		ierr = MatZeroEntries(K.ff);CHKERRQ(ierr);
		ierr = MatZeroEntries(K.fp);CHKERRQ(ierr);
		ierr = MatZeroEntries(K.pf);CHKERRQ(ierr);
		ierr = MatZeroEntries(K.pp);CHKERRQ(ierr);
	}
	if (assembleFunc)
	{
		ierr = VecZeroEntries(P_int.f);CHKERRQ(ierr);
		ierr = VecZeroEntries(P_int.p);CHKERRQ(ierr);
	}

	// Make sure the internal variables we are gonna update are the correct ones that correspond to the begining of the load step
	_prob->get_internal_vars()->copy_data(update_ISVs);

	// Assemble the stiffness matrix
	DenseMatrix<double> K_el;
	std::vector<double> P_el_int;
	std::vector<id_type> elem_dofs;
	std::vector<double> elem_sol;
	// double max = 0.0;


	for(Mesh::element_iterator it=mesh->active_elements_begin(), end=mesh->active_elements_end(); it!=end; ++it)
	{
		// Get the local element number
		id_type l_elem = mesh->global_to_local_elem((*it)->get_id());

		// Get the global dofs for this element
		dofs->get_elem_global_dofs((*it), elem_dofs); // elem_dofs will now contain a vector with all the global dofs associated with the element





		// Assemble the current displacement vector
		ProblemUtilities::formElementalVector(elem_sol, (*it), solution);
		
		// Assemble the elemental stiffness matrix
		assemble_elem_nonlinear((*it), K_el, P_el_int, delta_t, elem_sol, update_ISVs[l_elem], assembleFunc, assembleJac);

		if (assembleJac)
			set_matrix(K, K_el, elem_dofs);
		if (assembleFunc)
			set_vec(P_int, P_el_int, elem_dofs);

		// for(id_type i=0; i<K_el.n_rows(); ++i)
		// 	for(id_type j=0; j<K_el.n_cols(); ++j)
		// 		max = std::max(max, std::abs(K_el(i,j)));
	}

	// // Apply the hanging node constraints using the penalty method
	// if (_prob->get_enforcement_method() == PENALTY)
	// 	apply_penalty_constraint(Kff, Kfp, Kpf, Kpp, max);
	
	if (assembleJac)
	{
		ierr = MatAssemblyBegin(K.ff, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyBegin(K.fp, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyBegin(K.pf, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyBegin(K.pp, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}
	if (assembleFunc)
	{
		ierr = VecAssemblyBegin(P_int.f);CHKERRQ(ierr);
		ierr = VecAssemblyBegin(P_int.p);CHKERRQ(ierr);
	}
	if (assembleJac)
	{
		ierr = MatAssemblyEnd(K.ff, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(K.fp, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(K.pf, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(K.pp, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	}
	if (assembleFunc)
	{
		ierr = VecAssemblyEnd(P_int.f);CHKERRQ(ierr);
		ierr = VecAssemblyEnd(P_int.p);CHKERRQ(ierr);
	}
	
	return ierr;
}

void Assembler::assemble_elem_nonlinear(Elem* el,
									 	DenseMatrix<double>& K_el, std::vector<double>& P_el_int,
									 	double delta_t, std::vector<double>& elem_U_curr, std::vector<std::vector<double> >& init_elem_internal_vars,
									 	bool assembleFunc, bool assembleJac)
{
	// Get the mesh and dof objects
	Mesh* mesh = _prob->get_mesh();
	id_type nndof = _prob->get_dofs()->nndof();

	// Get the local element number
	id_type l_elem = mesh->global_to_local_elem(el->get_id());


	// Define the initial input parameters
	Material::input_params input;
	input.dim = mesh->dim();
	input.delta_t = delta_t;
	if (_prob->get_classification() == STRUCTURAL)
	{
		std::vector<int> switch_dim = {1, 3, 6};
		input.strain.resize(switch_dim[input.dim - 1]); // I don't know if I need this here. Just a placeholder in case I don't wanna compute the actual strain
		if (_prob->get_parameter("plane_strain") == 0.0)
			input.plane_strain = false;
		else
			input.plane_strain = true;
	}



	// If this is a non-intersected element
	if (!el->is_intersected())
	{
		// Initialize the Stiffness matrix and the internal force vector
		id_type nn = el->n_nodes();
		id_type ndof = nn*nndof;
		K_el.clear();
		P_el_int.clear();
		K_el.resize(ndof, ndof);
		P_el_int.resize(ndof);
		DenseMatrix<double> K_temp = K_el;
		std::vector<double> P_int_temp = P_el_int;
		
		// Get a pointer to the material object associated with this non-intersected element
		Material* curr_mat = mesh->get_element_material_global(el->get_id());

		// Loop over all of the quadrature points for this element
		id_type nqp = el->n_q_points();
		for(id_type qp=0; qp<nqp; ++qp)
		{
			// Get the internal variable list associated with this quadrature point
			input.internal_vars = &(init_elem_internal_vars[qp]);

			if (storedBmats())
			{
				KernelNonlinear (K_temp, P_int_temp,
								 mesh->get_shape(l_elem,qp), getBmat(l_elem,qp),
								 curr_mat, input,
								 elem_U_curr, assembleFunc, assembleJac);
			}
			else
			{
				KernelNonlinear (K_temp, P_int_temp,
								 mesh->get_shape(l_elem,qp), mesh->get_shape_grad(l_elem,qp),
								 curr_mat, input,
								 elem_U_curr, assembleFunc, assembleJac);
			}
			double& J = mesh->get_J(l_elem,qp);
			double& W = mesh->get_W(l_elem,qp);
			double wgt = J * W;
			if (assembleJac)
				K_el.AXPY(wgt, K_temp);
			if (assembleFunc)
				cblas_daxpy(P_el_int.size(), wgt, P_int_temp.data(), 1, P_el_int.data(), 1);
		}
	}
	
	// Otherwise, this element is intersected somehow
	else
	{

		// Initialize the Stiffness matrix and the internal and external force vectors
		id_type nn = el->n_nodes() + el->n_enrich_nodes();
		id_type ndof = nn*nndof;
		K_el.clear();
		P_el_int.clear();
		K_el.resize(ndof, ndof);
		P_el_int.resize(ndof);
		DenseMatrix<double> K_temp = K_el;
		std::vector<double> P_int_temp = P_el_int;

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
				// Get the internal variable list associated with this quadrature point
				input.internal_vars = &(init_elem_internal_vars[curr_qp]);

				if (storedBmats())
				{
					KernelNonlinear (K_temp, P_int_temp,
									 mesh->get_shape(l_elem,qp), getBmat(l_elem,curr_qp),
									 curr_mat, input,
									 elem_U_curr, assembleFunc, assembleJac);
				}
				else
				{
					KernelNonlinear (K_temp, P_int_temp,
									 mesh->get_shape(l_elem,qp), mesh->get_shape_grad(l_elem,curr_qp),
									 curr_mat, input,
									 elem_U_curr, assembleFunc, assembleJac);
				}
				double& J = mesh->get_J(l_elem,curr_qp);
				double& W = mesh->get_W(l_elem,curr_qp);
				double wgt = J * W;
				if (assembleJac)
					K_el.AXPY(wgt, K_temp);
				if (assembleFunc)
					cblas_daxpy(P_el_int.size(), wgt, P_int_temp.data(), 1, P_el_int.data(), 1);
				// Update the current qp
				curr_qp++;
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
				CohesiveElem* coh_el = el->get_cohesive_elem(ce);//------------------------------I
				// Initialize the cohesive contributions
				id_type n_dof_coh = coh_el->n_nodes() * nndof;
				DenseMatrix<double> K_coh(n_dof_coh, n_dof_coh);
				std::vector<double> P_int_coh(n_dof_coh);
				DenseMatrix<double> K_coh_temp = K_coh;
				std::vector<double> P_int_coh_temp = P_int_coh;

				// Get the nodal solution object for the current cohesive element
				std::vector<double> coh_U_curr(coh_el->n_nodes() * nndof); // The vectorized form of the previous structure

				for(id_type n=0; n<coh_el->n_nodes(); ++n) //------------------------------------I
				{
					id_type node = coh_elem_struct[ce][n];
					for(id_type d=0; d<nndof; ++d)
						coh_U_curr[n*nndof + d] = elem_U_curr[node*nndof + d];
				}

				// Get the material of the current cohesive element directly from the element
				Material* curr_mat = el->get_cohesive_material(ce);

				// Loop over the quadrature points
				id_type nqp = coh_el->n_q_points();
				for(id_type qp=0; qp<nqp; ++qp)
				{
					input.internal_vars = &(init_elem_internal_vars[curr_qp]);

					KernelCohesive (K_coh_temp, P_int_coh_temp,
									mesh->get_shape(l_elem,curr_qp), mesh->get_rot_matrix(l_elem,curr_qp),
									curr_mat, input,
									coh_U_curr, assembleFunc, assembleJac);

					double& J = mesh->get_J(l_elem,curr_qp);
					double& W = mesh->get_W(l_elem,curr_qp);
					double wgt = J * W;
					// K_coh.print();
					if (assembleJac)
						K_coh.AXPY(wgt, K_coh_temp);
					if (assembleFunc)
						cblas_daxpy(P_int_coh.size(), wgt, P_int_coh_temp.data(), 1, P_int_coh.data(), 1);
					// K_coh.print();
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
						if (assembleJac)
						{
							for(id_type n2=0; n2<coh_elem_struct[ce].size(); ++n2)
							{
								id_type col_node = coh_elem_struct[ce][n2];
								for(id_type d2=0; d2<nndof; ++d2)
								{
									id_type col_dof = col_node*nndof + d2;
									id_type col_dof_local = n2*nndof + d2;
									
									// Stiffness matrix contribution
									K_el(row_dof, col_dof) += K_coh(row_dof_local, col_dof_local);
								}
							}
						}

						// Internal force contribution
						if (assembleFunc)
							P_el_int[row_dof] += P_int_coh[row_dof_local];
					}
				} // End cohesive contribution update
			} // End cohesive elem loop
		} // End is coheisve if-statement

	} // End if intersected
}






















// Function used to update the prescribed displacement and external force vector
// This function will be called at the initialization of each new load stop in a nonlinear analysis
// Up is the prescribed dispacement vector, Pf_external is the external load applied to free dofs, and current_time_frac is (t_curr/t_end)
PetscErrorCode Assembler::assemble_new_load_step(Vec& Up, pvector& F_ext, double current_time_frac)
{
	// Some initiailizations
	PetscErrorCode ierr;
	BoundaryObject* boundary = _prob->get_boundary();
	DofObject* dofs = _prob->get_dofs();
	Mesh* mesh = _prob->get_mesh();

	// Assemble the prescribed displacement vector

	// Monotonic Loading
	std::vector<PetscInt> rows;
	std::vector<PetscScalar> vals;
	for(BoundaryObject::dirichlet_iterator it=boundary->dirichlet_begin(), end=boundary->dirichlet_end(); it!=end; ++it)
	{
		id_type g_node = (*it).first;
		for(std::map<id_type, double>::iterator it2=it->second.begin(), end2=it->second.end(); it2!=end2; ++it2)
		{
			id_type dof = (*it2).first;
			PetscScalar val = current_time_frac * (*it2).second; // Only apply part of this displacements
			PetscInt row = dofs->get_global_dof(g_node, dof) - dofs->n_global_free_dofs();  // get the local row number since this a constrained dof
			rows.push_back(row);
			vals.push_back(val);
		}
	}
	ierr = VecSetValues(Up, rows.size(), &rows[0], &vals[0], INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(Up);CHKERRQ(ierr);

	// //  Non-Monotonic Loading
	// std::vector<PetscInt> rows;
	// std::vector<PetscScalar> vals;
	// for(BoundaryObject::dirichlet_iterator it=boundary->dirichlet_begin(), end=boundary->dirichlet_end(); it!=end; ++it)
	// {
	// 	id_type g_node = (*it).first;
	// 	for(std::map<id_type, double>::iterator it2=it->second.begin(), end2=it->second.end(); it2!=end2; ++it2)
	// 	{
	// 		id_type dof = (*it2).first;
	// 		PetscScalar val;
	// 		if (current_time_frac < (1.0/3.0)) // Non-Monotonic Loading
	// 			val = current_time_frac * (*it2).second; // Only apply part of this displacements
	// 		else if (current_time_frac < (2.0/3.0))
	// 			val = 2*(0.5-current_time_frac) * (*it2).second;
	// 		else
	// 			val = 4 * (current_time_frac - 0.75) * (*it2).second;
	// 		PetscInt row = dofs->get_global_dof(g_node, dof) - dofs->n_global_free_dofs();  // get the local row number since this a constrained dof
	// 		rows.push_back(row);
	// 		vals.push_back(val);
	// 	}
	// }
	// ierr = VecSetValues(Up, rows.size(), &rows[0], &vals[0], INSERT_VALUES);CHKERRQ(ierr);
	// ierr = VecAssemblyBegin(Up);CHKERRQ(ierr);

	// //  Non-Monotonic Loading
	// std::vector<PetscInt> rows;
	// std::vector<PetscScalar> vals;
	// for(BoundaryObject::dirichlet_iterator it=boundary->dirichlet_begin(), end=boundary->dirichlet_end(); it!=end; ++it)
	// {
	// 	id_type g_node = (*it).first;
	// 	for(std::map<id_type, double>::iterator it2=it->second.begin(), end2=it->second.end(); it2!=end2; ++it2)
	// 	{
	// 		id_type dof = (*it2).first;
	// 		PetscScalar val;
	// 		if (current_time_frac < (1.0/3.0)) // Non-Monotonic Loading
	// 			val = 1.5 * current_time_frac * (*it2).second; // Only apply part of this displacements
	// 		else if (current_time_frac < (2.0/3.0))
	// 			val = (1.0-1.5*current_time_frac) * (*it2).second;
	// 		else
	// 			val = (3.0*current_time_frac - 2.0)* (*it2).second;
	// 		PetscInt row = dofs->get_global_dof(g_node, dof) - dofs->n_global_free_dofs();  // get the local row number since this a constrained dof
	// 		rows.push_back(row);
	// 		vals.push_back(val);
	// 	}
	// }
	// ierr = VecSetValues(Up, rows.size(), &rows[0], &vals[0], INSERT_VALUES);CHKERRQ(ierr);
	// ierr = VecAssemblyBegin(Up);CHKERRQ(ierr);



	// Loop through body forces and assemble the body force terms
	std::vector<double> P_el_ext;
	std::vector<id_type> elem_dofs;
	for (Problem::body_load_iterator it=_prob->body_loads_begin(), end=_prob->body_loads_end(); it!=end; ++it)
	{
		// Loop through all elements and calculate the contributions to the body force
		for (Mesh::element_iterator el_it=mesh->active_elements_begin(), el_end=mesh->active_elements_end(); el_it!=el_end; ++el_it)
		{
			// Get the global dofs for this element
			std::vector<id_type> elem_dofs;
			dofs->get_elem_global_dofs((*el_it), elem_dofs);

			// Do the elemental assembly
			assemble_elem_body_load(*el_it, *it, P_el_ext);

			// Set the time fraction
			for (id_type i=0; i<P_el_ext.size(); ++i)
				P_el_ext[i] *= current_time_frac;

			// Add it to PETSc
			set_vec(F_ext, P_el_ext, elem_dofs);
		}
	}


	// Assemble the external force (I wonder what will happen if I just don't asemble Pp_ext...)
	ierr = VecAssemblyBegin(F_ext.f);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(Up);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(F_ext.f);CHKERRQ(ierr);
	
	return ierr;
}



void Assembler::assemble_elem_body_load(Elem* el, BodyLoad* load, std::vector<double>& P_el_ext)
{
	Mesh* mesh = _prob->get_mesh();

	// Define the initial input parameters
	Material::input_params input;
	input.dim = mesh->dim();
	input.delta_t = 0.0;
	if (_prob->get_classification() == STRUCTURAL)
	{
		std::vector<int> switch_dim = {1, 3, 6};
		input.strain.resize(switch_dim[input.dim - 1]); // I don't know if I need this here. Just a placeholder in case I don't wanna compute the actual strain
		if (_prob->get_parameter("plane_strain") == 0.0)
			input.plane_strain = false;
		else
			input.plane_strain = true;
	}

	// Get the local element number
	id_type l_elem = mesh->global_to_local_elem(el->get_id());

	// Set up the load vector
	P_el_ext.clear();
	P_el_ext.resize((el->n_nodes() + el->n_enrich_nodes()) * _prob->nndof());
	std::vector<double> P_el_temp = P_el_ext;

	// If this is a nonintersected element
	if (!el->is_intersected())
	{	
		// Get a pointer to the material object associated with this non-intersected element
		Material* curr_mat = mesh->get_element_material_global(el->get_id());

		// Loop over all of the quadrature points for this element and update the stiffness matrix using the Kernel function for each point
		id_type nqp = el->n_q_points();
		for(id_type qp=0; qp<nqp; ++qp)
		{
			std::vector<double> gcoords = mesh->get_global_coords_undeformed(l_elem, qp);
			load->BF_Kernel(P_el_temp, gcoords,
							mesh->get_shape(l_elem, qp), mesh->get_shape_grad(l_elem, qp),
							mesh->get_J(l_elem, qp), mesh->get_W(l_elem, qp),
							curr_mat, input);
			for (id_type i=0; i<P_el_temp.size(); ++i)
				P_el_ext[i] += P_el_temp[i];
		}
	}
	
	// Otherwise, this element is intersected somehow
	else
	{
		id_type local_qp = 0;
		for(id_type ie=0; ie<el->n_integration_elem(); ++ie)
		{
			// Get the material of the current integration element using a mesh function
			Material* curr_mat = mesh->get_element_material_global(el->get_id(), ie); // Get the material pointer for the current integration element

			// Loop over the quadrature points
			id_type nqp = el->get_integration_elem(ie)->n_q_points();
			for(id_type qp=0; qp<nqp; ++qp)
			{
				std::vector<double> gcoords = mesh->get_global_coords_undeformed(l_elem, local_qp);
				load->BF_Kernel(P_el_temp, gcoords,
								mesh->get_shape(l_elem, local_qp), mesh->get_shape_grad(l_elem, local_qp),
								mesh->get_J(l_elem, local_qp), mesh->get_W(l_elem, local_qp),
								curr_mat, input);
				for (id_type i=0; i<P_el_temp.size(); ++i)
					P_el_ext[i] += P_el_temp[i];
				local_qp++;
			}
		}
	}
}













// ===========================================================================================================================
// STUFF HAVING TO DO WITH STORING THE B MATRICIES FOR SPECIFIC CLASSES OF PROBLEMS FOR SPEED
// ===========================================================================================================================
void Assembler::storeBmats()
{
	_Bmats.clear();
	if (storedBmats())
	{
		Mesh* mesh = _prob->get_mesh();
		_Bmats.resize( mesh->n_local_active_elem() );
		for (Mesh::element_iterator it=mesh->active_elements_begin(), end=mesh->active_elements_end(); it!=end; ++it)
		{
			id_type l_elem = mesh->global_to_local_elem((*it)->get_id());
			id_type nqp = mesh->n_volumetric_qps(l_elem);
			_Bmats[l_elem].resize( nqp );
			for (id_type qp=0; qp<nqp; ++qp)
				fillBmat(_Bmats[l_elem][qp], mesh->get_shape_grad(l_elem, qp));

			// Could store cohesive "Bmats" (Really Nmats) here too
		}
	}
	// Otherwise do nothing
}
DenseMatrix<double>& Assembler::getBmat(id_type local_e, id_type qp)
{
	if (storedBmats())
	{
		if (local_e < _Bmats.size())
		{
			if (qp < _Bmats[local_e].size())
				return _Bmats[local_e][qp];
			else
				err_message("Invalid Quadrature point.");
		}
		else
			err_message("Invalid element number.");
	}
	else
		err_message("Called getBmat for an assembler that does not support B matrix storage");
}















// ===========================================================================================================================
// STUFF HAVING TO DO WITH PARTITIONING THE ELEMNTAL VECTORS INTO FREE AND PRESCRIBED PORTIONS
// ===========================================================================================================================

// Sets the values for a matrix
PetscErrorCode Assembler::set_matrix(pmatrix& K, DenseMatrix<double>& mat, std::vector<id_type>& dofs)
{
	id_type nd = dofs.size();
	id_type ngfd = _prob->get_dofs()->n_global_free_dofs();
	PetscErrorCode ierr(0);
	for (id_type i=0; i<nd; ++i)
	{
		for (id_type j=0; j<nd; ++j)
		{
			if (dofs[i] < ngfd && dofs[j] < ngfd)
			{ierr = MatSetValue(K.ff, dofs[i], dofs[j], mat(i,j), ADD_VALUES);CHKERRQ(ierr);}
			else if (dofs[i] < ngfd && dofs[j] >= ngfd)
			{ierr = MatSetValue(K.fp, dofs[i], (dofs[j]-ngfd), mat(i,j), ADD_VALUES);CHKERRQ(ierr);}
			else if (dofs[i] >= ngfd && dofs[j] < ngfd)
			{ierr = MatSetValue(K.pf, (dofs[i]-ngfd), dofs[j], mat(i,j), ADD_VALUES);CHKERRQ(ierr);}
			else
			{ierr = MatSetValue(K.pp, (dofs[i]-ngfd), (dofs[j]-ngfd), mat(i,j), ADD_VALUES);CHKERRQ(ierr);}
		}
	}
	return ierr;
}
// Sets the values for a vector
PetscErrorCode Assembler::set_vec(pvector& P, std::vector<double>& vec, std::vector<id_type>& dofs)
{
	id_type nd = dofs.size();
	id_type ngfd = _prob->get_dofs()->n_global_free_dofs();
	PetscErrorCode ierr(0);
	for (id_type i=0; i<nd; ++i)
	{
		if (vec[i]!=0.0)
		{
			// if (dofs[i] < ngfd)
			// 	{ierr = VecSetValue(*P.f, dofs[i], vec[i], ADD_VALUES);CHKERRQ(ierr);}
			// else
			// 	{ierr = VecSetValue(*P.p, dofs[i]-ngfd, vec[i], ADD_VALUES);CHKERRQ(ierr);}
			if (dofs[i] < ngfd)
				{ierr = VecSetValue(P.f, dofs[i], vec[i], ADD_VALUES);CHKERRQ(ierr);}
			else
				{ierr = VecSetValue(P.p, dofs[i]-ngfd, vec[i], ADD_VALUES);CHKERRQ(ierr);}
		}
	}
	return ierr;
}




















// ===========================================================================================================================
// STUFF HAVING TO DO WITH HANGING NODE CONSTRAINTS
// ===========================================================================================================================


// Assembles the matrix used for the enforcement of hanging node constraints based on the penalty method
void Assembler::assemble_penalty_matrix(DenseMatrix<double>& P, std::vector<double> constraint_vals)
{
	id_type nn = constraint_vals.size();
	if (nn <= 2)
		err_message("A hanging node must have at least 2 neighbors");

	P.resize(nn, nn);
	for (id_type i=0; i<nn; ++i)
		for (id_type j=0; j<nn; ++j)
			P(i,j) = constraint_vals[i] * constraint_vals[j];
}
void Assembler::apply_penalty_constraint(pmatrix& K, double max)
{
	Mesh* mesh = _prob->get_mesh();
	DofObject* dofs = _prob->get_dofs();
	DenseMatrix<double> P;
	double weight = pow(10, 8+std::floor((std::log(max)/std::log(10.0)))); // penalty weight (follows "sqrt rule" for error reduction)
	for (Mesh::hanging_node_iterator it=mesh->hanging_nodes_begin(), end=mesh->hanging_nodes_end(); it!=end; ++it)
	{
		id_type node = it->first;
		std::vector<id_type> constraint_nodes(it->second.size());
		std::vector<double> constraint_vals(it->second.size());
		for (id_type n=0; n<constraint_nodes.size(); ++n)
		{
			constraint_nodes[n] = it->second[n].first;
			constraint_vals[n] = it->second[n].second;
		}
		constraint_nodes.push_back( node );
		constraint_vals.push_back( -1.0 );
		assemble_penalty_matrix(P, constraint_vals);
		P = weight * P;

		for (id_type d=0; d<_prob->nndof(); ++d)
		{
			// Get the vector of global dofs
			std::vector<id_type> dofs_vec(constraint_nodes.size());
			for (id_type n=0; n<constraint_nodes.size(); ++n)
				dofs_vec[n] = dofs->get_global_dof(constraint_nodes[n], d);

			set_matrix(K, P, dofs_vec);
		}
	}
}

