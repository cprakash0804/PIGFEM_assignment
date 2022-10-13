/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#include "Solver.h"
#include "Mesh.h"
#include "DofObject.h"
#include "Problem.h"
#include "NodalData.h"
#include "node.h"
#include "elem.h"
#include "Utilities.h"
#include "Writer.h"
#include "SubscaleModel.h"
#include "Assembler.h"
#include "Options.h"
#include <algorithm>



Solver::Solver()
	: _rel_tol(1e-8), _abs_tol(1e-8), _T_final(1.0),
	  _max_time_step(1e-2), _min_time_step(1e-5), _max_iter(25), _setup(false), _nonzero_start(false)
{
}

Solver::~Solver()
{
	// Destroy the index sets
	ISDestroy(&_global_f);
	ISDestroy(&_to_f);
	ISDestroy(&_global_p);
	ISDestroy(&_to_p);
}



PetscErrorCode Solver::setup()
{
	PetscErrorCode ierr;
	ierr = preallocate_matrix(_K);CHKERRQ(ierr);
	ierr = preallocate_vectors();CHKERRQ(ierr);
	ierr = preallocate_index_sets();CHKERRQ(ierr);
	_setup = true; // This FE system has been set up
	return ierr;
}




PetscErrorCode Solver::initKSP(KSP& ksp, bool create)
{
	PetscErrorCode ierr;
	Mesh* mesh = _prob->get_mesh();

	// Create the KSP context and set operators (If I'm making this for the first time)
	if (create)
	{
		ierr = KSPCreate(mesh->get_comm(),&ksp);CHKERRQ(ierr);
		ierr = KSPSetOperators(ksp, _K.ff, _K.ff);CHKERRQ(ierr); // Set the preconditioner to be the same as the matrix. This can probably change to make it run faster
	}
	
	// Pick a preconditioner
	PC pc;
	ierr = KSPGetPC(ksp, &pc);CHKERRQ(ierr);
	{	// Direct Solver
		// ierr = PCSetType(pc, PCLU);CHKERRQ(ierr);
		ierr = PCSetType(pc, PCCHOLESKY);CHKERRQ(ierr);
		// ierr = PCFactorSetMatSolverPackage(pc, MATSOLVERSUPERLU_DIST);CHKERRQ(ierr);
		ierr = PCFactorSetMatSolverPackage(pc, MATSOLVERMUMPS);CHKERRQ(ierr); // MUMPS seems to perform about 4-5x better than SuperLU

		// Pick a Krylov method (Don't need a Krylov method for LU)
		ierr = KSPSetType(ksp, KSPPREONLY);
	}
	{
		// // ierr = PCSetType(pc, PCJACOBI);CHKERRQ(ierr);
		// ierr = PCSetType(pc, PCSOR);CHKERRQ(ierr);

		// // Pick a Krylov method
		// // ierr = KSPSetType(ksp, KSPGMRES); CHKERRQ(ierr);
		// ierr = KSPSetType(ksp, KSPCG);CHKERRQ(ierr);
		// // ierr = KSPCGSetType(ksp, KSP_CG_SYMMETRIC);CHKERRQ(ierr);

		// // CG tolerences
		// ierr = KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);CHKERRQ(ierr);
	}
	{
		// ierr = PCSetType(pc, PCGAMG);CHKERRQ(ierr);
		// id_type count = 0;
		// id_type dim = mesh->dim();
		// id_type n_nodes = mesh->n_local_owned_nodes() + mesh->n_local_owned_enrich_nodes();
		// std::vector<PetscReal> node_coords(n_nodes * dim);
		// for (Mesh::node_iterator it=mesh->nodes_begin(), end=mesh->nodes_end(); it!=end; ++it)
		// {
		// 	if (mesh->own_node_global((*it)->get_id()))
		// 	{
		// 		std::vector<double> node = (*it)->get_coords();
		// 		for (id_type d=0; d<dim; ++d)
		// 			node_coords[count*dim+d] = node[d];
		// 		count++;
		// 	}
		// }
		// for (Mesh::enrich_node_iterator it=mesh->enrich_nodes_begin(), end=mesh->enrich_nodes_end(); it!=end; ++it)
		// {
		// 	if (mesh->own_node_global((*it)->get_id()))
		// 	{
		// 		std::vector<double> node = (*it)->get_coords();
		// 		for (id_type d=0; d<dim; ++d)
		// 			node_coords[count*dim+d] = node[d];
		// 		count++;
		// 	}
		// }
		// ierr = PCSetCoordinates(pc, dim, n_nodes, node_coords.data());

		// // Pick a Krylov method
		// ierr = KSPSetType(ksp, KSPGMRES); CHKERRQ(ierr);
		// ierr = KSPSetType(ksp, KSPCG);CHKERRQ(ierr);
		// ierr = KSPCGSetType(ksp, KSP_CG_SYMMETRIC);CHKERRQ(ierr);

		// // CG tolerences
		// ierr = KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);CHKERRQ(ierr);
	}
	{
		// ierr = PCSetType(pc, PCHYPRE);CHKERRQ(ierr);
		// ierr = PCHYPRESetType(pc, "boomeramg");

		// // Pick a Krylov method
		// ierr = KSPSetType(ksp, KSPGMRES); CHKERRQ(ierr);
		// ierr = KSPSetType(ksp, KSPCG);CHKERRQ(ierr);
		// ierr = KSPCGSetType(ksp, KSP_CG_SYMMETRIC);CHKERRQ(ierr);

		// // CG tolerences
		// ierr = KSPSetTolerances(ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);CHKERRQ(ierr);
	}

	return ierr;
}




void Solver::Output(id_type step_iter, double current_t)
{
	Utilities::timer timer;
	timer.start();

	// Convert the step iteration to a string
	std::ostringstream convert;
	convert << step_iter;
	std::string step_out = convert.str();

	std::string out_folder = _prob->getOptions()->getOption("-OUT_FOLDER");
	// Writing VTK files
	if (_prob->getOptions()->getBoolOption("-WRITE_VTK"))
	{
		std::string history_file = out_folder + _prob->getOptions()->getOption("-VTK_FILE");
		std::string filename = history_file + "." + step_out + ".vtk";
		_prob->get_writer()->write(filename, current_t);

		// If the mesh is cohesive and we're supposed to write cohesive VTKs, write them
		if (_prob->get_mesh()->is_cohesive() && _prob->getOptions()->getBoolOption("-OUTPUT_COH"))
		{
			filename = history_file + "_coh." + step_out + ".vtk";
			_prob->get_cohesive_writer()->write(filename, current_t);
		}
	}

	// If this is a subscale model then I want to do the homogenization
	if (_prob->subscale())
	{
		SubscaleModel* model = _prob->get_subscale_model();
		std::vector<double> current_macro_strain = _prob->get_assembler()->get_vec_parameter("current strain");
		std::vector<double> stress, perturbation_strain;
		model->homogenize_stress_strain(current_macro_strain, stress, perturbation_strain);
		std::string filename = out_folder + _prob->getOptions()->getOption("-SS_FILE");
		model->write_stress_strain(filename, step_iter, current_macro_strain, stress, perturbation_strain);
	}
	
	// If this problem involves sensitivity, write the current sensitivities
	if (_prob->sensitivity())
	{
		std::string filename = out_folder + _prob->getOptions()->getOption("-SENS_FILE");
		if (!_prob->getOptions()->getBoolOption("-OUTPUT_FINAL_ONLY")) // output like normal
		{
			if (!_prob->getOptions()->getBoolOption("-SENSITIVITY_FINAL_ONLY"))
				_prob->getSensitivityWriter()->writeConsecutive(filename, current_t, step_iter);
			else if (current_t >= _T_final)
				_prob->getSensitivityWriter()->writeConsecutive(filename, current_t, 0);
		}
		else // This is the first time this file is being written to
			_prob->getSensitivityWriter()->writeConsecutive(filename, current_t, 0);
	}

	// If this problem is cohesive and the flag to output the cohesive failure lengths was set, output them
	if (_prob->get_mesh()->is_cohesive() && _prob->getOptions()->getBoolOption("-OUTPUT_COH_FAIL"))
	{
		std::string filename = out_folder + _prob->getOptions()->getOption("-COH_FAIL_FILE");
		_prob->getCohesiveFailureWriter()->writeConsecutive(filename, current_t, step_iter);
	}

	// If we are outputting the solution along a path, write it here
	if (_prob->getOptions()->getBoolOption("-WRITE_PATH"))
	{
		std::string filename = out_folder + _prob->getOptions()->getOption("-PATH_OUT") + "." + step_out + ".csv";
		_prob->getPathWriter()->write(filename, current_t);
	}

	// If we are outputting the cohesive displacement jumps
	if (_prob->getOptions()->getBoolOption("-WRITE_COH_OPENINGS"))
	{
		std::string filename = out_folder + _prob->getOptions()->getOption("-COH_OPENINGS_OUT") + "." + step_out + ".csv";
		_prob->getCohOpeningsWriter()->write(filename, current_t);
	}

	_prob->_solver_write_time += timer.time_elapsed();
}






void Solver::set_dparameter(std::string name, double val)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if (name=="REL_TOL")
		_rel_tol = val;
	else if (name=="ABS_TOL")
		_abs_tol = val;
	else if (name=="MAX_TIME_STEP")
		_max_time_step = val;
	else if (name=="MIN_TIME_STEP")
		_min_time_step = val;
	else if (name=="T_FINAL")
		_T_final = val;
	else
		err_message("Invalid solver parameter!");
}
double Solver::get_dparameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if (name=="REL_TOL")
		return _rel_tol;
	else if (name=="ABS_TOL")
		return _abs_tol;
	else if (name=="MAX_TIME_STEP")
		return _max_time_step;
	else if (name=="MIN_TIME_STEP")
		return _min_time_step;
	else if (name=="T_FINAL")
		return _T_final;
	else
		err_message("Invalid solver parameter!");
}
void Solver::set_iparameter(std::string name, id_type val)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if (name=="MAX_ITER" || name=="MAX_ITERATIONS")
		_max_iter = val;
	else
		err_message("Invalid solver parameter!");
}
id_type Solver::get_iparameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if (name=="MAX_ITER" || name=="MAX_ITERATIONS")
		return _max_iter;
	else
		err_message("Invalid solver parameter!");
}




// Creates the index sets necessary for the VecScatters to store the solution on each local process
// These index sets are used in every call to store_solution and this function is called in the solve() function
PetscErrorCode Solver::preallocate_index_sets()
{
	PetscErrorCode ierr;
	if (_setup)
	{
		ierr = ISDestroy(&_global_f);CHKERRQ(ierr);
		ierr = ISDestroy(&_to_f);CHKERRQ(ierr);
		ierr = ISDestroy(&_global_p);CHKERRQ(ierr);
		ierr = ISDestroy(&_to_p);CHKERRQ(ierr);
	}

	// Get the mesh and dofs
	Mesh* mesh = _prob->get_mesh();
	DofObject* dofs = _prob->get_dofs();

	// Get the bounds of the to_vector
	id_type n_local_dofs = dofs->n_local_dofs();
	MPI_Scan(&n_local_dofs, &_r_start, 1, MPI_ID, MPI_SUM, mesh->get_comm());
	_r_start -= n_local_dofs;

	// Create the index sets
	PetscInt n_f=dofs->n_local_free_dofs(), n_p=dofs->n_local_const_dofs();
	PetscInt IS_gf[n_f], IS_tf[n_f], IS_gp[n_p], IS_tp[n_p];
	PetscInt curr_f=0, curr_p=0;
	for(Mesh::node_iterator it=mesh->nodes_begin(), end=mesh->nodes_end(); it!=end; ++it)
	{
		id_type l_node = mesh->global_to_local_node((*it)->get_id());
		std::vector<id_type>& gdofs = dofs->get_nodal_global_dofs_local(l_node);
		std::vector<id_type>& ldofs = dofs->get_nodal_local_dofs_local(l_node);
		for(id_type d=0; d<gdofs.size(); ++d)
		{
			if(gdofs[d] < dofs->n_global_free_dofs()) // This is a free dof
			{
				IS_tf[curr_f] = _r_start + ldofs[d];
				IS_gf[curr_f] = gdofs[d];
				curr_f++;
			}
			else // This is a constrained dof
			{
				IS_tp[curr_p] = _r_start + ldofs[d];
				IS_gp[curr_p] = gdofs[d] - dofs->n_global_free_dofs();
				curr_p++;
			}
		}
	}
	for(Mesh::enrich_node_iterator it=mesh->enrich_nodes_begin(), end=mesh->enrich_nodes_end(); it!=end; ++it)
	{
		id_type l_node = mesh->global_to_local_node((*it)->get_id());
		std::vector<id_type>& gdofs = dofs->get_nodal_global_dofs_local(l_node);
		std::vector<id_type>& ldofs = dofs->get_nodal_local_dofs_local(l_node);
		for(id_type d=0; d<gdofs.size(); ++d)
		{
			if(gdofs[d] < dofs->n_global_free_dofs()) // This is a free dof
			{
				IS_tf[curr_f] = _r_start + ldofs[d];
				IS_gf[curr_f] = gdofs[d];
				curr_f++;
			}
			else // This is a constrained dof ( Should never happen as all enriched dofs are free )
			{
				IS_tp[curr_p] = _r_start + ldofs[d];
				IS_gp[curr_p] = gdofs[d] - dofs->n_global_free_dofs();
				curr_p++;
			}
		}
	}
	ierr = ISCreateGeneral(mesh->get_comm(), n_f, IS_gf, PETSC_COPY_VALUES, &_global_f);
	ierr = ISCreateGeneral(mesh->get_comm(), n_f, IS_tf, PETSC_COPY_VALUES, &_to_f);
	ierr = ISCreateGeneral(mesh->get_comm(), n_p, IS_gp, PETSC_COPY_VALUES, &_global_p);
	ierr = ISCreateGeneral(mesh->get_comm(), n_p, IS_tp, PETSC_COPY_VALUES, &_to_p);

	return ierr;
}







// This function determines the local and global sizes of all matricies used in the global solution procedure
// This function also determines the sparsity pattern of all matricies from the geometry of the mesh
// Knowledge of the sparsity pattern of the matricies allows to memory preallocation which can increase assembly performance up to 50x
PetscErrorCode Solver::preallocate_matrix(pmatrix& K)
{
	// If this solver has already been used to solve a FE system, then I need to
	// destroy the matrix that's already in the system
	if (_setup)
	{
		_K.destroy();
	}
	Mesh* mesh = _prob->get_mesh();
	DofObject* dofs = _prob->get_dofs();

	// Can do the serial case, much simpler
	if(mesh->serial())
	{
		std::vector<PetscInt> Kff_nnz( dofs->n_global_free_dofs() );
		std::vector<PetscInt> Kfp_nnz( dofs->n_global_free_dofs() );
		std::vector<PetscInt> Kpf_nnz( dofs->n_global_const_dofs() );
		std::vector<PetscInt> Kpp_nnz( dofs->n_global_const_dofs() );	

		 // Determine the dof neighbor lists for every node
		std::vector<std::set<id_type> > dof_neighbors(mesh->n_local_nodes() + mesh->n_local_enrich_nodes());
		populate_local_neighbors(dof_neighbors);
		

		// Iterate through all neighbors lists and determine how many of the neighbors belong to the diagonal and off diagonal portions of each matrix
		combine_periodic_contributions(dof_neighbors);
		for(id_type n=0; n<dof_neighbors.size(); ++n)
		{
			id_type id;
			if (n<mesh->n_local_nodes())
				id = n;
			else
				id = n - mesh->n_local_nodes() + ENRICH_START;
			if (_prob->get_mesh()->check_node_responsibility(id)) // Own this node (meaningless in serial) and a primary periodic node (if in a periodic nodeset)
			{
				// Grad a reference to this nodes global dofs and loop over them (the rows)
				std::vector<id_type> dof_vec = dofs->get_nodal_global_dofs_local(id);

				for(id_type d=0; d<dof_vec.size(); ++d)
				{
					id_type row_dof = dof_vec[d];
					// Iterate through all the dof neighbors in this nodes list
					for(auto it=dof_neighbors[n].begin(), end=dof_neighbors[n].end(); it!=end; it++)
					{
						id_type col_dof = *it;
						// If the row is a free dof then it belongs to either K_ff or K_fp
						if(row_dof < dofs->n_global_free_dofs() )
						{
							if(col_dof < dofs->n_global_free_dofs() ) // If the column is a free dof it belongs to K_ff
								Kff_nnz[row_dof]++;
							else 					// Otherwise it belongs to K_fp
								Kfp_nnz[row_dof]++;
						}
						
						// Otherwise it belongs to K_pf or K_pp
						else
						{
							if(col_dof < dofs->n_global_free_dofs() )   // If the column is a free dof it belongs to K_ff
								Kpf_nnz[row_dof - dofs->n_global_free_dofs()]++;
							else 								// Otherwise it belongs to Kpp
								Kpp_nnz[row_dof - dofs->n_global_free_dofs()]++;
						}
					}
				}
			}
		}

		// ACtually create the matrices
		PetscErrorCode ierr;
		ierr = MatCreateSeqAIJ(mesh->get_comm(), dofs->n_global_free_dofs(), dofs->n_global_free_dofs(), 0, &Kff_nnz[0], &K.ff);CHKERRQ(ierr);
		ierr = MatCreateSeqAIJ(mesh->get_comm(), dofs->n_global_free_dofs(), dofs->n_global_const_dofs(), 0, &Kfp_nnz[0], &K.fp);CHKERRQ(ierr);
		ierr = MatCreateSeqAIJ(mesh->get_comm(), dofs->n_global_const_dofs(), dofs->n_global_free_dofs(), 0, &Kpf_nnz[0], &K.pf);CHKERRQ(ierr);
		ierr = MatCreateSeqAIJ(mesh->get_comm(), dofs->n_global_const_dofs(), dofs->n_global_const_dofs(), 0, &Kpp_nnz[0], &K.pp);CHKERRQ(ierr);
		// ierr = MatSetBlockSize(Kff, _prob->nndof());CHKERRQ(ierr);
		ierr = MatSetUp(K.ff);CHKERRQ(ierr);
		return ierr;
	}


	PetscErrorCode ierr;
	PetscInt n, N, n_const, N_const;
	PetscInt Kff_rstart, Kfp_rstart, Kpf_rstart, Kpp_rstart, Kff_rend, Kfp_rend, Kpf_rend, Kpp_rend;
	PetscInt Kff_cstart, Kfp_cstart, Kpf_cstart, Kpp_cstart, Kff_cend, Kfp_cend, Kpf_cend, Kpp_cend;
	n = dofs->n_local_owned_free_dofs();
	N = dofs->n_global_free_dofs();
	n_const = dofs->n_local_owned_const_dofs();
	N_const = dofs->n_global_const_dofs();
	
	// Create our matricies, defining types and sizes. Preallocation will be done later
	ierr = MatCreate(mesh->get_comm(), &K.ff);CHKERRQ(ierr);
	ierr = MatCreate(mesh->get_comm(), &K.fp);CHKERRQ(ierr);
	ierr = MatCreate(mesh->get_comm(), &K.pf);CHKERRQ(ierr); // Could technically be done with a subcommunicator and only involve processes that actually have constrained dofs
	ierr = MatCreate(mesh->get_comm(), &K.pp);CHKERRQ(ierr); // Because really these are just arbitrary partitions right now. Distributing the solutions back to the processes that own the nodes will be annoying

	ierr = MatSetType(K.ff, MATAIJ);CHKERRQ(ierr);
	ierr = MatSetType(K.fp, MATAIJ);CHKERRQ(ierr);
	ierr = MatSetType(K.pf, MATAIJ);CHKERRQ(ierr);
	ierr = MatSetType(K.pp, MATAIJ);CHKERRQ(ierr);

	ierr = MatSetSizes(K.ff, n, n, N, N);CHKERRQ(ierr);
	ierr = MatSetSizes(K.fp, n, n_const, N, N_const);CHKERRQ(ierr);  // Not really sure if I should be using PETSC_DECIDE or  PETSC_DETERMINE here
	ierr = MatSetSizes(K.pf, n_const, n, N_const, N);CHKERRQ(ierr);  // Ditto. Also not sure about the local column count
	ierr = MatSetSizes(K.pp, n_const, n_const, N_const, N_const);CHKERRQ(ierr);

	// ierr = MatSetBlockSize(Kff, _prob->nndof());CHKERRQ(ierr);
	ierr = MatSetUp(K.fp);CHKERRQ(ierr);   // Need to do this to get column ownership???? Hoping I can change it later
	ierr = MatSetUp(K.ff);CHKERRQ(ierr); 
	ierr = MatSetUp(K.pf);CHKERRQ(ierr);
	ierr = MatSetUp(K.pp);CHKERRQ(ierr);
	
	ierr = MatGetOwnershipRange(K.ff, &Kff_rstart, &Kff_rend);CHKERRQ(ierr);
	ierr = MatGetOwnershipRange(K.fp, &Kfp_rstart, &Kfp_rend);CHKERRQ(ierr);
	ierr = MatGetOwnershipRange(K.pf, &Kpf_rstart, &Kpf_rend);CHKERRQ(ierr);
	ierr = MatGetOwnershipRange(K.pp, &Kpp_rstart, &Kpp_rend);CHKERRQ(ierr);

	ierr = MatGetOwnershipRangeColumn(K.ff, &Kff_cstart, &Kff_cend);CHKERRQ(ierr);
	ierr = MatGetOwnershipRangeColumn(K.fp, &Kfp_cstart, &Kfp_cend);CHKERRQ(ierr);
	ierr = MatGetOwnershipRangeColumn(K.pf, &Kpf_cstart, &Kpf_cend);CHKERRQ(ierr);
	ierr = MatGetOwnershipRangeColumn(K.pp, &Kpp_cstart, &Kpp_cend);CHKERRQ(ierr);

	// Because the stiffness matrix is symmetric
	// ierr = MatSetOption(Kff, MAT_SYMMETRIC, PETSC_TRUE);CHKERRQ(ierr);
	// ierr = MatSetOption(Kpp, MAT_SYMMETRIC, PETSC_TRUE);CHKERRQ(ierr);
	
	// Create the vectors that will be passed to PETSc
	// I tried doing this with normal arrays but I overflowed the stack lol
	std::vector<PetscInt> Kff_d_nnz(Kff_rend-Kff_rstart);
	std::vector<PetscInt> Kff_o_nnz(Kff_rend-Kff_rstart);
	std::vector<PetscInt> Kfp_d_nnz(Kfp_rend-Kfp_rstart);
	std::vector<PetscInt> Kfp_o_nnz(Kfp_rend-Kfp_rstart);
	std::vector<PetscInt> Kpf_d_nnz(Kpf_rend-Kpf_rstart);
	std::vector<PetscInt> Kpf_o_nnz(Kpf_rend-Kpf_rstart);
	std::vector<PetscInt> Kpp_d_nnz(Kpp_rend-Kpp_rstart);
	std::vector<PetscInt> Kpp_o_nnz(Kpp_rend-Kpp_rstart);	

	// Determine the dof neighbor lists for every node
	// Normal nodes are numbered for the first n_local_nodes entries and the next n_local_enrich_nodes are for the enrichment nodes
	std::vector<std::set<id_type> > dof_neighbors(mesh->n_local_nodes() + mesh->n_local_enrich_nodes()); // Used to store the dof neighbors for each node
	std::set<id_type>::iterator set_it;
	std::set<id_type>::iterator set_end;
	populate_local_neighbors(dof_neighbors);
	

	 // Serialize information about nodes I don't own on an interface to send to other partitions and also see how many I should be recieving 
	std::map<int, int> n_n_sets_expected;
	std::map<int, std::vector<id_type> > dof_sets_ptr_send_lists;
	std::map<int, std::vector<id_type> > dof_sets_ind_send_lists;
	std::map<int, id_type> curr_inds1;
	for(Mesh::partition_iterator it=mesh->partition_interface_begin(), end=mesh->partition_interface_end(); it!=end; ++it)
	{
		id_type g_node = it->first;
		id_type l_node = mesh->global_to_local_node(g_node);
		int owner = mesh->get_node_owner_local(l_node);
		if(owner != mesh->get_rank()) // I don't own this node. Need to serialize list
		{
			if( dof_sets_ptr_send_lists.find(owner) == dof_sets_ptr_send_lists.end() ) // Never found this partition before so I have to have empty vectors to insert into
			{
				dof_sets_ptr_send_lists.insert(std::pair<int, std::vector<id_type> >(owner, std::vector<id_type>()));
				dof_sets_ind_send_lists.insert(std::pair<int, std::vector<id_type> >(owner, std::vector<id_type>()));
			}
			dof_sets_ptr_send_lists[owner].push_back(g_node);
			dof_sets_ptr_send_lists[owner].push_back(curr_inds1[owner]);
			set_it = dof_neighbors[l_node].begin();
			set_end = dof_neighbors[l_node].end();
			// Iterate through the nodes neighbor lists and add all dofs to the send list
			for(; set_it!=set_end; ++set_it)
			{
				dof_sets_ind_send_lists[owner].push_back(*set_it);
				curr_inds1[owner]++;
			}
			dof_sets_ptr_send_lists[owner].push_back(curr_inds1[owner]);
		}
		else // Otherwise I own the node So i need to expect information for it
		{
			for(id_type i=0; i<(*it).second.size(); ++i)
				if((*it).second[i] != mesh->get_rank())
					n_n_sets_expected[(*it).second[i]]++;
		}
	}
	if(mesh->IGFEM())
	{
		// Loop over enrichment nodes on a partition interface and add them to the serialized lists to send
		// NOTE: THIS ASSUMES THAT ONLY ONE INCLUSION CUTS THE ELEMENT. IF WE LATER ALLOW MORE THAN ONE INCLUSION TO CUT AN ELEMENT THEN 
		// ENRICHMENT NODES NOT ARISING FROM THE SAME INCLUSION DON'T NEIGHTBOR EACH OTHER (I think...)
		for(Mesh::enrich_partition_iterator it=mesh->enrich_partition_interface_begin(), end=mesh->enrich_partition_interface_end(); it!=end; ++it)
		{
			id_type g_node = (*it).first;
			id_type l_node = mesh->global_to_local_node(g_node);
			id_type l_e_node = l_node - ENRICH_START + mesh->n_local_nodes();
			int owner = mesh->get_node_owner_local(l_node);
			if(owner != mesh->get_rank()) // I don't own this node. Need to serialize list
			{
				if( dof_sets_ptr_send_lists.find(owner) == dof_sets_ptr_send_lists.end() ) // Never found this partition before so I have to have empty vectors to insert into
				{
					dof_sets_ptr_send_lists.insert(std::pair<int, std::vector<id_type> >(owner, std::vector<id_type>()));
					dof_sets_ind_send_lists.insert(std::pair<int, std::vector<id_type> >(owner, std::vector<id_type>()));
				}
				dof_sets_ptr_send_lists[owner].push_back(g_node);
				dof_sets_ptr_send_lists[owner].push_back(curr_inds1[owner]);
				set_it = dof_neighbors[l_e_node].begin();
				set_end = dof_neighbors[l_e_node].end();
				// Iterate through the nodes neighbor lists and add all dofs to the send list
				for(; set_it!=set_end; ++set_it)
				{
					dof_sets_ind_send_lists[owner].push_back(*set_it);
					curr_inds1[owner]++;
				}
				dof_sets_ptr_send_lists[owner].push_back(curr_inds1[owner]);
			}
			else // Otherwise I own the node So i need to expect information for it
			{
				for(id_type i=0; i<it->second.size(); ++i)
				{
					int part = it->second[i];
					if(part != mesh->get_rank())
						n_n_sets_expected[part]++;
				}
			}
		}
	}
	
	
	// Send the information to the other processors
	int n_sends = dof_sets_ptr_send_lists.size()*2;
	MPI_Request reqs1[n_sends];
	id_type n_sent = 0;
	for(auto it=dof_sets_ptr_send_lists.begin(), end=dof_sets_ptr_send_lists.end(); it!=end; ++it)
	{
		int rank = (*it).first;
		MPI_Isend(&(*it).second[0], (*it).second.size(), MPI_ID, rank, 0, mesh->get_comm(), &reqs1[n_sent*2]);
		MPI_Isend(&dof_sets_ind_send_lists[rank][0], dof_sets_ind_send_lists[rank].size(), MPI_ID, rank, 1, mesh->get_comm(), &reqs1[n_sent*2+1]);
		n_sent++;
	}
	
	// Recieve the nodal neighbor information from other processors
	MPI_Status status[2];
	for(auto it=n_n_sets_expected.begin(), end=n_n_sets_expected.end(); it!=end; ++it)
	{
		int rank = (*it).first;
		int count0, count1;
		std::vector<id_type> neighbor_ptr_recv_vec;
		std::vector<id_type> neighbor_ind_recv_vec;
		MPI_Probe(rank, 0, MPI_COMM_WORLD, &status[0]);
		MPI_Probe(rank, 1, MPI_COMM_WORLD, &status[1]);
		MPI_Get_count(&status[0], MPI_ID, &count0);
		MPI_Get_count(&status[1], MPI_ID, &count1);
		if((*it).second*3 != count0)  // Safety check which hopefully shouldn't happen. Mostly for debugging
			err_message("During dof neighbor check, the number of dof neighbor lists recieved did not match the number of dofs neighbor lists expected");
		neighbor_ptr_recv_vec.resize(count0);
		neighbor_ind_recv_vec.resize(count1);
		MPI_Recv(&neighbor_ptr_recv_vec[0], count0, MPI_ID, rank, 0, mesh->get_comm(), &status[0]);
		MPI_Recv(&neighbor_ind_recv_vec[0], count1, MPI_ID, rank, 1, mesh->get_comm(), &status[1]);
		
		// Parse the information into the node neighbor sets
		for(int i=0; i<count0; i=i+3)
		{
			id_type g_node = neighbor_ptr_recv_vec[i];
			id_type l_node = mesh->global_to_local_node(g_node);
			if(l_node >= ENRICH_START)
				l_node = l_node - ENRICH_START + mesh->n_local_nodes();
			id_type start = neighbor_ptr_recv_vec[i+1];
			id_type end = neighbor_ptr_recv_vec[i+2];
			dof_neighbors[l_node].insert(neighbor_ind_recv_vec.begin()+start, neighbor_ind_recv_vec.begin()+end);
		}
	}
	
	// Iterate through all neighbors lists and determine how many of the neighbors belong to the diagonal and off diagonal portions of each matrix
	combine_periodic_contributions(dof_neighbors);
	for (id_type n=0; n<dof_neighbors.size(); ++n)
	{
		// get the actual local node number
		id_type l_node = n;
		if (n >= mesh->n_local_nodes())
			l_node = n - mesh->n_local_nodes() + ENRICH_START;
		id_type node_id = mesh->get_node_local(l_node)->get_id();
		
		// I only need to do this if I actually own the node
		if (mesh->check_node_responsibility(node_id))
		{
			// Grad a reference to the global dofs of the local node and loop over all of the rows
			std::vector<id_type>& dof_vec = dofs->get_nodal_global_dofs_local(l_node);;

			for (id_type d=0; d<dof_vec.size(); ++d)
			{
				id_type row_dof = dof_vec[d];
				
				// Iterate through all the dof neighbors in this nodes list
				for (auto it=dof_neighbors[n].begin(), end=dof_neighbors[n].end(); it!=end; it++)
				{
					id_type col_dof = *it;
					// If the row is a free dof then it belongs to either K_ff or K_fp
					if(row_dof < dofs->n_global_free_dofs() )
					{
						id_type local_row = row_dof - Kff_rstart; // This is row that I own in the local chunk of rows
						// If the column is a free dof it belongs to K_ff
						if(col_dof < dofs->n_global_free_dofs() )
						{
							if(col_dof>=(id_type)Kff_cstart && col_dof<(id_type)Kff_cend)
								Kff_d_nnz[local_row]++;
							else
								Kff_o_nnz[local_row]++;
						}
						
						// Otherwise it belongs to K_fp
						else
						{
							id_type col_dof_p = col_dof - dofs->n_global_free_dofs();
							if(col_dof_p>=(id_type)Kfp_cstart && col_dof_p<(id_type)Kfp_cend)
								Kfp_d_nnz[local_row]++;
							else
								Kfp_o_nnz[local_row]++;
						}
					}
					
					
					
					// Otherwise it belongs to K_pf or K_pp
					else
					{
						id_type row_dof_p = row_dof - dofs->n_global_free_dofs();
						id_type local_row = row_dof_p - Kpp_rstart; // This is row that I own in the local chunk of rows
						// If the column is a free dof it belongs to K_ff
						if(col_dof < dofs->n_global_free_dofs() )
						{
							if(col_dof>=(id_type)Kpf_cstart && col_dof<(id_type)Kpf_cend)
								Kpf_d_nnz[local_row]++;
							else
								Kpf_o_nnz[local_row]++;
						}
						
						// Otherwise it belongs to K_fp
						else
						{
							id_type col_dof_p = col_dof - dofs->n_global_free_dofs();
							if(col_dof_p>=(id_type)Kpp_cstart && col_dof_p<(id_type)Kpp_cend)
								Kpp_d_nnz[local_row]++;
							else
								Kpp_o_nnz[local_row]++;
						}
					}
				}
			}
		}
	}
	
	// Actually set the preallocations
	ierr = MatMPIAIJSetPreallocation(K.ff, 0, &Kff_d_nnz[0], 0, &Kff_o_nnz[0]);CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(K.fp, 0, &Kfp_d_nnz[0], 0, &Kfp_o_nnz[0]);CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(K.pf, 0, &Kpf_d_nnz[0], 0, &Kpf_o_nnz[0]);CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(K.pp, 0, &Kpp_d_nnz[0], 0, &Kpp_o_nnz[0]);CHKERRQ(ierr);
	
	// Wait for the non-blocking sends to complete just in case
	MPI_Waitall(n_sends, reqs1, MPI_STATUSES_IGNORE);

	return ierr;
}

// Helper function for the preallocate matrices function
void Solver::populate_local_neighbors(std::vector<std::set<id_type> >& dof_neighbors)
{
	Mesh* mesh = _prob->get_mesh();
	DofObject* dofs = _prob->get_dofs();

	// Loop over all of the nodes of the mesh to add al the dof neighbors for each node
	for(Mesh::node_iterator it=mesh->nodes_begin(), end=mesh->nodes_end(); it!=end; ++it) // Loop over all nodes in the local mesh
	{
		// Get the current local node number (this should be equal to n... but just in case I change anything)
		id_type l_node = mesh->global_to_local_node((*it)->get_id());
		// Grab a reference to the node_elem vector for this node
		std::vector<id_type>& node_elem = mesh->get_node_elem_local(l_node);

		// Loop over all elements that node is connected to on the local partition
		for(id_type e=0; e<node_elem.size(); ++e)
		{
			// Figure out which element this is and grab a reference to its local elem_node table
			id_type el = mesh->global_to_local_elem(node_elem[e]); // local element number of e'th neighbor
			std::vector<id_type>& elem_node = mesh->get_elem_node_local(el);

			// Loop over each neighboring element's nodes
			for(id_type n2=0; n2<elem_node.size(); ++n2)
			{
				id_type l_node2 = elem_node[n2];
				// Grad a reference to the global dofs for this node and insert them all into the node's set
				std::vector<id_type>& gdofs = dofs->get_nodal_global_dofs_local(l_node2);
				dof_neighbors[l_node].insert(gdofs.begin(), gdofs.end());
			}

			// If the element has any enriched dofs, I need to add those to every nodes neighbors as well
			if(mesh->IGFEM())
			{
				Elem* elem = mesh->get_elem_local(el);
				for(id_type n2=0; n2<elem->n_enrich_nodes(); ++n2)
				{
					id_type l_node2 = mesh->global_to_local_node(elem->get_enrich_node(n2)->get_id());
					std::vector<id_type>& gdofs = dofs->get_nodal_global_dofs_local(l_node2);
					dof_neighbors[l_node].insert(gdofs.begin(), gdofs.end());
				}
			}
		}
	}
	// Loop over enrichment nodes and set their neighbors
	if(mesh->IGFEM())
	{
		for(Mesh::enrich_node_iterator it=mesh->enrich_nodes_begin(), end=mesh->enrich_nodes_end(); it!=end; ++it)
		{
			id_type l_node = mesh->global_to_local_node((*it)->get_id());
			id_type l_e_node = l_node - ENRICH_START + mesh->n_local_nodes();

			// Grab a reference to the node_elem table for this node and loop over all of this enrichment node's elements
			std::vector<id_type>& node_elem = mesh->get_node_elem_local(l_node);
			for(id_type e=0; e<node_elem.size(); ++e)
			{
				// Figure out which element this is and grab a reference to its local elem_node table
				id_type el = mesh->global_to_local_elem(node_elem[e]);
				std::vector<id_type>& elem_node = mesh->get_elem_node_local(el);

				// Loop over all of the normal nodes of each element
				for(id_type n2=0; n2<elem_node.size(); ++n2)
				{
					id_type l_node2 = elem_node[n2];
					// Grad a reference to the global dofs for this node and insert them all into the node's set
					std::vector<id_type>& gdofs = dofs->get_nodal_global_dofs_local(l_node2);
					dof_neighbors[l_e_node].insert(gdofs.begin(), gdofs.end());
				}

				// Loop over the other enriched dofs on this element
				// If the element has any enriched dofs, I need to add those to every nodes neighbors as well
				Elem* elem = mesh->get_elem_local(el);
				for(id_type n2=0; n2<elem->n_enrich_nodes(); ++n2)
				{
					id_type l_node2 = mesh->global_to_local_node(elem->get_enrich_node(n2)->get_id());
					std::vector<id_type>& gdofs = dofs->get_nodal_global_dofs_local(l_node2);
					dof_neighbors[l_e_node].insert(gdofs.begin(), gdofs.end());
				}
			}
		}

		// Loop over any hanging node constraints (Assuming we are using the penalty method)
		// Any hanging node constraints that involve enrichments will cause additional allocs
		for (Mesh::hanging_node_iterator it=mesh->hanging_nodes_begin(), end=mesh->hanging_nodes_end(); it!=end; ++it)
		{
			id_type hang_node = it->first;
			std::vector<id_type> constraint_nodes(it->second.size() + 1);
			for (id_type n=0; n<constraint_nodes.size(); ++n)
				constraint_nodes[n] = it->second[n].first;
			constraint_nodes[constraint_nodes.size() - 1] = hang_node;

			// Loop over all of the nodes in the constraint and set its dof neighbors
			for (id_type n=0; n< constraint_nodes.size(); ++n)
			{
				id_type l_node = mesh->global_to_local_node(constraint_nodes[n]);
				id_type l_e_node = l_node;
				if (l_node > ENRICH_START)
					l_e_node = l_node - ENRICH_START + mesh->n_local_nodes();

				for (id_type n2=0; n2<constraint_nodes.size(); ++n2)
				{
					id_type l_node2 = mesh->global_to_local_node(constraint_nodes[n2]);
					std::vector<id_type>& gdofs = dofs->get_nodal_global_dofs_local(l_node2);
					dof_neighbors[l_e_node].insert(gdofs.begin(), gdofs.end());
				}
			}
		}
	}
}


void Solver::combine_periodic_contributions(std::vector<std::set<id_type> >& dof_neighbors)
{
	Mesh* mesh = _prob->get_mesh();
	// Loop over all periodic nodesets and add the non-primary contributions to the primary contributions
	for (id_type nset = 0; nset<mesh->n_periodic_sets(); ++nset)
	{
		std::vector<id_type> set = mesh->get_periodic_nodeset(nset);
		id_type l_node_primary = mesh->global_to_local_node(set[0]);

		for (id_type n=1; n<set.size(); ++n)
		{
			id_type l_node = mesh->global_to_local_node(set[n]);
			dof_neighbors[l_node_primary].insert(dof_neighbors[l_node].begin(), dof_neighbors[l_node].end());
		}
	}
}











// Function to store scatter the PETSc vectors and store them in my custom storage
PetscErrorCode Solver::store_solution(pvector& U, NodalData* data)
{
	return store_solution(U.f, U.p, data);
}








PetscErrorCode Solver::store_solution(Vec& Uf, Vec& Up, NodalData* data)
{
	PetscErrorCode ierr;
	Mesh* mesh = _prob->get_mesh();
	DofObject* dofs = _prob->get_dofs();

	// Simpler serial case
	if(mesh->serial())
	{
		// Get read only pointers to the underlying vectors themselves
		const PetscScalar *Ufree, *Uconst;
		ierr = VecGetArrayRead(Uf, &Ufree);CHKERRQ(ierr);
		ierr = VecGetArrayRead(Up, &Uconst);CHKERRQ(ierr);

		// Store them in my structures
		for(Mesh::node_iterator it=mesh->nodes_begin(), end=mesh->nodes_end(); it!=end; ++it)
		{
			id_type l_node = mesh->global_to_local_node((*it)->get_id());
			std::vector<id_type>& gdofs = dofs->get_nodal_global_dofs_local(l_node);
			for(id_type d=0; d<gdofs.size(); ++d)
			{
				id_type dof = gdofs[d];
				if(dof < dofs->n_global_free_dofs()) // This is a free dof
					data->get_value_local(l_node, d) = Ufree[dof];
				else // This is a constrained dof
					data->get_value_local(l_node, d) = Uconst[dof-dofs->n_global_free_dofs()];
			}
		}
		if(mesh->IGFEM())
		{
			for(Mesh::enrich_node_iterator it=mesh->enrich_nodes_begin(), end=mesh->enrich_nodes_end(); it!=end; ++it)
			{
				id_type l_node = mesh->global_to_local_node((*it)->get_id());
				std::vector<id_type>& gdofs = dofs->get_nodal_global_dofs_local(l_node);
				for(id_type d=0; d<gdofs.size(); ++d)
				{
					id_type dof = gdofs[d];
					if(dof < dofs->n_global_free_dofs()) // This is a free dof
						data->get_value_local(l_node, d) = Ufree[dof];
					else // This is a constrained dof (NOTE: This should actually never happen since all enrich nodes are supposed to be free but just in case)
						data->get_value_local(l_node, d) = Uconst[dof-dofs->n_global_free_dofs()];
				}
			}
		}

		ierr = VecRestoreArrayRead(Uf, &Ufree);CHKERRQ(ierr);
		ierr = VecRestoreArrayRead(Up, &Uconst);CHKERRQ(ierr);
	}



	// Now I can handle the parallel case
	// -------------------------------------------------------------------------------------------------------------
	else
	{
		VecScatter scat_Uf, scat_Up;
		Vec U_sol;

		// Define recieve vectors
		ierr = VecCreateMPI(mesh->get_comm(), dofs->n_local_dofs(), PETSC_DETERMINE, &U_sol);CHKERRQ(ierr);

		// Determine the first row f the recieve vector that I own

		// Create the Scatter Contexts (Using previously allocated index sets (preallocate_index_sets)) and perform the scatter
		ierr = VecScatterCreate(Uf, _global_f, U_sol, _to_f, &scat_Uf);CHKERRQ(ierr);
		ierr = VecScatterCreate(Up, _global_p, U_sol, _to_p, &scat_Up);CHKERRQ(ierr);
		ierr = VecScatterBegin(scat_Uf, Uf, U_sol, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
		ierr = VecScatterBegin(scat_Up, Up, U_sol, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
		ierr = VecScatterEnd(scat_Uf, Uf, U_sol, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
		ierr = VecScatterEnd(scat_Up, Up, U_sol, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);

		// Retrieve the solution values from the scattered vectors
		PetscInt ix[dofs->n_local_dofs()];
		PetscScalar y_U[dofs->n_local_dofs()];
		for(id_type i=0; i<dofs->n_local_dofs(); ++i)
			ix[i] = i + _r_start;
		ierr = VecGetValues(U_sol, dofs->n_local_dofs(), ix, y_U);CHKERRQ(ierr);

		// Store the solution values in the storage structures passed in
		for(Mesh::node_iterator it=mesh->nodes_begin(), end=mesh->nodes_end(); it!=end; ++it)
		{
			id_type l_node = mesh->global_to_local_node((*it)->get_id());
			std::vector<id_type>& ldofs = dofs->get_nodal_local_dofs_local(l_node);
			for(id_type d=0; d<ldofs.size(); ++d)
				data->get_value_local(l_node, d) = y_U[ldofs[d]];
		}
		for(Mesh::enrich_node_iterator it=mesh->enrich_nodes_begin(), end=mesh->enrich_nodes_end(); it!=end; ++it)
		{
			id_type l_node = mesh->global_to_local_node((*it)->get_id());
			std::vector<id_type>& ldofs = dofs->get_nodal_local_dofs_local(l_node);
			for(id_type d=0; d<ldofs.size(); ++d)
				data->get_value_local(l_node, d) = y_U[ldofs[d]];
		}

		// Clean Up memory
		ierr = VecScatterDestroy(&scat_Uf);CHKERRQ(ierr);
		ierr = VecScatterDestroy(&scat_Up);CHKERRQ(ierr);
		ierr = VecDestroy(&U_sol);CHKERRQ(ierr);
	}

	return ierr;
}




// FIXME: If this is ever used for anything other than the residual strain problem
// then something needs to be done about fixing the ISV handing
PetscErrorCode Solver::initialSolution(NodalData& starting_solution)
{
	PetscErrorCode ierr(0);
	DofObject* dofs = _prob->get_dofs();
	Mesh* mesh = _prob->get_mesh();
	ierr = VecDuplicate(_U.p, &_Up_initial); // Create the right structure
	id_type ngfd = _prob->get_dofs()->n_global_free_dofs();

	for(Mesh::node_iterator it=mesh->nodes_begin(), end=mesh->nodes_end(); it!=end; ++it)
	{
		id_type id = (*it)->get_id();
		if (mesh->check_node_responsibility(id))
		{
			std::vector<id_type>& gdofs = dofs->get_nodal_global_dofs_global( id );
			for (id_type d=0; d<gdofs.size(); ++d)
			{
				double val = starting_solution.get_value_global(id, d);
				if (gdofs[d] < ngfd) // This is a free dof
					{ierr = VecSetValue(_U.f, gdofs[d], val, INSERT_VALUES);CHKERRQ(ierr);}
				else
					{ierr = VecSetValue(_Up_initial, gdofs[d]-ngfd, val, INSERT_VALUES);CHKERRQ(ierr);}
			}
		}
	}
	for (Mesh::enrich_node_iterator it=mesh->enrich_nodes_begin(), end=mesh->enrich_nodes_end(); it!=end; ++it)
	{
		id_type id = (*it)->get_id();
		if (mesh->check_node_responsibility(id))
		{
			std::vector<id_type>& gdofs = dofs->get_nodal_global_dofs_global( id );
			for (id_type d=0; d<gdofs.size(); ++d)
			{
				double val = starting_solution.get_value_global(id, d);
				if (gdofs[d] < ngfd) // This is a free dof
					{ierr = VecSetValue(_U.f, gdofs[d], val, INSERT_VALUES);CHKERRQ(ierr);}
				else
					{ierr = VecSetValue(_Up_initial, gdofs[d]-ngfd, val, INSERT_VALUES);CHKERRQ(ierr);}
			}
		}
	}

	// Assemble the solution
	ierr = VecAssemblyBegin(_U.f);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(_Up_initial);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(_U.f);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(_Up_initial);CHKERRQ(ierr);

	// Store the solution in the current solution object
	ierr = VecCopy(_Up_initial, _U.p);CHKERRQ(ierr);
	store_solution(_U.f, _Up_initial, _prob->get_solution()); // Note thatvthis line actually just copies the data into itself for the residual strain problem but I mught as well keep it general

	_nonzero_start = true;

	return ierr;
}
