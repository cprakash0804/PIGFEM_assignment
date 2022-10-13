/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated July 2017

##################################################################################
*/
#include "SolverNonlinear.h"
#include "Problem.h"
#include "Assembler.h"
#include "Mesh.h"
#include "DofObject.h"
#include "NodalData.h"
#include "InternalVars.h"
#include "Writer.h"
#include "SubscaleModel.h"
#include "SensitivitySolver.h"
#include "elem.h"
#include "Options.h"
#include "Utilities.h"
#include <sstream>

// PETScNonlinearKSPSolver
SolverNonlinear::SolverNonlinear()
{
}

SolverNonlinear::~SolverNonlinear()
{
}


// This function determines the local and global sizes of all vectors used in the global solution procedure
PetscErrorCode SolverNonlinear::preallocate_vectors()
{
	PetscErrorCode ierr;

	// If this FE system has already been set up before, I need to clear ou the old vectors
	if (_setup)
	{
		_U.destroy();
		_F_ext.destroy();
		_P_int.destroy();
	}

	PetscInt n, N, n_const, N_const;
	DofObject* dofs = _prob->get_dofs();
	Mesh* mesh = _prob->get_mesh();
	n = dofs->n_local_owned_free_dofs();
	N = dofs->n_global_free_dofs();
	n_const = dofs->n_local_owned_const_dofs();
	N_const = dofs->n_global_const_dofs();

	if(mesh->serial())
	{
		ierr = VecCreateSeq(mesh->get_comm(), N, &_U.f);CHKERRQ(ierr);
		ierr = VecCreateSeq(mesh->get_comm(), N, &_P_int.f);CHKERRQ(ierr);
		ierr = VecCreateSeq(mesh->get_comm(), N, &_F_ext.f);CHKERRQ(ierr);
		ierr = VecCreateSeq(mesh->get_comm(), N_const, &_U.p);CHKERRQ(ierr);
		ierr = VecCreateSeq(mesh->get_comm(), N_const, &_P_int.p);CHKERRQ(ierr);
		ierr = VecCreateSeq(mesh->get_comm(), N_const, &_F_ext.p);CHKERRQ(ierr);
	}
	else
	{
		ierr = VecCreateMPI(mesh->get_comm(), n, N, &_U.f);CHKERRQ(ierr);
		ierr = VecCreateMPI(mesh->get_comm(), n, N, &_P_int.f);CHKERRQ(ierr);
		ierr = VecCreateMPI(mesh->get_comm(), n, N, &_F_ext.f);CHKERRQ(ierr);
		ierr = VecCreateMPI(mesh->get_comm(), n_const, N_const, &_U.p);CHKERRQ(ierr);
		ierr = VecCreateMPI(mesh->get_comm(), n_const, N_const, &_P_int.p);CHKERRQ(ierr);
		ierr = VecCreateMPI(mesh->get_comm(), n_const, N_const, &_F_ext.p);CHKERRQ(ierr);
	}

	return ierr;
}



PetscErrorCode SolverNonlinear::setup()
{
	PetscErrorCode ierr;

	if (_setup)
	{
		ierr = VecDestroy(&_Res_f);CHKERRQ(ierr);
	}

	ierr = Solver::setup();CHKERRQ(ierr);

	ierr = VecDuplicate(_U.f, &_Res_f);CHKERRQ(ierr);

	return ierr;
}














// Solves the global system using PETSc.
// Adaptively increments the applied loads
PetscErrorCode SolverNonlinear::solve()
{
	// Some initializations that shouldn't be hardcoded but I can pass in later I guess
	_delta_t = _max_time_step;
	id_type max_time_step_repeat = 20;

	// Get pointers to the main objects
	Mesh* mesh = _prob->get_mesh();
	Assembler* assem = _prob->get_assembler();
	InternalVars* int_vars = _prob->get_internal_vars();

	// Error code initialization
	PetscErrorCode ierr;

	// Local vector initializations
	// ierr = VecZeroEntries(*_U.f);CHKERRQ(ierr);
	// ierr = VecZeroEntries(*_U.p);CHKERRQ(ierr);
	// ierr = VecZeroEntries(*_F_ext.f);CHKERRQ(ierr);
	ierr = VecZeroEntries(_U.f);CHKERRQ(ierr);
	ierr = VecZeroEntries(_U.p);CHKERRQ(ierr);
	ierr = VecZeroEntries(_F_ext.f);CHKERRQ(ierr);
	
	// Load step initializations
	double Res_norm;
	unsigned int current_time_step_repeat = 0;
	solver_state previous;
	std::vector<std::vector<std::vector<double> > > update_ISVs;
	int_vars->preallocate_internal_vars_object(update_ISVs);
	ierr = VecDuplicate(_U.f, &previous.Uf);CHKERRQ(ierr);
	ierr = VecCopy(_U.f, previous.Uf);CHKERRQ(ierr);
	unsigned int step_iter = 0;
	double current_t = 0.0;
	bool sens_final_only = _prob->getOptions()->getBoolOption("-SENSITIVITY_FINAL_ONLY");
	bool out_final_only = _prob->getOptions()->getBoolOption("-OUTPUT_FINAL_ONLY");
	Utilities::timer timer;
	int convergence = 0;   // variable used for output just to show what kinda of convergence we had
	bool problem_controlled_dt = false;
	double prev_dt = _delta_t;


	// Write some things to screen and to file
	PIGFEMPrint("\n\nStarting time step algorithm.\n");
	PIGFEMPrint("Cnv cheatsheet:\n  1: Converged\n -1: NR-diverged\n -2: Large damage rate\n");
	PIGFEMPrint("LStp      T(s)      TStep(s)    Rpt   Cnv    Res\n");
	PIGFEMPrint("----   ----------   --------    ---   ---  --------\n");
	if (!out_final_only)
		Output( step_iter, current_t );


	// LOAD STEP LOOP
	//---------------------------
	while(current_t < _T_final)
	{
		// LOAD STEP REPITITION LOOP
		//---------------------------
		current_time_step_repeat = 0;
		//can_increase_time_step = true;
		bool converged = false;
		bool exceeds_rates = false;
		bool problem_repeat = false;
		bool last_implicit_step_failed = false;

//----------------------------------------------------------------------------------------------------------------
//
//if (current_t<2*_delta_t)
//{
		//InternalVars* a =_prob->get_internal_vars();
		Utilities::penetration_store(1);
		ierr = VecZeroEntries(_U.p);CHKERRQ(ierr);
		ierr = VecZeroEntries(_F_ext.f);CHKERRQ(ierr);
		timer.start();
		ierr = assem->assemble_new_load_step(_U.p, _F_ext, (current_t+_delta_t)/_T_final);CHKERRQ(ierr);

		if (getNonzeroStart())
			{ierr = VecAXPY(_U.p, 1.0, _Up_initial);CHKERRQ(ierr);}
		_prob->_assemble_time += timer.time_elapsed();
		store_solution(_U, _prob->get_solution());
		Assembler* assem = _prob->get_assembler();
		assem->assemble_nonlinear(_K, _P_int, _delta_t,
								  update_ISVs, _prob->get_solution(), false, true);
		int_vars->copy_data_in(update_ISVs); 
		Utilities::penetration_store(1); 
//}
//------------------------------------------------------------------------------------------------------------------





		while(!converged || exceeds_rates || problem_repeat)
		{
			// Check to see if the number of times through this load step is sufficient for divergence
			if (!problem_controlled_dt) // I want to be able to do special things in a repetition cause by the problem (like a 0 dt)
			{
				if(current_time_step_repeat==max_time_step_repeat ||
				   _delta_t < _min_time_step)
				{
					// The last time I was here was the very last implicit step so I guess I'll just give up here
					if (last_implicit_step_failed)
					{
						// Failed to converge so I'll output my statistics file and exit here
						_prob->writeStatistics();
						if (mesh->get_rank() == 0)
							err_message("Sorry, the Newton-Raphson solver did not converge!");
					}

					// I'll try the explicit solver
					PIGFEMPrint("ENTERING THE EXPLICIT SOLVER. I HOPE THIS WORKS.\n");
					last_implicit_step_failed = true;
					bool cont = true;
					KSP explicit_ksp;
					ierr = initKSP(explicit_ksp, true);CHKERRQ(ierr);
					ierr = explicit_solve(100, current_t, update_ISVs, explicit_ksp, cont);CHKERRQ(ierr);
					ierr = KSPDestroy(&explicit_ksp);CHKERRQ(ierr);

					// Write Output
					ierr = VecCopy(_P_int.p, _F_ext.p);
					store_solution(_F_ext, _prob->get_external_load());
					step_iter++;
					if (out_final_only && current_t >= _T_final)
						Output( 0, current_t );
					else
						Output( step_iter, current_t );

					// If the explicit solver didn't get all the way to the end of execution
					if (cont)
					{
						_delta_t *= 4; // Update the time step back to the time step just before failure
						double time_remaining = _T_final - current_t;
						if (time_remaining < _delta_t)
							_delta_t = time_remaining;
						converged = false;
					}
					else
						continue;
				}
			}

			// Increment the number of times I've gone through this time step
			current_time_step_repeat++;


			// Compute the new free external force and prescribed solution vectors
			ierr = VecZeroEntries(_U.p);CHKERRQ(ierr);
			ierr = VecZeroEntries(_F_ext.f);CHKERRQ(ierr);
			timer.start();
			ierr = assem->assemble_new_load_step(_U.p, _F_ext, (current_t+_delta_t)/_T_final);CHKERRQ(ierr);
			if (getNonzeroStart())
				{ierr = VecAXPY(_U.p, 1.0, _Up_initial);CHKERRQ(ierr);}
			_prob->_assemble_time += timer.time_elapsed();
			store_solution(_U, _prob->get_solution());




			// Do the actual nonlinear solve step
			ierr = nonlinearStep(update_ISVs, converged, Res_norm);CHKERRQ(ierr);



			// If NR converged and it didn't exceed the rate limits, then update the current time
			if (converged && !exceeds_rates)
			{
				// Check if the problem says we can continue with the next load step
				problem_repeat = _prob->repeat_time_step();

				// If I have to repeat it mark the convergence appropriately and restart
				if (problem_repeat)
				{
					convergence = -3;
					VecCopy(previous.Uf, _U.f);
				}

				// Otherwise, I can actually continue with the next load step
				else
				{
					convergence = 1;
					current_t += _delta_t;
				}
			}
			// Otherwise, restore variables to their previous state and restart the load step
			else
			{
				PIGFEMPrint("REPEATING LOADSTEP" << std::endl);

				// If things didn't converge and the problem is controlling the load steps then currently we're out of luck
				if (problem_controlled_dt)
					err_message("Sorry, the Newton-Raphson solver did not converge with the problem-controlled load step!");

				if(exceeds_rates)
					convergence = -2;
				if(!converged)	// Obviously not converging takes precendence over exceeding the rates
					convergence = -1;

				// Rest the solution field (Don't have to do internal variables because they were reset in the residual norm check)
				VecCopy(previous.Uf, _U.f);
			}

			// Print out some stuff to the screen or a file
			char buf[200];
			sprintf(buf, " %d     %1.5e  %1.2e     %d     %d   %1.2e",
            		step_iter+1, current_t, _delta_t,
            		current_time_step_repeat, convergence,
            		Res_norm);
			PIGFEMPrint(buf << std::endl);	

			// Update the time step. If NR did not converge or if I exceed the internal variable rates, decrease the dt,
			// otherwise I might increase it depending on certain things
			prev_dt = _delta_t;
			update_time_step(_delta_t, current_t, _T_final,
							 _min_time_step, _max_time_step, converged,
							 exceeds_rates, problem_repeat,
							 problem_controlled_dt, _rel_tol, _abs_tol,
							 current_time_step_repeat);

		}	// End load step repeat loop --------------------------











		// Finished a load step, store/write some stuff
		// Compute the reaction forces (Pp_ext) at this load step
		ierr = VecCopy(_P_int.p, _F_ext.p);
		store_solution(_F_ext, _prob->get_external_load());

		// Compute sensitivities
		if (_prob->sensitivity())
		{
			if (!out_final_only || int_vars->haveISVs())
			{
				// briefly change the delta_dt back to what it was
				std::swap(_delta_t, prev_dt);
				timer.start();
				_prob->get_sensitivity_solver()->solveProblem();
				_prob->_sensitivity_time += timer.time_elapsed();
				std::swap(_delta_t, prev_dt);
			}
			else
			{
				if (!sens_final_only || (sens_final_only && current_t>=_T_final))
				{
					// briefly change the delta_dt back to what it was
					std::swap(_delta_t, prev_dt);
					timer.start();
					_prob->get_sensitivity_solver()->solveProblem();
					_prob->_sensitivity_time += timer.time_elapsed();
					std::swap(_delta_t, prev_dt);
				}
			}
		}

		// Update internal variables
		int_vars->copy_data_in(update_ISVs);

		// Store the converged solution vector as the new previous
		VecCopy(_U.f, previous.Uf);

		// Increment the step iteration counter
		step_iter++;
		
		// Write output
		if (out_final_only)
		{
			if (current_t >= _T_final)
				Output( 0, current_t );
		}
		else
			Output( step_iter, current_t ); // step_iter is 0 to initialize the files that get written to every step


	}	// End Load step loop --------------------------

	// Clean up memory
	ierr = VecDestroy(&previous.Uf);CHKERRQ(ierr);

	// Set Boolean
	_prob->get_solved() = true;

	return ierr;
}























































// Function to modify the time step
// This function can be changed to more smartly adapt the timestep
void SolverNonlinear::update_time_step(double& delta_t, double& curr_t, double final_t,
											   double min_dt, double max_dt, bool converged,
											   bool exceeds_rates, bool problem_repeat,
											   bool& problem_controlled_dt, double& rel_tol, double& abs_tol,
											   unsigned int curr_time_step_repeat)
{
	// Constants that control how much I'll increase or decease the time step
	double decrease = 0.25;
	double increase = 1.5;

	// If I didn't converge, or if my convergence exceed the internal variable rate limits, decrease the time step
	if (!converged || exceeds_rates)
		delta_t = decrease*delta_t;

	// Otherwise, If I converged without decreasing the time step then I'll try increasing it a bit as long as I'm not supposed to repeat the same load step
	else
	{
		if (!problem_controlled_dt)
			if(curr_time_step_repeat==1)
				delta_t = increase*delta_t;

		// Now, after everthing else, I'll see if the problem wants to control the time this next load step
		problem_controlled_dt = _prob->update_time_step(delta_t, rel_tol, abs_tol);
	}

	// Limit the time step just in case
	double time_remaining = final_t - curr_t;
	if(delta_t > max_dt)
		delta_t = max_dt;
	if(delta_t > time_remaining)		// If the time remaining in the simulation is small enough that the chosen dt would cause the simulation to go over time, we'll limit it
	{
		if(time_remaining >= min_dt)	// If the time remaining is larger than the minimum dt, then set the dt to the time remaining		
			delta_t = time_remaining;
		else							// Otherwise the time remaining is so small that even the minimum time step would cause us to go overtime, so basically we're done now
			curr_t = final_t;
	}
}






/*
 * Perform a series of explicit solves as a last resort to maybe get past a difficult part in the implicit solver
 */
PetscErrorCode SolverNonlinear::explicit_solve(id_type n_iterations, double& current_time, std::vector<std::vector<std::vector<double> > >& update_ISVs, KSP& ksp, bool& cont)
{
	// Objects
	PetscErrorCode ierr;
	Assembler* assem = _prob->get_assembler();
	InternalVars* int_vars = _prob->get_internal_vars();
	Utilities::timer timer;

	// Declare delta_U
	Vec delta_Uf;
	ierr = VecDuplicate(_U.f, &delta_Uf);CHKERRQ(ierr);
	double delta_t = _min_time_step;	// Set our explicit time step as the absolute minimum
	double Res_norm;

	// id_type n_iter_per_out = _max_time_step/_min_time_step;
	id_type n_iter_till_end = static_cast<id_type>((_T_final - current_time) / _min_time_step + 0.5); // Round up the number of time steps left
	if (n_iter_till_end < n_iterations)
	{
		n_iterations = n_iter_till_end;
		cont = false;
	}
	else
		cont = true;

	for (id_type step=0; step<n_iterations; ++step)
	{
		// Output for purely explicit solver
		// if (step%n_iter_per_out == 0)
		// 	Output(step/n_iter_per_out, current_time);

		// Assemble the new time step
		timer.start();
		ierr = assem->assemble_new_load_step(_U.p, _F_ext, (current_time + delta_t)/_T_final);CHKERRQ(ierr);
		if (getNonzeroStart())
			{ierr = VecAXPY(_U.p, 1.0, _Up_initial);CHKERRQ(ierr);}
		_prob->_assemble_time += timer.time_elapsed();

		// Assemble the new stiffness matrix and inernal force vectors
		timer.start();
		ierr = assem->assemble_nonlinear(_K, _P_int, delta_t, update_ISVs, _prob->get_solution(), true, true);CHKERRQ(ierr);
		_prob->_assemble_time += timer.time_elapsed();

		// Compute the new residual vector
		ierr = VecCopy(_F_ext.f, _Res_f);CHKERRQ(ierr);
		ierr = VecAXPY(_Res_f, -1.0, _P_int.f);CHKERRQ(ierr);
		ierr = VecNorm(_Res_f, NORM_2, &Res_norm);CHKERRQ(ierr);

		// Some output
		PIGFEMPrint("\tEXPLICIT SOLVE STEP " << step << " Current time: " << current_time << " Initial Residual: " << Res_norm << std::endl);

		if (!std::isfinite(Res_norm))
			err_message("Invalid residual norm found during an explicit solve step!");

		// Solve the system once
		timer.start();
		ierr = KSPSolve(ksp, _Res_f, delta_Uf);CHKERRQ(ierr);	// Solve the system
		_prob->_solve_time += timer.time_elapsed();

		// Update the solution values
		ierr = VecAXPY(_U.f, 1.0, delta_Uf);CHKERRQ(ierr);
		store_solution(_U, _prob->get_solution());;

		// Copy over stuff from the finished time step
		int_vars->copy_data_in(update_ISVs);

		// Update the current time
		current_time += delta_t;
	}

	return ierr;
}

















































// Function that calculates the average rate of change of all of the internal variables
void SolverNonlinear::compute_int_var_avgs(InternalVars* vars,
												   std::map<std::string, std::vector<double> >& avg,
												   std::map<std::string, std::vector<double> >& maxes)
{
	// Declare the local variables that will total the local contributions
	std::map<std::string, std::vector<double> > avg_local, maxes_local;

	// Clear the current avg and maxes
	for(InternalVars::internal_var_iterator it = vars->internal_var_begin(), end = vars->internal_var_end(); it!=end; it++)
		avg_local.insert(std::pair<std::string, std::vector<double> >(it->first, std::vector<double>(it->second, 0.0)));
	avg = avg_local;
	maxes = avg_local;
	maxes_local = avg_local;

	// Loop over all of the elements
	Mesh* mesh = _prob->get_mesh();
	std::map<std::string, id_type> nqp_local; // Count of the local number of qudrature points for each material
	for(Mesh::element_iterator it=mesh->active_elements_begin(), end=mesh->active_elements_end(); it!=end; ++it)
	{
		// Get local element number
		id_type local_e = mesh->global_to_local_elem((*it)->get_id());

		Material* mat;
		if (!(*it)->is_intersected()) // non-intersected
		{
			// Get the material associated witht this element
			mat = mesh->get_element_material_global((*it)->get_id());
			std::string name = mat->get_name();

			// Loop over the quadrature points and pull the values
			id_type nqp = (*it)->n_q_points();
			for(id_type qp=0; qp<nqp; ++qp)
			{
				for(id_type iv=0; iv<mat->n_internal_vars(); ++iv)
				{
					avg_local[name][iv] += vars->get_internal_var_local(local_e, qp, iv);
					maxes_local[name][iv] = std::max(maxes_local[name][iv], vars->get_internal_var_local(local_e, qp, iv));
				}
			}

			// increment the number of local quadrature points
			nqp_local[name] += nqp;
		}



		else // intersected
		{
			id_type elem_qp = 0;
			for(id_type ie=0; ie<(*it)->n_integration_elem(); ++ie)
			{
				Elem* int_el = (*it)->get_integration_elem(ie);
				mat = mesh->get_element_material_global((*it)->get_id(), ie);
				std::string name = mat->get_name();

				// Loop over the quadrature points of
				id_type nqp = int_el->n_q_points();
				for(id_type qp=0; qp<nqp; ++qp)
				{
					for(id_type iv=0; iv<mat->n_internal_vars(); ++iv)
					{
						avg_local[name][iv] += vars->get_internal_var_local(local_e, elem_qp+qp, iv);
						maxes_local[name][iv] = std::max(maxes_local[name][iv], vars->get_internal_var_local(local_e, elem_qp+qp, iv));
					}
				}

				// FIXME: Need to do cohesive elements here as well

				// increment the number of local quadrature points
				elem_qp += nqp;
				nqp_local[name] += nqp;
			}
		}
	}



	// Perform global reductions to determine the global rates and maxes
	for(auto it=nqp_local.begin(), end=nqp_local.end(); it!=end; ++it)
	{
		std::string name = it->first;

		// If there are actually any internal variables for this material, then do the reductions
		std::map<std::string, std::vector<double> >::iterator it2 = avg_local.find(name);
		if(it2->second.size() != 0)
		{
			// Only sum the nuimber of quadrature points on the first time calling this function. Saves one global reduction every call except the first
			if(_nqp_global.size()==0)
			{
				_nqp_global.insert(std::pair<std::string, id_type>(name, 0));
				MPI_Allreduce(&it->second, &_nqp_global[name], 1, MPI_ID, MPI_SUM, mesh->get_comm()); // Global reduction on the number of quadrature points
			}

			MPI_Allreduce(&it2->second[0], &avg[name][0], it2->second.size(), MPI_DOUBLE, MPI_SUM, mesh->get_comm()); // Global reduction of average value of internal variable
			MPI_Allreduce(&maxes_local[name][0], &maxes[name][0], maxes[name].size(), MPI_DOUBLE, MPI_MAX, mesh->get_comm()); // Global reductino on the maximum value of each internal variable

			// Actually do the global averaging
			for(id_type iv=0; iv<avg[name].size(); ++iv)
				avg[name][iv] = avg[name][iv]/_nqp_global[name];
		}
	}
}





bool SolverNonlinear::check_rates(std::map<std::string, std::vector<double> >& avgs,
								  std::map<std::string, std::vector<double> >& avg_prev,
								  std::map<std::string, std::vector<double> >& maxes,
								  std::map<std::string, std::vector<double> >& rate_limits,
								  std::map<std::string, std::vector<double> >& check_limit)
{
	// Loop over all of the materials
	for(auto it=avgs.begin(), end=avgs.end(); it!=end; ++it)
	{
		std::string name = it->first;

		// Loop over all of the internal variables for that materials
		for(id_type iv=0; iv<it->second.size(); ++iv)
		{

			// If i'm supposed to check it
			if(maxes[name][iv] > check_limit[name][iv])
			{
				double rate = avgs[name][iv] - avg_prev[name][iv];

				// If the rate exceeds the prescibed rate, return true
				if(rate > rate_limits[name][iv])
					return true;
			}
		}
	}

	// If I haven't returned true yet, hen return false
	return false;
}

