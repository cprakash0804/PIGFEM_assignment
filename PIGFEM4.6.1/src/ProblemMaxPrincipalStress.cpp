/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "ProblemMaxPrincipalStress.h"
#include "AssemblerSmallStrainNonlinear.h"
#include "Mesh.h"
#include "InternalVars.h"
#include "elem.h"
#include "Utilities.h"
#include "mpi.h"
#include <cmath>
#include <iostream>


ProblemMaxPrincipalStress::ProblemMaxPrincipalStress()
	: _prev_step_had_damage(false), _plane_strain(true), _prob_rel_tol(1e-6), _prob_abs_tol(1e-6)
{}

// This function will be overridden in derived classes to set the actual assembler and solver that will be used
void ProblemMaxPrincipalStress::generate_assembler()
{
	_assembler = new AssemblerSmallStrainNonlinearStructural;
}

// This function return the number of degrees of freedom present on each node
id_type ProblemMaxPrincipalStress::nndof()
{
	return _mesh->dim();
}

// Returns the problem type enum for this problem
problem_type ProblemMaxPrincipalStress::get_type()
{
	return PROBLEM_MAX_PRINCIPAL_STRESS;
}

// Set any possible problem specific parameters
void ProblemMaxPrincipalStress::set_parameter(std::string name, double val)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if (name=="PLANE_STRAIN" || name=="PLANE STRAIN")
	{
		if (val==0.0) // stupid hack because virtual functions can't be templated
			_plane_strain = false;
		else
			_plane_strain = true;
	}
	else if (name=="SMAX" || name=="SIGMAMAX" || name=="S_MAX" || name=="SIGMA_MAX")
		_sigma_max = val;
	else if (name=="MAGNITUDE" || name=="MAX_DAMAGE_ORDER")
		_max_order_decrease = static_cast<unsigned int>(val + 0.5);
	else if (name=="DAMAGE_STEP")
		_order_decrease_step = static_cast<unsigned int>(val + 0.5);
	else if (name=="REL_TOL")
		_prob_rel_tol = val;
	else if (name=="ABS_TOL")
		_prob_abs_tol = val;
	else
		err_message("Please input a valid parameter name.");
}

// Get whatever possible problem specific parameter
double ProblemMaxPrincipalStress::get_parameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if (name=="PLANE_STRAIN" || name=="PLANE STRAIN")
	{
		if (_plane_strain)
			return 1.0;
		else
			return 0.0;
	}
	else if (name=="SMAX" || name=="SIGMAMAX" || name=="S_MAX" || name=="SIGMA_MAX")
		return _sigma_max;
	else if (name=="MAGNITUDE" || name=="MAX_DAMAGE_ORDER")
		return (double)_max_order_decrease;
	else if (name=="DAMAGE_STEP")
		return (double)_order_decrease_step;
	else if (name=="REL_TOL")
		return _prob_rel_tol;
	else if (name=="ABS_TOL")
		return _prob_abs_tol;
	else
		err_message("Please input a valid parameter name.");
}



void ProblemMaxPrincipalStress::init()
{
	Problem::init();

	_curr_order_decrease.resize(_mesh->n_local_elem());
	for (Mesh::element_iterator it=_mesh->active_elements_begin(), end=_mesh->active_elements_end(); it!=end; ++it)
	{
		Elem* el = (*it);
		id_type l_elem = _mesh->global_to_local_elem(el->get_id());
		if (!el->is_intersected())
			_curr_order_decrease[l_elem].resize(1); // Only one element
		else
			_curr_order_decrease[l_elem].resize(el->n_integration_elem()); // One for each integration element
	}
}


void ProblemMaxPrincipalStress::fillBmat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& grad_x)
{
	ProblemUtilities::Assemble_Small_Strain_B_Mat(B, grad_x);
}




bool ProblemMaxPrincipalStress::repeat_time_step()
{
	// Loop over all active elem, checking the maximum principal stresses in each element
	// If an element is damaged then I begin progressively adding damage to it
	bool local_damage_this_step = false;
	for (Mesh::element_iterator it=_mesh->active_elements_begin(), end=_mesh->active_elements_end(); it!=end; ++it)
	{
		Elem* el = (*it);
		id_type l_elem = _mesh->global_to_local_elem(el->get_id());

		
		if (!el->is_intersected())
		{
			// If this element has already been damaged, either update the damage or move on
			if (_curr_order_decrease[l_elem][0] > 0)
			{
				if (_curr_order_decrease[l_elem][0] >= _max_order_decrease) // As damaged as we're gonna allow
					continue;
				else
				{
					local_damage_this_step = true;
					id_type orders_remaining = _max_order_decrease - _curr_order_decrease[l_elem][0];
					if (orders_remaining > _order_decrease_step) // I'll decrease the order by the step amount
						_curr_order_decrease[l_elem][0] += _order_decrease_step;
					else 	// Otherwise, go to the maximum amount
						_curr_order_decrease[l_elem][0] += orders_remaining;

					// Loop over all of the quadrature points and set the new damage state
					id_type nqp = el->n_q_points();
					for (id_type qp=0; qp<nqp; ++qp)
						get_internal_vars()->get_internal_var_local(l_elem, qp, 0) = 1.0 - 1.0/pow(10, _curr_order_decrease[l_elem][0]);
				}
			}
			else
			{
				std::vector<std::vector<double> > strains, stresses;
				compute_stress_strain(strains, stresses, el);
				bool elem_damaged = false;
				if (_mesh->get_element_material_local(l_elem)->get_type() == LINEAR_ELASTIC_ISOTROPIC_PROBLEM_CONTROLLED_DAMAGE)
				{
					for (id_type qp=0; qp<stresses.size(); ++qp)
					{
						std::vector<double> principal = Utilities::principal_stress_strain(stresses[qp]);
						if (principal[0] > _sigma_max)
						{
							elem_damaged = true;
							break;
						}
					}
				}

				if (elem_damaged)
				{
					_curr_order_decrease[l_elem][0] += _order_decrease_step; // Decrease the order of damage by the step amount
					local_damage_this_step = true;
					for (id_type qp=0; qp<stresses.size(); ++qp)
						get_internal_vars()->get_internal_var_local(l_elem, qp, 0) = 1.0 - 1.0/pow(10, _curr_order_decrease[l_elem][0]);
				}
			}
		}

		// This element is intersected
		else
		{
			// Compute the stresses and strains
			std::vector<std::vector<double> > strains, stresses;
			compute_stress_strain(strains, stresses, el);

			id_type initial_qp = 0;
			for (id_type ie=0; ie<el->n_integration_elem(); ++ie)
			{
				id_type nqp = el->get_integration_elem(ie)->n_q_points();

				// If this element has already been damaged, either update the damage or move on
				if (_curr_order_decrease[l_elem][ie] > 0)
				{
					if (_curr_order_decrease[l_elem][ie] >= _max_order_decrease) // As damaged as we're gonna allow
						continue;
					else
					{
						local_damage_this_step = true;
						id_type orders_remaining = _max_order_decrease - _curr_order_decrease[l_elem][ie];
						if (orders_remaining > _order_decrease_step) // I'll decrease the order by the step amount
							_curr_order_decrease[l_elem][ie] += _order_decrease_step;
						else 	// Otherwise, go to the maximum amount
							_curr_order_decrease[l_elem][ie] += orders_remaining;

						// Loop over all of the quadrature points and set the new damage state
						for (id_type qp=0; qp<nqp; ++qp)
							get_internal_vars()->get_internal_var_local(l_elem, initial_qp+qp, 0) = 1.0 - 1.0/pow(10, _curr_order_decrease[l_elem][ie]);
					}
				}
				else
				{
					bool elem_damaged = false;
					if (_mesh->get_element_material_local(l_elem, ie)->get_type() == LINEAR_ELASTIC_ISOTROPIC_PROBLEM_CONTROLLED_DAMAGE)
					{
						id_type nqp = el->get_integration_elem(ie)->n_q_points();
						for (id_type qp=0; qp<nqp; ++qp)
						{
							std::vector<double> principal = Utilities::principal_stress_strain(stresses[initial_qp+qp]);
							if (principal[0] > _sigma_max)
							{
								elem_damaged = true;
								break;
							}
						}
					}

					if (elem_damaged)
					{
						_curr_order_decrease[l_elem][ie] += _order_decrease_step; // Decrease the order of damage by the step amount
						local_damage_this_step = true;
						for (id_type qp=0; qp<nqp; ++qp)
							get_internal_vars()->get_internal_var_local(l_elem, initial_qp+qp, 0) = 1.0 - 1.0/pow(10, _curr_order_decrease[l_elem][ie]);
					}
				}

				initial_qp += nqp;
			}
		}
	}

	// Perform global reduction to find if there was damage anywhere
	bool damage_this_step;
	MPI_Allreduce(&local_damage_this_step, &damage_this_step, 1, MPI_C_BOOL, MPI_LOR, _mesh->get_comm());

	// Logic here:
	//	If the previous load step had damage, then I am the middle of so called "pseudo
	//		load steps", with a time step of 0. Its really just reachieving equilibrium.
	//		If this is the case then if there was also damage this step I need to continue
	//		repeating this zero-time load step. If there was no damage this step then I need
	//		to set the _prev_step_had_damage to false and revert back to the original solver dt
	//		and continue with normal load steps
	//	If the previous step had no damage and this step had no damage, continue like normal.
	//		If this step had damage though, I need to store the solver dt passed in (will revert
	//		back to this later) and start the so called "pseudo load step"
	if (_prev_step_had_damage)
	{
		if (damage_this_step)
		{
			_convergence_case = 0;
			return true;
		}
		else
		{
			if (_mesh->get_rank()==0) std::cout << "EQUILIBRIUM RE-ACHIEVED\n";
			_convergence_case = 1;
			_prev_step_had_damage = false;
			return false;
		}
	}
	else
	{
		if (damage_this_step)
		{
			if (_mesh->get_rank()==0) std::cout << "RE-ACHIEVING EQUILIBRIUM\n";
			_convergence_case = 2;
			_prev_step_had_damage = true;
			return false; // Don't actually repeat the load step because I need to maintain the current time
		}
		else
		{
			_convergence_case = 3;
			return false;
		}
	}
}


// Function to update the time step and return whether or not the problem is controlling the time step
bool ProblemMaxPrincipalStress::update_time_step(double& dt, double& rel_tol, double& abs_tol)
{
	double min_dt = 0.0;

	switch (_convergence_case)
	{
		case 0:	// Had damage last step and I still have damage happening
			dt = min_dt;
			rel_tol = _prob_rel_tol;
			abs_tol = _prob_abs_tol;
			return true;
		case 1:	// Had damage last step but no more happened this step. Leaving problem controlled load steps
			dt = _solver_dt;
			rel_tol = _solver_rel_tol;
			abs_tol = _solver_abs_tol;
			return false;
		case 2:	// Didn't have any damage before but now I do. Entering a problem controlled load step
			_solver_dt = dt; // Store the dt used to get to this point
			_solver_rel_tol = rel_tol;
			_solver_abs_tol = abs_tol;
			dt = min_dt;
			rel_tol = _prob_rel_tol;
			abs_tol = _prob_abs_tol;
			return true;
		case 3:	// No damage before, no damage now
			return false;
		default:
			err_message("Invalid convergence case for the Maximum Principal Stress problem.");
	}
}
