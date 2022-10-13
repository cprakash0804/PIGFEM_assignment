/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "ProblemNonlinearStructural_Shrinkage.h"
#include "AssemblerSmallStrainNonlinear_Shrinkage.h"
#include "Mesh.h"
#include "BoundaryObject.h"
#include "Utilities.h"
#include "Solver.h"

ProblemNonlinearStructural_Shrinkage::ProblemNonlinearStructural_Shrinkage()
	: _plane_strain(false), _delta_T(0.0), _delta_T_set(false), _apply_loads(true)
{
}

// This function will be overridden in derived classes to set the actual assembler and solver that will be used
void ProblemNonlinearStructural_Shrinkage::generate_assembler()
{
	_assembler = new AssemblerSmallStrainNonlinear_Shrinkage;
	if (_delta_T_set)
		_assembler->set_parameter("DELTA_T", _delta_T);
}

// This function return the number of degrees of freedom present on each node
id_type ProblemNonlinearStructural_Shrinkage::nndof()
{
	return _mesh->dim();
}

// Returns the problem type enum for this problem
problem_type ProblemNonlinearStructural_Shrinkage::get_type()
{
	return NONLINEAR_STRUCTURAL_SHRINKAGE;
}

// Set any possible problem specific parameters
void ProblemNonlinearStructural_Shrinkage::set_parameter(std::string name, double val)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if (name=="PLANE_STRAIN" || name=="PLANE STRAIN")
	{
		if (val==0.0) // stupid hack because virtual functions can't be templated
			_plane_strain = false;
		else
			_plane_strain = true;
	}
	else if (name=="APPLY_LOADS" || name=="APPLY LOADS" || name=="LOAD")
	{
		if (val==0.0) // stupid hack because virtual functions can't be templated
			_apply_loads = false;
		else
			_apply_loads = true;
	}
	else if (name=="TEMP" || name=="DELTA_T" || name=="DELTA_TEMP")
	{
		if (_assembler!=NULL)
			_assembler->set_parameter(name, val);
		else
		{
			_delta_T = val;
			_delta_T_set = true;
		}
	}
	else
		err_message(name << " is not a valid parameter name");
}
double ProblemNonlinearStructural_Shrinkage::get_parameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if (name=="PLANE_STRAIN" || name=="PLANE STRAIN")
	{
		if (_plane_strain)
			return 1.0;
		else
			return 0.0;
	}
	else if (name=="APPLY_LOADS" || name=="APPLY LOADS" || name=="LOAD")
	{
		if (_apply_loads) // stupid hack because virtual functions can't be templated
			return 1.0;
		else
			return 0.0;
	}
	else if (name=="TEMP" || name=="DELTA_T" || name=="DELTA_TEMP")
		return _assembler->get_parameter(name);
	else if (name=="APPLY_TIME" || name=="APPLICATION_TIME")
		return _assembler->get_parameter(name);
	else
		err_message(name << " is not a valid parameter name");
}



void ProblemNonlinearStructural_Shrinkage::fillBmat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& grad_x)
{
	ProblemUtilities::Assemble_Small_Strain_B_Mat(B, grad_x);
}



// This problem needs a special solve sequence as it first needs to relax the structure and then load is applied
void ProblemNonlinearStructural_Shrinkage::solve_problem()
{
	// Call the initialization routine
	if (!_init)
		err_message("Please initialize the problem before calling solve!");

	// Actually solve the problem!
	Utilities::timer timer;
	timer.start();
	PIGFEMPrint("Running the relaxation phase of the residual strain problem");

	// Apply the traction-free boundary conditions for the relaxation phase
	BoundaryObject* current_boundary = get_boundary();
	BoundaryObject* traction_free_boundary = new BoundaryObject();
	traction_free_boundary->attachProblem(this);
	if ( !_mesh->generated() )
		err_message("Can't apply constrained traction free boundary conditions to a mesh that was not generated");
	else
	{
		// Add the tractoin free BCs
		id_type top_left_node = (_mesh->nElemInDim(1)+1) * _mesh->nElemInDim(2);
		for(unsigned int d=0; d<_mesh->dim(); ++d)
			traction_free_boundary->set_dirichlet_bc(0, d, 0);
		if (_mesh->dim()>=2)
			traction_free_boundary->set_dirichlet_bc(top_left_node, 0, 0.0);
		if (_mesh->dim()==3)
			traction_free_boundary->set_dirichlet_bc(top_left_node, 2, 0.0);

		// temporarily replace the boundarObject in the problme with the traction free one
		_boundary = traction_free_boundary;
	}

	// Apply solver preferences for the relaxaion solve
	double T_final = get_solver()->get_dparameter("T_FINAL");
	double max_dt = get_solver()->get_dparameter("MAX_TIME_STEP");
	get_solver()->set_dparameter("MAX_TIME_STEP", T_final); // One load step

	// Actually do the solve
	setup();
	get_solver()->solve();
	_assembler->set_parameter("COMPLETE", 1.0);
	








	// Restore boundary conditions and solver to settings before the traction free solution
	_boundary = current_boundary;
	delete traction_free_boundary;
	get_solver()->set_dparameter("MAX_TIME_STEP", max_dt); // One load step




	










	// Actually solve the loading problem!
	if (_apply_loads)
	{
		// Transfer the solution
		setup(); // Rebuild the finite element system for the new boundary conditions
		get_solver()->initialSolution(*get_solution());

		// Apply the loading boundary conditions
		PIGFEMPrint("Running the loading phase of the residual strain problem");
		
		get_solver()->solve();
	}
	_solved = true;

	double total_solve_time = timer.time_elapsed();
	_misc_solve_time = total_solve_time - _assemble_time - _solve_time - _solver_write_time - _sensitivity_time;


	writeStatistics();
}