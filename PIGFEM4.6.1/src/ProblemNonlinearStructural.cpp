/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated August 2017

##################################################################################
*/
#include "ProblemNonlinearStructural.h"
#include "AssemblerSmallStrainNonlinear.h"
#include "Mesh.h"
#include "Utilities.h"

ProblemNonlinearStructural::ProblemNonlinearStructural()
	: _plane_strain(true)
{}

// This function will be overridden in derived classes to set the actual assembler and solver that will be used
void ProblemNonlinearStructural::generate_assembler()
{
	_assembler = new AssemblerSmallStrainNonlinearStructural;
}

// This function return the number of degrees of freedom present on each node
id_type ProblemNonlinearStructural::nndof()
{
	return _mesh->dim();
}

// Returns the problem type enum for this problem
problem_type ProblemNonlinearStructural::get_type()
{
	return NONLINEAR_STRUCTURAL;
}

// Set any possible problem specific parameters
void ProblemNonlinearStructural::set_parameter(std::string name, double val)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if (name=="PLANE_STRAIN" || name=="PLANE STRAIN")
	{
		if (val==0.0) // stupid hack because virtual functions can't be templated
			_plane_strain = false;
		else
			_plane_strain = true;
	}
	else
		err_message("Please input a valid parameter name.");
}

// Get whatever possible problem specific parameter
double ProblemNonlinearStructural::get_parameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if (name=="PLANE_STRAIN" || name=="PLANE STRAIN")
	{
		if (_plane_strain)
			return 1.0;
		else
			return 0.0;
	}
	else
		err_message("Please input a valid parameter name.");
}


void ProblemNonlinearStructural::fillBmat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& grad_x)
{
	ProblemUtilities::Assemble_Small_Strain_B_Mat(B, grad_x);
}