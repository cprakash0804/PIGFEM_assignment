/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "ProblemLinearElasticity.h"
#include "AssemblerLinearElasticity.h"
#include "Mesh.h"
#include "Utilities.h"

ProblemLinearElasticity::ProblemLinearElasticity()
	: _plane_strain(true)
{}

// This function will be overridden in derived classes to set the actual assembler and solver that will be used
void ProblemLinearElasticity::generate_assembler()
{
	_assembler = new AssemblerLinearElasticity;
}

// This function return the number of degrees of freedom present on each node
id_type ProblemLinearElasticity::nndof()
{
	return _mesh->dim();
}

// Returns the problem type enum for this problem
problem_type ProblemLinearElasticity::get_type()
{
	return LINEAR_ELASTICITY;
}

// Set any possible problem specific parameters
void ProblemLinearElasticity::set_parameter(std::string name, double val)
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
double ProblemLinearElasticity::get_parameter(std::string name)
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



void ProblemLinearElasticity::fillBmat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& grad_x)
{
	ProblemUtilities::Assemble_Small_Strain_B_Mat(B, grad_x);
}