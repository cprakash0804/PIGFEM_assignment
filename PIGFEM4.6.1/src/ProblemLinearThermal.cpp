/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "ProblemLinearThermal.h"
#include "AssemblerLinearThermal.h"
#include "Utilities.h"

ProblemLinearThermal::ProblemLinearThermal()
{}

// This function will be overridden in derived classes to set the actual assembler and solver that will be used
void ProblemLinearThermal::generate_assembler()
{
	_assembler = new AssemblerLinearThermal;
}

// This function return the number of degrees of freedom present on each node
id_type ProblemLinearThermal::nndof()
{
	return 1;
}

// Returns the problem type enum for this problem
problem_type ProblemLinearThermal::get_type()
{
	return LINEAR_HEAT;
}

void ProblemLinearThermal::fillBmat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& grad_x)
{
	ProblemUtilities::Assemble_Thermal_B_Mat(B, grad_x);
}