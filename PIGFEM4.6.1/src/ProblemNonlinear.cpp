/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "ProblemNonlinear.h"
#include "SolverNonlinearKSP.h"
#include "SolverNonlinearSNES.h"




/*
 * Constructor, doesn't do anything
 */
ProblemNonlinear::ProblemNonlinear()
{}




/*
 * Set the solver to the nonlinear solver
 */
void ProblemNonlinear::generate_solver()
{
	_solver = new SolverNonlinearSNES;
}




/*
 * Nothing special for the init function
 */
