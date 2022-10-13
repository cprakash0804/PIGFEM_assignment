/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "ProblemLinear.h"
#include "SolverLinear.h"
#include "Mesh.h"
#include <iostream>
#include <fstream>




/*
 * Constructor, doesn't do anything
 */
ProblemLinear::ProblemLinear()
{}




/*
 * Set the solver to the linear solver
 */
void ProblemLinear::generate_solver()
{
	_solver = new PETScLinearKSPSolver;
}




/*
 * Initializer. Makes sure all materials are linear
 * Then calls basic Problem initializer
 */
void ProblemLinear::init()
{
	// Make sure that all materials are linear
	std::vector<Material*> mats = _mesh->get_materials();
	for(id_type m=0; m<mats.size(); ++m)
		if(!mats[m]->linear())
			err_message("Your mesh contains a nonlinear material in a linear analysis!");

	// Otherwise, we're fine so run the rest of the initialization process
	Problem::init();
}
