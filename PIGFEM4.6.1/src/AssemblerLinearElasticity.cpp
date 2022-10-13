/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "AssemblerLinearElasticity.h"
#include "Problem.h"




bool AssemblerLinearElasticity::storedBmats()
{
	return false;
}


// Kernel function that actually does the computes the physics behind the specific problem
// Computes the local contribution to the stiffness matrix and internal load vector
// Also updates the internal variable storage in the input variable
void AssemblerLinearElasticity::KernelLinear(DenseMatrix<double>& K_el,
										  const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
										  Material* mat, Material::input_params& input)
{
	// If I wanted to actually get the stress I would compute the strain here. It doesn't really matter though
	// I also don't have the current solution so whatever

	// Compute material constitutive matrix
	Material::output_params* output = mat->Constitutive(input);   // Computes the D-matrix for the current material
	DenseMatrix<double>& D = output->Dmat;

	// ProblemUtilities::SmallStrainKFast(K_el, B, D);
	ProblemUtilities::SmallStrainKFast(K_el, shape_grad, D);
}
