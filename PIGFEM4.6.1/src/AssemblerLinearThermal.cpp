/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "AssemblerLinearThermal.h"



bool AssemblerLinearThermal::storedBmats()
{
	return false;
}

// Kernel function that actually does the computes the physics behind the specific problem
// Computes the local contribution to the stiffness matrix and internal load vector
// Also updates the internal variable storage in the input variable
void AssemblerLinearThermal::KernelLinear(DenseMatrix<double>& K_el,
										  const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
										  Material* mat, Material::input_params& input)
{
	// Declare some matrices
	DenseMatrix<double> B, Bt;

	// Define the strain-displacement matrix
	ProblemUtilities::Assemble_Thermal_B_Mat(B, shape_grad);
	B.transpose(Bt);

	// Compute material constitutive matrix
	// Don't actually need to do this here since the problem is linear but it helps with
	//	determining whether or not the problme is plane strain or not for structural problems
	Material::output_params* output = mat->Constitutive(input);   // Computes the D-matrix for the current material
	DenseMatrix<double>& D = output->Dmat;

	// Compute Bt*D*B
	K_el = Bt*D*B;
}
