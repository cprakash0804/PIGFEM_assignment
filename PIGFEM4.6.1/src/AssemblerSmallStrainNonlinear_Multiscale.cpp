/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "AssemblerSmallStrainNonlinear_Multiscale.h"
#include "Problem.h"
#include <algorithm>



AssemblerSmallStrainNonlinear_Multiscale::AssemblerSmallStrainNonlinear_Multiscale()
	: _strain_set(false)
{
}

AssemblerSmallStrainNonlinear_Multiscale::~AssemblerSmallStrainNonlinear_Multiscale()
{
	_max_macro_strain.clear();
	_curr_macro_strain.clear();
	_strain_set = false;
}
void AssemblerSmallStrainNonlinear_Multiscale::clear()
{
	_max_macro_strain.clear();
	_curr_macro_strain.clear();
	_strain_set = false;
	Assembler::clear();
}


bool AssemblerSmallStrainNonlinear_Multiscale::storedBmats()
{
	return true;
}


// Kernel function that actually does the computes the physics behind the specific problem
// Computes the local contribution to the stiffness matrix and internal load vector
// Also updates the internal variable storage in the input variable
//	B-matrix variant
void AssemblerSmallStrainNonlinear_Multiscale::KernelNonlinear(DenseMatrix<double>& K_el, std::vector<double>& P_el_int,
													const std::vector<double>& shape, const DenseMatrix<double>& B,
													Material* mat, Material::input_params& input,
													std::vector<double>& elem_U_curr, bool assembleFunc, bool assembleJac)
{
	// Compute the current strain
	std::vector<double>& strain = input.strain;
	ProblemUtilities::SmallStrainFast(strain, B, elem_U_curr);

	// Add the current macroscopic strain to it
	Utilities::_VecAXPY(1.0, _curr_macro_strain, strain); // strain = 1.0*thermal_strain + strain

	// Compute the constitutive relations
	Material::output_params* output = mat->Constitutive(input);

	// Compute the stiffness matrix
	if (assembleJac)
		ProblemUtilities::SmallStrainKFast(K_el, B, output->Dmat);

	// Compute the internal force
	if (assembleFunc)
		ProblemUtilities::SmallStrainInternalForceFast(P_el_int, B, output->stress);
}

void AssemblerSmallStrainNonlinear_Multiscale::KernelNonlinear(DenseMatrix<double>& K_el, std::vector<double>& P_el_int,
													const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
													Material* mat, Material::input_params& input,
													std::vector<double>& elem_U_curr, bool assembleFunc, bool assembleJac)
{
	// Compute the current strain
	std::vector<double>& strain = input.strain;
	ProblemUtilities::SmallStrainFast(strain, shape_grad, elem_U_curr);

	// Add the current macroscopic strain to it
	Utilities::_VecAXPY(1.0, _curr_macro_strain, strain); // strain = 1.0*thermal_strain + strain

	// Compute the constitutive relations
	Material::output_params* output = mat->Constitutive(input);

	// Compute the stiffness matrix
	if (assembleJac)
		ProblemUtilities::SmallStrainKFast(K_el, shape_grad, output->Dmat);

	// Compute the internal force
	if (assembleFunc)
		ProblemUtilities::SmallStrainInternalForceFast(P_el_int, shape_grad, output->stress);
}




void AssemblerSmallStrainNonlinear_Multiscale::KernelCohesive(DenseMatrix<double>& K_coh, std::vector<double>& P_int_coh,
															 const std::vector<double>& shape, const DenseMatrix<double>& rotation_matrix,
															 Material* mat, Material::input_params& input,
															 std::vector<double>& coh_U_curr, bool assembleFunc, bool assembleJac)
{
	// Compute the current opening vector in the ntt coordinate system
	ProblemUtilities::cohesiveDeltaFast(input.delta, shape, coh_U_curr, rotation_matrix);

	// Compute the constitutive relations
	Material::output_params* output = mat->Constitutive(input);

	// Compute the stiffness matrix
	if (assembleJac)
		ProblemUtilities::cohesiveKFast(K_coh, shape, output->Dmat, rotation_matrix);

	// Compute the internal force
	if (assembleFunc)
		ProblemUtilities::cohesiveInternalForceFast(P_int_coh, shape, output->traction, rotation_matrix);
}



// Select the appropriate B-matrix fill function
void AssemblerSmallStrainNonlinear_Multiscale::fillBmat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& shape_grad)
{
	ProblemUtilities::Assemble_Small_Strain_B_Mat(B, shape_grad);
}







// Function to assemble the presribed displacements and free loads at the beginning of each load step
// (Or at the beginning of the analysis for linear)
// This is made virtual so that calling it can update the assembler state at every load step is desired
PetscErrorCode AssemblerSmallStrainNonlinear_Multiscale::assemble_new_load_step(Vec& Up, pvector& F_ext, double current_time_frac)
{
	// Make sure I've assigned a strain value
	if (!_strain_set)
		err_message("Please set a macroscopic strain value before attemping to solve a multiscale problem.");

	// Update the current time here
	_curr_macro_strain.resize(_max_macro_strain.size());
	for (id_type i=0; i<_max_macro_strain.size(); ++i)
		_curr_macro_strain[i] = current_time_frac*_max_macro_strain[i];

	// Call the base class function to do all of the same things
	return Assembler::assemble_new_load_step(Up, F_ext, current_time_frac);
}



/*
 * Function used to easily set a vector-valued parameter
 */
void AssemblerSmallStrainNonlinear_Multiscale::set_vec_parameter(std::string name, std::vector<double> val)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="STRAIN" || name=="MACRO_STRAIN" || name=="MACRO STRAIN")
	{
		_max_macro_strain = val;
		_curr_macro_strain.clear(); _curr_macro_strain.resize(_max_macro_strain.size()); // Set current strain to 0
		_strain_set = true;
	}
	else
		err_message(name << " is not a valid parameter name");
}
std::vector<double> AssemblerSmallStrainNonlinear_Multiscale::get_vec_parameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="STRAIN" || name=="MACRO_STRAIN" || name=="MACRO STRAIN")
	{
		if (_strain_set)
			return _max_macro_strain;
		else
			err_message("Please set a strain value before attempting to get it.");
	}
	else if (name=="CURRENT STRAIN" || name=="CURRENT_STRAIN")
		return _curr_macro_strain;
	else
		err_message(name << " is not a valid parameter name");
}