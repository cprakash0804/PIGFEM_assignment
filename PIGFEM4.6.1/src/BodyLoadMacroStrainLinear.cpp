#include "BodyLoadMacroStrainLinear.h"
#include "Utilities.h"



/*
 * Main function that each Boy force will have to implement
 * Takes the (x, y, z) coordinates and returns the 
 */
void BodyLoadMacroStrainLinear::BF_Kernel(std::vector<double>& P_el, const std::vector<double>& points,
										  const std::vector<double>& N, const std::vector<std::vector<double> >& dN,
										  const double& J, const double& W,
										  Material* mat, Material::input_params& input)
{
	// Get the local D matrix
	Material::output_params* output = mat->Constitutive(input);   // Computes the D-matrix for the current material
	DenseMatrix<double>& D = output->Dmat;

	// Compute -Bt*D*\epsilon_bar
	std::vector<double> ret1 = D * _macro_strain;
	double wgt = W*J;
	for (id_type i=0; i<ret1.size(); ++i)
		ret1[i] *= -1.0 * wgt;
	std::vector<double> ret;
	ProblemUtilities::SmallStrainInternalForceFast(P_el, dN, ret1); // Computes Bt * ret1
}


/*
 * Function used to easily set a vector-valued parameter
 */
void BodyLoadMacroStrainLinear::set_vec_parameter(std::string name, std::vector<double> val)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="STRAIN" || name=="MACRO_STRAIN" || name=="MACRO STRAIN")
		_macro_strain = val;
	else
		err_message("Please input a valid parameter name.");
}
std::vector<double> BodyLoadMacroStrainLinear::get_vec_parameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="STRAIN" || name=="MACRO_STRAIN" || name=="MACRO STRAIN")
		return _macro_strain;
	else
		err_message("Please input a valid parameter name.");
}
