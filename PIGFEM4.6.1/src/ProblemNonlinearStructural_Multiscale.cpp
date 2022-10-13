/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "ProblemNonlinearStructural_Multiscale.h"
#include "AssemblerSmallStrainNonlinear_Multiscale.h"
#include "Mesh.h"
#include "Utilities.h"

ProblemNonlinearStructural_Multiscale::ProblemNonlinearStructural_Multiscale()
	: _plane_strain(false)
{
	_subscale = true;
}

// This function will be overridden in derived classes to set the actual assembler and solver that will be used
void ProblemNonlinearStructural_Multiscale::generate_assembler()
{
	_assembler = new AssemblerSmallStrainNonlinear_Multiscale;
	if (_strain_set)
		_assembler->set_vec_parameter("macro strain", _max_macro_strain);
}

// This function return the number of degrees of freedom present on each node
id_type ProblemNonlinearStructural_Multiscale::nndof()
{
	return _mesh->dim();
}

// Returns the problem type enum for this problem
problem_type ProblemNonlinearStructural_Multiscale::get_type()
{
	return NONLINEAR_STRUCTURAL_MULTISCALE;
}

// Set any possible problem specific parameters
void ProblemNonlinearStructural_Multiscale::set_parameter(std::string name, double val)
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
		err_message(name << " is not a valid parameter name");
}
double ProblemNonlinearStructural_Multiscale::get_parameter(std::string name)
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
		err_message(name << " is not a valid parameter name");
}


/*
 * Function used to easily set a vector-valued parameter
 */
void ProblemNonlinearStructural_Multiscale::set_vec_parameter(std::string name, std::vector<double> val)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="STRAIN" || name=="MACRO_STRAIN" || name=="MACRO STRAIN")
	{
		if (_assembler!=NULL)
			_assembler->set_vec_parameter(name, val);
		else
		{
			_max_macro_strain = val;
			_strain_set = true;
		}
	}		
	else
		err_message(name << " is not a valid parameter name");
}
std::vector<double> ProblemNonlinearStructural_Multiscale::get_vec_parameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="STRAIN" || name=="MACRO_STRAIN" || name=="MACRO STRAIN")
		return _assembler->get_vec_parameter(name);
	else
		err_message(name << " is not a valid parameter name");
}



void ProblemNonlinearStructural_Multiscale::fillBmat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& grad_x)
{
	ProblemUtilities::Assemble_Small_Strain_B_Mat(B, grad_x);
}