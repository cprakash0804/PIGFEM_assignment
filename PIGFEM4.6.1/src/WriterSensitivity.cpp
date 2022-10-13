/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#include "WriterSensitivity.h"
#include "Problem.h"
#include "Mesh.h"
#include "SensitivitySolver.h"

/*
 * The Big 3
 * Constructor need the Problem
 * Destructor doesn't do anything
 * Copy Constructor copies the problem pointer
 */
WriterSensitivity::WriterSensitivity(Problem* prob)
	: Writer(prob)
{}
WriterSensitivity::WriterSensitivity()
{}
WriterSensitivity::~WriterSensitivity()
{}
WriterSensitivity::WriterSensitivity(const WriterSensitivity& other)
{
	_prob = other.get_prob();
	store_mesh();
}



/*
 * Main function used to write problem to any inherited file type
 */
void WriterSensitivity::writeConsecutive(std::string filename, double curr_t, id_type step)
{
	if (_prob->sensitivity())
	{
		if (_prob->get_mesh()->get_rank() == 0)
		{
			// Open file for writing or appending base on the load step
			std::ofstream myfile;
			if (step==0)
				myfile.open(filename.c_str(), std::ofstream::out);
			else
				myfile.open(filename.c_str(), std::ofstream::app);
			if (!myfile.good())
			{
				std::string out = "Error opening the output file " + filename;
				err_message( out.data() );
			}
			myfile.precision(16);

			writeFromStream(myfile, curr_t);

			myfile.close();
		}
	}
}



/*
 * The actual function that will do all of the writing (maybe by calling other functions)
 */
void WriterSensitivity::writeFromStream(std::ofstream& myfile, double curr_t)
{
	myfile << curr_t << ",";

	for (id_type f=0; f<_prob->n_SensitivityFunctions(); ++f)
	{
		if (f!=0)
			myfile << ",";
		myfile << _prob->get_sensitivity_solver()->get_function_val(f);
		for (id_type p=0; p<_prob->n_SensitivityParameters(); ++p)
			myfile << "," << _prob->get_sensitivity_solver()->get_sensitivity(f, p);	
	}
	myfile << "\n";
}