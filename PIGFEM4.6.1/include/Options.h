/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#ifndef _OPTIONS_H_
#define _OPTIONS_H_
#include "common.h"
#include <string>
#include <map>


/*
 * Class that just contains a bunch of options for running a PCIGFEM simulation
 */
class Options
{
	private:

		std::map<std::string, std::string> _string_options;
		std::map<std::string, bool> _bool_options;

	public:
		Options();

		// Function to read the command line input
		void parseInput(int argc, char* argv[]);

		std::string getOption(std::string flag);
		bool getBoolOption(std::string flag);

		// Functions to set options
		// These functions have checks to make sure that they are valid options
		void setBoolOption(std::string flag, std::string val);
		void setOutputFolder(std::string folder);
		void setVTKPrefix(std::string prefix);
		void setSensFile(std::string sens_file);
		void setSSFile(std::string ss_file);
		void setStatsFile(std::string stats_file);
		void setInputFile(std::string input);
		void setMeshFile(std::string mesh_file);
		void setLogFile(std::string mesh_file);
		void setCohFailureFile(std::string coh_fail_file);
		void setPathInput(std::string path_in);
		void setPathOutput(std::string path_out);
		void setCohOpeningsOutput(std::string path_out);
};

#endif