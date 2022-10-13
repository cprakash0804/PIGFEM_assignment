#include "Options.h"
#include "mpi.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>

bool writeToScreen = false; // Initial value that will never be used hopefully
std::string logFile = "Output/wrong.txt";


// Set all default options
// The only required one is _in_file
Options::Options()
{
	_string_options["-IN"] = "none";
	_string_options["-OUT_FOLDER"] = "Output/";
	_string_options["-VTK_FILE"] = "History";
	_string_options["-SENS_FILE"] = "sensitivities.csv";
	_string_options["-STATS_FILE"] = "statistics.txt";
	_string_options["-SS_FILE"] = "stress-strain.csv";
	_string_options["-MESH"] = "";
	_string_options["-LOG"] = "log.txt";
	_string_options["-COH_FAIL_FILE"] = "cohesive_failure_lengths.csv";
	_string_options["-PATH_IN"] = "";
	_string_options["-PATH_OUT"] = "path";
	_string_options["-COH_OPENINGS_OUT"] = "coh_openings";

	_bool_options["-WRITE_VTK"] = true;
	_bool_options["-OUTPUT_COH"] = false;
	_bool_options["-OUTPUT_COH_FAIL"] = false;
	_bool_options["-OUTPUT_ALL_COH_FAIL"] = false;
	_bool_options["-WRITE_PATH"] = false;
	_bool_options["-WRITE_COH_OPENINGS"] = false;
	_bool_options["-SUBSCALE"] = false;
	_bool_options["-SENSITIVITY"] = false;
	_bool_options["-SENSITIVITY_FINAL_ONLY"] = false;
	_bool_options["-OUTPUT_FINAL_ONLY"] = false;
	_bool_options["-WRITE_TO_SCREEN"] = true;

	// Global control variable defintions
	writeToScreen = true;
	logFile = "Output/log.txt";
}






// Function to read the command line input
void Options::parseInput(int argc, char* argv[])
{
	if (argc%2 == 0)
		err_message("Every input parameter must have a matching flag");

	// First entry is the name of the executable
	bool out_specified = false;
	for (int i=1; i<argc; i+=2)
	{
		std::string flag(argv[i]);
		std::string val(argv[i+1]);
		std::transform(flag.begin(), flag.end(), flag.begin(), ::toupper);

		if (flag=="-OUT_FOLDER")
		{
			setOutputFolder(val);
			out_specified = true;
		}
		else if (flag=="-IN")
			setInputFile(val);
		else if (flag=="-VTK_FILE")
			setVTKPrefix(val);
		else if (flag=="-SENS_FILE")
			setSensFile(val);
		else if (flag=="-STATS_FILE")
			setStatsFile(val);
		else if (flag=="-SS_FILE")
			setSSFile(val);
		else if (flag=="-MESH")
			setMeshFile(val);
		else if (flag=="-LOG")
			setLogFile(val);
		else if (flag=="-COH_FAIL_FILE")
			setCohFailureFile(val);
		else if (flag=="-PATH_IN")
			setPathInput(val);
		else if (flag=="-PATH_OUT")
			setPathOutput(val);
		else if (flag=="-COH_OPENINGS_OUT")
			setCohOpeningsOutput(val);
		else if (flag=="-WRITE_VTK")
			setBoolOption(flag, val);
		else if (flag=="-OUTPUT_COH")
			setBoolOption(flag, val);
		else if (flag=="-OUTPUT_COH_FAIL")
			setBoolOption(flag, val);
		else if (flag=="-OUTPUT_ALL_COH_FAIL")
			setBoolOption(flag, val);
		else if (flag=="-WRITE_PATH")
			setBoolOption(flag, val);
		else if (flag=="-WRITE_COH_OPENINGS")
			setBoolOption(flag, val);
		else if (flag=="-SUBSCALE")
			setBoolOption(flag, val);
		else if (flag=="-SENSITIVITY")
			setBoolOption(flag, val);
		else if (flag=="-SENSITIVITY_FINAL_ONLY")
			setBoolOption(flag, val);
		else if (flag=="-OUTPUT_FINAL_ONLY")
			setBoolOption(flag, val);
		else if (flag=="-WRITE_TO_SCREEN")
			setBoolOption(flag, val);
		else // Unknown option, assume it's a string
		{
			_string_options[flag] = val;
		}
	}

	// If the output folder wasn't specified, I need to call it just to make sure the default exists
	if (!out_specified)
		setOutputFolder(_string_options["-OUT_FOLDER"]);

	// If we're suposed to be writing the path, we need to have specified an input file
	if (getBoolOption("-WRITE_PATH"))
		if (getOption("-PATH_IN")=="")
			err_message("To write a path solution output, an input file of points must be specified");
}






std::string Options::getOption(std::string flag)
{
	std::map<std::string, std::string>::iterator it = _string_options.find(flag);
	if (it != _string_options.end())
		return it->second;
	else
	{
		std::string out = "Command line parameter " + flag + " was not set";
		err_message(out.data());
	}
}
bool Options::getBoolOption(std::string flag)
{
	std::map<std::string, bool>::iterator it = _bool_options.find(flag);
	if (it != _bool_options.end())
		return it->second;
	else
	{
		std::string out = "Command line parameter " + flag + " was not set";
		err_message(out.data());
	}
}
























// Functions to set options
void Options::setBoolOption(std::string flag, std::string val)
{
	std::transform(val.begin(), val.end(), val.begin(), ::toupper);
	if (val=="FALSE" || val=="F" || val=="NO" || val=="N")
		_bool_options[flag] = false;
	else if (val=="TRUE" || val=="T" || val=="YES" || val=="Y")
		_bool_options[flag] = true;
	else
	{
		std::string out = "Invalid value " + val + " for a boolean option.";
		err_message(out.data());
	}
	if (flag=="-WRITE_TO_SCREEN")
	{
		if (val=="FALSE" || val=="F" || val=="NO" || val=="N")
			writeToScreen = false;
		else if (val=="TRUE" || val=="T" || val=="YES" || val=="Y")
			writeToScreen = true;
		else
		{
			std::string out = "Invalid value " + val + " for a boolean option.";
			err_message(out.data());
		}
	}
}
// These functions have checks to make sure that they are valid options
void Options::setOutputFolder(std::string folder)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	// Make it a subfolder of Output
	if (folder.find("Output/", 0) == std::string::npos)
		folder = "Output/" + folder;
	if (rank==0)
	{
		// Check if this folder exists already and create it if it doesn't
		// NOTE: This is not portable at all... This will only work on some Linux systems
		struct stat st;
		if (stat(folder.data(), &st) == 0) // found something with this file
		{
			if (!S_ISDIR(st.st_mode)) // What we found isn't a directory...
				err_message("Invalid output directory! Output directory is already a file");
		}
		else // Didn't find anything, have to create it
		{
			if (mkdir(folder.data(), 0777) != 0) // Failed to create the folder
				err_message("Error making the output directory!");
		}
	}

	// Add a "/" to the end if it isn't there
	if (folder.find("/", folder.length()-1) == std::string::npos)
		folder = folder + "/";

	logFile = folder + _string_options["-LOG"];

	_string_options["-OUT_FOLDER"] = folder;
}
void Options::setVTKPrefix(std::string prefix)
{
	// Nothing really special to check here that I know of
	_string_options["-VTK_FILE"] = prefix;
}
void Options::setSensFile(std::string sens_file)
{
	// Check if the the end of the filefile is ".csv"
	if (sens_file.find(".csv", sens_file.length()-4) == std::string::npos)
		sens_file = sens_file + ".csv";

	_string_options["-SENS_FILE"] = sens_file;
}
void Options::setSSFile(std::string ss_file)
{
	// Check if the the end of the filefile is ".csv"
	if (ss_file.find(".csv", ss_file.length()-4) == std::string::npos)
		ss_file = ss_file + ".csv";

	_string_options["-SS_FILE"] = ss_file;
}
void Options::setStatsFile(std::string stats_file)
{
	// Check if the the end of the filefile is ".txt"
	if (stats_file.find(".txt", stats_file.length()-4) == std::string::npos)
		stats_file = stats_file + ".txt";

	_string_options["-STATS_FILE"] = stats_file;
}
void Options::setInputFile(std::string input)
{
	// Add the Input folder to the beginning of the filename
	if (input.find("Input/") != 0)
		input = "Input/" + input;
	// I don't want to restrict the kind of file extensions i CAN HAVE HERE
	_string_options["-IN"] = input;
}
void Options::setMeshFile(std::string mesh_file)
{
	// Again, nothing special. Might be able to read things beside ".inp" files later
	_string_options["-MESH"] = mesh_file;
}
void Options::setLogFile(std::string log_file)
{
	_string_options["-LOG"] = log_file;
	logFile = _string_options["-OUT_FOLDER"] + log_file;
}
void Options::setCohFailureFile(std::string coh_fail_file)
{
	// Check if the the end of the filefile is ".csv"
	if (coh_fail_file.find(".csv", coh_fail_file.length()-4) == std::string::npos)
		coh_fail_file = coh_fail_file + ".csv";

	_string_options["-COH_FAIL_FILE"] = coh_fail_file;
}
void Options::setPathInput(std::string path_in)
{
	// Set the Input part of the input file
	if (path_in.find("Input/") != 0)
		path_in = "Input/" + path_in;

	// Should probably actually check that this file exists
	_string_options["-PATH_IN"] = path_in;
}
void Options::setPathOutput(std::string path_out)
{
	if (path_out.find(".csv", path_out.length()-4) != std::string::npos)
		path_out = path_out.substr(0, path_out.length()-4);

	_string_options["-PATH_OUT"] = path_out;
}
void Options::setCohOpeningsOutput(std::string coh_openings_out)
{
	if (coh_openings_out.find(".csv", coh_openings_out.length()-4) != std::string::npos)
		coh_openings_out = coh_openings_out.substr(0, coh_openings_out.length()-4);

	_string_options["-COH_OPENINGS_OUT"] = coh_openings_out;
}