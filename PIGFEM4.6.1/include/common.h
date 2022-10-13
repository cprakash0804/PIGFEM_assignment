/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _COMMON_H_
#define _COMMON_H_

#include <stdio.h>
#include <execinfo.h>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <climits>
#include <vector>
#include <sstream>
#include "petscmat.h"
#include "petscvec.h"
#include "mpi.h"



// Define structs that will be used in the solution of the partitioned fintite element system
// =====================================================================================
struct pmatrix
{
	Mat ff;
	Mat fp;
	Mat pf;
	Mat pp;

	pmatrix()
	{
		ff = PETSC_NULL;
		fp = PETSC_NULL;
		pf = PETSC_NULL;
		pp = PETSC_NULL;
	}
	void destroy()
	{
		if (ff != PETSC_NULL)
		{
			MatDestroy(&ff);
			ff = PETSC_NULL;
		}
		if (fp != PETSC_NULL)
		{
			MatDestroy(&fp);
			fp = PETSC_NULL;
		}
		if (pf != PETSC_NULL)
		{
			MatDestroy(&pf);
			pf = PETSC_NULL;
		}
		if (pp != PETSC_NULL)
		{
			MatDestroy(&pp);
			pp = PETSC_NULL;
		}
	}
	~pmatrix()
	{
		destroy();
	}
};
struct pvector
{
	Vec f;
	Vec p;

	pvector()
	{
		f = PETSC_NULL;
		p = PETSC_NULL;
	}
	void destroy()
	{
		if (f != PETSC_NULL)
		{
			VecDestroy(&f);
			f = PETSC_NULL;
		}
		if (p != PETSC_NULL)
		{
			VecDestroy(&p);
			p = PETSC_NULL;
		}
	}
	~pvector()
	{
		destroy();
	}
};



// Define any types that are used globally throughout the code
// =====================================================================================

enum elem_type {
	INVALID_ELEM,
	POINT1,
	EDGE2,
	EDGE3,
	TRI3,
	TRI6,
	QUAD4,
	QUAD8,
	TET4, 
	TET10,
	HEX8,
	HEX20,
	HEX27,
	PRISM6,
	PRISM15,
	PRISM18,
	POLYGON};

enum coh_elem_type {
	INVALID_COHELEM,
	COHPOINT1,
	COHEDGE2,
	COHEDGE3,
	COHTRI3,
	COHTRI6,
	COHQUAD4,
	COHQUAD8};
	
enum material_type {
	INVALID_MATERIAL,
	LINEAR_ELASTIC_ISOTROPIC,
	LINEAR_THERMAL,
	CONTINUUM_DAMAGE,
	OP_COHESIVE,
	OP_COHESIVE_NO_UNLOADING,
	LINEAR_ELASTIC_ISOTROPIC_MAX_PRINCIPAL_STRESS,
	LINEAR_ELASTIC_TRANSVERSELY_ISOTROPIC,
	LINEAR_ELASTIC_ISOTROPIC_PROBLEM_CONTROLLED_DAMAGE,
	XN_COHESIVE,
	XN_COHESIVE_NO_UNLOADING};

enum classification {
	INVALID_CLASSIFICATION,
	STRUCTURAL,
	THERMAL};

enum inclusion_type {
	INVALID_INCLUSION,
	PLANE,
	ELLIPSE,
	ELLIPSOID,
	CIRCLE,
	SPHERE};

enum problem_type {
	LINEAR_ELASTICITY,
	NONLINEAR_STRUCTURAL,
	NONLINEAR_DAMAGE,
	LINEAR_HEAT,
	NONLINEAR_HEAT,
	PROBLEM_MAX_PRINCIPAL_STRESS,
	NONLINEAR_STRUCTURAL_MULTISCALE,
	NONLINEAR_STRUCTURAL_SHRINKAGE};

enum enforcement_method {
	PENALTY,
	PI_MATRIX};

enum body_load_type {
	INVALID_BODY_LOAD,
	MULTISCALE_LINEAR};

enum sensitivity_parameter_type {
	MATRIX_ASSEMBLER,
	MATERIAL,
	LOAD,
	SHAPE};


// Declare some global variables
extern std::vector<std::string> elem_type_names;
extern std::vector<std::string> coh_elem_type_names;
extern std::vector<std::string> material_type_names;
extern std::vector<std::string> problem_type_names;


// Definition of id_type type.
typedef unsigned int id_type; // unsigned 32 bit
#define MPI_ID MPI_UNSIGNED
#define SPEC "u"

/*
typedef unsigned long int id_type; // unsigned 32 bit
#define MPI_ID MPI_UNSIGNED_LONG
#define SPEC "lu"
*/

// Definition of the short id type (For memory purposes when I known numbers won't be that big)
typedef unsigned char short_id_type;

// This variable is the start of indexing for local and global enricment node numbers
// NOTE: if the id_type is changed above, the definition of this variable must be changed at Mesh.cpp:14
// extern id_type ENRICH_START;
const id_type ENRICH_START = UINT_MAX/2+1;


#define PI 3.14159265358979323846264338327950

// Define some colors for printing things out nicely
// =====================================================================================
#define RESET "\033[0m"
#define BLACK "\033[30m"		/* Black */
#define RED "\033[31m"			/* Red */
#define GREEN "\033[32m"	 	/* Green */
#define YELLOW "\033[33m"		/* Yellow */
#define BLUE "\033[34m"			/* Blue */
#define MAGENTA "\033[35m"		/* Magenta */
#define CYAN "\033[36m"			/* Cyan */
#define WHITE "\033[37m"		/* White */



// Define a few macros that will help me throughout the code
// =====================================================================================

#ifndef EXEC_FUNC
#define EXEC_FUNC
inline
std::string exec(const char* cmd) {
    char buffer[128];
    std::string result = "";
    FILE* pipe = popen(cmd, "r");
    if (!pipe) throw std::runtime_error("popen() failed!");
    try {
        while (!feof(pipe)) {
            if (fgets(buffer, 128, pipe) != NULL)
                result += buffer;
        }
    } catch (...) {
        pclose(pipe);
        throw;
    }
    pclose(pipe);
    return result;
}
#endif


#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

#define print_stack()															\
	do{																			\
		void *trace[16];															\
		char **messages = (char **)NULL;											\
		int i, trace_size = 0;														\
																					\
		trace_size = backtrace(trace, 16);											\
																					\
		messages = backtrace_symbols(trace, trace_size);							\
																					\
		for (i=0; i<trace_size; ++i)												\
		{																			\
			PIGFEMPrint("[bt] Execution path: ");									\
																					\
			int p = 0;																\
			while(messages[i][p] != '(' && messages[i][p] != ' '					\
				&& messages[i][p] != 0)												\
				 ++p;																\
																					\
			char syscom[256];														\
			sprintf(syscom,"addr2line %p -e %.*s", trace[i], p, messages[i]);		\
			std::string line = exec(syscom);										\
			PIGFEMPrint(line << "\n");												\
		}																			\
	}while(0)
/*
															\
*/
	// system(syscom);
	// printf("[bt] Execution path:\n");
	// printf("[bt] #%d %s\n", i, messages[i]);																

#define err_message(message)																\
	do{																						\
		PIGFEMPrint("\n\n" << RED << "ERROR: " << message << RESET << "\n");					\
		print_stack();																		\
		exit(0);																			\
	} while(0)

// fprintf(stderr, RED "\n\nERROR: %s" RESET "\n", message);


/* 
 * Global print function. Handles output to screen based on state of
 * writeToScreen global variable set in the Options structure.
 * Also writes to a log file if given
 */
extern bool writeToScreen;
extern std::string logFile;
#define PIGFEMPrint(message)										\
	do{																\
		int rank;													\
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);						\
		if (rank==0)												\
		{															\
			if (writeToScreen)										\
				std::cout << message;								\
																	\
			std::ofstream myfile;									\
			myfile.open(logFile.c_str(), std::ofstream::app);		\
			myfile << message;										\
			myfile.close();											\
		}															\
	} while(0)


#endif //_COMMON_H_
