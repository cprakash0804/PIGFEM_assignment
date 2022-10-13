/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated April 2017

##################################################################################
*/
#ifndef _Writer_Coh_Fail_Lengths_H_
#define _Writer_Coh_Fail_Lengths_H_

#include "Writer.h"
#include "Utilities.h"
#include <map>




/*
 * This class is used to write how much of the cohesive surface
 *  in the problem have passed the point where it woud fail
 */
class WriterCohesiveFailureLengths : public Writer
{
	public:


		/*
		 * The Big 3
		 * Constructor need the Problem
		 * Destructor doesn't do anything
		 * Copy Constructor copies the problem pointer and stores the mesh
		 */
		WriterCohesiveFailureLengths(Problem* prob);
		WriterCohesiveFailureLengths();
		virtual ~WriterCohesiveFailureLengths();
		WriterCohesiveFailureLengths(const WriterCohesiveFailureLengths& other);

		/*
		 * Main function used to write problem to any inherited file type
		 */
		virtual void writeConsecutive(std::string filename, double curr_t, id_type step);


		/*
		 * Stores the mesh in this Writer's internal storage
		 * Can be recalled if the mesh is adapted
		 */
		virtual void store_mesh();


		virtual void setParameter(std::string param, bool val);
		virtual bool getBoolParameter(std::string param);

	protected:


		/*
         * The actual function that will do all of the writing (maybe by calling other functions)
         */
        virtual void writeFromStream(std::ofstream& myfile, double curr_t);
        void writeLengths(std::ofstream& myfile);


        // Parameter that governs whether or not to output ALL of the cohesive data and not just a summary of it
        bool _individual;


        // vector of the local lengths of each cohesive zone 
		std::vector<double> _local_lengths;
};





#endif