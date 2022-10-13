/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated April 2017

##################################################################################
*/
#ifndef _Writer_COH_OPENINGS_H_
#define _Writer_COH_OPENINGS_H_

#include "Writer.h"
#include "Utilities.h"
#include <map>




/*
 * This class is used to write the Problem sensitivity data to
 *  a single file over the course of the simulation
 */
class WriterCohesiveOpenings : public Writer
{
	public:


		/*
		 * The Big 3
		 * Constructor need the Problem
		 * Destructor doesn't do anything
		 * Copy Constructor copies the problem pointer and stores the mesh
		 */
		WriterCohesiveOpenings(Problem* prob);
		WriterCohesiveOpenings();
		virtual ~WriterCohesiveOpenings();
		WriterCohesiveOpenings(const WriterCohesiveOpenings& other);


		/*
		 * Main function used to write problem to any inherited file type
		 */
		virtual void write(std::string filename, double curr_t);


		/*
		 * Stores the mesh in this Writer's internal storage
		 * Can be recalled if the mesh is adapted
		 */
		virtual void store_mesh() {};


	protected:


		/*
         * The actual function that will do all of the writing (maybe by calling other functions)
         */
        virtual void writeFromStream(std::ofstream& myfile, double curr_t);

};





#endif