/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#ifndef _WRITER_H_
#define _WRITER_H_

#include <string>
#include <vector>
#include "common.h"

// Predeclarations
class Problem;
class CohesiveElem;
class Material;


/*
 * This class is a base class used to output the
 * Problem to a file
 */
class Writer
{
	protected:


		/*
		 * The Problem that this writer outputs
		 */
		Problem* _prob;


        /*
         * The actual function that will do all of the writing (maybe by calling other functions)
         */
        virtual void writeFromStream(std::ofstream& myfile, double curr_t) = 0;


	public:


		/*
		 * The Big 3
		 * Constructor need the Problem
		 * Destructor doesn't do anything
		 * Copy Constructor copies the problem pointer
		 */
		Writer(Problem* prob);
		Writer();
		virtual ~Writer();
		Writer(const Writer& other);


		/*
		 * Main function used to write problem to any inherited file type
		 */
		virtual void write(std::string filename, double curr_t) {};
		virtual void writeConsecutive(std::string filename, double curr_t, id_type step) {};


		/*
		 * Stores the mesh in this Writer's internal storage
		 * Can be recalled if the mesh is adapted
		 */
		virtual void store_mesh() {};


		/*
		 * Returns the Problem pointer
		 */
		Problem* get_prob() const {return _prob;}


		virtual void setParameter(std::string param, bool val) {};
		virtual bool getBoolParameter(std::string param) {return false;};
		virtual void setParameter(std::string param, std::vector<double> val) {};
		virtual std::vector<double> getVecParameter(std::string param) {return std::vector<double>();};
};



#endif
