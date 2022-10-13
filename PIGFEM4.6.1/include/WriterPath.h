/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated April 2017

##################################################################################
*/
#ifndef _Writer_PATH_H_
#define _Writer_PATH_H_

#include "Writer.h"
#include "Utilities.h"
#include <map>


class Material;

/*
 * This class is used to write the Problem sensitivity data to
 *  a single file over the course of the simulation
 */
class WriterPath : public Writer
{
	public:


		/*
		 * The Big 3
		 * Constructor need the Problem
		 * Destructor doesn't do anything
		 * Copy Constructor copies the problem pointer and stores the mesh
		 */
		WriterPath(Problem* prob);
		WriterPath();
		virtual ~WriterPath();
		WriterPath(const WriterPath& other);


		/*
		 * Main function used to write problem to any inherited file type
		 */
		virtual void write(std::string filename, double curr_t);


		/*
		 * Stores the mesh in this Writer's internal storage
		 * Can be recalled if the mesh is adapted
		 */
		virtual void store_mesh();


		virtual void setParameter(std::string param, bool val);
		virtual bool getBoolParameter(std::string param);
		virtual void setParameter(std::string param, std::vector<double> val);
		virtual std::vector<double> getVecParameter(std::string param);


	protected:


		/*
         * The actual function that will do all of the writing (maybe by calling other functions)
         */
        virtual void writeFromStream(std::ofstream& myfile, double curr_t);


        // The list of all points that I am finding a solutions for
        std::vector<std::vector<double> > _points;
        std::vector<int> _on_partition;

        // The list of the nodal shape functions and their gradients at each point
        std::vector<id_type> _elements; // global element numbers
        std::vector<std::vector<double> > _N;
        std::vector<std::vector<std::vector<double> > > _dN;
        std::vector<Material*> _mats;

        // Whether or not to output the stresses with the solution
        bool _output_stress;

};





#endif