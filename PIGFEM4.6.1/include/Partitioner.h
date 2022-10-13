/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _PARTITIONER_H_
#define _PARTITIONER_H_

#include "common.h"

// Predeclarations
class Mesh;


/*
 * Base class that is used to derive mesh partitioners
 */
class Partitioner
{
	protected:


		/*
		 * The mesh object that this object partitions
		 */
		Mesh* _mesh;


	public:


		/*
		 * The Big 3
		 * The constructor takes the mesh that will be partitioned
		 * Destructor doesn't do anything
		 * Copy constructor copies over the mesh pointer
		 */
		Partitioner(Mesh* mesh);
		Partitioner();
		~Partitioner();
		Partitioner(const Partitioner& other);


		/*
		 * Get the mesh pointer associated with this partitioner
		 */
		Mesh* get_mesh() const {return _mesh;};


		/*
		 * The main function that is called to partition the mesh
		 */
		virtual void partition() = 0;
};



#endif