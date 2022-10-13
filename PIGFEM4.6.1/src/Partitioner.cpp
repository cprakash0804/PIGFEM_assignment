/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "Partitioner.h"


/*
 * The Big 3
 * The constructor takes the mesh that will be partitioned
 * Destructor doesn't do anything
 * Copy constructor copies over the mesh pointer
 */
Partitioner::Partitioner(Mesh* mesh)
	: _mesh(mesh)
{}
Partitioner::Partitioner()
{}
Partitioner::~Partitioner()
{}
Partitioner::Partitioner(const Partitioner& other)
{
	_mesh = other.get_mesh();
}