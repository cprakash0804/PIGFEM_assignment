/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#include "Writer.h"
#include "Problem.h"

/*
 * The Big 3
 * Constructor need the Problem
 * Destructor doesn't do anything
 * Copy Constructor copies the problem pointer
 */
Writer::Writer(Problem* prob)
	: _prob(prob)
{}
Writer::Writer()
	: _prob(NULL)
{}
Writer::~Writer()
{}
Writer::Writer(const Writer& other)
{
	_prob = other.get_prob();
}
