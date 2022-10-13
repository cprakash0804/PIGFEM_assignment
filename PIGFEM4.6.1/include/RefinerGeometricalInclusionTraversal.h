/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _REFINER_INCLUSION_TRAVERSE_H_
#define _REFINER_INCLUSION_TRAVERSE_H_

#include "RefinerGeometricalInclusion.h"


class RefinerGeometricalInclusionTraversal : public RefinerGeometricalInclusion
{
	public:


		/*
		 * The Big 3
		 * The constructor takes the problem that will be partitioned
		 * Destructor doesn't do anything
		 * Copy constructor copies over the problem pointer
		 */
		RefinerGeometricalInclusionTraversal(Problem* prob);
		~RefinerGeometricalInclusionTraversal();
		RefinerGeometricalInclusionTraversal(const RefinerGeometricalInclusionTraversal& other);


		/*
		 * The main interface with this class. Will be called to refine the mesh.
		 * Can be overridden for specific implementation in derived classes
		 */
		virtual void refine();


	protected:


		/*
		 * Traverse the exterior of inclusions from an intial refinement point
		 * This provides a better representation of the actual inclusion object
		 */
		void traverse_inclusion_exterior(std::set<id_type>& refined_elements,
										 const std::vector<unsigned char>& elem_code);


		/*
		 * Adds to the set refined_elements all those elements that are neighbors of the new_elements.
		 * Also adds lists of new remote elements that should be refined in parallel
		 */
		void find_local_traversal_additions(std::set<id_type>& refined_elements, std::set<id_type>& new_elements,
											std::map<int, std::vector<id_type> >& remote_additions,
											const std::vector<unsigned char>& elem_code);
};




#endif