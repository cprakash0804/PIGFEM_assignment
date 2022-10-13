/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _REFINER_INCLUSION_H_
#define _REFINER_INCLUSION_H_

#include "Refiner.h"


class RefinerGeometricalInclusion : public Refiner
{
	public:


		/*
		 * The Big 3
		 * The constructor takes the problem that will be partitioned
		 * Destructor doesn't do anything
		 * Copy constructor copies over the problem pointer
		 */
		RefinerGeometricalInclusion(Problem* prob);
		RefinerGeometricalInclusion();
		~RefinerGeometricalInclusion();
		RefinerGeometricalInclusion(const RefinerGeometricalInclusion& other);


		/*
		 * The main interface with this class. Will be called to refine the mesh.
		 * Can be overridden for specific implementation in derived classes
		 */
		virtual void refine();


	protected:


		/*
		 * This is the function that must be implemented to generate a set of elements to refine
		 * (Stubbed out here)
		 */
		virtual void generate_refine_set(std::set<id_type>& elements);


		/*
		 * Count the number of unique inclusions that intersect each local element
		 */
		void generate_element_code(std::vector<unsigned char>& elem_code);


		/*
		 * Takes the element code and adds all the elements that satisfy the criterion of
		 * being intersected by more then one inclusion
		 */
		void add_multiple_intersected_elements(std::set<id_type>& new_elements, const std::vector<unsigned char>& elem_code);
};




#endif