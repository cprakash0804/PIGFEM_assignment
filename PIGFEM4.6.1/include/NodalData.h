/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _NODAL_DATA_H_
#define _NODAL_DATA_H_

#include "common.h"
#include "DenseMatrix.h"
#include <vector>

// Predeclarations
class Mesh;





/*
 * This class encapsulates the storage of the data associated
 * with the nodes from a problem on the local mesh partition.
 * This data could be the displacement, temperature, external force, etc
 */
 class NodalData
 {
 	private:

 		/*
 		 * These structures store the solution data corresponding
 		 * to all of the nodes on the local mesh partition.
 		 * _solution stores the data on the normal nodes
 		 * _enrich_solution stores the enrichment solution on the
 		 * enrichment nodes
 		 */
 		DenseMatrix<double> _data;
		DenseMatrix<double> _enrich_data;

		id_type _nndof;
		bool _preallocated;


		/*
		 * This is the problem for which this solution is based on
		 */
		Mesh* _mesh;

	public:

		/*
		 * The Big 3.
		 * Constructor requires a Problem pointer to be built
		 * Destructor and Copy Constructor don't do anything
		 */
		NodalData(Mesh* mesh);
		~NodalData();
		NodalData(const NodalData& other);
		id_type nndof() const {return _nndof;};


		/*
		 * Function preallocates storage based on the current state
		 * of the Problem. May be called after mesh refinement to
		 * add storage
		 */
		void preallocate_storage(id_type nndof);
		bool preallocated() const {return _preallocated;};


		/*
		 * Return writeable reference to the nodal data
		 */
		// const double& get_value_local(id_type idx, id_type dof) const;
		// const double& get_value_global(id_type id, id_type dof) const;
		// const double& get_value(id_type id, id_type dof) const {return get_value_global(id, dof);};
		// double& get_value_local(id_type idx, id_type dof)
		// {
		// 	return const_cast<double&>(static_cast<const NodalData*>(this)->get_value_local(idx, dof));
		// }
		// double& get_value_global(id_type id, id_type dof)
		// {
		// 	return const_cast<double&>(static_cast<const NodalData*>(this)->get_value_global(id, dof));
		// }
		// double& get_value(id_type id, id_type dof)
		// {
		// 	return const_cast<double&>(static_cast<const NodalData*>(this)->get_value(id, dof));
		// }
		double& get_value_local(id_type idx, id_type dof);
		double& get_value_global(id_type id, id_type dof);
		double& get_value(id_type id, id_type dof) {return get_value_global(id, dof);};
		double get_value_local(id_type idx, id_type dof) const;
		double get_value_global(id_type id, id_type dof) const;
		double get_value(id_type id, id_type dof) const {return get_value_global(id, dof);};

		/*
		 * Returns the pointer to the problem object
		 */
		Mesh* get_mesh() const {return _mesh;};


 };






#endif