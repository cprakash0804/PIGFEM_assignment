/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _LT_ASSEMBLER_H_
#define _LT_ASSEMBLER_H_
#include "Assembler.h"
#include "common.h"


class AssemblerLinearThermal : public Assembler
{
	public:
		virtual ~AssemblerLinearThermal() {};

		virtual bool storedBmats();

	protected:

		// Kernel function that actually does the computes the physics behind the specific problem
		// Computes the local contribution to the stiffness matrix and internal load vector
		// Also updates the internal variable storage in the input variable
		virtual void KernelLinear(DenseMatrix<double>& K_el,
								  const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
								  Material* mat, Material::input_params& input);
};


#endif
