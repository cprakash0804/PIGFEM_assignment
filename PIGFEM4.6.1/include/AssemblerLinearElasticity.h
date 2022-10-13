#ifndef _LE_ASSEMBLER_H_
#define _LE_ASSEMBLER_H_
#include "Assembler.h"
#include "common.h"


class AssemblerLinearElasticity : public Assembler
{
	public:
		virtual ~AssemblerLinearElasticity() {};

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
