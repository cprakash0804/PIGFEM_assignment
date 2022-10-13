#include "SensitivityRightLoad.h"
#include "Problem.h"
#include "SensitivityParameter.h"
#include "SensitivitySolver.h"
#include "Mesh.h"
#include "DofObject.h"
#include "BoundaryObject.h"
#include <set>
#include <vector>




SensitivityRightLoad::SensitivityRightLoad()
	: _setup(false), _area(0.0)
{}
SensitivityRightLoad::~SensitivityRightLoad()
{
	if (_setup)
		VecDestroy(&_L);
}

void SensitivityRightLoad::setup()
{
	assembleLVector(_prob);
	_setup = true;
}




PetscErrorCode SensitivityRightLoad::assembleLVector(Problem* prob)
{
	// Assemble the L vector
	PetscErrorCode ierr;
	if (_setup)
		VecDestroy(&_L);
	
	DofObject* dofs = _prob->get_dofs();
	Mesh* mesh = _prob->get_mesh();
	BoundaryObject* boundary = _prob->get_boundary();
	id_type N = dofs->n_global_const_dofs();
	id_type Nf = dofs->n_global_free_dofs();
	if (_prob->get_mesh()->serial())
	{
		ierr = VecCreateSeq(mesh->get_comm(), N, &_L);CHKERRQ(ierr);
	}
	else
	{
		id_type n = dofs->n_local_owned_const_dofs();
		ierr = VecCreateMPI(mesh->get_comm(), n, N, &_L);CHKERRQ(ierr);
	}
	ierr = VecZeroEntries(_L);CHKERRQ(ierr);

	// Loop over all of the dirichlet bcs and if a node is on the right ahnd side and has an x-dirichlet bc add it to the L vector
	std::set<id_type> right = mesh->get_nodeset("right");
	std::vector<PetscInt> rows;
	std::vector<PetscScalar> vals;
	for(BoundaryObject::dirichlet_iterator it=boundary->dirichlet_begin(), end=boundary->dirichlet_end(); it!=end; ++it)
	{
		id_type g_node = (*it).first;
		if ( !mesh->own_node_global(g_node) )	// This mesh partition doesn't own this node
			continue;
		if (right.find(g_node) == right.end())	// This node isn't on the right side
			continue;

		bool has_x_dof = false;
		for(std::map<id_type, double>::iterator it2=it->second.begin(), end2=it->second.end(); it2!=end2; ++it2)
		{
			id_type dof = (*it2).first;
			if (dof == 0)
			{
				has_x_dof = true;
				break;
			}
		}
		if (has_x_dof)
		{
			rows.push_back( dofs->get_global_dof(g_node, 0) - Nf );
			vals.push_back( 1.0 );
		}
	}
	ierr = VecSetValues(_L, rows.size(), rows.data(), vals.data(), INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(_L);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(_L);CHKERRQ(ierr);

	// Compute the area of the right hand side
	std::vector<double> mins, maxes;
	mesh->getDomainLimits(mins, maxes);
	if (mesh->dim() == 1)
		_area = 1.0;
	else if (mesh->dim() == 2)
		_area = maxes[1]-mins[1];
	else if (mesh->dim() == 3)
		_area = (maxes[1]-mins[1]) * (maxes[2]-mins[2]);

	return ierr;
}








double SensitivityRightLoad::evaluate(pmatrix* K, pvector* U, pvector* F_ext)
{
	PetscScalar ret;
	VecDot(_L, F_ext->p, &ret);
	ret /= _area;
	return ret;
}

void SensitivityRightLoad::assemble_dfdUf(Vec* dfdUf,
										  pmatrix* K, pvector* U, pvector* F_ext)
{
	MatMultTranspose((K->pf), _L, *dfdUf);
	VecScale(*dfdUf, 1.0/_area);
}

void SensitivityRightLoad::assemble_dfdUp(Vec* dfdUp,
										  pmatrix* K, pvector* U, pvector* F_ext)
{
	MatMultTranspose((K->pp), _L, *dfdUp);
	VecScale(*dfdUp, 1.0/_area);
}

double SensitivityRightLoad::get_parameter_partial(SensitivityParameter* param, SensitivitySolver* solver,
												   pmatrix* K, pvector* U, pvector* F_ext)
{
	PetscScalar ret;
	VecDot(_L, *solver->get_dPp(param->get_id()), &ret);
	ret /= _area;
	return ret;
}