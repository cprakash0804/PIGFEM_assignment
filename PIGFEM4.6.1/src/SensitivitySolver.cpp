#include "SensitivitySolver.h"
#include "Problem.h"
#include "Solver.h"
#include "SensitivityParameter.h"
#include "SensitivityFunction.h"
#include "DofObject.h"
#include "Mesh.h"
#include "InternalVars.h"


SensitivitySolver::SensitivitySolver()
	: _internal_vars(NULL), _setup(false)
{}

SensitivitySolver::~SensitivitySolver()
{
	// Destroy all of the PETSc vectors
	for (id_type i=0; i<_dUp_dd.size(); ++i)
		VecDestroy(&_dUp_dd[i]);
	for (id_type i=0; i<_dFf_dd.size(); ++i)
		VecDestroy(&_dFf_dd[i]);

	if (_internal_vars != NULL)
		delete _internal_vars;
}



// Initilize everything and make sure all structures are allocated
PetscErrorCode SensitivitySolver::init()
{
	_function_vals.resize(n_functions());
	_sensitivities.resize(n_functions());
	for (id_type func=0; func<n_functions(); ++func)
	{
		_functions[func]->attachProblem(_prob);
		_functions[func]->init();
		_sensitivities[func].resize(n_parameters());
	}

	for (id_type param=0; param<n_parameters(); ++param)
	{
		_parameters[param]->attachProblem(_prob);
		_parameters[param]->init();
	}

	_sensitivity_solve_time = 0.0;
	_sensitivity_substitution_time = 0.0;
	_sensitivity_assemble_time = 0.0;

	if (_internal_vars != NULL)
		delete _internal_vars;
	_internal_vars = new InternalVars(_prob, true);

	return 0;
}

PetscErrorCode SensitivitySolver::setup()
{
	PetscErrorCode ierr;
	if (_setup)
	{
		for (id_type func=0; func<_delf_delU.size(); ++func)
			_delf_delU[func].destroy();
		for (id_type param=0; param<_delP_deld.size(); ++param)
			_delP_deld[param].destroy();
	}

	_delf_delU.resize(n_functions());
	for (id_type func=0; func<n_functions(); ++func)
	{
		ierr = preallocateVector(&_delf_delU[func].f, true);CHKERRQ(ierr);
		ierr = preallocateVector(&_delf_delU[func].p, false);CHKERRQ(ierr);
		_functions[func]->setup();
	}
	_delP_deld.resize(n_parameters());
	for (id_type param=0; param<n_parameters(); ++param)
	{
		ierr = preallocateVector(&_delP_deld[param].f, true);CHKERRQ(ierr);
		ierr = preallocateVector(&_delP_deld[param].p, false);CHKERRQ(ierr);
		_parameters[param]->setup();
	}

	// Get the pointers to the PETSc solver objetcs so we can use them here
	Solver* solver = _prob->get_solver();
	_K = &solver->_K;
	_U = &solver->_U;
	_F_ext = &solver->_F_ext;

	_setup = true;
	return ierr;
}



// Function to preallocate the appropriate vectors
PetscErrorCode SensitivitySolver::preallocateVector(Vec* vec, bool free_portion)
{
	PetscErrorCode ierr;
	PetscInt n, N, n_pre, N_pre;
	DofObject* dofs = _prob->get_dofs();
	Mesh* mesh = _prob->get_mesh();
	n = dofs->n_local_owned_free_dofs();
	N = dofs->n_global_free_dofs();
	n_pre = dofs->n_local_owned_const_dofs();
	N_pre = dofs->n_global_const_dofs();

	if(mesh->serial())
	{
		if (free_portion)
		{
			ierr = VecCreateSeq(mesh->get_comm(), N, vec);CHKERRQ(ierr);
		}
		else
		{
			ierr = VecCreateSeq(mesh->get_comm(), N_pre, vec);CHKERRQ(ierr);
		}
	}
	else
	{
		if (free_portion)
		{
			ierr = VecCreateMPI(mesh->get_comm(), n, N, vec);CHKERRQ(ierr);
		}
		else
		{
			ierr = VecCreateMPI(mesh->get_comm(), n_pre, N_pre, vec);CHKERRQ(ierr);
		}
	}

	return ierr;
}











// Add a function to find the sensitivity of
void SensitivitySolver::addFunction(SensitivityFunction* func)
{
	_functions.push_back(func);
}



/*
 * Function to add parameters to the problem
 */
void SensitivitySolver::addParameter(SensitivityParameter* param)
{
	if (param->get_type() == MATERIAL ||
		param->get_type() == SHAPE)
	{
		_parameters.push_back(param);
	}

	else if (param->get_type() == LOAD)
	{
		err_message("No framework for sensitivity wrt load implemented");
	}

	else
		err_message("Unknown sensitivty parameter type!");
}




double SensitivitySolver::get_function_val(id_type func)
{
	if (func < n_functions())
		return _function_vals[func];
	else
		err_message("Invalid function.");
}
double SensitivitySolver::get_sensitivity(id_type func, id_type param)
{
	if (func < n_functions())
	{
		if (param < _sensitivities[func].size())
			return _sensitivities[func][param];
		else
			err_message("Invalid parameter.");
	}
	else
		err_message("Invalid function.");
}







// Function to call all of the appropriate assemblers to assemble the current derivatives
PetscErrorCode SensitivitySolver::assembleDerivatives()
{
	PetscErrorCode ierr;
	bool have_load_params = false;
	for (id_type i=0; i<n_parameters(); ++i)
		if (_parameters[i]->get_type() == LOAD)
			have_load_params = true;

	// Assemble the objective derivatives
	for (id_type i=0; i<n_functions(); ++i)
	{
		ierr = VecZeroEntries(_delf_delU[i].f);CHKERRQ(ierr);
		_functions[i]->assemble_dfdUf(&_delf_delU[i].f,
									  _K, _U, _F_ext);
		if (have_load_params)
		{
			ierr = VecZeroEntries(_delf_delU[i].p);CHKERRQ(ierr);
			_functions[i]->assemble_dfdUp(&_delf_delU[i].p,
										  _K, _U, _F_ext);
		}
	}

	// Assemble the partials wrt the parameters
	for (id_type i=0; i<n_parameters(); ++i)
	{
		ierr = VecZeroEntries(_delP_deld[i].f);CHKERRQ(ierr);
		ierr = VecZeroEntries(_delP_deld[i].p);CHKERRQ(ierr);
		_parameters[i]->assemble_dPdd(_delP_deld[i], _prob);
		if (_parameters[i]->get_type() == LOAD)
		{
			ierr = VecZeroEntries(_dUp_dd[i]);CHKERRQ(ierr);
			_parameters[i]->assemble_load_derivatives(&_dUp_dd[i], &_dFf_dd[i], _prob);
		}
	}

	// The hard part. get the partials of the objective functions wrt all of the parameters
	_delf_deld.resize(n_functions());
	for (id_type i=0; i<n_functions(); ++i)
	{
		_delf_deld[i].resize(n_parameters());
		for (id_type j=0; j<n_parameters(); ++j)
		{
			_delf_deld[i][j] = _functions[i]->get_parameter_partial(_parameters[j], this,
																	_K, _U, _F_ext);
		}
	}

	return ierr;
}