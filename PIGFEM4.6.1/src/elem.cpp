/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#include "elem.h"
#include "point1.h"
#include "edge2.h"
#include "tri3.h"
#include "quad4.h"
#include "hex8.h"
#include "tet4.h"
#include "prism6.h"
#include "Polygon.h"
#include "Utilities.h"
#include <algorithm>


double Elem::_inside_tol = 1e-10;


Elem::Elem()
	: _global_id(0), _is_intersected(false), _detected(false)
{
}
Elem::Elem(std::vector<Node*> nodes, id_type id)
{
	_nodes = nodes;
	_global_id = id;
	_detected = false;
	_is_intersected = false;
}
Elem::~Elem()
{
	// Since the integration elements are only stored heere they need to be deleted here
	for(id_type ie=0; ie<n_integration_elem(); ++ie)
		delete _integrationElem[ie];
	for(id_type ce=0; ce<n_cohesive_elem(); ++ce)
		delete _cohesiveElem[ce];
}

void Elem::clear()
{
	_nodes.clear();
	_global_id = 0;
	//_elemDetector = 0;
	_detected = false;
	_enrichment_on_inclusion.clear();
	_cutEdge.clear();
	_int_elem_struct.clear();
	_int_elem_mat.clear();
	_coh_mat.clear(); // Pointers will be stored and deleted in th mesh
	for(id_type ie=0; ie<n_integration_elem(); ++ie)
		delete _integrationElem[ie];
	_integrationElem.clear();
	_enrichment_nodes.clear();
	_is_intersected = false;
}
void Elem::clearEnrichments()
{
	_is_intersected = false;
	_detected = false;
	_cutEdge.clear();
	_enrichment_on_inclusion.clear();
	_int_elem_struct.clear();
	_coh_elem_struct.clear();
	_int_elem_mat.clear(); // Don't need to delete because they exist in the mesh
	_coh_mat.clear();
	for (id_type i=0; i<_integrationElem.size(); ++i)
		delete _integrationElem[i];
	for (id_type i=0; i<_cohesiveElem.size(); ++i)
		delete _cohesiveElem[i];
	_integrationElem.clear();
	_cohesiveElem.clear();
	_enrichment_nodes.clear(); // These are also handled by the mesh
}

void Elem::copy(Elem & other_elem)
{
	other_elem.set_id(_global_id);
	other_elem.set_nodes(_nodes);
	//other_elem.getElemDetector() = _elemDetector;
	other_elem.get_detected() = _detected;
	other_elem.getIntersectedEdges() = _cutEdge;
	other_elem.getInclusionNumbers() =_enrichment_on_inclusion;
	other_elem.getIntegrationElemStructure() = _int_elem_struct;
	other_elem.getIntegrationElemMaterials() = _int_elem_mat;
	other_elem.getIntegrationElem() = _integrationElem;
	other_elem.getEnrichmentNodes() = _enrichment_nodes;
	other_elem.get_intersected() = _is_intersected;
}

Elem::Elem(Elem & other_elem)
{
	copy(other_elem);
}

Node& Elem::operator()(id_type i)
{
	if(i < n_nodes())
		return *_nodes[i];
	else if(i < (n_nodes()+n_enrich_nodes()))
		return *_enrichment_nodes[i-n_nodes()];
	else
		err_message("Node index must be less than the number of nodes.");
}

void Elem::ShapeFunctions(std::vector<double> rcoords, std::vector<double>& N, std::vector<std::vector<double> >& dNdx, double& J)
{
	// Fill in the shape functions and local gradients
	N = compute_shape(rcoords);

	std::vector<std::vector<double> > dNdr = compute_shape_grad(rcoords);

	// Transform the local gradients to the global gradients (and get the jacobian)
	transform_gradients(dNdr, dNdx, J);
}
void Elem::transform_gradients(std::vector<std::vector<double> >& dNdr, std::vector<std::vector<double> >& dNdx, double& J)
{
	DenseMatrix<double> Jac(dim(), dim());	// The child mapping jacobian
	for (id_type alpha=0; alpha<dim(); ++alpha)
		for (id_type i=0; i<dim(); ++i)
			for (id_type n=0; n<n_nodes(); ++n)
				Jac(alpha,i) += dNdr[n][alpha] * ((*_nodes[n])(i));

	J = Utilities::det(Jac);
	DenseMatrix<double> Jinv = Utilities::inverse(Jac);

	dNdx.resize(n_nodes());
	for (id_type n=0; n<n_nodes(); ++n)
		dNdx[n] = Jinv * dNdr[n];
}




void Elem::ShapeFunctionPartials(std::vector<double> rcoords, int child, const std::vector<double>& nodal_velocity,
								 std::vector<double>& dN_dd, std::vector<std::vector<double> >& d2N_dxdd,
								 std::vector<double>& v, double& div_v)
{
	if (!_is_intersected)
		err_message("Shape Function Partials shouldn't be called on a non-intersected element!");

	// Do the child element computations
	Elem* int_el = _integrationElem[child];
	id_type n_child = int_el->n_nodes();
	std::vector<double> dN_dd_child;
	std::vector<std::vector<double> > d2N_dxdd_child;

	// Build the local structure of the integration elem
	std::vector<id_type> child_struct(n_child);
	for (id_type n=0; n<n_child; ++n)
	{
		id_type n_id = int_el->get_node(n)->get_id();
		for (id_type n1=0; n1<n_nodes(); ++n1)
			if (_nodes[n1]->get_id() == n_id)
				child_struct[n] = n1;
		for (id_type n1=0; n1<n_enrich_nodes(); ++n1)
			if (_enrichment_nodes[n1]->get_id() == n_id)
				child_struct[n] = n1 + n_nodes();
	}

	// Assemble the child element velocity nodal vector
	std::vector<double> child_nodal_velocity(dim() * n_child);
	for (id_type n=0; n<n_child; ++n)
		for (id_type d=0; d<dim(); ++d)
			child_nodal_velocity[n*dim() + d] = nodal_velocity[child_struct[n]*dim() + d];

	std::vector<double> x; // The global coordinate values
	int_el->ChildShapeFunctionPartials(rcoords, child_nodal_velocity, dN_dd_child, d2N_dxdd_child, x, v, div_v);

	// Now need to find the parent r coordinates
	bool inside;
	std::vector<double> prcoords = inverse_map(x, inside);
	if (!inside)
		err_message("Invalid inverse map in shape function partials computation");

	// Compute the parent jacobian and jacobian partial
	std::vector<std::vector<double> > dN_dr = compute_shape_grad(prcoords);
	std::vector<DenseMatrix<double> > d2N_drdr = compute_shape_grad_grad(prcoords);

	DenseMatrix<double> Jac(dim(), dim());	// The parent mapping jacobian
	DenseMatrix<double> dJac(dim(), dim());	// The parent mapping jacobian partial
	for (id_type alpha=0; alpha<dim(); ++alpha)
		for (id_type i=0; i<dim(); ++i)
			for (id_type n=0; n<n_nodes(); ++n)
			{
				Jac(alpha,i) += dN_dr[n][alpha] * ((*_nodes[n])(i));
				dJac(alpha,i) += dN_dr[n][alpha] * nodal_velocity[n*dim() + i]; // This should always be 0 actually
			}
	for (id_type n=0; n<n_nodes(); ++n)
	{
		std::vector<double> ret = d2N_drdr[n] * (Jac * v);
		for (id_type alpha=0; alpha<dim(); ++alpha)
			for (id_type i=0; i<dim(); ++i)
				dJac(alpha,i) = ret[alpha] * ((*_nodes[n])(i));
	}

	DenseMatrix<double> Jinv = Utilities::inverse(Jac);

	// Compute the derivatives of the parent shape functions and shape function gradients
	dN_dd.resize(n_nodes() + n_enrich_nodes());
	d2N_dxdd.resize(n_nodes() + n_enrich_nodes());
	std::fill(dN_dd.begin(), dN_dd.end(), 0.0);
	std::fill(d2N_dxdd.begin(), d2N_dxdd.end(), std::vector<double>(dim(), 0.0));
	for (id_type n=0; n<n_nodes(); ++n)
	{
		std::vector<double> dNdx = Jinv * dN_dr[n];
		dN_dd[n] = Utilities::dot(dNdx, v);

		std::vector<double> temp1 = dJac * (Jinv * dN_dr[n]);
		std::vector<double> temp2 = d2N_drdr[n] * (Jinv * v);
		std::vector<double> temp3 = Utilities::minus(temp2, temp1);
		d2N_dxdd[n] = Jinv * temp3;
		// std::vector<double> temp2 = d2N_drdr[n] * (Jinv * v);
		// d2N_dxdd[n] = Jinv * temp2;
	}

	// Add the contributions from the enrichment nodes
	for (id_type n=0; n<child_struct.size(); ++n)
	{
		if (child_struct[n] >= n_nodes())
		{
			dN_dd[child_struct[n]] = dN_dd_child[n];
			d2N_dxdd[child_struct[n]] = d2N_dxdd_child[n];
		}
	}

	// // Add in the second part of the contribution to the divergence term
	// for (id_type n=0; n<n_nodes(); ++n)
	// 	for (id_type d=0; d<dim(); ++d)
	// 		div_v += d2N_dxdd[n][d] * ((*_nodes[n])(d));
	// for (id_type n=0; n<n_enrich_nodes(); ++n)
	// 	for (id_type d=0; d<dim(); ++d)
	// 		div_v += d2N_dxdd[n+n_nodes()][d] *((*_enrichment_nodes[n])(d));
}



void Elem::ChildShapeFunctionPartials(const std::vector<double>& rcoords, const std::vector<double>& child_nodal_velocity,
									  std::vector<double>& dN_dd_child, std::vector<std::vector<double> >& d2N_dxdd_child,
									  std::vector<double>& x, std::vector<double>& v, double& div_v)
{
	std::vector<double> N = compute_shape(rcoords);
	std::vector<std::vector<double> > dN_dr = compute_shape_grad(rcoords);
	dN_dd_child.resize(n_nodes());
	for (id_type n=0; n<n_nodes(); ++n)
		dN_dd_child[n] = 0.0;

	// Compute the velocity at the current integration point
	v.resize(dim());
	x.resize(dim());
	for (id_type d=0; d<dim(); ++d)
	{
		v[d] = 0.0;
		for (id_type n=0; n<n_nodes(); ++n)
		{
			x[d] += N[n] * ((*_nodes[n])(d));
			v[d] += N[n] * child_nodal_velocity[n*dim() + d];
		}
	}

	DenseMatrix<double> Jac(dim(), dim());	// The child mapping jacobian
	DenseMatrix<double> dJac(dim(), dim());	// The child mapping jacobian partial
	for (id_type alpha=0; alpha<dim(); ++alpha)
		for (id_type i=0; i<dim(); ++i)
			for (id_type n=0; n<n_nodes(); ++n)
			{
				Jac(alpha,i) += dN_dr[n][alpha] * ((*_nodes[n])(i));
				dJac(alpha,i) += dN_dr[n][alpha] * child_nodal_velocity[n*dim() + i];
			}

	// Copmute the inverse of the jacobian
	DenseMatrix<double> Jinv = Utilities::inverse(Jac);

	// Compute the partial of the shape function gradients and the divergence of the velocity
	d2N_dxdd_child.resize(n_nodes());
	div_v = 0.0;
	for (id_type n=0; n<n_nodes(); ++n)
	{
		d2N_dxdd_child[n] = (-1.0 * Jinv) * (dJac * (Jinv * dN_dr[n]) );
		std::vector<double> dNdx = Jinv * dN_dr[n];
		for (id_type d=0; d<dim(); ++d)
			div_v += dNdx[d] * child_nodal_velocity[n*dim() + d];
	}
}




void Elem::print(std::ostream & os) const
{
	os << "Element id: " << _global_id << " (";
	switch(get_type())
	{
		case INVALID_ELEM:
			os << "Invalid Element";
			break;
		case POINT1:
			os << "Point1";
			break;
		case EDGE2:
			os << "Edge2";
			break;
		case EDGE3:
			os << "Edge3";
			break;
		case TRI3:
			os << "Tri3";
			break;
		case TRI6:
			os << "Tri6";
			break;
		case QUAD4:
			os << "Quad4";
			break;
		case QUAD8:
			os << "Quad8";
			break;
		case TET4:
			os << "Tet4";
			break;
		case TET10:
			os << "Tet10";
			break;
		case HEX8:
			os << "Hex8";
			break;
		case HEX20:
			os << "Hex20";
			break;
		case HEX27:
			os << "Hex27";
			break;
		case PRISM6:
			os << "Prism6";
			break;
		case PRISM15:
			os << "Prism15";
			break;
		case PRISM18:
			os << "Prism18";
			break;
		default:
			os << "Non-supported type";
			break;
	}
	os << ")" << std::endl;
	for(id_type i=0; i<n_nodes(); ++i)
		os << "\t" << *_nodes[i];
	for(id_type i=0; i<n_enrich_nodes(); ++i)
		os << "\tEnrichment: " << *_enrichment_nodes[i];
}


Elem* Elem::build(elem_type etype)
{
	switch(etype)
	{
		case POINT1:
			return new Point1;
		case EDGE2:
			return new Edge2;
		case TRI3:
			return new Tri3;
		case QUAD4:
			return new Quad4;
		case TET4:
			return new Tet4;
		case HEX8:
			return new Hex8;
		case PRISM6:
			return new Prism6;
		case EDGE3:
		case TRI6:
		case QUAD8:
		case TET10:
		case HEX20:
		case HEX27:
		case PRISM15:
		case PRISM18:
		default:
			char buf[100];
			sprintf(buf, "Attempting to build a non-supported elemnet type of type %i!", etype);
			err_message( buf );
	}
}



id_type Elem::get_cut_edge(id_type e)
{
	if (e < _cutEdge.size())
		return _cutEdge[e];
	else
		err_message("Attempted to access an invalid cut edge.");
}
id_type Elem::get_inclusion_from_enrich(id_type e)
{
	if (e < _enrichment_nodes.size())
	{
		if (e < _enrichment_on_inclusion.size()) // Normal node
			return _enrichment_on_inclusion[e];
		else // Mirror node
			return _enrichment_on_inclusion[e - _enrichment_on_inclusion.size()];
	}
	else
		err_message("Attempted to access and invalid cut edge.");
}

void Elem::generate_integration_elem()
{
	if(!_detected)
		err_message("Must call element detection before attempting to create the elements inegration subelements.");
	
	// Build all the integration elements
	for(id_type i=0; i<_int_elem_struct.size(); ++i)
	{
		Elem* el;
		switch(dim())
		{
			case 1:
				if(_int_elem_struct[i].size()==2) el = Elem::build(EDGE2);
				else err_message("Attempted to build an unknown 1D integration subelement.");
				break;
			case 2:
				if(_int_elem_struct[i].size()==3) el = Elem::build(TRI3);
				else if(_int_elem_struct[i].size()==4) el = Elem::build(QUAD4);
				else if(_int_elem_struct[i].size() > 4) el = new Polygon(_int_elem_struct[i].size()); // Manually create a polygon because they'll only be created here and no where else
				else err_message("Attempted to build an unknown 2D integration subelement.");
				break;
			case 3:
				if(_int_elem_struct[i].size()==4) el = Elem::build(TET4);
				else if(_int_elem_struct[i].size()==6) el = Elem::build(PRISM6);
				else if(_int_elem_struct[i].size()==8) el = Elem::build(HEX8);
				else err_message("Attempted to build an unknown 3D integration subelement.");
				break;
			default:
				err_message("Integration Element Generation: Somehow the dimension of an element is not 1, 2, or 3...");
		}
		for(id_type j=0; j<_int_elem_struct[i].size(); ++j)
		{
			if(_int_elem_struct[i][j]<n_nodes())
				el->set_node(j) = _nodes[_int_elem_struct[i][j]];
			else
			{
				id_type enode = _int_elem_struct[i][j] - n_nodes();
				if (enode<0 || enode>=_enrichment_nodes.size())
					err_message("Error impending here.");
				el->set_node(j) = _enrichment_nodes[_int_elem_struct[i][j] - n_nodes()];
			}
		}
		el->set_id(i);
		_integrationElem.push_back(el);
	}

	// Build all of the cohesive elements
	for(id_type i=0; i<_coh_elem_struct.size(); ++i)
	{
		CohesiveElem* el;
		switch(dim())
		{
			case 1:
				if(_coh_elem_struct[i].size()==2) el = CohesiveElem::build(COHPOINT1);
				else err_message("Attempted to build an unknown 0D cohesive subelement.");
				break;
			case 2:
				if(_coh_elem_struct[i].size()==4) el = CohesiveElem::build(COHEDGE2);
				else err_message("Attempted to build an unknown 1D cohesive subelement.");
				break;
			case 3:
				if(_coh_elem_struct[i].size()==6) el = CohesiveElem::build(COHTRI3);
				else if(_coh_elem_struct[i].size()==8) el = CohesiveElem::build(COHQUAD4);
				else err_message("Attempted to build an unknown 2D cohesive subelement.");
				break;
			default:
				err_message("Cohesive Element Generation: Somehow the dimension of an element is not 1, 2, or 3...");
		}
		for(id_type j=0; j<_coh_elem_struct[i].size(); ++j)
		{
			if(_coh_elem_struct[i][j]<n_nodes())
				el->set_node(j) = _nodes[_coh_elem_struct[i][j]];
			else
				el->set_node(j) = _enrichment_nodes[_coh_elem_struct[i][j] - n_nodes()];
		}
		el->set_id(i);
		_cohesiveElem.push_back(el);
	}

	// Release the useless memory in the struct variables
	// _int_elem_struct.swap( std::vector<std::vector<short_id_type> >() );
	//_coh_elem_struct.swap( std::vector<std::vector<short_id_type> >() );
}



// Compute the objectve fucntion that we are trying to zero
int inv_map_f(const gsl_vector * x, void *params, gsl_vector * f)
{
	// Get the goal global position
	const std::vector<double> goal = ((struct inv_map_params *) params)->goal;
	Elem* elem = ((struct inv_map_params *) params)->elem;
	id_type dim = elem->dim();
	if(goal.size() != dim)
		err_message("The dimensions of the goal vector must match the number of dimensions of the element.");

	// Get the current parent coordinates
	std::vector<double> pcoords(dim);
	for(id_type d=0; d<dim; ++d)
		pcoords[d] = gsl_vector_get(x, d);

	// Compute the shape function values at this point
	std::vector<double> N = elem->compute_shape(pcoords);

	// Translate these to global coordinates
	double gcoords[dim];
	for(id_type d=0; d<dim; ++d)
		gcoords[d] = 0.0;
	for(id_type n=0; n<elem->n_nodes(); ++n)
		for(id_type d=0; d<dim; ++d)
			gcoords[d] += (*elem->get_node(n))(d) * N[n];

	// Compute the difference between this x and the goal
	for(id_type d=0; d<dim; ++d)
		gsl_vector_set(f, d, (gcoords[d]-goal[d]));

	return GSL_SUCCESS;
}

// Compute the jacobian of the function that we are trying to zero
int inv_map_df(const gsl_vector * x, void *params, gsl_matrix * J)
{
	// Get the goal global position
	const std::vector<double> goal = ((struct inv_map_params *) params)->goal;
	Elem* elem = ((struct inv_map_params *) params)->elem;
	id_type dim = elem->dim();
	if(goal.size() != dim)
		err_message("The dimensions of the goal vector must match the number of dimensions of the element.");

	// Get the current parent coordinates
	std::vector<double> pcoords(dim);
	for(id_type d=0; d<dim; ++d)
		pcoords[d] = gsl_vector_get(x, d);

	// Compute the shape function gradients at this point
	std::vector<std::vector<double> > dN = elem->compute_shape_grad(pcoords);

	// Compute the Jacobian fo the transformation
	double jac[dim][dim];
	for(id_type d1=0; d1<dim; ++d1) // Loop ver the global coordinates
		for(id_type d2=0; d2<dim; ++d2) // Lopp over the local coordinates
			for(id_type n=0; n<elem->n_nodes(); ++n) // Loop over the nodes
				jac[d1][d2] += (*elem->get_node(n))(d1) * dN[n][d2];

	// Store the jacobian
	for(id_type d1=0; d1<dim; ++d1)
		for(id_type d2=0; d2<dim; ++d2)
			gsl_matrix_set(J, d1, d2, jac[d1][d2]);

	return GSL_SUCCESS;
}

// Compute the objective and the Jacobian
int inv_map_fdf(const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J)
{
	// Get the goal global position
	const std::vector<double> goal = ((struct inv_map_params *) params)->goal;
	Elem* elem = ((struct inv_map_params *) params)->elem;
	id_type dim = elem->dim();
	if(goal.size() != dim)
		err_message("The dimensions of the goal vector must match the number of dimensions of the element.");

	// Get the current parent coordinates
	std::vector<double> pcoords(dim);
	for(id_type d=0; d<dim; ++d)
		pcoords[d] = gsl_vector_get(x, d);

	// Compute the shape function values at this point
	std::vector<double> N = elem->compute_shape(pcoords);

	// Translate these to global coordinates
	double gcoords[dim];
	for(id_type d=0; d<dim; ++d)
		gcoords[d] = 0.0;
	for(id_type n=0; n<elem->n_nodes(); ++n)
		for(id_type d=0; d<dim; ++d)
			gcoords[d] += (*elem->get_node(n))(d) * N[n];

	// Compute the difference between this x and the goal
	for(id_type d=0; d<dim; ++d)
		gsl_vector_set(f, d, (gcoords[d]-goal[d]));

	// Compute the shape function gradients at this point
	std::vector<std::vector<double> > dN = elem->compute_shape_grad(pcoords);

	// Compute the Jacobian fo the transformation
	double jac[dim][dim];
	for(id_type d1=0; d1<dim; ++d1) // Loop over the global coordinates
		for(id_type d2=0; d2<dim; ++d2) // Loop over the local coordinates
			for(id_type n=0; n<elem->n_nodes(); ++n) // Loop over the nodes
				jac[d1][d2] += (*elem->get_node(n))(d1) * dN[n][d2];


	// Store the jacobian
	for(id_type d1=0; d1<dim; ++d1)
		for(id_type d2=0; d2<dim; ++d2)
			gsl_matrix_set(J, d1, d2, jac[d1][d2]);

	return GSL_SUCCESS;
}



















// Function to compute the inverse map of the globalcoordinates
// This will be used for elements with nonlinear shape function.
// For simpler elements (triangle) explicit forms can be found so this will be overridden
std::vector<double> Elem::inverse_map(std::vector<double> gcoords, bool& inside)
{
	if(gcoords.size() != dim())
		err_message("The number of coordinates passed into the inverse map must match the dimension of the element.");
	// Should probably also have a check to see if the global coordinates are even within the bounds of the element...


	// Create pointers to the solver objects
	const gsl_multiroot_fdfsolver_type *T;
	gsl_multiroot_fdfsolver *s;

	// Define some parameters
	int status;
	size_t iter = 0;
	const size_t n = dim();
	struct inv_map_params p = {gcoords, this};

	// Attach the objective and Jacobian functions
	gsl_multiroot_function_fdf f = {&inv_map_f,
									&inv_map_df,
									&inv_map_fdf,
									n, &p};

	// Set the initial guess (all zeros is as good as anything I can think of I guess)
	double x_init[n];
	for(size_t nn=0; nn<n; ++nn) x_init[nn] = 0.0;
	gsl_vector * x = gsl_vector_alloc (n);
	for(size_t it=0; it<n; ++it)
		gsl_vector_set (x, it, x_init[it]);

	// Set the type of solver and allocate memory for the solver object
	T = gsl_multiroot_fdfsolver_hybridj; // Can play with the type to see what works best
	s = gsl_multiroot_fdfsolver_alloc (T, n);
	gsl_multiroot_fdfsolver_set (s, &f, x);

	// print_state (iter, s);
	do
	{
		iter++;

		// Iterate the solver
		status = gsl_multiroot_fdfsolver_iterate (s);
		// print_state (iter, s);
		if (status)
			break;

		// Check the residual magnitude
		status = gsl_multiroot_test_residual (s->f, 1e-12);
	}
	while (status == GSL_CONTINUE && iter < 1000);

	// Retrieve the solution
	std::vector<double> sol(dim());
	for(id_type d=0; d<dim(); ++d)
		sol[d] = gsl_vector_get(s->x, d);

	// Set the boolean for if this point is actually inside this element
	inside = coords_inside(sol);

	// Cleanup memory
	gsl_multiroot_fdfsolver_free (s);
	gsl_vector_free (x);
	return sol;
}






bool Elem::belongs_to_edge(id_type edge, id_type id)
{
	for (id_type n=0; n<n_nodes_on_edge(edge); ++n)
		if(_nodes[edge_nodes_map(edge, n)]->get_id() == id)
			return true;
	return false;
}

bool Elem::belongs_to_side(id_type side, id_type id)
{
	for (id_type n=0; n<n_nodes_on_side(side); ++n)
		if(_nodes[side_nodes_map(side, n)]->get_id() == id)
			return true;
	return false;
}


int Elem::get_edge_from_nodes(std::vector<id_type>& nodes)
{
	// Get the edges that the nodes belongs to
	std::vector<std::set<id_type> > possible_edges(nodes.size());
	for (id_type n=0; n<nodes.size(); ++n)
	{
		for (id_type e=0; e<n_edges(); ++e)
			if ( belongs_to_edge(e, nodes[n]) )
				possible_edges[n].insert(e);
	}

	std::set<id_type> intersect, temp_intersect;
	intersect = possible_edges[0];
	for (id_type n=1; n<nodes.size(); ++n)
	{
		std::set_intersection(intersect.begin(), intersect.end(),
							  possible_edges[n].begin(), possible_edges[n].end(),
							  std::inserter(temp_intersect, temp_intersect.begin()));
		intersect.swap( temp_intersect );
		temp_intersect.clear();
	}

	// If there's more than one match then something went wrong
	if(intersect.size() != 1)
		return -1;

	// Otherwise just return the first element
	else
		return *intersect.begin();
}


int Elem::get_side_from_nodes(std::vector<id_type>& nodes)
{
	// Get the sides that the nodes belongs to
	std::vector<std::set<id_type> > possible_sides(nodes.size());
	for (id_type n=0; n<nodes.size(); ++n)
	{
		for (id_type e=0; e<n_sides(); ++e)
			if ( belongs_to_side(e, nodes[n]) )
				possible_sides[n].insert(e);
	}

	std::set<id_type> intersect, temp_intersect;
	intersect = possible_sides[0];
	for (id_type n=1; n<nodes.size(); ++n)
	{
		std::set_intersection(intersect.begin(), intersect.end(),
							  possible_sides[n].begin(), possible_sides[n].end(),
							  std::inserter(temp_intersect, temp_intersect.begin()));
		intersect.swap( temp_intersect );
	}

	// If there's more than one match then something went wrong
	if(intersect.size() != 1)
		return -1;

	// Otherwise just return the first element
	else
		return *intersect.begin();
}






