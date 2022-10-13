/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#include "CohesiveElem.h"
#include "CohPoint1.h"
#include "CohEdge2.h"
#include "CohTri3.h"
#include "CohQuad4.h"
#include "Utilities.h"
#include <cmath>

CohesiveElem::CohesiveElem()
	: _global_id(0)
{
}
CohesiveElem::CohesiveElem(std::vector<Node*> nodes, id_type id)
{
	_nodes = nodes;
	_global_id = id;
}
CohesiveElem::~CohesiveElem()
{
}

void CohesiveElem::clear()
{
	_nodes.clear();
	_global_id = 0;
	
}

void CohesiveElem::copy(CohesiveElem & other_elem)
{
	other_elem.set_id(_global_id);
	other_elem.set_nodes(_nodes);
}

CohesiveElem::CohesiveElem(CohesiveElem & other_elem)
{
	copy(other_elem);
}

Node& CohesiveElem::operator()(id_type i)
{
	if(i < n_nodes())
		return *_nodes[i];
	else
		err_message("Node index must be less than the number of nodes.");
}

void CohesiveElem::print(std::ostream & os) const
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
		default:
			os << "Non-supported type";
			break;
	}
	os << ")" << std::endl;
	for(id_type i=0; i<n_nodes(); ++i)
		os << "\t" << *_nodes[i];
}



// Output helper function
std::vector<id_type> CohesiveElem::plot_elem_ids()
{
	std::vector<id_type> v(n_nodes());
	for (id_type n=0; n<_nodes.size(); ++n)
		v[n] = _nodes[n]->get_id();
	return v;
}




CohesiveElem* CohesiveElem::build(coh_elem_type etype)
{
	switch(etype)
	{
		case POINT1:
			return new CohPoint1;
		case EDGE2:
			return new CohEdge2;
		case TRI3:
			return new CohTri3;
		case QUAD4:
			return new CohQuad4;
		case EDGE3:
		case TRI6:
		case QUAD8:
		default:
			err_message("Attempting to build a non-supported elemnet type of type!");
	}
}

















void CohesiveElem::ShapeFunctions(std::vector<double> rcoords, std::vector<double>& N, DenseMatrix<double>& Rot, double& J)
{
	// Fill in the shape functions and local gradients (currently don't do the gradients for cohesive elements
	id_type nn = n_nodes() / 2;
	N.resize( nn );
	for (id_type n=0; n<nn; ++n)
		N[n] = compute_shape(rcoords, n);

	// Store the rotation matrix and the area jacobian
	compute_rotation(rcoords, Rot, J);
}





void CohesiveElem::compute_rotation(const std::vector<double>& coords, DenseMatrix<double>& mat, double& dA) const
{
	// Assemble the nodal coordinates of the undeformed surface
	std::vector<std::vector<double> > node_coords(n_nodes()/2);
	for(id_type n=0; n<node_coords.size(); ++n)
	{
		node_coords[n].resize(dim()+1);
		for(id_type d=0; d<=dim(); ++d)
			node_coords[n][d] = (*_nodes[n])(d);
	}

	// The actual call to the inherited function
	compute_rotation_general(coords, node_coords, mat, dA);
}






// This is not used at the moment but in case I feel like changing the code late I'll leave it
void CohesiveElem::compute_deformed_rotation(const std::vector<double>& coords, const std::vector<double>& UEL,
											 DenseMatrix<double>& mat, double& dA) const
{
	// Assemble the nodal coordinates of the deformed surface
	id_type half_nodes = n_nodes()/2;
	std::vector<std::vector<double> > node_coords(n_nodes()/2);
	for(id_type n=0; n<node_coords.size(); ++n)
	{
		node_coords[n].resize(dim()+1);
		for(id_type d=0; d<(dim()+1); ++d)
		{
			double node_minus = (*_nodes[n])(d) + UEL[n*dim() + d];
			double node_plus = (*_nodes[n+half_nodes])(d) + UEL[(n+half_nodes)*dim() + d];
			//double node_minus = (*_nodes[n])(d) + UEL[n][d];
			//double node_plus = (*_nodes[n+half_nodes])(d) + UEL[n+half_nodes][d];
			node_coords[n][d] = (node_minus+node_plus)/2;
		}
	}

	// The actual call to the inherited function
	compute_rotation_general(coords, node_coords, mat, dA);
}



void CohesiveElem::compute_rotation_general(const std::vector<double>& coords, const std::vector<std::vector<double> >& node_coords,
										DenseMatrix<double>& mat, double& dA) const
{	// Compute dx/dr
	std::vector<std::vector<double> > Jac(dim(), std::vector<double>(dim()+1, 0.0));
	for(id_type n=0; n<node_coords.size(); ++n) // Loop over the nodes
	{
		std::vector<double> parent_grad = compute_shape_grad(coords, n);
		for (id_type alpha=0; alpha<dim(); ++alpha)
			for (id_type i=0; i<(dim()+1); ++i)
				Jac[alpha][i] += parent_grad[alpha] * node_coords[n][i];
	}

	// Get the vector normal to the other vectors (cross product)
	std::vector<double> norm;
	if (dim() == 1) // 2D simulation
		norm = {-Jac[0][1], Jac[0][0]};
	else if (dim() == 2) // 3D simulation
		norm = Utilities::cross(Jac[0], Jac[1]);

	std::vector<double> magnitudes(dim());
	for (id_type alpha=0; alpha<dim(); ++alpha)
		magnitudes[alpha] = sqrt(Utilities::dot(Jac[alpha], Jac[alpha]));

	dA = sqrt(Utilities::dot(norm, norm));

	mat.resize(dim()+1, dim()+1);
	for (id_type d=0; d<(dim()+1); ++d)
		mat(dim(), d) = norm[d] / dA;
	for (id_type d=0; d<(dim()+1); ++d)
		for (id_type alpha=0; alpha<dim(); ++alpha)
			mat(alpha, d) = Jac[alpha][d] / magnitudes[alpha];
}





void CohesiveElem::RotationSensitivity(const std::vector<double>& rcoords, const std::vector<double>& coh_velocity, 
									   DenseMatrix<double>& dRot, std::vector<double>& v, double& div_v_surf)
{
	// Assume undeformed here...
	std::vector<std::vector<double> > node_coords(n_nodes()/2, std::vector<double>(dim()+1));
	for(id_type n=0; n<node_coords.size(); ++n)
		for(id_type d=0; d<=dim(); ++d)
			node_coords[n][d] = (*_nodes[n])(d);

	// Get all of the nodal shape function gradients (in local coords)
	std::vector<std::vector<double> > parent_grad(node_coords.size());
	for (id_type n=0; n< node_coords.size(); ++n)
		parent_grad[n] = compute_shape_grad(rcoords, n);

	// Compute dx/dr
	std::vector<std::vector<double> > Jac(dim(), std::vector<double>(dim()+1, 0.0));
	std::vector<std::vector<double> > dJac(dim(), std::vector<double>(dim()+1, 0.0));
	for(id_type n=0; n<node_coords.size(); ++n) // Loop over the nodes
	{
		for (id_type alpha=0; alpha<dim(); ++alpha)
			for (id_type i=0; i<(dim()+1); ++i)
			{
				Jac[alpha][i] += parent_grad[n][alpha] * node_coords[n][i];
				dJac[alpha][i] += parent_grad[n][alpha] * coh_velocity[n*(dim()+1) + i];
			}	
	}

	// Compute the velocity here
	v.resize(dim()+1);
	for (id_type d=0; d<v.size(); ++d)
	{
		v[d] = 0.0;
		for (id_type n=0; n<node_coords.size(); ++n)
		{
			double N = compute_shape(rcoords, n);
			v[d] += N * coh_velocity[n*(dim()+1) + d];
		}
	}

	// Compute the derivative of the rotation matrix
	dRot.resize(dim()+1, dim()+1);
	if (dim() == 1)
		dRot(0,0) = 0.0;
	if (dim() == 1) // 2D problem
	{
		double theta = atan2(Jac[0][1], Jac[0][0]);
		double A = sqrt(pow(Jac[0][0],2) + pow(Jac[0][1],2));
		double dA = (Jac[0][0]*dJac[0][0]+Jac[0][1]*dJac[0][1])/A;
		double dtheta;
		if (std::fabs(cos(theta)) < 1e-3) // Default to using the cosine formulation unless it is ill conditioned here
			dtheta = (A*dJac[0][0]-Jac[0][0]*dA)/(-1.0*sin(theta)*A*A);
		else
			dtheta = (A*dJac[0][1]-Jac[0][1]*dA)/(cos(theta)*A*A);
		double cn = cos(theta+PI/2.0);
		double sn = sin(theta+PI/2.0);
		std::vector<double> tangent = {cos(theta), sin(theta)}; // really just the scale Jacobian
		std::vector<double> norm = {cn, sn}; // again, should be the same thing as in compute_rotation_general
		std::vector<double> dtangent = {-sin(theta)*dtheta, cos(theta)*dtheta};
		std::vector<double> dnorm = {-sn*dtheta, cn*dtheta};
		dRot(0,0) = dtangent[0];
		dRot(0,1) = dtangent[1];
		dRot(1,0) = dnorm[0];
		dRot(1,1) = dnorm[1];

		// Useful in the divergence calculation
		double val = 1.0 / (norm[1]*Jac[0][0] - norm[0]*Jac[0][1]);

		// Compute the divergence of the velocity
		div_v_surf = 0.0;
		for(id_type n=0; n<node_coords.size(); ++n) // Loop over the nodes
		{
			std::vector<double> grad_x = {parent_grad[n][0] * norm[1] * val, -1.0 * parent_grad[n][0] * norm[0] * val};
			for (id_type i=0; i<=dim(); ++i)
				div_v_surf += coh_velocity[n*(dim()+1) + i] * grad_x[i];
		}
	}
	else if (dim() == 2) // 3D problem
	{
		std::vector<double>& t0 = Jac[0];
		std::vector<double>& dt0_dd = dJac[0];
		double t0_mag = sqrt(Utilities::dot(t0, t0));
		std::vector<double>& t1 = Jac[1];
		std::vector<double>& dt1_dd = dJac[1];
		double t1_mag = sqrt(Utilities::dot(t1, t1));
		std::vector<double> norm = Utilities::cross(Jac[0], Jac[1]);
		std::vector<double> dnorm_dd = Utilities::_plus(Utilities::cross(dJac[0], Jac[1]), Utilities::cross(Jac[0], dJac[1])); // Product rule
		double norm_mag = sqrt(Utilities::dot(norm, norm));

		std::vector<std::vector<double> > vec_dRot(3);
		vec_dRot[0] = Utilities::_scale(Utilities::_minus( dt0_dd, Utilities::_scale( t0, (Utilities::dot( dt0_dd, t0 ))*(1./pow(t0_mag,2)) ) ),
								  (1./t0_mag));
		vec_dRot[1] = Utilities::_scale(Utilities::_minus( dt1_dd, Utilities::_scale( t1, (Utilities::dot( dt1_dd, t1 ))*(1./pow(t1_mag,2)) ) ),
								  (1./t1_mag));
		vec_dRot[2] = Utilities::_scale(Utilities::_minus( dnorm_dd, Utilities::_scale( norm, (Utilities::dot( dnorm_dd, norm ))*(1./pow(norm_mag,2)) ) ),
								  (1./norm_mag));
		for (id_type alpha=0; alpha<3; ++alpha)
			for (id_type i=0; i<3; ++i)
				dRot(alpha, i) = vec_dRot[alpha][i];

		// Useful in the divergence calculation
		double val = 1.0 /
					 (norm[0]*(Jac[0][1]*Jac[1][2]-Jac[0][2]*Jac[1][1]) -
					  norm[1]*(Jac[0][0]*Jac[1][2]-Jac[0][2]*Jac[1][0]) +
					  norm[2]*(Jac[0][0]*Jac[1][1]-Jac[0][1]*Jac[1][0]));

		// Compute the divergence of the velocity
		div_v_surf = 0.0;
		for(id_type n=0; n<node_coords.size(); ++n) // Loop over the nodes
		{
			double r = parent_grad[n][0];
			double s = parent_grad[n][1];
			std::vector<double> grad_x = {val * (norm[1]*(s*Jac[0][2]-r*Jac[1][2]) + norm[2]*(r*Jac[1][1]-s*Jac[0][1])),
								   -1.0 * val * (norm[0]*(s*Jac[0][2]-r*Jac[1][2]) + norm[2]*(r*Jac[1][0]-s*Jac[0][0])),
										  val * (norm[0]*(s*Jac[0][1]-r*Jac[1][1]) + norm[1]*(r*Jac[1][0]-s*Jac[0][0]))};
			for (id_type i=0; i<=dim(); ++i)
				div_v_surf += coh_velocity[n*(dim()+1) + i] * grad_x[i];
		}
	}
}