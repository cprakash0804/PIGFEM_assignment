/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated October 2016

##################################################################################
*/
#include "Utilities.h"
#include "elem.h"
#include "NodalData.h"
#include "gsl_cblas.h"  // For matrix computations
#include "gsl_linalg.h"
#include <limits> // std::numeric_limits





// Define some of my global variables here
// --------------------------------------------------------------------
std::vector<std::string> elem_type_names = {"INVALID_ELEM",
											"POINT1",
											"EDGE2",
											"EDGE3",
											"TRI3",
											"TRI6",
											"QUAD4",
											"QUAD8",
											"TET4", 
											"TET10",
											"HEX8",
											"HEX20",
											"HEX27",
											"PRISM6",
											"PRISM15",
											"PRISM18",
											"POLYGON"};
std::vector<std::string> coh_elem_type_names = {"INVALID_COHELEM",
												"COHPOINT1",
												"COHEDGE2",
												"COHEDGE3",
												"COHTRI3",
												"COHTRI6",
												"COHQUAD4",
												"COHQUAD8"};
std::vector<std::string> material_type_names = {"INVALID_MATERIAL",
												"LINEAR_ELASTIC_ISOTROPIC",
												"LINEAR_THERMAL",
												"CONTINUUM_DAMAGE",
												"OP_COHESIVE",
												"OP_COHESIVE_NO_UNLOADING",
												"LINEAR_ELASTIC_ISOTROPIC_MAX_PRINCIPAL_STRESS",
												"LINEAR_ELASTIC_TRANSVERSELY_ISOTROPIC",
												"LINEAR_ELASTIC_ISOTROPIC_PROBLEM_CONTROLLED_DAMAGE"};
std::vector<std::string> problem_type_names = {"LINEAR ELASTICITY",
											   "NONLINEAR STRUCTURAL",
											   "NONLINEAR DAMAGE",
											   "LINEAR HEAT",
											   "NONLINEAR HEAT",
											   "PROBLEM MAX PRINCIPAL STRESS",
											   "NONLINEAR STRUCTURAL MULTISCALE",
											   "NONLINEAR STRUCTURAL THERMAL EXPANSION"};


// Function to assemble the "B" matrix used in assembly of the small strain elasticity stiffness matrix.
void ProblemUtilities::Assemble_Small_Strain_B_Mat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& grad_x)
{
	id_type n_nodes = grad_x.size();
	id_type dim;
	if(n_nodes!=0)
		dim = grad_x[0].size();
	else
		err_message("In the assemble B matrix function, the grad_x variable has to have something in it!");

	int nrows;
	if(dim==1)
		nrows=1;
	else if(dim==2)
		nrows=3;
	else if(dim==3)
		nrows=6;
	else
		err_message("Dimension for the B matrix must be less than or equal to 3.");

	B.resize(nrows, dim*n_nodes);
	for(id_type i=0; i<B.n_rows(); ++i)
		for(id_type j=0; j<B.n_cols(); ++j)
			B(i,j) = 0.0;
	
	// Assemble the axial strain components
	for(id_type n=0; n<n_nodes; ++n)
	{
		for(id_type i=0; i<dim; ++i)
		{
			int col = dim*n + i;
			B(i, col) = grad_x[n][i];
		}
	}
	
	// Assemble the shear strain components
	switch(dim)
	{
		case 2: // Two dimensional
			for(id_type n=0; n<n_nodes; ++n)
			{
				// gamma_xy
				B(2, 2*n) = grad_x[n][1];
				B(2, 2*n+1) = grad_x[n][0];
			}
			break;
		case 3:
			for(id_type n=0; n<n_nodes; ++n)
			{
				B(3, dim*n+1) = grad_x[n][2];	// dN_n/dz
				B(3, dim*n+2) = grad_x[n][1];	// dN_n/dy
				B(4, dim*n)   = grad_x[n][2];	// dN_n/dz
				B(4, dim*n+2) = grad_x[n][0];	// dN_n/dx
				B(5, dim*n)   = grad_x[n][1];	// dN_n/dy
				B(5, dim*n+1) = grad_x[n][0];	// dN_n/dx
			}
			break;
		default:
			break;
	}
}








// Function to assemble the "B" matrix used in assembly of the finite strain elasticity stiffness matrix.
void ProblemUtilities::Assemble_Finite_Strain_B_Mat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& grad_x)
{
	// Need to write this
}





// Function to assemble the "B" matrix used in assembly of the elasticity stiffness matrix.
void ProblemUtilities::Assemble_Thermal_B_Mat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& grad_x)
{
	id_type n_nodes = grad_x.size();
	id_type dim;
	if(n_nodes!=0)
		dim = grad_x[0].size();
	else
		err_message("In the assemble B matrix function, the grad_x variable has to have something in it!");

	B.resize(dim, n_nodes);
	for(id_type i=0; i<n_nodes; ++i)
		for(id_type j=0; j<dim; ++j)
			B(j, i) = grad_x[i][j];
}



// Function to assemble the "B" matrix used in assembly of the elasticity stiffness matrix.
void ProblemUtilities::Assemble_Cohesive_N_Mat(DenseMatrix<double>& N, const std::vector<double>& shapes, id_type dim)
{
	id_type n_nodes = shapes.size(); // Really there are twice as many nodes as there are shape functinos because of the mirror nodes

	N.clear();
	N.resize(dim, 2*dim*n_nodes);
	
	for(id_type n=0; n<n_nodes; ++n)
	{
		for(id_type i=0; i<dim; ++i)
		{
			N(i, dim*n+i) = -1.0*shapes[n]; // Shape function multipliers for the lower surface (nodes come first)
			N(i, dim*(n+n_nodes)+i) = shapes[n]; // Shape function multipliers for the upper surface (mirror nodes come second)
		}
	}
}














// This function essentially does Bt*D*B taking into account all of the zeros in the B matrix so its faster
void ProblemUtilities::SmallStrainKFast(DenseMatrix<double>& K, const std::vector<std::vector<double> >& dN, const DenseMatrix<double>& D)
{
	unsigned char n = dN.size();
	unsigned char dim = dN[0].size(); // Assumes there are nodes
	id_type ndof = n*dim;
	K.resize(ndof, ndof);

	switch (dim)
	{
		case 2: // Probably the most common case
			for (id_type row=0; row<ndof; ++row)
			{
				unsigned char i_star = row%2;
				unsigned char i_node = row/2;
				double part0, part1, part2;
				if (i_star == 0)
				{
					part0 = dN[i_node][0]*D(0,0) + dN[i_node][1]*D(2,0);
					part1 = dN[i_node][0]*D(0,1) + dN[i_node][1]*D(2,1);
					part2 = dN[i_node][0]*D(0,2) + dN[i_node][1]*D(2,2);
				}
				else
				{
					part0 = dN[i_node][1]*D(1,0) + dN[i_node][0]*D(2,0);
					part1 = dN[i_node][1]*D(1,1) + dN[i_node][0]*D(2,1);
					part2 = dN[i_node][1]*D(1,2) + dN[i_node][0]*D(2,2);
				}
				for (id_type col=0; col<ndof; ++col)
				{
					unsigned char j_star = col%2;
					unsigned char j_node = col/2;
					if (j_star==0)
						K(row, col) = part0*dN[j_node][0] + part2*dN[j_node][1];
					else
						K(row, col) = part1*dN[j_node][1] + part2*dN[j_node][0];
				}
			}
			break;




		case 3:
			for (id_type row=0; row<ndof; ++row)
			{
				unsigned char i_star = row%3;
				unsigned char i_node = row/3;
				double part0, part1, part2, part3, part4, part5;
				if (i_star == 0)
				{
					part0 = dN[i_node][0]*D(0,0) + dN[i_node][2]*D(4,0) + dN[i_node][1]*D(5,0);
					part1 = dN[i_node][0]*D(0,1) + dN[i_node][2]*D(4,1) + dN[i_node][1]*D(5,1);
					part2 = dN[i_node][0]*D(0,2) + dN[i_node][2]*D(4,2) + dN[i_node][1]*D(5,2);
					part3 = dN[i_node][0]*D(0,3) + dN[i_node][2]*D(4,3) + dN[i_node][1]*D(5,3);
					part4 = dN[i_node][0]*D(0,4) + dN[i_node][2]*D(4,4) + dN[i_node][1]*D(5,4);
					part5 = dN[i_node][0]*D(0,5) + dN[i_node][2]*D(4,5) + dN[i_node][1]*D(5,5);
				}
				else if (i_star == 1)
				{
					part0 = dN[i_node][1]*D(1,0) + dN[i_node][2]*D(3,0) + dN[i_node][0]*D(5,0);
					part1 = dN[i_node][1]*D(1,1) + dN[i_node][2]*D(3,1) + dN[i_node][0]*D(5,1);
					part2 = dN[i_node][1]*D(1,2) + dN[i_node][2]*D(3,2) + dN[i_node][0]*D(5,2);
					part3 = dN[i_node][1]*D(1,3) + dN[i_node][2]*D(3,3) + dN[i_node][0]*D(5,3);
					part4 = dN[i_node][1]*D(1,4) + dN[i_node][2]*D(3,4) + dN[i_node][0]*D(5,4);
					part5 = dN[i_node][1]*D(1,5) + dN[i_node][2]*D(3,5) + dN[i_node][0]*D(5,5);
				}
				else
				{
					part0 = dN[i_node][2]*D(2,0) + dN[i_node][1]*D(3,0) + dN[i_node][0]*D(4,0);
					part1 = dN[i_node][2]*D(2,1) + dN[i_node][1]*D(3,1) + dN[i_node][0]*D(4,1);
					part2 = dN[i_node][2]*D(2,2) + dN[i_node][1]*D(3,2) + dN[i_node][0]*D(4,2);
					part3 = dN[i_node][2]*D(2,3) + dN[i_node][1]*D(3,3) + dN[i_node][0]*D(4,3);
					part4 = dN[i_node][2]*D(2,4) + dN[i_node][1]*D(3,4) + dN[i_node][0]*D(4,4);
					part5 = dN[i_node][2]*D(2,5) + dN[i_node][1]*D(3,5) + dN[i_node][0]*D(4,5);
				}
				for (id_type col=0; col<ndof; ++col)
				{
					unsigned char j_star = col%3;
					unsigned char j_node = col/3;
					if (j_star == 0)
						K(row, col) = part0*dN[j_node][0] + part4*dN[j_node][2] + part5*dN[j_node][1];
					else if (j_star == 1)
						K(row, col) = part1*dN[j_node][1] + part3*dN[j_node][2] + part5*dN[j_node][0];
					else
						K(row, col) = part2*dN[j_node][2] + part3*dN[j_node][1] + part4*dN[j_node][0];
				}
			}
			break;




		case 1: // No difference in this case
			for (id_type row=0; row<ndof; ++row)
				for (id_type col=0; col<ndof; ++col)
					K(row, col) = dN[row][0] * D(0,0) * dN[col][0];
			break;
		default:
			err_message("Invalid dimension in fast small strain stiffness assembly");
	}
}


// This function essentially does B*U taking into account all of the zeros in the B matrix so its faster
void ProblemUtilities::SmallStrainFast(std::vector<double>& strain, const std::vector<std::vector<double> >& dN, const std::vector<double>& U)
{
	unsigned char n = dN.size();
	unsigned char dim = dN[0].size(); // Assumes there are nodes
	id_type ndof = n*dim;

	unsigned char ncomp;
	if (dim == 2)
		ncomp = 3;
	else if (dim == 3)
		ncomp = 6;
	else if (dim == 1)
		ncomp = 1;
	else
		err_message("Invalid dimension in fast small strain assembly");

	strain.resize(ncomp);
	std::fill(strain.begin(), strain.end(), 0.0);

	switch (dim)
	{
		case 2: // Most common case probably
			for (id_type node=0; node<n; ++node)
				strain[0] += dN[node][0] * U[node*2];
			for (id_type node=0; node<n; ++node)
				strain[1] += dN[node][1] * U[node*2+1];
			for (id_type node=0; node<n; ++node)
				strain[2] += dN[node][1]*U[node*2] + dN[node][0]*U[node*2+1];
			break;




		case 3:
			for (id_type node=0; node<n; ++node)
				strain[0] += dN[node][0] * U[node*3];
			for (id_type node=0; node<n; ++node)
				strain[1] += dN[node][1] * U[node*3+1];
			for (id_type node=0; node<n; ++node)
				strain[2] += dN[node][2] * U[node*3+2];
			for (id_type node=0; node<n; ++node)
				strain[3] += dN[node][2]*U[node*3+1] + dN[node][1]*U[node*3+2];
			for (id_type node=0; node<n; ++node)
				strain[4] += dN[node][2]*U[node*3] + dN[node][0]*U[node*3+2];
			for (id_type node=0; node<n; ++node)
				strain[5] += dN[node][1]*U[node*3] + dN[node][0]*U[node*3+1];
			break;




		case 1:
			for (id_type dof=0; dof<ndof; ++dof)
				strain[0] += dN[dof][0] * U[dof];
			break;



		default:
			err_message("Invalid dimension in fast small strain assembly");
	}
}

// This function performs Bt*sigma taking into account the zero structure of Bt
void ProblemUtilities::SmallStrainInternalForceFast(std::vector<double>& P, const std::vector<std::vector<double> >& dN, const std::vector<double>& sigma)
{
	unsigned char n = dN.size();
	unsigned char dim = dN[0].size(); // Assumes there are nodes
	id_type ndof = n*dim;

	P.resize(ndof);

	switch (dim)
	{
		case 2:
			for (id_type node=0; node<n; ++node)
			{
				P[node*2] = dN[node][0]*sigma[0] + dN[node][1]*sigma[2];
				P[node*2+1] = dN[node][1]*sigma[1] + dN[node][0]*sigma[2];
			}
			break;




		case 3:
			for (id_type node=0; node<n; ++node)
			{
				P[node*3] = dN[node][0]*sigma[0] + dN[node][2]*sigma[4] + dN[node][1]*sigma[5];
				P[node*3+1] = dN[node][1]*sigma[1] + dN[node][2]*sigma[3] + dN[node][0]*sigma[5];
				P[node*3+2] = dN[node][2]*sigma[2] + dN[node][1]*sigma[3] + dN[node][0]*sigma[4];
			}
			break;




		case 1:
			for (id_type node=0; node<n; ++node)
				P[node] = dN[n][0] * sigma[0];
			break;




		default:
			err_message("Invalid dimension in fast small strain internal force assembly");
	}
}


















// This function does Bt*D*B  (B-matrix variant)
void ProblemUtilities::SmallStrainKFast(DenseMatrix<double>& K, const DenseMatrix<double>& B, const DenseMatrix<double>& D)
{
	int in = D.n_rows();
	int ex = B.n_cols();
	DenseMatrix<double> BtD(ex, in);
	K.resize(ex, ex);
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, ex, in, in, 1.0, B.data(), ex, D.data(), in, 0.0, BtD.data(), in);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, ex, ex, in, 1.0, BtD.data(), in, B.data(), ex, 0.0, K.data(), ex);
}

// This function does B*U (B-matrix variant)
void ProblemUtilities::SmallStrainFast(std::vector<double>& strain, const DenseMatrix<double>& B, const std::vector<double>& U)
{
	strain.resize(B.n_rows());
	cblas_dgemv(CblasRowMajor, CblasNoTrans, B.n_rows(), B.n_cols(), 1.0, B.data(), B.n_cols(), U.data(), 1, 0.0, strain.data(), 1);
}

// This function performs Bt*sigma (B-matrix variant)
void ProblemUtilities::SmallStrainInternalForceFast(std::vector<double>& P, const DenseMatrix<double>& B, const std::vector<double>& sigma)
{
	P.resize(B.n_cols());
	cblas_dgemv(CblasRowMajor, CblasTrans, B.n_rows(), B.n_cols(), 1.0, B.data(), B.n_cols(), sigma.data(), 1, 0.0, P.data(), 1);
}













// This function computes the opening accross a cohesive surface from the shape functions, solution vector, and rotatino matrix
void ProblemUtilities::cohesiveDeltaFast(std::vector<double>& delta, const std::vector<double>& N, const std::vector<double>& U, const DenseMatrix<double>& R)
{
	id_type dim = R.n_rows();
	id_type nnodes = N.size();
	delta.clear();
	delta.resize(dim);

	for (id_type d=0; d<dim; ++d)
		for (id_type n=0; n<nnodes; ++n)
			delta[d] += N[n] * (U[(n+nnodes)*dim + d] - U[n*dim + d]);

	delta = R*delta;
}

// This function compute the stiffness matrix contribution from a cohesive opening
void ProblemUtilities::cohesiveKFast(DenseMatrix<double>& K, const std::vector<double>& N, const DenseMatrix<double>& D, const DenseMatrix<double>& R)
{
	DenseMatrix<double> Rt, RDR;
	R.transpose(Rt);
	RDR = Rt*D*R;

	id_type dim = R.n_rows();
	id_type nnodes = N.size();
	id_type ndof = 2 * nnodes * dim; // 2 because of mirror nodes
	K.resize(ndof, ndof);

	for (id_type n1=0; n1<nnodes; ++n1)
	{
		for (id_type n2=0; n2<nnodes; ++n2)
		{
			for (id_type d1=0; d1<dim; ++d1)
			{
				for (id_type d2=0; d2<dim; ++d2)
				{
					double val = N[n1] * N[n2] * RDR(d1, d2);
					K(n1*dim + d1, n2*dim + d2) = val;
					K((n1+nnodes)*dim + d1, n2*dim + d2) = -val;
					K(n1*dim + d1, (n2+nnodes)*dim + d2) = -val;
					K((n1+nnodes)*dim + d1, (n2+nnodes)*dim + d2) = val;
				}
			}
		}
	}
}

// This function compute the internal force contribution contribution from a cohesive opening
void ProblemUtilities::cohesiveInternalForceFast(std::vector<double>& P, const std::vector<double>& N, const std::vector<double>& trac, const DenseMatrix<double>& R)
{
	DenseMatrix<double> Rt;
	R.transpose(Rt);
	std::vector<double> tracR = Rt*trac;

	id_type dim = R.n_rows();
	id_type nnodes = N.size();
	id_type ndof = 2 * nnodes * dim; // 2 because of mirror nodes
	P.resize(ndof);

	for (id_type n=0; n<nnodes; ++n)
	{
		for (id_type d=0; d<dim; ++d)
		{
			double val = N[n] * tracR[d];
			P[n*dim + d] = -val;
			P[(n+nnodes)*dim + d] = val;
		}
	}
}



// This function forms the elemental solution vector
void ProblemUtilities::formElementalVector(std::vector<double>& elem_vec, Elem* el, NodalData* data)
{
	// Assemble the current displacement vector
	id_type nn = el->n_nodes() + el->n_enrich_nodes();
	id_type nndof = data->nndof();
	id_type ndof = nn*nndof;
	elem_vec.resize(ndof);
	for(id_type n=0; n<nn; ++n)
	{
		id_type g_node;
		if (n < el->n_nodes())
			g_node = el->get_node(n)->get_id();
		else
			g_node = el->get_enrich_node(n-el->n_nodes())->get_id();

		for(id_type d=0; d<nndof; ++d)
			elem_vec[n*nndof + d] = data->get_value_global(g_node, d);
	}
}












bool Utilities::penetration_store(bool a)
{
	static bool p = 0;
	if (a==0)
	{
		return p;
	}else{
		if(p==0)
		{
			p=1;
		}else{
			p=0;
		}
		return p;
	}
}

bool Utilities::penetration_step(bool a)
{
	static bool p = 0;
	if (a==0)
	{
		return p;
	}else{
		if(p==0)
		{
			p=1;
		}else{
			p=0;
		}
		return p;
	}
}




double Utilities::constrain_angle(double x)
{
	x = fmod(x + PI,2.0*PI);
	if (x < 0)
		x += 2.0*PI;
	return x - PI;
}


// Kroeneker Delta
int Utilities::kron(int i, int j)
{
	if(i==j)
		return 1;
	else
		return 0;
}


// Converts stress or strain to the principal values (2D requires plane strain boolean, defaults to false)
std::vector<double> Utilities::principal_stress_strain(const std::vector<double>& vec/*, bool plane_strain*/)
{
	id_type dim = 0;
	if (vec.size() == 1)
		dim = 1;
	else if (vec.size() == 3)
		dim = 2;
	else if (vec.size() == 6)
		dim = 3;
	else
		err_message("Invalid size of the stress vector.");

	std::vector<double> principal(dim);

	// 1 Dimensional
	if (dim == 1)
		principal[0] = vec[0];

	// 2 Dimensional
	else if (dim==2)
	{
		double part1 = (vec[0] + vec[1]) / 2;
		double part2 = sqrt(pow((vec[0]-vec[1])/2, 2) + pow(vec[2], 2));
		principal[0] = part1 + part2;
		principal[1] = part1 - part2;
	}
/*
	{
		// Plane strain case
		if (plane_strain)
		{
			double part1 = (vec[0] + vec[1]) / 2;
			double part2 = sqrt(pow((vec[0]-vec[1])/2, 2) + pow(vec[2], 2));
			principal[0] = part1 + part2;
			principal[1] = part1 - part2;
			principal.push_back( _nu * (vec[0]+vec[1]) ); // Add in out of plane stress as well and we'll just sort it in descending order
			std::sort(principal.begin(), principal.end(), std::greater<double>());
		}

		// Plane stress
		else
		{
			double part1 = (vec[0] + vec[1]) / 2;
			double part2 = sqrt(pow((vec[0]-vec[1])/2, 2) + pow(vec[2], 2));
			principal[0] = part1 + part2;
			principal[1] = part1 - part2;
		}
	}
*/

	// 3 Dimensional
	// Principal stresses are defined as the eigenvalues of the stress tensor, so find those
	else
	{
		double data[] = {vec[0], vec[5], vec[4],
						 vec[5], vec[1], vec[3],
						 vec[4], vec[3], vec[2]};
		gsl_matrix_view m = gsl_matrix_view_array(data, 3, 3);
		gsl_vector *eval = gsl_vector_alloc(3);
		gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(3);
		gsl_eigen_symm(&m.matrix, eval, w); // Actually calculate the eigen values of the stress matrix
		gsl_eigen_symm_free(w); // free the eigen system
		for (int i=0; i<3; ++i)
			principal[i] = gsl_vector_get(eval, i);
		std::sort(principal.begin(), principal.end(), std::greater<double>()); // Sort in descending order
		gsl_vector_free(eval);
	}

	return principal;
}






double Utilities::von_mises(const std::vector<double>& stress)
{
	if(stress.size() == 1)
		return stress[0];
	else if (stress.size() == 3)
		return sqrt(stress[0]*stress[0] - stress[0]*stress[1] + stress[1]*stress[1] + 3.0*stress[2]*stress[2]);
	else if (stress.size() == 6)
		return sqrt( 0.5*( pow(stress[0]-stress[1], 2) + pow(stress[1]-stress[2], 2) + pow(stress[2]-stress[0], 2) + 6.0*(pow(stress[3], 2) + pow(stress[4], 2) + pow(stress[5], 2))) );
	else
		err_message("Invalid stress vector for von Mises calculation");
}


float Utilities::von_mises(const std::vector<float>& stress)
{
	if(stress.size() == 1)
		return stress[0];
	else if (stress.size() == 3)
		return sqrt(stress[0]*stress[0] - stress[0]*stress[1] + stress[1]*stress[1] + 3.0*stress[2]*stress[2]);
	else if (stress.size() == 6)
		return sqrt( 0.5*( pow(stress[0]-stress[1], 2) + pow(stress[1]-stress[2], 2) + pow(stress[2]-stress[0], 2) + 6.0*(pow(stress[3], 2) + pow(stress[4], 2) + pow(stress[5], 2))) );
	else
		err_message("Invalid stress vector for von Mises calculation");
}




// Some basic vector operations
double Utilities::dot(const std::vector<double>& v1, const std::vector<double>& v2)
{
	if (v1.size() != v2.size())
		err_message("vectors must be the same size for a dot product");

	return cblas_ddot(v1.size(), v1.data(), 1, v2.data(), 1);
	// double ret = 0.0;
	// for (id_type i=0; i<v1.size(); i++)
	// 	ret += v1[i]*v2[i];

	// return ret;
}
double Utilities::perp(const std::vector<double>& v1, const std::vector<double>& v2)
{
	if (v1.size()!=v2.size() || v1.size()!=2)
		err_message("vectors must be the same size for a perp product");

	return v1[0]*v2[1] - v1[1]*v2[0];
}
std::vector<double> Utilities::cross(const std::vector<double>& v1, const std::vector<double>& v2)
{
	if (v1.size() != v2.size())
		err_message("To cross vectors they must be the same size");

	std::vector<double> ret(v1.size());
	if (v1.size() == 2)
		ret[2] = v1[0]*v2[1] - v1[1]*v2[0];
	else if (v1.size() == 3)
	{
		ret[0] = v1[1]*v2[2] - v1[2]*v2[1];
		ret[1] = -(v1[0]*v2[2] - v1[2]*v2[0]);
		ret[2] = v1[0]*v2[1] - v1[1]*v2[0];
	}
	else
		err_message("Unknown definition of a cross product for vectors that aren't 2D or 3D");
	return ret;
}
std::vector<double> Utilities::plus(const std::vector<double>& v1, const std::vector<double>& v2)
{
	if (v1.size() != v2.size())
		err_message("vectors must be the same size to add themt");

	std::vector<double> ret(v1.size());
	for (id_type i=0; i<v1.size(); ++i)
		ret[i] = v1[i] + v2[i];
	return ret;
}
std::vector<double> Utilities::minus(const std::vector<double>& v1, const std::vector<double>& v2)
{
	if (v1.size() != v2.size())
		err_message("vectors must be the same size to subtract themt");

	std::vector<double> ret(v1.size());
	for (id_type i=0; i<v1.size(); ++i)
		ret[i] = v1[i] - v2[i];
	return ret;
}
std::vector<double> Utilities::VecAXPY(double alpha, const std::vector<double>& x, const std::vector<double>& y)
{
	if (x.size() != y.size())
		err_message("vectors must be the same size to perform an AXPY operation");

	std::vector<double> ret(x.size());
	for (id_type i=0; i<x.size(); ++i)
		ret[i] = alpha*x[i] + y[i];
	return ret;
}
void Utilities::_VecAXPY(double alpha, std::vector<double>& x, std::vector<double>& y)
{
	if (x.size() != y.size())
		err_message("vectors must be the same size to perform an AXPY operation");

	cblas_daxpy(x.size(), alpha, x.data(), 1, y.data(), 1);
}

std::vector<double> Utilities::scale(const std::vector<double>& vec, double alpha)
{
	std::vector<double> ret(vec.size());
	for (id_type i=0; i<vec.size(); ++i)
		ret[i] = alpha*vec[i];
	return ret;
}


/*
 * PASS BY VALUE VERSIONS
 */
double Utilities::_dot(std::vector<double> v1, std::vector<double> v2)
{
	if (v1.size() != v2.size())
		err_message("vectors must be the same size for a dot product");

	return cblas_ddot(v1.size(), v1.data(), 1, v2.data(), 1);
}
double Utilities::_perp(std::vector<double> v1, std::vector<double> v2)
{
	if (v1.size()!=v2.size() || v1.size()!=2)
		err_message("vectors must be the same size for a perp product");

	return v1[0]*v2[1] - v1[1]*v2[0];
}
std::vector<double> Utilities::_cross(std::vector<double> v1, std::vector<double> v2)
{
	if (v1.size() != v2.size())
		err_message("To cross vectors they must be the same size");

	std::vector<double> ret(v1.size());
	if (v1.size() == 2)
		ret[2] = v1[0]*v2[1] - v1[1]*v2[0];
	else if (v1.size() == 3)
	{
		ret[0] = v1[1]*v2[2] - v1[2]*v2[1];
		ret[1] = -(v1[0]*v2[2] - v1[2]*v2[0]);
		ret[2] = v1[0]*v2[1] - v1[1]*v2[0];
	}
	else
		err_message("Unknown definition of a cross product for vectors that aren't 2D or 3D");
	return ret;
}
std::vector<double> Utilities::_plus(std::vector<double> v1, std::vector<double> v2)
{
	_VecAXPY(1.0, v1, v2);
	return v2;
}
std::vector<double> Utilities::_minus(std::vector<double> v1, std::vector<double> v2)
{
	_VecAXPY(-1.0, v2, v1); // Does v1=-v2+v1
	return v1;
}
std::vector<double> Utilities::_scale(std::vector<double> vec, double alpha)
{
	cblas_dscal(vec.size(), alpha, vec.data(), 1);
	return vec;
}


// Really should use GSL for this but I don't feel like figuring out how to use that right now
double Utilities::det(DenseMatrix<double>& J)
{
	if (J.n_rows() != J.n_cols())
		err_message("For computation of a determinant, the matrix must be square!");

	if (J.n_rows() == 1)
		return J(0,0);
	else if (J.n_rows() == 2)
		return J(0,0)*J(1,1) - J(0,1)*J(1,0);
	else if (J.n_rows() == 3)
		return J(0,0)*(J(1,1)*J(2,2)-J(1,2)*J(2,1)) - J(0,1)*(J(1,0)*J(2,2)-J(1,2)*J(2,0)) + J(0,2)*(J(1,0)*J(2,1)-J(1,1)*J(2,0));
	else
		err_message("Determinant is only supported for up to 3D matricies");

}
DenseMatrix<double> Utilities::inverse(DenseMatrix<double>& J)
{
	if (J.n_rows() != J.n_cols())
		err_message("For computation of an inverse, the matrix must be square!");

	id_type dim = J.n_rows();
	DenseMatrix<double> ret(dim, dim);
	if (J.n_rows() == 1)
		ret(0,0) = 1.0;
	else if (J.n_rows() == 2)
	{
		ret(0,0) = J(1,1);
		ret(1,1) = J(0,0);
		ret(0,1) = -1.0*J(0,1);
		ret(1,0) = -1.0*J(1,0);
	}
	else if (J.n_rows() == 3)
	{
		ret(0,0) = J(1,1)*J(2,2) - J(2,1)*J(1,2);
		ret(0,1) = -1.0*(J(0,1)*J(2,2) - J(2,1)*J(0,2));
		ret(0,2) = J(0,1)*J(1,2) - J(1,1)*J(0,2);
		
		ret(1,0) = -1.0*(J(1,0)*J(2,2) - J(2,0)*J(1,2));
		ret(1,1) = J(0,0)*J(2,2) - J(2,0)*J(0,2);
		ret(1,2) = -1.0*(J(0,0)*J(1,2) - J(1,0)*J(0,2));
		
		ret(2,0) = J(1,0)*J(2,1) - J(2,0)*J(1,1);
		ret(2,1) = -1.0*(J(0,0)*J(2,1) - J(2,0)*J(0,1));
		ret(2,2) = J(0,0)*J(1,1) - J(1,0)*J(0,1);
	}
	else
		err_message("Inverse is only supported for up to 3D matricies");

	double d = Utilities::det(J);
	ret *= (1.0 / d);
	return ret;
}






// Finds the intersection of Line P0-P1 and line V0-V1
// Returns the multiple of the distance from P0 to P1 that the intersection occurs at
double Utilities::LineLineIntersection(const std::vector<double>& P0, const std::vector<double>& dS,
										 const std::vector<double>& V0, const std::vector<double>& e)
{
	double D = -Utilities::perp(e, dS);		// dot(ne, dS)
	if (fabs(D) < 1e-14) // S is nearly parallel to this edge
		return std::numeric_limits<double>::max();
	else
	{
		std::vector<double> temp1 = Utilities::minus(P0, V0);
		double N = Utilities::perp(e, temp1);	// -dot(ne, P0-V0)
		return N / D;
	}
}


// Find the intersection of a line starting at P0 in the dirction dS and the plane with normal n and point V0
// Returns the number of multiples of dS that occur bitween Po and the intersection
double Utilities::LinePlaneIntersection(const std::vector<double>& P0, const std::vector<double>& dS,
										const std::vector<double>& V0, const std::vector<double>& n)
{
	double D = Utilities::dot(n, dS); // The amount of the dS vector normal to the plane
	if (fabs(D) < 1e-14) // dS is essential parallel to the plane
		return std::numeric_limits<double>::max(); // Hopefuly avoid a divide by zero situation
	else
	{
		std::vector<double> temp1 = Utilities::minus(V0, P0);
		double N = Utilities::dot(n, temp1);	// -dot(ne, P0-V0)
		return N / D;
	}
}




// Returns whether or not the given point P is inside of the projected triangle defined by a, b, c
// Algorithms from: http://blackpawn.com/texts/pointinpoly/
bool Utilities::pointInTriangleSS(const std::vector<double>& P, const std::vector<double>& a,
								  const std::vector<double>& b, const std::vector<double>& c)
{
	if (sameSide(P, a, b, c) && sameSide(P, b, a, c) && sameSide(P, c, a, b))
		return true;
	else
		return false;
}
bool Utilities::pointInTriangleBary(const std::vector<double>& P, const std::vector<double>& a,
									const std::vector<double>& b, const std::vector<double>& c)
{
	std::vector<double> v0 = Utilities::minus(c, a);
	std::vector<double> v1 = Utilities::minus(b, a);
	std::vector<double> v2 = Utilities::minus(P, a);

	double dot00 = Utilities::dot(v0, v0);
	double dot01 = Utilities::dot(v0, v1);
	double dot02 = Utilities::dot(v0, v2);
	double dot11 = Utilities::dot(v1, v1);
	double dot12 = Utilities::dot(v1, v2);

	// Compute barycentric coordinates
	double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
	double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	if (u >= 0)
	{
		double v = (dot00 * dot12 - dot01 * dot02) * invDenom;
		return (v >= 0.0) && (u + v < 1.0);
	}
	else
		return false;
}
// Helper function for pointInTriangleSS. Deermines is P1 and P2 are in the same half plane defined by AB
bool Utilities::sameSide(const std::vector<double>& P1, const std::vector<double>& P2,
						 const std::vector<double>& a, const std::vector<double>& b)
{
	std::vector<double> AB = Utilities::minus(b, a);
	std::vector<double> AP1 = Utilities::minus(P1, a);
	std::vector<double> AP2 = Utilities::minus(P2, a);
	std::vector<double> cp1 = Utilities::cross(AB, AP1);
	std::vector<double> cp2 = Utilities::cross(AB, AP2);
	if (Utilities::dot(cp1, cp2) >= 0)
		return true;
	else
		return false;
}



// Splits a string into pieces using the delimiter specified
std::vector<std::string> Utilities::splitString(std::string obj, std::string del)
{
	std::vector<std::string> res;
	char * pch;
	pch = strtok (&obj[0], del.data());
	while (pch != NULL)
	{
		std::string token(pch);
		res.push_back(token);
		pch = strtok (NULL, del.data());
	}
	return res;
}
