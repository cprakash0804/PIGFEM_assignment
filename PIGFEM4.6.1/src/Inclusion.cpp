/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "Inclusion.h"
#include "InclusionPlane.h"
#include "InclusionEllipse.h"
#include "InclusionEllipsoid.h"


Inclusion::Inclusion()
	: _id(0), _mat(NULL)
{
}
Inclusion::Inclusion(Inclusion & other_inc)
{
	other_inc.set_id(_id);
	other_inc.set_material(_mat);
}
Inclusion::~Inclusion()
{
}




std::vector<double> Inclusion::find_intersection(std::vector<Node>& nodes)
{
	if(nodes.size() != 2)
		err_message("Currently intersection finding with more than 2 nodes is not supported.");
	// Should check to make sure that the sign of the objective function is different at the boundaries
	struct int_params params = {nodes, this};
	if((int_f(0.0, &params)*int_f(1.0, &params)) >= 0)
	{
		char buf[150], buf2[150];
		strcpy(buf, "No valid intersecion to be found between nodes %");
		strcat(buf, SPEC);
		strcat(buf, " and %");
		strcat(buf, SPEC);
		sprintf(buf2, buf, nodes[0].get_id(), nodes[1].get_id());
		err_message( buf2 );
	}
	
	// Define some parameters
	int status;
	int iter = 0, max_iter = 25;
	const gsl_root_fdfsolver_type *T;
	gsl_root_fdfsolver *s;
	double alpha0, alpha = 0.5; // Initial guess
	
	gsl_function_fdf FDF;
	FDF.f = &int_f_wrapper;
	FDF.df = &int_df_wrapper;
	FDF.fdf = &int_fdf_wrapper;
	FDF.params = &params;

	T = gsl_root_fdfsolver_steffenson; // This is an accelerated Newton method
	s = gsl_root_fdfsolver_alloc (T);
	gsl_root_fdfsolver_set (s, &FDF, alpha);

	do
	{
		// Perform a newton-like iteration
		iter++;
		status = gsl_root_fdfsolver_iterate (s);

		// Check the change in the alpha value from that iteration
		alpha0 = alpha;
		alpha = gsl_root_fdfsolver_root (s);
		status = gsl_root_test_delta (alpha, alpha0, 1e-7, 0); // Perform the actual check
	}
	while (status == GSL_CONTINUE && iter < max_iter);

	// For some reason the Newton solution didnt coonverge to the correct value
	// To be sure to get an accurate value we'll just try a bisection search
	if(iter==max_iter || alpha<0.0 || alpha>1.0)
	{
		const gsl_root_fsolver_type *T2;
		gsl_root_fsolver *s2;

		//inclusion* ptr2 = this;
		//auto ptr = [=](double x)->double{return ptr2->};



		gsl_function F;
		F.function = &int_f_wrapper;
		F.params = &params;

		T2 = gsl_root_fsolver_brent; // This is an accelerated bracketing method
		s2 = gsl_root_fsolver_alloc (T2);
		double a_lo = 0.0, a_hi = 1.0;
		gsl_root_fsolver_set (s2, &F, a_lo, a_hi);

		iter = 0;
		max_iter = 100;
		do
		{
			// Perform a bracketing iteration
			iter++;
			gsl_root_fsolver_iterate (s2);

			// Check the bracketing bounds
			alpha = gsl_root_fsolver_root (s2);
			a_lo = gsl_root_fsolver_x_lower (s2);
			a_hi = gsl_root_fsolver_x_upper (s2);
			status = gsl_root_test_delta(a_lo, a_hi, 1e-7, 0);
		}
		while (status==GSL_CONTINUE && iter<max_iter);

		gsl_root_fsolver_free (s2);
	}


	gsl_root_fdfsolver_free (s);
	
	// Generate the vector of coordinates
	std::vector<double> ret;
	ret.push_back( nodes[0](0) + alpha*(nodes[1](0)-nodes[0](0)) );
	ret.push_back( nodes[0](1) + alpha*(nodes[1](1)-nodes[0](1)) );
	ret.push_back( nodes[0](2) + alpha*(nodes[1](2)-nodes[0](2)) );

	return ret;
}
