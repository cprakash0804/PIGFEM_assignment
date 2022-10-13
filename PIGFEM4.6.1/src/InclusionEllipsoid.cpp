/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#include "InclusionEllipsoid.h"
#include "SensitivityShapeParameterStructural.h"
#include "SensitivityShapeParameterThermal.h"
#include <limits> // std::numeric_limits
#include <cmath> // sqrt





Ellipsoid_Inclusion::Ellipsoid_Inclusion()
	: _alpha(0.0), _ca(1.0), _sa(0.0), _beta(0.0), _cb(1.0), _sb(0.0), _gamma(0.0), _cg(1.0), _sg(0.0), _a(0.0), _b(0.0), _c(0.0)
{
}
Ellipsoid_Inclusion::Ellipsoid_Inclusion(Ellipsoid_Inclusion & other_inc)
{
	copy(&other_inc);
}
Ellipsoid_Inclusion::~Ellipsoid_Inclusion()
{
	_center.clear();
}


Inclusion* Ellipsoid_Inclusion::allocate_and_copy()
{
	Inclusion* new_inc = allocate();   // Hopefully this works to call the appropriate allocate function based on what instance calls this function
	copy(new_inc);
	return new_inc;
}
Inclusion* Ellipsoid_Inclusion::allocate()
{
	return new Ellipsoid_Inclusion;
}
// This function assumes the current material is set up correctly. Which it should't be possible for it not to be
void Ellipsoid_Inclusion::copy(Inclusion* other_inc)
{
	other_inc->set_material(_mat);
	other_inc->set_id(_id);
	other_inc->set_vec_parameter("Center", _center);
	other_inc->set_parameter("a", _a);
	other_inc->set_parameter("b", _b);
	other_inc->set_parameter("c", _c);
	other_inc->set_parameter("alpha", _alpha);
	other_inc->set_parameter("beta", _beta);
	other_inc->set_parameter("gamma", _gamma);
}





void Ellipsoid_Inclusion::set_parameter(std::string name, double val)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if (name=="X")
	{
		_center.resize(3);
		_center[0] = val;
	}
	else if (name=="Y")
	{
		_center.resize(3);
		_center[1] = val;
	}
	else if (name=="Z")
	{
		_center.resize(3);
		_center[2] = val;
	}
	else if(name=="A" || name=="MAJOR")
		_a = val;
	else if(name=="B" || name=="SEMIMINOR" || name=="SEMI MINOR" || name=="SEMI_MINOR" || name=="SEMI-MINOR")
		_b = val;
	else if(name=="C" || name=="MINOR")
		_c = val;
	else if(name=="ALPHA" || name=="FIRST EULER ANGLE" || name=="FIRST_EULER_ANGLE")
	{
		_alpha = val;
		_ca = cos(-_alpha); // Negative signs are because when transforming to local coords we do the inverse rotations
		_sa = sin(-_alpha);
	}
	else if(name=="BETA" || name=="SECOND EULER ANGLE" || name=="SECOND_EULER_ANGLE")
	{
		_beta = val;
		_cb = cos(-_beta);
		_sb = sin(-_beta);
	}
	else if(name=="GAMMA" || name=="THIRD EULER ANGLE" || name=="THIRD_EULER_ANGLE")
	{
		_gamma = val;
		_cg = cos(-_gamma);
		_sg = sin(-_gamma);
	}
	else
	{
		char buf[50];
		sprintf(buf, "%s is not a valid paramter name.", name.c_str());
		err_message( buf );
	}
}
double Ellipsoid_Inclusion::get_parameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if (name=="X")
	{
		if (_center.size()==3)
			return _center[0];
		else
			err_message("Invalid parameter call in a Ellipsoid");
	}
	else if (name=="Y")
	{
		if (_center.size()==3)
			return _center[1];
		else
			err_message("Invalid parameter call in a Ellipsoid");
	}
	else if (name=="Z")
	{
		if (_center.size()==3)
			return _center[2];
		else
			err_message("Invalid parameter call in a Ellipsoid");
	}
	else if(name=="A" || name=="MAJOR")
		return _a;
	else if(name=="B" || name=="SEMIMINOR" || name=="SEMI MINOR" || name=="SEMI_MINOR" || name=="SEMI-MINOR")
		return _b;
	else if(name=="C" || name=="MINOR")
		return _c;
	else if(name=="ALPHA" || name=="FIRST EULER ANGLE" || name=="FIRST_EULER_ANGLE")
		return _alpha;
	else if(name=="BETA" || name=="SECOND EULER ANGLE" || name=="SECOND_EULER_ANGLE")
		return _beta;
	else if(name=="GAMMA" || name=="THIRD EULER ANGLE" || name=="THIRD_EULER_ANGLE")
		return _gamma;
	else
	{
		char buf[50];
		sprintf(buf, "%s is not a valid paramter name.", name.c_str());
		err_message( buf );
	}
}




void Ellipsoid_Inclusion::set_vec_parameter(std::string name, std::vector<double>& val)
{
	if(val.size()>3 || val.size()==0)
	{
		char buf [100];
		sprintf(buf, "The number of values in the vector must be less than 4 and greater than 0 for the %s parameter.", name.c_str());
		err_message( buf );
	}
	if(val.size() != 3)
		val.resize(3);

	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="CENTER" || name=="P" || name=="ORIGIN" || name=="POINT")
		_center = val;
	else
	{
		char buf[50];
		sprintf(buf, "%s is not a valid paramter name.", name.c_str());
		err_message( buf );
	}
}
std::vector<double> Ellipsoid_Inclusion::get_vec_parameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="CENTER" || name=="P" || name=="ORIGIN" || name=="POINT")
		return _center;
	else
	{
		char buf[50];
		sprintf(buf, "%s is not a valid paramter name.", name.c_str());
		err_message( buf );
	}
}



std::vector<double> Ellipsoid_Inclusion::get_domain_lims()
{
	double max = _a;
	if (_b > max)
		max = _b;
	if (_c > max)
		max = _c;
	std::vector<double> ret = {_center[0]-max, _center[1]-max, _center[2]-max, _center[0]+max, _center[1]+max, _center[1]+max};
	return ret;
}

id_type Ellipsoid_Inclusion::detect_node(Node* node)
{
	double global[3] = {(*node)(0), (*node)(1), (*node)(2)};
	double local[3];
	convert_to_local(global, local);

	double val = pow(local[0]/_a, 2) + pow(local[1]/_b, 2) + pow(local[2]/_c, 2) - 1.0;

	double eps = 1e-12;
	if(val > eps)
		return 0;
	else if(val < -1.0*eps)
		return 2;
	else
		return 1;
}

// Don't feel like acually figuring out how to get the surface normal so we'll just return a vector from he center to the node. Should be close enough
std::vector<double> Ellipsoid_Inclusion::get_surface_normal(Node* node)
{
	std::vector<double> ret = {(*node)(0)-_center[0], (*node)(1)-_center[1], (*node)(2)-_center[2]};
	double norm = sqrt(ret[0]*ret[0] + ret[1]*ret[1] + ret[2]*ret[2]);
	ret[0] = ret[0]/norm;
	ret[1] = ret[1]/norm;
	ret[2] = ret[2]/norm;
	return ret;
}




std::vector<double> Ellipsoid_Inclusion::getIntersectionVelocity(std::vector<Node*>& nodes, std::vector<double> intersect, std::string param_name)
{
	std::transform(param_name.begin(), param_name.end(), param_name.begin(), ::toupper); // Capitilize the name
	double al = (*nodes[0])(0) - (*nodes[1])(0);
	double bl = (*nodes[0])(1) - (*nodes[1])(1);
	double cl = (*nodes[0])(2) - (*nodes[1])(2);
	double global[3] = {intersect[0], intersect[1], intersect[2]};
	double local[3];
	convert_to_local(global, local); // local is the x,y,z primes in my notes
	/*global = {intersect[0] - _center[0],
			  intersect[1] - _center[1],
			  intersect[2] - _center[2]}; */
	global[0] = intersect[0] - _center[0];
	global[1] = intersect[1] - _center[1];
	global[2] = intersect[2] - _center[2];

	// Some temporary variables
	double m = (2.0*local[0]/(_a*_a))*(_cb*_cg) + (2.0*local[1]/(_b*_b))*(_sa*_sb*_cg+_ca*_sg) + (2.0*local[2]/(_c*_c))*(-_ca*_sb*_cg+_sa*_sg);
	double n = (2.0*local[0]/(_a*_a))*(_cb*_sg) + (2.0*local[1]/(_b*_b))*(-_sa*_sb*_sg+_ca*_cg) + (2.0*local[2]/(_c*_c))*(_ca*_sb*_sg+_sa*_cg);
	double o = (2.0*local[0]/(_a*_a))*(_sb) + (2.0*local[1]/(_b*_b))*(-_sa*_cb) + (2.0*local[2]/(_c*_c))*(_ca*_cb);

	double val = 1.0 / (al*m + bl*n + cl*o);
	if (param_name=="X" || param_name=="XC" || param_name=="X_CENTER" || param_name=="X CENTER")
		val *= (2.0*local[0]/(_a*_a))*(_cb*_cg);
	else if (param_name=="Y" || param_name=="YC" || param_name=="Y_CENTER" || param_name=="Y CENTER")
		val *= (2.0*local[1]/(_b*_b))*(-_sa*_sb*_sg+_ca*_cg);
	else if (param_name=="Z" || param_name=="ZC" || param_name=="Z_CENTER" || param_name=="Z CENTER")
		val *= (2.0*local[2]/(_c*_c))*(_ca*_cb);
	else if(param_name=="A" || param_name=="MAJOR")
		val *= (2.0*pow(local[0],2)/pow(_a,3));
	else if(param_name=="B" || param_name=="SEMIMINOR" || param_name=="SEMI MINOR" || param_name=="SEMI_MINOR" || param_name=="SEMI-MINOR")
		val *= (2.0*pow(local[1],2)/pow(_b,3));
	else if(param_name=="C" || param_name=="MINOR")
		val *= (2.0*pow(local[2],2)/pow(_b,3));
	else if(param_name=="ALPHA" || param_name=="FIRST EULER ANGLE" || param_name=="FIRST_EULER_ANGLE")
		val *= ((2.0*local[0]/(_a*_a)) * ( 0 ) +
				(2.0*local[1]/(_b*_b)) * ( (-_ca*_sb*_cg+_sa*_sg)*global[0] + (_ca*_sb*_sg+_sa*_cg)*global[1] + (_ca*_cb)*global[2]) + 
				(2.0*local[2]/(_c*_c)) * ( (-_sa*_sb*_cg-_ca*_sg)*global[0] + (_sa*_sb*_sg-_ca*_cg)*global[1] + (_sa*_cb)*global[2]) );
	else if(param_name=="BETA" || param_name=="SECOND EULER ANGLE" || param_name=="SECOND_EULER_ANGLE")
		val *= ((2.0*local[0]/(_a*_a)) * ( (_sb*_cg)*global[0] - (_sb*_sg)*global[1] - _cb*global[2]) +
				(2.0*local[1]/(_b*_b)) * ( (-_sa*_cb*_cg)*global[0] + (_sa*_cb*_sg)*global[1] + (-_sa*_sb)*global[2]) + 
				(2.0*local[2]/(_c*_c)) * ( (_ca*_cb*_cg)*global[0] + (-_ca*_cb*_sg)*global[1] + (_ca*_sb)*global[2]) );
	else if(param_name=="GAMMA" || param_name=="THIRD EULER ANGLE" || param_name=="THIRD_EULER_ANGLE")
		val *= ((2.0*local[0]/(_a*_a)) * ( (_cb*_sg)*global[0] - (-_cb*_cg)*global[1] ) +
				(2.0*local[1]/(_b*_b)) * ( (_sa*_sb*_sg-_ca*_cg)*global[0] + (_sa*_sb*_cg+_ca*_sg)*global[1]) + 
				(2.0*local[2]/(_c*_c)) * ( (-_ca*_sb*_sg-_sa*_cg)*global[0] + (-_ca*_sb*_cg+_sa*_sg)*global[1]) );
	else
		err_message("Invalid shape parameter for velocity computation.");

	std::vector<double> ret = {al*val, bl*val, cl*val};
	return ret;
}


SensitivityShapeParameter* Ellipsoid_Inclusion::getSensitivityParameter(std::string param_name, classification type)
{
	SensitivityShapeParameter* ret;
	if (type==STRUCTURAL)
		ret = new SensitivityShapeParameterStructural;
	else if (type==THERMAL)
		ret = new SensitivityShapeParameterThermal;
	else
		err_message("Invalid Problem classification for a shape sensitivity parameter");

	ret->set_inc_id(_id);

	std::transform(param_name.begin(), param_name.end(), param_name.begin(), ::toupper); // Capitilize the name

	if (param_name=="X" || param_name=="XC" || param_name=="X_CENTER" || param_name=="X CENTER")
		ret->set_param_name("X");
	else if (param_name=="Y" || param_name=="YC" || param_name=="Y_CENTER" || param_name=="Y CENTER")
		ret->set_param_name("Y");
	else if (param_name=="Z" || param_name=="ZC" || param_name=="Z_CENTER" || param_name=="Z CENTER")
		ret->set_param_name("Z");
	else if(param_name=="A" || param_name=="MAJOR")
		ret->set_param_name("A");
	else if(param_name=="B" || param_name=="SEMIMINOR" || param_name=="SEMI MINOR" || param_name=="SEMI_MINOR" || param_name=="SEMI-MINOR")
		ret->set_param_name("B");
	else if(param_name=="C" || param_name=="MINOR")
		ret->set_param_name("C");
	else if(param_name=="ALPHA" || param_name=="FIRST EULER ANGLE" || param_name=="FIRST_EULER_ANGLE")
		ret->set_param_name("ALPHA");
	else if(param_name=="BETA" || param_name=="SECOND EULER ANGLE" || param_name=="SECOND_EULER_ANGLE")
		ret->set_param_name("BETA");
	else if(param_name=="GAMMA" || param_name=="THIRD EULER ANGLE" || param_name=="THIRD_EULER_ANGLE")
		ret->set_param_name("GAMMA");
	else
	{
		delete ret;
		err_message("Please input a valid parameter name.");
	}

	return ret;
}

std::vector<std::string> Ellipsoid_Inclusion::getSensitivityParameterNames()
{
	std::vector<std::string> names = {"X", "Y", "Z", "A", "B", "C", "ALPHA", "BETA", "GAMMA"};
	return names;
}



// NOTE: THIS ALPHA IS NOT THE SAME AS THE ANGLE ALPHA, IT IS A LINE SEARCH PARAMETER
// I should probably just change one of them...
double Ellipsoid_Inclusion::int_f(double alpha, void *params)
{
	struct int_params *p = (struct int_params*) params;
	std::vector<Node> nodes = p->nodes;
	double global[3] = {nodes[0](0)+alpha*(nodes[1](0)-nodes[0](0)),
						nodes[0](1)+alpha*(nodes[1](1)-nodes[0](1)),
						nodes[0](2)+alpha*(nodes[1](2)-nodes[0](2))};
	double local[3];

	// Convert from the global coords to the local coords
	convert_to_local(global, local);

	double ret = pow(local[0]/_a, 2) + pow(local[1]/_b, 2) + pow(local[2]/_c, 2) - 1.0;
	return ret;
}
double Ellipsoid_Inclusion::int_df(double alpha, void *params)
{
	struct int_params *p = (struct int_params*) params;
	std::vector<Node> nodes = p->nodes;
	double global[3] = {nodes[0](0)+alpha*(nodes[1](0)-nodes[0](0)),
						nodes[0](1)+alpha*(nodes[1](1)-nodes[0](1)),
						nodes[0](2)+alpha*(nodes[1](2)-nodes[0](2))};
	double local[3];

	// Convert from the global coords to the local coords
	convert_to_local(global, local);

	// Compute the partial derivative of the coords wrt alpha
	double dxg_dalpha[3] = {nodes[1](0)-nodes[0](0),
							nodes[1](1)-nodes[0](1),
							nodes[1](2)-nodes[0](2)};
	double dxl_dalpha[3] = {(_cb*_cg)*dxg_dalpha[0] - (_cb*_sg)*dxg_dalpha[1] + _sb*dxg_dalpha[2],
							(_sa*_sb*_cg+_ca*_sg)*dxg_dalpha[0] + (-_sa*_sb*_sg+_ca*_cg)*dxg_dalpha[1] + (-_sa*_cb)*dxg_dalpha[2],
							(-_ca*_sb*_cg+_sa*_sg)*dxg_dalpha[0] + (_ca*_sb*_sg+_sa*_cg)*dxg_dalpha[1] + (_ca*_cb)*dxg_dalpha[2]};

	double ret = (1/(_a*_a))*2.0*local[0]*dxl_dalpha[0] +
				 (1/(_b*_b))*2.0*local[1]*dxl_dalpha[1] +
				 (1/(_c*_c))*2.0*local[2]*dxl_dalpha[2];
	return ret;
}
void Ellipsoid_Inclusion::int_fdf(double alpha, void *params, double *f, double *df)
{
	struct int_params *p = (struct int_params*) params;
	std::vector<Node> nodes = p->nodes;
	double global[3] = {nodes[0](0)+alpha*(nodes[1](0)-nodes[0](0)),
						nodes[0](1)+alpha*(nodes[1](1)-nodes[0](1)),
						nodes[0](2)+alpha*(nodes[1](2)-nodes[0](2))};
	double local[3];

	// Convert from the global coords to the local coords
	convert_to_local(global, local);

	// Evaluate the objective
	*f = pow(local[0]/_a, 2) + pow(local[1]/_b, 2) + pow(local[2]/_c, 2) - 1.0;

	// Compute the partial derivative of the global coords wrt alpha
	double dxg_dalpha[3] = {nodes[1](0)-nodes[0](0),
							nodes[1](1)-nodes[0](1),
							nodes[1](2)-nodes[0](2)};
	double dxl_dalpha[3] = {(_cb*_cg)*dxg_dalpha[0] - (_cb*_sg)*dxg_dalpha[1] + _sb*dxg_dalpha[2],
							(_sa*_sb*_cg+_ca*_sg)*dxg_dalpha[0] + (-_sa*_sb*_sg+_ca*_cg)*dxg_dalpha[1] + (-_sa*_cb)*dxg_dalpha[2],
							(-_ca*_sb*_cg+_sa*_sg)*dxg_dalpha[0] + (_ca*_sb*_sg+_sa*_cg)*dxg_dalpha[1] + (_ca*_cb)*dxg_dalpha[2]};

	*df = (1/(_a*_a))*2.0*local[0]*dxl_dalpha[0] +
		  (1/(_b*_b))*2.0*local[1]*dxl_dalpha[1] +
		  (1/(_c*_c))*2.0*local[2]*dxl_dalpha[2];
}

void Ellipsoid_Inclusion::convert_to_local(double* global, double* local)
{
	// Subtract off the center coordinates
	global[0] -= _center[0];
	global[1] -= _center[1];
	global[2] -= _center[2];

	// Perform the rotations in reverse order
	local[0] = (_cb*_cg)*global[0] - (_cb*_sg)*global[1] + _sb*global[2];
	local[1] = (_sa*_sb*_cg+_ca*_sg)*global[0] + (-_sa*_sb*_sg+_ca*_cg)*global[1] + (-_sa*_cb)*global[2];
	local[2] = (-_ca*_sb*_cg+_sa*_sg)*global[0] + (_ca*_sb*_sg+_sa*_cg)*global[1] + (_ca*_cb)*global[2];
}
