/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#include "InclusionSphere.h"
#include "SensitivityShapeParameterStructural.h"
#include "SensitivityShapeParameterThermal.h"
#include <limits> // std::numeric_limits
#include <cmath> // sqrt





Sphere_Inclusion::Sphere_Inclusion()
	: _r(0.0)
{
}
Sphere_Inclusion::Sphere_Inclusion(Sphere_Inclusion & other_inc)
{
	copy(&other_inc);
}
Sphere_Inclusion::~Sphere_Inclusion()
{
	_center.clear();
}


Inclusion* Sphere_Inclusion::allocate_and_copy()
{
	Inclusion* new_inc = allocate();   // Hopefully this works to call the appropriate allocate function based on what instance calls this function
	copy(new_inc);
	return new_inc;
}
Inclusion* Sphere_Inclusion::allocate()
{
	return new Sphere_Inclusion;
}
// This function assumes the current material is set up correctly. Which it should't be possible for it not to be
void Sphere_Inclusion::copy(Inclusion* other_inc)
{
	other_inc->set_material(_mat);
	other_inc->set_id(_id);
	other_inc->set_vec_parameter("Center", _center);
	other_inc->set_parameter("r", _r);
}





void Sphere_Inclusion::set_parameter(std::string name, double val)
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
	else if(name=="R" || name=="RADIUS")
		_r = val;
	else
	{
		char buf[50];
		sprintf(buf, "%s is not a valid paramter name.", name.c_str());
		err_message( buf );
	}
}
double Sphere_Inclusion::get_parameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if (name=="X")
	{
		if (_center.size()==3)
			return _center[0];
		else
			err_message("Invalid parameter call in a Sphere");
	}
	else if (name=="Y")
	{
		if (_center.size()==3)
			return _center[1];
		else
			err_message("Invalid parameter call in a Sphere");
	}
	else if (name=="Z")
	{
		if (_center.size()==3)
			return _center[2];
		else
			err_message("Invalid parameter call in a Sphere");
	}
	else if(name=="R" || name=="RADIUS")
		return _r;
	else
	{
		char buf[50];
		sprintf(buf, "%s is not a valid paramter name.", name.c_str());
		err_message( buf );
	}
}




void Sphere_Inclusion::set_vec_parameter(std::string name, std::vector<double>& val)
{
	if(val.size()>3 || val.size()==0)
	{
		char buf [100];
		sprintf(buf, "The number of values in the vector must be less than 3 and greater than 0 for the %s parameter.", name.c_str());
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
std::vector<double> Sphere_Inclusion::get_vec_parameter(std::string name)
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



std::vector<double> Sphere_Inclusion::get_domain_lims()
{
	std::vector<double> ret = {_center[0]-_r, _center[1]-_r, _center[2]-_r, _center[0]+_r, _center[1]+_r, _center[1]+_r};
	return ret;
}

id_type Sphere_Inclusion::detect_node(Node* node)
{
	double global[3] = {(*node)(0) - _center[0],
						(*node)(1) - _center[1],
						(*node)(2) - _center[2]};

	double val = pow(global[0]/_r, 2) + pow(global[1]/_r, 2) + pow(global[2]/_r, 2) - 1.0;

	double eps = 1e-12;
	if(val > eps)
		return 0;
	else if(val < -1.0*eps)
		return 2;
	else
		return 1;
}

// Don't feel like acually figuring out how to get the surface normal so we'll just return a vector from he center to the node. Should be close enough
std::vector<double> Sphere_Inclusion::get_surface_normal(Node* node)
{
	std::vector<double> ret = {(*node)(0)-_center[0], (*node)(1)-_center[1], (*node)(2)-_center[2]};
	double norm = sqrt(ret[0]*ret[0] + ret[1]*ret[1] + ret[2]*ret[2]);
	ret[0] = ret[0]/norm;
	ret[1] = ret[1]/norm;
	ret[2] = ret[2]/norm;
	return ret;
}

std::vector<double> Sphere_Inclusion::getIntersectionVelocity(std::vector<Node*>& nodes, std::vector<double> intersect, std::string param_name)
{
	std::transform(param_name.begin(), param_name.end(), param_name.begin(), ::toupper); // Capitilize the name
	double al = (*nodes[1])(0) - (*nodes[0])(0);
	double bl = (*nodes[1])(1) - (*nodes[0])(1);
	double cl = (*nodes[1])(2) - (*nodes[0])(2);
	double x_bar = intersect[0] - _center[0];
	double y_bar = intersect[1] - _center[1];
	double z_bar = intersect[2] - _center[2];

	double val = 0.5 / (al*x_bar + bl*y_bar + cl*z_bar);
	if (param_name=="X" || param_name=="XC" || param_name=="XS")
		val *= (2. * x_bar);
	else if (param_name=="Y" || param_name=="YC" || param_name=="YS")
		val *= (2. * y_bar);
	else if (param_name=="Z" || param_name=="ZC" || param_name=="ZS")
		val *= (2. * z_bar);
	else if (param_name=="R" || param_name=="RC" || param_name=="RS" || param_name=="RADIUS")
		val *= (2. * _r);
	else
		err_message("Invalid shape parameter for velocity computation.");

	std::vector<double> ret = {al*val, bl*val, cl*val};
	return ret;
}




SensitivityShapeParameter* Sphere_Inclusion::getSensitivityParameter(std::string param_name, classification type)
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
	else if (param_name=="R" || param_name=="RC" || param_name=="RS" || param_name=="RADIUS")
		ret->set_param_name("R");
	else
	{
		delete ret;
		err_message("Please input a valid parameter name.");
	}

	return ret;
}

std::vector<std::string> Sphere_Inclusion::getSensitivityParameterNames()
{
	std::vector<std::string> names = {"X", "Y", "Z", "R"};
	return names;
}












// NOTE: THIS ALPHA IS NOT THE SAME AS THE ANGLE ALPHA, IT IS A LINE SEARCH PARAMETER
// I should probably just change one of them...
double Sphere_Inclusion::int_f(double alpha, void *params)
{
	struct int_params *p = (struct int_params*) params;
	std::vector<Node> nodes = p->nodes;
	double global[3] = {nodes[0](0)+alpha*(nodes[1](0)-nodes[0](0)),
						nodes[0](1)+alpha*(nodes[1](1)-nodes[0](1)),
						nodes[0](2)+alpha*(nodes[1](2)-nodes[0](2))};
	// Subtract off the center coordinates
	global[0] -= _center[0];
	global[1] -= _center[1];
	global[2] -= _center[2];

	double ret = pow(global[0]/_r,2) + pow(global[1]/_r,2) + pow(global[2]/_r, 2) - 1.0;
	return ret;
}
double Sphere_Inclusion::int_df(double alpha, void *params)
{
	struct int_params *p = (struct int_params*) params;
	std::vector<Node> nodes = p->nodes;
	double global[3] = {nodes[0](0)+alpha*(nodes[1](0)-nodes[0](0)),
						nodes[0](1)+alpha*(nodes[1](1)-nodes[0](1)),
						nodes[0](2)+alpha*(nodes[1](2)-nodes[0](2))};
	// Subtract off the center coordinates
	global[0] -= _center[0];
	global[1] -= _center[1];
	global[2] -= _center[2];

	// Compute the partial derivative of the coords wrt alpha
	double dx_dalpha[3] = {nodes[1](0)-nodes[0](0),
							nodes[1](1)-nodes[0](1),
							nodes[1](2)-nodes[0](2)};

	double ret = (1.0/(_r*_r)) * 2.0 * (global[0]*dx_dalpha[0] + global[1]*dx_dalpha[1] + global[2]*dx_dalpha[2]);
	return ret;
}
void Sphere_Inclusion::int_fdf(double alpha, void *params, double *f, double *df)
{
	struct int_params *p = (struct int_params*) params;
	std::vector<Node> nodes = p->nodes;
	double global[3] = {nodes[0](0)+alpha*(nodes[1](0)-nodes[0](0)),
						nodes[0](1)+alpha*(nodes[1](1)-nodes[0](1)),
						nodes[0](2)+alpha*(nodes[1](2)-nodes[0](2))};
	// Subtract off the center coordinates
	global[0] -= _center[0];
	global[1] -= _center[1];
	global[2] -= _center[2];
	
	// Evaluate the objective
	*f = pow(global[0]/_r,2) + pow(global[1]/_r,2) + pow(global[2]/_r, 2) - 1.0;

	// Compute the partial derivative of the global coords wrt alpha
	double dx_dalpha[3] = {nodes[1](0)-nodes[0](0),
						   nodes[1](1)-nodes[0](1),
						   nodes[1](2)-nodes[0](2)};

	*df = (1.0/(_r*_r)) * 2.0 * (global[0]*dx_dalpha[0] + global[1]*dx_dalpha[1] + global[2]*dx_dalpha[2]);
}
