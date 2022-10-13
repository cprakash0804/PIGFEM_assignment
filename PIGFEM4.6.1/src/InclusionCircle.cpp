/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#include "InclusionCircle.h"
#include "SensitivityShapeParameterStructural.h"
#include "SensitivityShapeParameterThermal.h"
#include <limits> // std::numeric_limits
#include <cmath> // sqrt



Circle_Inclusion::Circle_Inclusion()
	: _r(0.0)
{
}
Circle_Inclusion::Circle_Inclusion(Circle_Inclusion & other_inc)
{
	copy(&other_inc);
}
Circle_Inclusion::~Circle_Inclusion()
{
	_center.clear();
}



Inclusion* Circle_Inclusion::allocate_and_copy()
{
	Inclusion* new_inc = allocate();   // Hopefully this works to call the appropriate allocate function based on what instance calls this function
	copy(new_inc);
	return new_inc;
}
Inclusion* Circle_Inclusion::allocate()
{
	return new Circle_Inclusion;
}
// This function assumes the current material is set up correctly. Which it should't be possible for it not to be
void Circle_Inclusion::copy(Inclusion* other_inc)
{
	other_inc->set_material(_mat);
	other_inc->set_id(_id);
	other_inc->set_vec_parameter("Center", _center);
	other_inc->set_parameter("r", _r);
}



void Circle_Inclusion::set_parameter(std::string name, double val)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if (name=="X")
	{
		_center.resize(2);
		_center[0] = val;
	}
	else if (name=="Y")
	{
		_center.resize(2);
		_center[1] = val;
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
double Circle_Inclusion::get_parameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if (name=="X")
	{
		if (_center.size()==2)
			return _center[0];
		else
			err_message("Invalid parameter call in a Circle");
	}
	else if (name=="Y")
	{
		if (_center.size()==2)
			return _center[1];
		else
			err_message("Invalid parameter call in a Circle");
	}
	else if (name=="R" || name=="RC" || name=="RADIUS")
		return _r;
	else
	{
		char buf[50];
		sprintf(buf, "%s is not a valid paramter name.", name.c_str());
		err_message( buf );
	}
}




void Circle_Inclusion::set_vec_parameter(std::string name, std::vector<double>& val)
{
	if(val.size()>2 || val.size()==0)
	{
		char buf [100];
		sprintf(buf, "The number of values in the vector must be less than 3 and greater than 0 for the %s parameter.", name.c_str());
		err_message( buf );
	}
	if(val.size() != 2)
		val.resize(2);

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
std::vector<double> Circle_Inclusion::get_vec_parameter(std::string name)
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



std::vector<double> Circle_Inclusion::get_domain_lims()
{
	double min = -1.0*std::numeric_limits<double>::max();
	double max = std::numeric_limits<double>::max();
	std::vector<double> ret = {_center[0]-_r, _center[1]-_r, min, _center[0]+_r, _center[1]+_r, max};
	return ret;
}

id_type Circle_Inclusion::detect_node(Node* node)
{
	double x = (*node)(0) - _center[0];
	double y = (*node)(1) - _center[1];
	double val = (sqrt(x*x + y*y) - _r) / _r;

	double eps = 1e-12;
	if(val > eps)
		return 0;
	else if(val < -1.0*eps)
		return 2;
	else
		return 1;
}

// Don't feel like acually figuring out how t get the surface normal so we'll just return a vector from he center to the node. Should be close enough
std::vector<double> Circle_Inclusion::get_surface_normal(Node* node)
{
	std::vector<double> ret = {(*node)(0)-_center[0], (*node)(1)-_center[1], 0.0};
	double norm = sqrt(ret[0]*ret[0] + ret[1]*ret[1]);
	ret[0] = ret[0]/norm;
	ret[1] = ret[1]/norm;
	return ret;
}


std::vector<double> Circle_Inclusion::getIntersectionVelocity(std::vector<Node*>& nodes, std::vector<double> intersect, std::string param_name)
{
	std::transform(param_name.begin(), param_name.end(), param_name.begin(), ::toupper); // Capitilize the name
	double al = (*nodes[0])(0) - (*nodes[1])(0);
	double bl = (*nodes[0])(1) - (*nodes[1])(1);
	double x = intersect[0] - _center[0];
	double y = intersect[1] - _center[1];

	double val = 1.0 / (al*x + bl*y);
	if (param_name=="X" || param_name=="XC" || param_name=="X_CENTER" || param_name=="X CENTER")
		val *= x;
	else if (param_name=="Y" || param_name=="YC" || param_name=="Y_CENTER" || param_name=="Y CENTER")
		val *= y;
	else if (param_name=="R" || param_name=="RC" || param_name=="RADIUS")
		val *= _r;
	else
		err_message("Invalid shape parameter for velocity computation.");

	std::vector<double> ret = {al*val, bl*val};
	return ret;
}



SensitivityShapeParameter* Circle_Inclusion::getSensitivityParameter(std::string param_name, classification type)
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
	else if (param_name=="R" || param_name=="RC" || param_name=="RADIUS")
		ret->set_param_name("R");
	else
	{
		delete ret;
		err_message("Please input a valid parameter name.");
	}

	return ret;
}


std::vector<std::string> Circle_Inclusion::getSensitivityParameterNames()
{
	std::vector<std::string> names = {"X", "Y", "R"};
	return names;
}














double Circle_Inclusion::int_f(double alpha, void *params)
{
	struct int_params *p = (struct int_params*) params;
	std::vector<Node> nodes = p->nodes;
	double gcoords[2] = {nodes[0](0)+alpha*(nodes[1](0)-nodes[0](0)),
						 nodes[0](1)+alpha*(nodes[1](1)-nodes[0](1))};
	gcoords[0] -= _center[0];
	gcoords[1] -= _center[1];

	double ret = pow(gcoords[0]/_r, 2) + pow(gcoords[1]/_r, 2) - 1;
	return ret;
}
double Circle_Inclusion::int_df(double alpha, void *params)
{
	struct int_params *p = (struct int_params*) params;
	std::vector<Node> nodes = p->nodes;
	double gcoords[2] = {nodes[0](0)+alpha*(nodes[1](0)-nodes[0](0)),
						 nodes[0](1)+alpha*(nodes[1](1)-nodes[0](1))};
	gcoords[0] -= _center[0];
	gcoords[1] -= _center[1];
	double d_dalpha[2] = {nodes[1](0)-nodes[0](0),
						  nodes[1](1)-nodes[0](1)};

	double ret = (1/(_r*_r)) * 2.0 * (gcoords[0] * d_dalpha[0] + gcoords[1] * d_dalpha[1]);
	return ret;
}
void Circle_Inclusion::int_fdf(double alpha, void *params, double *f, double *df)
{
	struct int_params *p = (struct int_params*) params;
	std::vector<Node> nodes = p->nodes;
	double gcoords[2] = {nodes[0](0)+alpha*(nodes[1](0)-nodes[0](0)),
						 nodes[0](1)+alpha*(nodes[1](1)-nodes[0](1))};
	gcoords[0] -= _center[0];
	gcoords[1] -= _center[1];
	double d_dalpha[2] = {nodes[1](0)-nodes[0](0),
						  nodes[1](1)-nodes[0](1)};
						  
	*f = pow(gcoords[0]/_r, 2) + pow(gcoords[1]/_r, 2) - 1;
	*df = (1/(_r*_r)) * 2.0 * (gcoords[0] * d_dalpha[0] + gcoords[1] * d_dalpha[1]);
}
