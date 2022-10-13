/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "InclusionEllipse.h"
#include "SensitivityShapeParameterStructural.h"
#include "SensitivityShapeParameterThermal.h"
#include <limits> // std::numeric_limits
#include <cmath> // sqrt



Ellipse_Inclusion::Ellipse_Inclusion()
	: _alpha(0.0), _ca(1.0), _sa(0.0), _a(0.0), _b(0.0)
{
}
Ellipse_Inclusion::Ellipse_Inclusion(Ellipse_Inclusion & other_inc)
{
	copy(&other_inc);
}
Ellipse_Inclusion::~Ellipse_Inclusion()
{
	_center.clear();
}



Inclusion* Ellipse_Inclusion::allocate_and_copy()
{
	Inclusion* new_inc = allocate();   // Hopefully this works to call the appropriate allocate function based on what instance calls this function
	copy(new_inc);
	return new_inc;
}
Inclusion* Ellipse_Inclusion::allocate()
{
	return new Ellipse_Inclusion;
}
// This function assumes the current material is set up correctly. Which it should't be possible for it not to be
void Ellipse_Inclusion::copy(Inclusion* other_inc)
{
	other_inc->set_material(_mat);
	other_inc->set_id(_id);
	other_inc->set_vec_parameter("Center", _center);
	other_inc->set_parameter("a", _a);
	other_inc->set_parameter("b", _b);
	other_inc->set_parameter("alpha", _alpha);
}



void Ellipse_Inclusion::set_parameter(std::string name, double val)
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
	else if(name=="A" || name=="MAJOR")
		_a = val;
	else if(name=="B" || name=="MINOR")
		_b = val;
	else if(name=="ALPHA" || name=="ANGLE" || name=="AXIS")
	{
		_alpha = val;
		_ca = cos(_alpha);
		_sa = sin(_alpha);
	}
	else
	{
		char buf[50];
		sprintf(buf, "%s is not a valid paramter name.", name.c_str());
		err_message( buf );
	}
}
double Ellipse_Inclusion::get_parameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if (name=="X")
	{
		if (_center.size()==2)
			return _center[0];
		else
			err_message("Invalid parameter call in a Ellipse");
	}
	else if (name=="Y")
	{
		if (_center.size()==2)
			return _center[1];
		else
			err_message("Invalid parameter call in a Ellipse");
	}
	else if(name=="A" || name=="MAJOR")
		return _a;
	else if(name=="B" || name=="MINOR")
		return _b;
	else if(name=="ALPHA" || name=="ANGLE" || name=="AXIS")
		return _alpha;
	else
	{
		char buf[50];
		sprintf(buf, "%s is not a valid paramter name.", name.c_str());
		err_message( buf );
	}
}




void Ellipse_Inclusion::set_vec_parameter(std::string name, std::vector<double>& val)
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
	if(name=="MAJOR" || name=="MAJOR AXIS" || name=="MAJOR_AXIS" || name=="AXIS")
	{
		_alpha = atan2(val[1], val[0]);
		_ca = cos(_alpha);
		_sa = sin(_alpha);
	}
	else if(name=="CENTER" || name=="P" || name=="ORIGIN" || name=="POINT")
		_center = val;
	else
	{
		char buf[50];
		sprintf(buf, "%s is not a valid paramter name.", name.c_str());
		err_message( buf );
	}
}
std::vector<double> Ellipse_Inclusion::get_vec_parameter(std::string name)
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



std::vector<double> Ellipse_Inclusion::get_domain_lims()
{
	double min = -1.0*std::numeric_limits<double>::min();
	double max = std::numeric_limits<double>::max();
	double max_rad = _a;
	if (_b > max_rad)
		max_rad = _b;;
	std::vector<double> ret = {_center[0]-max_rad, _center[1]-max_rad, min, _center[0]+max_rad, _center[1]+max_rad, max};
	return ret;
}

id_type Ellipse_Inclusion::detect_node(Node* node)
{
	double gcoords[2] = {(*node)(0), (*node)(1)};
	gcoords[0] -= _center[0];
	gcoords[1] -= _center[1];
	double lcoords[2] = {gcoords[0]*_ca+gcoords[1]*_sa,
						 gcoords[0]*_sa-gcoords[1]*_ca};

	double val = pow(lcoords[0]/_a, 2) + pow(lcoords[1]/_b, 2) - 1;

	double eps = 1e-12;
	if(val > eps)
		return 0;
	else if(val < -1.0*eps)
		return 2;
	else
		return 1;
}

// Don't feel like acually figuring out how t get the surface normal so we'll just return a vector from he center to the node. Should be close enough
std::vector<double> Ellipse_Inclusion::get_surface_normal(Node* node)
{
	std::vector<double> ret = {(*node)(0)-_center[0], (*node)(1)-_center[1], 0.0};
	double norm = sqrt(ret[0]*ret[0] + ret[1]*ret[1]);
	ret[0] = ret[0]/norm;
	ret[1] = ret[1]/norm;
	return ret;
}




std::vector<double> Ellipse_Inclusion::getIntersectionVelocity(std::vector<Node*>& nodes, std::vector<double> intersect, std::string param_name)
{
	std::transform(param_name.begin(), param_name.end(), param_name.begin(), ::toupper); // Capitilize the name
	double A, B, C;
	A = pow(_ca/_a, 2) + pow(_sa/_b, 2);
	B = 2.0*_ca*_sa*(1.0/(_a*_a) - 1.0/(_b*_b));
	C = pow(_sa/_a, 2) + pow(_ca/_b, 2);

	double al = (*nodes[0])(0) - (*nodes[1])(0);
	double bl = (*nodes[0])(1) - (*nodes[1])(1);
	double x = intersect[0] - _center[0];
	double y = intersect[1] - _center[1];

	double val = 1.0 / (al*(2.0*A*x+B*y) + bl*(2.0*C*y+B*x));
	if (param_name=="X" || param_name=="XC")
		val *= (2.0*A*x + B*y);
	else if (param_name=="Y" || param_name=="YC")
		val *= (2.0*C*y + B*x);
	else if (param_name=="A" || param_name=="MAJOR")
	{
		double dA = (-2.0/_a)*pow(_ca/_a, 2);
		double dB = -4.0*_ca*_sa/pow(_a, 3);
		double dC = (-2.0/_a)*pow(_sa/_a, 2);
		val *= (-dA*x*x - dB*x*y - dC*y*y);
	}
	else if (param_name=="B" || param_name=="MINOR")
	{
		double dA = (-2.0/_b)*pow(_sa/_b, 2);
		double dB = 4.0*_ca*_sa/pow(_b, 3);
		double dC = (-2.0/_b)*pow(_ca/_b, 2);
		val *= (-dA*x*x - dB*x*y - dC*y*y);
	}
	else if (param_name=="ALPHA" || param_name=="ANGLE" || param_name=="AXIS")
	{
		double dA = -B;
		double dB = 2.0*(A-C);
		double dC = B;
		val *= (-dA*x*x - dB*x*y - dC*y*y);
	}
	else
		err_message("Invalid shape parameter for velocity computation.");

	std::vector<double> ret = {al*val, bl*val};
	return ret;
}

SensitivityShapeParameter* Ellipse_Inclusion::getSensitivityParameter(std::string param_name, classification type)
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
	else if (param_name=="A" || param_name=="MAJOR")
		ret->set_param_name("A");
	else if (param_name=="B" || param_name=="MINOR")
		ret->set_param_name("B");
	else if (param_name=="ALPHA" || param_name=="ANGLE" || param_name=="AXIS")
		ret->set_param_name("ALPHA");
	else
	{
		delete ret;
		err_message("Please input a valid parameter name.");
	}

	return ret;
}

std::vector<std::string> Ellipse_Inclusion::getSensitivityParameterNames()
{
	std::vector<std::string> names = {"X", "Y", "A", "B", "ALPHA"};
	return names;
}






// NOTE: THIS ALPHA IS NOT THE SAME AS THE ANGLE ALPHA, IT IS A LINE SEARCH PARAMETER
// I should probably just change one of them...
double Ellipse_Inclusion::int_f(double alpha, void *params)
{
	struct int_params *p = (struct int_params*) params;
	std::vector<Node> nodes = p->nodes;
	double gcoords[2] = {nodes[0](0)+alpha*(nodes[1](0)-nodes[0](0)),
						 nodes[0](1)+alpha*(nodes[1](1)-nodes[0](1))};
	gcoords[0] -= _center[0];
	gcoords[1] -= _center[1];
	double lcoords[2] = {gcoords[0]*_ca+gcoords[1]*_sa,
						 gcoords[0]*_sa-gcoords[1]*_ca};

	double ret = pow(lcoords[0]/_a, 2) + pow(lcoords[1]/_b, 2) - 1;
	return ret;
}
double Ellipse_Inclusion::int_df(double alpha, void *params)
{
	struct int_params *p = (struct int_params*) params;
	std::vector<Node> nodes = p->nodes;
	double gcoords[2] = {nodes[0](0)+alpha*(nodes[1](0)-nodes[0](0)),
						 nodes[0](1)+alpha*(nodes[1](1)-nodes[0](1))};
	gcoords[0] -= _center[0];
	gcoords[1] -= _center[1];
	double d_dalpha[2] = {nodes[1](0)-nodes[0](0),
						  nodes[1](1)-nodes[0](1)};
	double lcoords[2] = {gcoords[0]*_ca+gcoords[1]*_sa,
						 gcoords[0]*_sa-gcoords[1]*_ca};

	double ret = (1/(_a*_a))*2.0*lcoords[0]*(d_dalpha[0]*_ca+d_dalpha[1]*_sa) + 
		  		 (1/(_b*_b))*2.0*lcoords[1]*(d_dalpha[0]*_sa-d_dalpha[1]*_ca);
	return ret;
}
void Ellipse_Inclusion::int_fdf(double alpha, void *params, double *f, double *df)
{
	struct int_params *p = (struct int_params*) params;
	std::vector<Node> nodes = p->nodes;
	double gcoords[2] = {nodes[0](0)+alpha*(nodes[1](0)-nodes[0](0)),
						 nodes[0](1)+alpha*(nodes[1](1)-nodes[0](1))};
	gcoords[0] -= _center[0];
	gcoords[1] -= _center[1];
	double d_dalpha[2] = {nodes[1](0)-nodes[0](0),
						  nodes[1](1)-nodes[0](1)};
	double lcoords[2] = {gcoords[0]*_ca+gcoords[1]*_sa,
						 gcoords[0]*_sa-gcoords[1]*_ca};

	*f = pow(lcoords[0]/_a, 2) + pow(lcoords[1]/_b, 2) - 1;
	*df = (1/(_a*_a))*2.0*lcoords[0]*(d_dalpha[0]*_ca+d_dalpha[1]*_sa) + 
		  (1/(_b*_b))*2.0*lcoords[1]*(d_dalpha[0]*_sa-d_dalpha[1]*_ca);
}
