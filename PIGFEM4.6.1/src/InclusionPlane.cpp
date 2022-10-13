/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated February 2017

##################################################################################
*/
#include "InclusionPlane.h"
#include "Utilities.h"
#include "SensitivityShapeParameterStructural.h"
#include "SensitivityShapeParameterThermal.h"
#include <limits> // std::numeric_limits
#include <cmath> // sqrt


Plane_Inclusion::Plane_Inclusion()
{
}
Plane_Inclusion::Plane_Inclusion(Plane_Inclusion & other_inc)
{
	copy(&other_inc);
}
Plane_Inclusion::~Plane_Inclusion()
{
	_normal.clear();
	_point.clear();
}


Inclusion* Plane_Inclusion::allocate_and_copy()
{
	Inclusion* new_inc = allocate();   // Hopefully this works to call the appropriate allocate function based on what instance calls this function
	copy(new_inc);
	return new_inc;
}
Inclusion* Plane_Inclusion::allocate()
{
	return new Plane_Inclusion;
}
// This function assumes the current material is set up correctly. Which it should't be possible for it not to be
void Plane_Inclusion::copy(Inclusion* other_inc)
{
	other_inc->set_material(_mat);
	other_inc->set_id(_id);
	other_inc->set_vec_parameter("Normal", _normal);
	other_inc->set_vec_parameter("Point", _point);
}




void Plane_Inclusion::set_parameter(std::string name, double val)
{
	if (name=="X")
	{
		_point.resize(3);
		_point[0] = val;
	}
	else if (name=="Y")
	{
		_point.resize(3);
		_point[1] = val;
	}
	else if (name=="Z")
	{
		_point.resize(3);
		_point[2] = val;
	}
	else if (name=="THETA")
		_theta = val;
	else if (name=="PHI")
		_phi = val;
	else
	{
		std::string out = name + " is not a valid parameter name!";
		err_message( out.data() );
	}
}
double Plane_Inclusion::get_parameter(std::string name)
{
	if (name=="X")
	{
		if (_point.size()==3)
			return _point[0];
		else
			err_message("Invalid parameter call in a Plane");
	}
	else if (name=="Y")
	{
		if (_point.size()==3)
			return _point[1];
		else
			err_message("Invalid parameter call in a Plane");
	}
	else if (name=="Z")
	{
		if (_point.size()==3)
			return _point[2];
		else
			err_message("Invalid parameter call in a Plane");
	}
	else if (name=="THETA")
		return _theta;
	else if (name=="PHI")
		return _phi;
	else
	{
		std::string out = name + " is not a valid parameter name!";
		err_message( out.data() );
	}
}




void Plane_Inclusion::set_vec_parameter(std::string name, std::vector<double>& val)
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
	if(name=="NORMAL" || name=="NORM")
	{
		// Make sure that the vector is normalized
		double mag = sqrt(val[0]*val[0] + val[1]*val[1] + val[2]*val[2]);
		for(id_type i=0; i<val.size(); ++i)
			val[i] = val[i]/mag;
		_normal = val;
	}
	else if(name=="POINT" || name=="P")
	{
		_point = val;
		for (id_type i=_point.size(); i<3; ++i)
			_point.push_back(0.0); // Make sure that it is storing a 3D point (easier math later)
	}
	else
	{
		char buf[50];
		sprintf(buf, "%s is not a valid paramter name.", name.c_str());
		err_message( buf );
	}
}
std::vector<double> Plane_Inclusion::get_vec_parameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="NORMAL" || name=="NORM")
		return _normal;
	else if(name=="POINT" || name=="P")
		return _point;
	else
	{
		char buf[50];
		sprintf(buf, "%s is not a valid paramter name.", name.c_str());
		err_message( buf );
	}
}

std::vector<double> Plane_Inclusion::get_domain_lims()
{
	double max = std::numeric_limits<double>::max();
	std::vector<double> ret = {-max, -max, -max, max, max, max};
	for(int i=0; i<3; ++i)
	{
		if(_normal[i] == -1.0)
			ret[i] = _point[i];
		else if(_normal[i] == 1.0)
			ret[i+3] = _point[i];
	}
	return ret;
}

id_type Plane_Inclusion::detect_node(Node* node)
{
	// Build a vector from the node to the point in the plane
	std::vector<double> vec = {_point[0]-(*node)(0), _point[1]-(*node)(1), _point[2]-(*node)(2)};

	// Normalize this vector
	double mag = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	for(int i=0; i<3; ++i)
		vec[i] = vec[i]/mag;

	// Perform the dot product between the vector and the normal
	double dot = vec[0]*_normal[0] + vec[1]*_normal[1] + vec[2]*_normal[2];

	double eps = 1e-12;
	if(dot < -1.0*eps)
		return 0;
	else if(dot > eps)
		return 2;
	else
		return 1;
}

std::vector<double> Plane_Inclusion::get_surface_normal(Node* node)
{
	return _normal;
}



std::vector<double> Plane_Inclusion::getIntersectionVelocity(std::vector<Node*>& nodes, std::vector<double> intersect, std::string param_name)
{
	std::transform(param_name.begin(), param_name.end(), param_name.begin(), ::toupper); // Capitilize the name
	double al = (*nodes[0])(0) - (*nodes[1])(0);
	double bl = (*nodes[0])(1) - (*nodes[1])(1);
	double cl = (*nodes[0])(2) - (*nodes[1])(2);

	double val = 1.0 / (al*_normal[0] + bl*_normal[1] + cl*_normal[2]);
	if (param_name=="X" || param_name=="XC" || param_name=="XS")
		val *= _normal[0];
	else if (param_name=="Y" || param_name=="YC" || param_name=="YS")
		val *= _normal[1];
	else if (param_name=="Z" || param_name=="ZC" || param_name=="ZS")
		val *= _normal[2];
	else if (param_name=="THETA")
		err_message("Sensitiviy not derived for Plane THETA sensitivity");
	else if (param_name=="PHI")
		err_message("Sensitiviy not derived for Plane PHI sensitivity");
	else
		err_message("Invalid shape parameter for velocity computation.");

	std::vector<double> ret = {al*val, bl*val, cl*val};
	std::cout << "Placeholder\n";
	return ret;
}




SensitivityShapeParameter* Plane_Inclusion::getSensitivityParameter(std::string param_name, classification type)
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
	else if (param_name=="THETA")
		ret->set_param_name("THETA");
	else if (param_name=="PHI")
		ret->set_param_name("PHI");
	else
	{
		delete ret;
		err_message("Please input a valid parameter name.");
	}

	return ret;
}


std::vector<std::string> Plane_Inclusion::getSensitivityParameterNames()
{
	std::vector<std::string> names = {"X", "Y", "Z", "THETA", "PHI"};
	return names;
}




// Finding an intersection is easy for a plane so this function is overridden here
std::vector<double> Plane_Inclusion::find_intersection(std::vector<Node>& nodes)
{
	std::vector<double> P0 = {nodes[0](0), nodes[0](1), nodes[0](2)};
	std::vector<double> P1 = {nodes[1](0), nodes[1](1), nodes[1](2)};
	std::vector<double> dS = Utilities::minus(P1, P0);

	double alpha = Utilities::LinePlaneIntersection(P0, dS, _point, _normal);
	if (alpha < 0.0 || alpha > 1.0)
		err_message("Invalid intersection found between 2 points and a plane inclusion.");

	return Utilities::VecAXPY(alpha, dS, P0);
}

double Plane_Inclusion::int_f(double alpha, void *params)
{
	err_message("Plane_Inclusion::int_f is not implemented.");
}
double Plane_Inclusion::int_df(double alpha, void *params)
{
	err_message("Plane_Inclusion::int_df is not implemented.");
}
void Plane_Inclusion::int_fdf(double alpha, void *params, double *f, double *df)
{
	err_message("Plane_Inclusion::int_fdf is not implemented.");
}
