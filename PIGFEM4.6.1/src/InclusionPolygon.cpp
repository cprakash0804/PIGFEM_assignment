/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated October 2016

##################################################################################
*/
#include "InclusionPolygon.h"
#include "Utilities.h"
#include <limits> // std::numeric_limits
#include <cmath> // sqrt


Polygon_Inclusion::Polygon_Inclusion()
{
}
Polygon_Inclusion::Polygon_Inclusion(Polygon_Inclusion & other_inc)
{
	copy(&other_inc);
}
Polygon_Inclusion::~Polygon_Inclusion()
{
	_vertices.clear();
	_bounds.clear();
}



Inclusion* Polygon_Inclusion::allocate_and_copy()
{
	Inclusion* new_inc = allocate();   // Hopefully this works to call the appropriate allocate function based on what instance calls this function
	copy(new_inc);
	return new_inc;
}
Inclusion* Polygon_Inclusion::allocate()
{
	return new Polygon_Inclusion;
}
// This function assumes the current material is set up correctly. Which it should't be possible for it not to be
void Polygon_Inclusion::copy(Inclusion* other_inc)
{
	other_inc->set_material(_mat);
	other_inc->set_id(_id);
	other_inc->set_vec_parameter("vertices", _vertices);
}




void Polygon_Inclusion::set_parameter(std::string name, double val)
{
	err_message("No parameters to be set in the Polygon_Inclusion class.");
}
double Polygon_Inclusion::get_parameter(std::string name)
{
	err_message("No parameters to get in the Polygon_Inclusion class.");
}




void Polygon_Inclusion::set_vec_parameter(std::string name, std::vector<double>& val)
{
	if((val.size()%2)!=0 || val.size()<6)
	{
		char buf [100];
		sprintf(buf, "The number of values in the vector must be a multiple of 2 and greater than or equal to 6 for a Polygon_Inclusion for the %s parameter.", name.c_str());
		err_message( buf );
	}

	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="VERTICES" || name=="VERTEX")
		set_vertices(val);
	else
	{
		char buf[50];
		sprintf(buf, "%s is not a valid paramter name.", name.c_str());
		err_message( buf );
	}
}
std::vector<double> Polygon_Inclusion::get_vec_parameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="VERTICES" || name=="VERTEX")
		return _vertices;
	else
	{
		char buf[50];
		sprintf(buf, "%s is not a valid paramter name.", name.c_str());
		err_message( buf );
	}
}

void Polygon_Inclusion::set_vertices(const std::vector<double>& vert)
{
	_vertices = vert;
	double min = -1.0*std::numeric_limits<double>::max();
	double max = std::numeric_limits<double>::max();
	double minx = max;
	double maxx = min;
	double miny = max;
	double maxy = min;
	for (id_type v=0; v<n_vertices(); ++v)
	{
		minx = (vx(v) < minx) ? vx(v) : minx;
		miny = (vy(v) < miny) ? vy(v) : miny;
		maxx = (vx(v) > minx) ? vx(v) : maxx;
		maxy = (vy(v) > miny) ? vy(v) : maxy;
	}
	_bounds = {minx, miny, min, maxx, maxy, max};
}

/*
 * Convenient wrapper functions to access the vertex coordinates
 */
double& Polygon_Inclusion::vx(const id_type& vertex)
{
	return _vertices[vertex*2];
}
double& Polygon_Inclusion::vy(const id_type& vertex)
{
	return _vertices[vertex*2 + 1];
}

std::vector<double> Polygon_Inclusion::get_domain_lims()
{
	return _bounds;
}


/*
Copyright (c) 1970-2003, Wm. Randolph Franklin

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
	1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.
	2. Redistributions in binary form must reproduce the above copyright notice in the documentation and/or other materials provided with the distribution.
	3. The name of W. Randolph Franklin may not be used to endorse or promote products derived from this Software without specific prior written permission.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO
THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS
OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
id_type Polygon_Inclusion::detect_node(Node* node)
{
	int c = 0;
	id_type nvert = n_vertices();
	double testx = (*node)(0), testy = (*node)(1);
	for (id_type i=0, j=(nvert-1); i<nvert; j = i++)
	{
		if ( ((vy(i)>testy) != (vy(j)>testy)) && 
			 (testx < (vx(j)-vx(i)) * (testy-vy(i)) / (vy(j)-vy(i)) + vx(i)) )
			c = !c;
	}

	if (c==0)	// Outside
		return 0;
	else		// Currently don't support nodes lying on polygon so just say its inside
		return 2;
}

std::vector<double> Polygon_Inclusion::get_surface_normal(Node* node)
{
	err_message("Currently nodes should never be found to lie on a polygon so this function should never be called...");
}








std::vector<double> Polygon_Inclusion::find_intersection(std::vector<Node>& nodes)
{
	std::vector<double> P0 = {nodes[0](0), nodes[0](1)};
	std::vector<double> P1 = {nodes[1](0), nodes[1](1)};
	std::vector<double> dS = Utilities::minus(P1, P0);
	for (id_type i=0, j=(n_vertices()-1); i<n_vertices(); j=i++)
	{
		std::vector<double> V0 = {vx(j), vy(j)};
		std::vector<double> V1 = {vx(i), vy(i)};
		std::vector<double> e = Utilities::minus(V1, V0);
		double alpha = Utilities::LineLineIntersection(P0, dS, V0, e);

		if (alpha < 0.0 || alpha > 1.0)	// Intersection occurs outside of the bounds of P0-P1
			continue;
		else	// Check if the intersection occurs in the Vi-Vi+1 edge
		{
			double D = -Utilities::perp(dS, e);		// dot(nS, e)
			std::vector<double> temp1 = Utilities::minus(V0, P0);
			double N = Utilities::perp(dS, temp1);	// -dot(ne, P0-V0)
			double beta = N / D;
			if (beta < 0.0 || beta > 1.0)
				continue;
			else	// Found the right edge
				return Utilities::VecAXPY(alpha, dS, P0);
		}
	}

	err_message("Failed to find an intersection point for a Polygonal inclusion.");
}





















double Polygon_Inclusion::int_f(double alpha, void *params)
{
	err_message("Polygon_Inclusion::int_f is not implemented.");
}
double Polygon_Inclusion::int_df(double alpha, void *params)
{
	err_message("Polygon_Inclusion::int_df is not implemented.");
}
void Polygon_Inclusion::int_fdf(double alpha, void *params, double *f, double *df)
{
	err_message("Polygon_Inclusion::int_fdf is not implemented.");
}
