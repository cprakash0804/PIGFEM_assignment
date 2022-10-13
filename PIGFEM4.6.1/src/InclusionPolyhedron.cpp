/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated October 2016

##################################################################################
*/
#include "InclusionPolyhedron.h"
#include "Utilities.h"
#include <limits> // std::numeric_limits
#include <cmath> // sqrt


Polyhedron_Inclusion::Polyhedron_Inclusion()
{
}
Polyhedron_Inclusion::Polyhedron_Inclusion(Polyhedron_Inclusion & other_inc)
{
	copy(&other_inc);
}
Polyhedron_Inclusion::~Polyhedron_Inclusion()
{
	_facets.clear();
	_normals.clear();
	_bounds.clear();
}


Inclusion* Polyhedron_Inclusion::allocate_and_copy()
{
	Inclusion* new_inc = allocate();   // Hopefully this works to call the appropriate allocate function based on what instance calls this function
	copy(dynamic_cast<Polyhedron_Inclusion*>(new_inc));
	return new_inc;
}

Inclusion* Polyhedron_Inclusion::allocate()
{
	return new Polyhedron_Inclusion;
}
// This function assumes the current material is set up correctly. Which it should't be possible for it not to be
void Polyhedron_Inclusion::copy(Polyhedron_Inclusion* other_inc)
{
	other_inc->set_material(_mat);
	other_inc->set_id(_id);
	other_inc->set_facets(_facets);
	other_inc->set_normals(_normals);
}




void Polyhedron_Inclusion::set_parameter(std::string name, double val)
{
	err_message("No parameters to be set in the Polyhedron_Inclusion class.");
}
double Polyhedron_Inclusion::get_parameter(std::string name)
{
	err_message("No parameters to get in the Polyhedron_Inclusion class.");
}
void Polyhedron_Inclusion::set_vec_parameter(std::string name, std::vector<double>& val)
{
	err_message("No parameters to get in the Polyhedron_Inclusion class.");
/*
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="FACE" || name=="FACES" || name=="FACET" || name=="FACETS")
	{
		 if (val.size()%9 != 0)
			err_message("There must be 9 coordinate values for every facet!");
		set_facets(val);
	}
	else if(name=="NORM" || name=="NORMAL" || name=="NORMS" || name=="NORMALS")
		_normals = val;
	else
	{
		char buf[50];
		sprintf(buf, "%s is not a valid paramter name.", name.c_str());
		err_message( buf );
	}
*/
}
std::vector<double> Polyhedron_Inclusion::get_vec_parameter(std::string name)
{
	err_message("No parameters to get in the Polyhedron_Inclusion class.");
/*
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="FACE" || name=="FACES" || name=="FACET" || name=="FACETS")
		return _facets;
	else if(name=="NORM" || name=="NORMAL" || name=="NORMS" || name=="NORMALS")
		return _normals;
	else
	{
		char buf[50];
		sprintf(buf, "%s is not a valid paramter name.", name.c_str());
		err_message( buf );
	}
*/
}





void Polyhedron_Inclusion::set_facets(const std::vector<std::vector<double> >& facets)
{
	_facets = facets;
	double min = -1.0*std::numeric_limits<double>::max();
	double max = std::numeric_limits<double>::max();
	double minx = max;
	double maxx = min;
	double miny = max;
	double maxy = min;
	double minz = max;
	double maxz = min;
	for (id_type f=0; f<n_facets(); ++f)
	{
		minx = (coord(f, 0, 0) < minx) ? coord(f, 0, 0) : minx;	maxx = (coord(f, 0, 0) > maxx) ? coord(f, 0, 0) : maxx;
		minx = (coord(f, 1, 0) < minx) ? coord(f, 1, 0) : minx;	maxx = (coord(f, 1, 0) > maxx) ? coord(f, 1, 0) : maxx;
		minx = (coord(f, 2, 0) < minx) ? coord(f, 2, 0) : minx;	maxx = (coord(f, 2, 0) > maxx) ? coord(f, 2, 0) : maxx;
		
		miny = (coord(f, 0, 1) < miny) ? coord(f, 0, 1) : miny;	maxy = (coord(f, 0, 1) > maxy) ? coord(f, 0, 1) : maxy;
		miny = (coord(f, 1, 1) < miny) ? coord(f, 1, 1) : miny;	maxy = (coord(f, 1, 1) > maxy) ? coord(f, 1, 1) : maxy;
		miny = (coord(f, 2, 1) < miny) ? coord(f, 2, 1) : miny;	maxy = (coord(f, 2, 1) > maxy) ? coord(f, 2, 1) : maxy;

		minz = (coord(f, 0, 2) < minz) ? coord(f, 0, 2) : minz;	maxz = (coord(f, 0, 2) > maxz) ? coord(f, 0, 2) : maxz;
		minz = (coord(f, 1, 2) < minz) ? coord(f, 1, 2) : minz;	maxz = (coord(f, 1, 2) > maxz) ? coord(f, 1, 2) : maxz;
		minz = (coord(f, 2, 2) < minz) ? coord(f, 2, 2) : minz;	maxz = (coord(f, 2, 2) > maxz) ? coord(f, 2, 2) : maxz;
	}
	_bounds = {minx, miny, minz, maxx, maxy, maxz};
}

/*
 * Convenient wrapper functions to access the vertex coordinates
 */
double Polyhedron_Inclusion::coord(id_type facet, id_type vertex, id_type dir)
{
	return _facets[facet][vertex*3 + dir];
}

std::vector<double> Polyhedron_Inclusion::get_domain_lims()
{
	return _bounds;
}


std::vector<double> Polyhedron_Inclusion::get_surface_normal(Node* node)
{
	err_message("Currently nodes should never be found to lie on a polygon so this function should never be called...");
}


bool Polyhedron_Inclusion::in_y_z_space(const std::vector<double>& P, id_type face)
{
	std::vector<double> ys = {coord(face,0,1), coord(face,1,1), coord(face,2,1)};
	std::vector<double> zs = {coord(face,0,2), coord(face,1,2), coord(face,2,2)};
	double ymin = *std::min_element(ys.begin(), ys.end());
	double zmin = *std::min_element(zs.begin(), zs.end());

	if (P[1] >= ymin && P[2] >= zmin)
	{
		double ymax = *std::max_element(ys.begin(), ys.end());
		double zmax = *std::max_element(zs.begin(), zs.end());
		if (P[1] <= ymax && P[2] <= zmax)
			return true;
		else
			return false;
	}
	else
		return false;
}



/*
 * Test is a node is inside a polyhedral inclusion using a ray casting techinque
 */
id_type Polyhedron_Inclusion::detect_node(Node* node)
{
	int c = 0;
	std::vector<double> P = {(*node)(0), (*node)(1), (*node)(2)};
	std::vector<double> dS = {_bounds[3]-_bounds[0], 0.0, 0.0}; // cast in the +x diection (with a vector that's roughly the same length as the object to reduce truncation)

	// Loop over all of the faces
	for (id_type f=0; f<n_facets(); ++f)
	{
		if (in_y_z_space(P, f)) // If the ray cast in the +x direction will actually pass throught the same y-z bounding box as the current face
		{
			// Find how far in the +x direction the intersection with the plane occurs
			std::vector<double> V0(_facets[f].begin(), _facets[f].begin()+3); // Get the coords of the 0th vertex of the current face
			std::vector<double>& n_f = get_normal(f); // Get the normal vector to the face (won't be normalized)
			double t = Utilities::LinePlaneIntersection(P, dS, V0, n_f);
			// If the intersection occurs in the positive direction
			if (t > 0)
			{
				std::vector<double> I = Utilities::VecAXPY(t, dS, P);
				
				std::vector<double> V1(_facets[f].begin()+3, _facets[f].begin()+6);
				std::vector<double> V2(_facets[f].begin()+6, _facets[f].begin()+9);
				if (Utilities::pointInTriangleBary(I, V0, V1, V2)) // If the intersection point actually occurs within the triangle
					c = !c;
			}
		}
	}

	if (c==0)	// Outside
		return 0;
	else		// Currently don't support nodes lying on polygon so just say its inside
		return 2;
}






std::vector<double> Polyhedron_Inclusion::find_intersection(std::vector<Node>& nodes)
{
	std::vector<double> P0 = {nodes[0](0), nodes[0](1), nodes[0](2)};
	std::vector<double> P1 = {nodes[1](0), nodes[1](1), nodes[1](2)};
	std::vector<double> dS = Utilities::minus(P1, P0);
	for (id_type f=0; f<n_facets(); f++)
	{
		// Find how far in the +x direction the intersection with the plane occurs
		std::vector<double> V0(_facets[f].begin(), _facets[f].begin()+3); // Get the coords of the 0th vertex of the current face
		std::vector<double>& n_f = get_normal(f); // Get the normal vector to the face (won't be normalized)
		double t = Utilities::LinePlaneIntersection(P0, dS, V0, n_f);
		// If the intersection occurs in the positive direction
		if (t < 0.0 || t > 1.0)
			continue;
		else
		{
			std::vector<double> I = Utilities::VecAXPY(t, dS, P0);
			std::vector<double> V1(_facets[f].begin()+3, _facets[f].begin()+6);
			std::vector<double> V2(_facets[f].begin()+6, _facets[f].begin()+9);
			if (Utilities::pointInTriangleBary(I, V0, V1, V2)) // If the intersection point actually occurs within the triangle
				return I;
			else
				continue;
		}
	}

	err_message("Failed to find an intersection point for a Polyhedral inclusion.");
}





















double Polyhedron_Inclusion::int_f(double alpha, void *params)
{
	err_message("Polyhedron_Inclusion::int_f is not implemented.");
}
double Polyhedron_Inclusion::int_df(double alpha, void *params)
{
	err_message("Polyhedron_Inclusion::int_df is not implemented.");
}
void Polyhedron_Inclusion::int_fdf(double alpha, void *params, double *f, double *df)
{
	err_message("Polyhedron_Inclusion::int_fdf is not implemented.");
}
