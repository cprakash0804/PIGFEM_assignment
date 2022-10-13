/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "Mesh.h"
#include "metis.h"
#include "abaqus_io.h"
#include "BoundaryObject.h"
#include "PartitionerSerial.h"
#include "DofObject.h"
#include "node.h"
#include "elem.h"
#include "material.h"
#include "Inclusion.h"
#include "Utilities.h"
#include <iostream>
#include <sys/time.h>
#include <string>
#include <fstream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <iterator>
#include <sstream>


static char help[] = "MPI-based parallel IGFEM nonlinear finite element code using PETSc.\n\n";



// Default constructor
Mesh::Mesh(int* argc, char*** argv)
	: _n_global_nodes(0), _n_global_elem(0), _dim(1), _rank(-1),
	  _init(false), _partitioned(false), _node_elem_generated(false),
	  _materials_assigned(false), _serial(false), _IGFEM(false), _was_generated(false), _periodic(false),
	  _n_global_enrich_nodes(0), _n_local_owned_enrich_nodes(0),
	  _nodes_detected(false), _elements_detected(false), _is_cohesive(false)
{
	init(argc, argv);
}

// Copy constructor
Mesh::Mesh(Mesh & other_mesh)
{
	// FIXME: I should probably make a copy constructor
}

// Routine to initialize anything that needs to be initialized before running stuff
PetscErrorCode Mesh::init(int *argc, char*** argv)
{
	_argc = argc;
	_argv = argv;
	MPI_Init(argc, argv);
	PetscErrorCode ierr = PetscInitialize(_argc, _argv, (char*)0, help);CHKERRQ(ierr);
	_communicator = MPI_COMM_WORLD;
	MPI_Comm_rank(_communicator, &_rank);
	MPI_Comm_size(_communicator, &_nranks);
	if (_nranks==1)
	{
		_partitioned = true;
		_serial = true;
	}
	_init = true;
	return ierr;
}

// Destructor
Mesh::~Mesh()
{
	clear();
	MPI_Finalize();
}

void Mesh::clear()
{
	// Clear all dynamically allocated memory
	for(id_type i=0; i<n_local_elem(); ++i)
		delete _elem[i];
	_elem.clear();
	
	for(id_type i=0; i<n_local_nodes(); ++i)
		delete _nodes[i];
	_nodes.clear();
	
	for(id_type i=0; i<_materials.size(); ++i)
		delete _materials[i];
	_materials.clear();
	
	// Clear all other data structures
	_node_owners.clear();
	_global_to_local_node.clear();
	_global_to_local_elem.clear();
	_elem_node.clear();
	_node_elem.clear();
	_remote_node_elem.clear();
	_node_partition_interface.clear();
	_elemsets.clear();
	_elem_to_material.clear();
	_nodesets.clear();
	_sidesets.clear();
	_dim = 0;
	_rank = 0;
	_nranks = 0;
	_partitioned = false;
	_node_elem_generated = false;
	_materials_assigned = false;
	_serial = false;
	_IGFEM = false;

	// Periodic stuff
	_periodic_nodesets.clear();
	_periodic_id_to_set.clear();
	_n_global_non_periodic_nodes = 0;
	_n_local_non_periodic_nodes = 0;
	_n_local_owned_non_periodic_nodes = 0;
	

	// IGFEM stuff
	for(id_type n=0; n< _enrich_nodes.size(); ++n)
		delete _enrich_nodes[n];
	for(id_type i=0; i<_inclusions.size(); ++i)
		delete _inclusions[i];
	_enrich_nodes.clear();
	_inclusions.clear();
	_global_to_local_enrich_node.clear();
	_enrich_owners.clear();
	_node_detect.clear();
	_n_global_enrich_nodes = 0;
	_n_local_owned_enrich_nodes = 0;
	_enrich_nodes_interface_lists.clear();
	_enrich_node_elem.clear();
	_nodes_detected = false;
	_elements_detected = false;
	_periodic = false;
}


void Mesh::Finalize()
{
	PetscFinalize();
}

void Mesh::set_periodic(bool periodic)
{
	_periodic = periodic;

	if (_periodic)
	{
		if (n_local_nodes() != 0)
		{
			if (_was_generated)
			{
				if (serial())
					build_periodic_nodesets(_Nx, _Ny, _Nz);
				else
					err_message("Can't build the periodic nodesets of a parallel mesh that has already been partitioned.");
			}
			else
				err_message("Can't build the periodic nodesets of a mesh that was not generated inside PCIGFEM.");
		}
	}
}

id_type Mesh::nElemInDim(id_type dim) const
{
	switch (dim)
	{
		case 1:
			return _Nx;
			break;
		case 2:
			return _Ny;
			break;
		case 3:
			return _Nz;
			break;
		default:
			err_message("Invalid mesh dimension");
	}
}


// Generate a mesh from some parameters
// Currently generates a line/square/cube mesh with dimensions size (max_x x max_y x max_z)
// Number of elements is (Nx x Ny x Nz)
void Mesh::generate_mesh(elem_type e_type, double min_x, double max_x, id_type Nx, double min_y, double max_y, id_type Ny, double min_z, double max_z, id_type Nz)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling generate_mesh.");

	_Nx = Nx;
	_Ny = Ny;
	_Nz = Nz;

	id_type nel;
	if (e_type==EDGE2) nel=2;
	else if (e_type==TRI3) nel=3;
	else if (e_type==QUAD4) nel=4;
	else if (e_type==TET4) nel=4;
	else if (e_type==PRISM6) nel=6;
	else if (e_type==HEX8) nel=8;
	else err_message("Element type selected is currently not supported.");

	// Check to make sure that a valid element type was selected
	switch(e_type)
	{
		case EDGE2:
			if (Nx==0 || Ny!=0 || Nz!=0)
				err_message("Unsupported element type for the given mesh dimension.");
			break;
		case TRI3:
		case QUAD4:
			if (Nx==0 || Ny==0 || Nz!=0)
				err_message("Unsupported element type for the given mesh dimension.");
			break;
		case TET4:
		case PRISM6:
		case HEX8:
			if (Nx==0 || Ny==0 || Nz==0)
				err_message("Unsupported element type for the given mesh dimension.");
			break;
		case INVALID_ELEM:
		case TET10:
		case HEX20:
		case HEX27:
		case PRISM15:
		case PRISM18:
		default:
			err_message("Unsupported element type for mesh generation");
	}


	// Set the minimum and maximum boundaries
	_minimum_bounds.push_back(min_x);
	_maximum_bounds.push_back(max_x);
	if (Ny!=0)
	{
		_minimum_bounds.push_back(min_y);
		_maximum_bounds.push_back(max_y);
		if (Nz!=0) // 2D: this means this would be the top and bottom
		{
			_minimum_bounds.push_back(min_z);
			_maximum_bounds.push_back(max_z);
		}
	}
	
	if (_rank == 0) // Only Generate the mesh on rank 0. The mesh will be partitioned later
	{
		// Figure out how big each element will be
		double el_size_x = (max_x-min_x)/Nx;
		double el_size_y = 0;
		double el_size_z = 0;
		_nodesets.insert(std::pair<std::string, std::set<id_type> >("LEFT", std::set<id_type>()));
		_nodesets.insert(std::pair<std::string, std::set<id_type> >("RIGHT", std::set<id_type>()));
		_sidesets.insert(std::pair<std::string, std::set<std::pair<id_type, id_type> > >
							("LEFT", std::set<std::pair<id_type, id_type> >()));
		_sidesets.insert(std::pair<std::string, std::set<std::pair<id_type, id_type> > >
							("RIGHT", std::set<std::pair<id_type, id_type> >()));

		if (Ny!=0)
		{
			el_size_y = (max_y-min_y)/Ny;
			if (Nz==0) // 2D: this means this would be the top and bottom
			{
				_nodesets.insert(std::pair<std::string, std::set<id_type> >("TOP", std::set<id_type>()));
				_nodesets.insert(std::pair<std::string, std::set<id_type> >("BOTTOM", std::set<id_type>()));
				_sidesets.insert(std::pair<std::string, std::set<std::pair<id_type, id_type> > >
							("TOP", std::set<std::pair<id_type, id_type> >()));
				_sidesets.insert(std::pair<std::string, std::set<std::pair<id_type, id_type> > >
							("BOTTOM", std::set<std::pair<id_type, id_type> >()));
			}
			else
			{
				_nodesets.insert(std::pair<std::string, std::set<id_type> >("FRONT", std::set<id_type>()));
				_nodesets.insert(std::pair<std::string, std::set<id_type> >("BACK", std::set<id_type>()));
				_sidesets.insert(std::pair<std::string, std::set<std::pair<id_type, id_type> > >
							("FRONT", std::set<std::pair<id_type, id_type> >()));
				_sidesets.insert(std::pair<std::string, std::set<std::pair<id_type, id_type> > >
							("BACK", std::set<std::pair<id_type, id_type> >()));
				el_size_z = (max_z-min_z)/Nz;
				_nodesets.insert(std::pair<std::string, std::set<id_type> >("TOP", std::set<id_type>()));
				_nodesets.insert(std::pair<std::string, std::set<id_type> >("BOTTOM", std::set<id_type>()));
				_sidesets.insert(std::pair<std::string, std::set<std::pair<id_type, id_type> > >
							("TOP", std::set<std::pair<id_type, id_type> >()));
				_sidesets.insert(std::pair<std::string, std::set<std::pair<id_type, id_type> > >
							("BOTTOM", std::set<std::pair<id_type, id_type> >()));
			}
		}

		// Create the nodes of the mesh (independent of element type)
		//_global_to_local_node.reserve((Nx+1)*(Ny+1)*(Nz+1)); // Preallocate storage to avoid rehashing of table
		_global_to_local_node.rehash(std::ceil(((Nx+1)*(Ny+1)*(Nz+1)) / _global_to_local_node.max_load_factor()));
		for(id_type z_iter=0; z_iter<=Nz; ++z_iter)
		{
			for(id_type y_iter=0; y_iter<=Ny; ++y_iter)
			{
				for(id_type x_iter=0; x_iter<=Nx; ++x_iter)
				{
					add_node(min_x+x_iter*el_size_x, min_y+y_iter*el_size_y, min_z+z_iter*el_size_z, (x_iter+y_iter*(Nx+1)+z_iter*(Nx+1)*(Ny+1)), 0);

					if (x_iter==0) _nodesets["LEFT"].insert(_nodes.size()-1);
					if (x_iter==Nx) _nodesets["RIGHT"].insert(_nodes.size()-1);
					if (Ny != 0)
					{
						if (Nz==0)	// 2D
						{
							if (y_iter==0) _nodesets["BOTTOM"].insert(_nodes.size()-1);
							if (y_iter==Ny) _nodesets["TOP"].insert(_nodes.size()-1);
						}
						else 		// 3D
						{
							if (y_iter==0) _nodesets["FRONT"].insert(_nodes.size()-1);
							if (y_iter==Ny) _nodesets["BACK"].insert(_nodes.size()-1);
							if (z_iter==0) _nodesets["BOTTOM"].insert(_nodes.size()-1);
							if (z_iter==Nz) _nodesets["TOP"].insert(_nodes.size()-1);
						}
					}
				}
			}
		}
		
		std::vector<id_type> el_nodes;
		if (Ny == 0) // 1D, Using Edge2 elements
		{
			_dim = 1;
			el_nodes.resize(nel);
			//_global_to_local_elem.reserve(Nx); // Preallocate storage to avoid rehashing
			_global_to_local_elem.rehash(std::ceil(Nx / _global_to_local_elem.max_load_factor()));
			for(id_type x_iter=0; x_iter<Nx; ++x_iter)
			{
				el_nodes[0] = x_iter;
				el_nodes[1] = x_iter+1;
				Elem* el = add_elem(EDGE2, el_nodes, _elem.size());
				if (x_iter==0) _sidesets["LEFT"].insert(std::pair<id_type, id_type>(el->get_id(), 0));
				if (x_iter==(Nx-1)) _sidesets["RIGHT"].insert(std::pair<id_type, id_type>(el->get_id(), 1));
			}
		}
		else if (Nz == 0) // 2D
		{
			_dim = 2;
			el_nodes.resize(nel);
			// Using Quad4 Elements
			if (e_type==QUAD4)
			{
				//_global_to_local_elem.reserve(Nx*Ny); // Preallocate storage to avoid rehashing
				_global_to_local_elem.rehash(std::ceil((Nx*Ny) / _global_to_local_elem.max_load_factor()));
				for(id_type y_iter=0; y_iter<Ny; ++y_iter)
				{
					for(id_type x_iter=0; x_iter<Nx; ++x_iter)
					{
						el_nodes[0] = (x_iter)+(y_iter)*(Nx+1);
						el_nodes[1] = (x_iter+1)+(y_iter)*(Nx+1);
						el_nodes[2] = (x_iter+1)+(y_iter+1)*(Nx+1);
						el_nodes[3] = (x_iter)+(y_iter+1)*(Nx+1);
						Elem* el = add_elem(QUAD4, el_nodes, _elem.size());
						if (x_iter==0) _sidesets["LEFT"].insert(std::pair<id_type, id_type>(el->get_id(), 3));
						if (x_iter==(Nx-1)) _sidesets["RIGHT"].insert(std::pair<id_type, id_type>(el->get_id(), 1));
						if (y_iter==0) _sidesets["BOTTOM"].insert(std::pair<id_type, id_type>(el->get_id(), 0));
						if (y_iter==(Ny-1)) _sidesets["TOP"].insert(std::pair<id_type, id_type>(el->get_id(), 2));
					}
				}
			}
			// Using Tri3 Elements
			else if (e_type==TRI3)
			{
				//_global_to_local_elem.reserve(2*Nx*Ny); // Preallocate storage to avoid rehashing
				_global_to_local_elem.rehash(std::ceil((2*Nx*Ny) / _global_to_local_elem.max_load_factor()));
				for(id_type y_iter=0; y_iter<Ny; ++y_iter)
				{
					for(id_type x_iter=0; x_iter<Nx; ++x_iter)
					{
						if ((x_iter+y_iter)%2 == 0)
						{
							el_nodes[0] = (x_iter)+(y_iter)*(Nx+1);
							el_nodes[1] = (x_iter+1)+(y_iter)*(Nx+1);
							el_nodes[2] = (x_iter+1)+(y_iter+1)*(Nx+1);
							Elem* el1 = add_elem(TRI3, el_nodes, _elem.size());
							el_nodes[0] = (x_iter)+(y_iter)*(Nx+1);
							el_nodes[1] = (x_iter+1)+(y_iter+1)*(Nx+1);
							el_nodes[2] = (x_iter)+(y_iter+1)*(Nx+1);
							Elem* el2 = add_elem(TRI3, el_nodes, _elem.size());
							if (x_iter==0) _sidesets["LEFT"].insert(std::pair<id_type, id_type>(el2->get_id(), 2));
							if (x_iter==(Nx-1)) _sidesets["RIGHT"].insert(std::pair<id_type, id_type>(el1->get_id(), 1));
							if (y_iter==0) _sidesets["BOTTOM"].insert(std::pair<id_type, id_type>(el1->get_id(), 0));
							if (y_iter==(Ny-1)) _sidesets["TOP"].insert(std::pair<id_type, id_type>(el2->get_id(), 1));
						}
						else
						{
							el_nodes[0] = (x_iter)+(y_iter)*(Nx+1);
							el_nodes[1] = (x_iter+1)+(y_iter)*(Nx+1);
							el_nodes[2] = (x_iter)+(y_iter+1)*(Nx+1);
							Elem* el1 = add_elem(TRI3, el_nodes, _elem.size());
							el_nodes[0] = (x_iter+1)+(y_iter)*(Nx+1);
							el_nodes[1] = (x_iter+1)+(y_iter+1)*(Nx+1);
							el_nodes[2] = (x_iter)+(y_iter+1)*(Nx+1);
							Elem* el2 = add_elem(TRI3, el_nodes, _elem.size());
							if (x_iter==0) _sidesets["LEFT"].insert(std::pair<id_type, id_type>(el1->get_id(), 2));
							if (x_iter==(Nx-1)) _sidesets["RIGHT"].insert(std::pair<id_type, id_type>(el2->get_id(), 0));
							if (y_iter==0) _sidesets["BOTTOM"].insert(std::pair<id_type, id_type>(el1->get_id(), 0));
							if (y_iter==(Ny-1)) _sidesets["TOP"].insert(std::pair<id_type, id_type>(el2->get_id(), 1));
						}
					}
				}
			}
		}
		else // 3D, Using Hex8 elements
		{
			_dim = 3;
			el_nodes.resize(nel);
			if (e_type==HEX8)
			{
				//_global_to_local_elem.reserve(Nx*Ny*Nz); // Preallocate storage to avoid rehashing
				_global_to_local_elem.rehash(std::ceil((Nx*Ny*Nz) / _global_to_local_elem.max_load_factor()));
				for(id_type z_iter=0; z_iter<Nz; ++z_iter)
				{
					for(id_type y_iter=0; y_iter<Ny; ++y_iter)
					{
						for(id_type x_iter=0; x_iter<Nx; ++x_iter)
						{
							el_nodes[0] = ((x_iter)+(y_iter)*(Nx+1) + (z_iter)*(Nx+1)*(Ny+1));
							el_nodes[1] = ((x_iter+1)+(y_iter)*(Nx+1) + (z_iter)*(Nx+1)*(Ny+1));
							el_nodes[2] = ((x_iter+1)+(y_iter+1)*(Nx+1) + (z_iter)*(Nx+1)*(Ny+1));
							el_nodes[3] = ((x_iter)+(y_iter+1)*(Nx+1) + (z_iter)*(Nx+1)*(Ny+1));
							el_nodes[4] = ((x_iter)+(y_iter)*(Nx+1) + (z_iter+1)*(Nx+1)*(Ny+1));
							el_nodes[5] = ((x_iter+1)+(y_iter)*(Nx+1) + (z_iter+1)*(Nx+1)*(Ny+1));
							el_nodes[6] = ((x_iter+1)+(y_iter+1)*(Nx+1) + (z_iter+1)*(Nx+1)*(Ny+1));
							el_nodes[7] = ((x_iter)+(y_iter+1)*(Nx+1) + (z_iter+1)*(Nx+1)*(Ny+1));
							Elem* el = add_elem(HEX8, el_nodes, _elem.size());
							if (x_iter==0) _sidesets["LEFT"].insert(std::pair<id_type, id_type>(el->get_id(), 4));
							if (x_iter==(Nx-1)) _sidesets["RIGHT"].insert(std::pair<id_type, id_type>(el->get_id(), 2));
							if (y_iter==0) _sidesets["FRONT"].insert(std::pair<id_type, id_type>(el->get_id(), 1));
							if (y_iter==(Ny-1)) _sidesets["BACK"].insert(std::pair<id_type, id_type>(el->get_id(), 3));
							if (z_iter==0) _sidesets["BOTTOM"].insert(std::pair<id_type, id_type>(el->get_id(), 0));
							if (z_iter==(Nz-1)) _sidesets["TOP"].insert(std::pair<id_type, id_type>(el->get_id(), 5));
						}
					}
				}
			}
			else if (e_type==PRISM6)
			{
				//_global_to_local_elem.reserve(2*Nx*Ny*Nz); // Preallocate storage to avoid rehashing
				_global_to_local_elem.rehash(std::ceil((2*Nx*Ny*Nz) / _global_to_local_elem.max_load_factor()));
				for(id_type z_iter=0; z_iter<Nz; ++z_iter)
				{
					for(id_type y_iter=0; y_iter<Ny; ++y_iter)
					{
						for(id_type x_iter=0; x_iter<Nx; ++x_iter)
						{
							el_nodes[0] = ((x_iter)+(y_iter)*(Nx+1) + (z_iter)*(Nx+1)*(Ny+1));		// 0
							el_nodes[1] = ((x_iter+1)+(y_iter)*(Nx+1) + (z_iter)*(Nx+1)*(Ny+1));	// 1
							el_nodes[2] = ((x_iter)+(y_iter+1)*(Nx+1) + (z_iter)*(Nx+1)*(Ny+1));	// 3
							el_nodes[3] = ((x_iter)+(y_iter)*(Nx+1) + (z_iter+1)*(Nx+1)*(Ny+1));	// 4
							el_nodes[4] = ((x_iter+1)+(y_iter)*(Nx+1) + (z_iter+1)*(Nx+1)*(Ny+1));	// 5
							el_nodes[5] = ((x_iter)+(y_iter+1)*(Nx+1) + (z_iter+1)*(Nx+1)*(Ny+1));	// 7
							Elem* el1 = add_elem(PRISM6, el_nodes, _elem.size());
							el_nodes[0] = ((x_iter+1)+(y_iter)*(Nx+1) + (z_iter)*(Nx+1)*(Ny+1));	// 1
							el_nodes[1] = ((x_iter+1)+(y_iter+1)*(Nx+1) + (z_iter)*(Nx+1)*(Ny+1));	// 2
							el_nodes[2] = ((x_iter)+(y_iter+1)*(Nx+1) + (z_iter)*(Nx+1)*(Ny+1));	// 3
							el_nodes[3] = ((x_iter+1)+(y_iter)*(Nx+1) + (z_iter+1)*(Nx+1)*(Ny+1));	// 5
							el_nodes[4] = ((x_iter+1)+(y_iter+1)*(Nx+1) + (z_iter+1)*(Nx+1)*(Ny+1));// 6
							el_nodes[5] = ((x_iter)+(y_iter+1)*(Nx+1) + (z_iter+1)*(Nx+1)*(Ny+1));	// 7
							Elem* el2 = add_elem(PRISM6, el_nodes, _elem.size());
							if (x_iter==0) _sidesets["LEFT"].insert(std::pair<id_type, id_type>(el1->get_id(), 3));
							if (x_iter==(Nx-1)) _sidesets["RIGHT"].insert(std::pair<id_type, id_type>(el2->get_id(), 3));
							if (y_iter==0) _sidesets["FRONT"].insert(std::pair<id_type, id_type>(el1->get_id(), 1));
							if (y_iter==(Ny-1)) _sidesets["BACK"].insert(std::pair<id_type, id_type>(el2->get_id(), 1));
							if (z_iter==0)
							{
								_sidesets["BOTTOM"].insert(std::pair<id_type, id_type>(el1->get_id(), 0));
								_sidesets["BOTTOM"].insert(std::pair<id_type, id_type>(el2->get_id(), 0));
							}
							if (z_iter==(Nz-1))
							{
								_sidesets["TOP"].insert(std::pair<id_type, id_type>(el1->get_id(), 4));
								_sidesets["TOP"].insert(std::pair<id_type, id_type>(el2->get_id(), 4));
							}
						}
					}
				}
			}
			else if (e_type==TET4) // Divide this into 6 tets instead of 5 to maintain consistent edges across unit cells
			{
				//_global_to_local_elem.reserve(6*Nx*Ny*Nz); // Preallocate storage to avoid rehashing
				_global_to_local_elem.rehash(std::ceil((6*Nx*Ny*Nz) / _global_to_local_elem.max_load_factor()));
				for(id_type z_iter=0; z_iter<Nz; ++z_iter)
				{
					for(id_type y_iter=0; y_iter<Ny; ++y_iter)
					{
						for(id_type x_iter=0; x_iter<Nx; ++x_iter)
						{
							el_nodes[0] = ((x_iter)+(y_iter)*(Nx+1) + (z_iter)*(Nx+1)*(Ny+1));			// 0
							el_nodes[1] = ((x_iter+1)+(y_iter+1)*(Nx+1) + (z_iter)*(Nx+1)*(Ny+1));		// 2
							el_nodes[2] = ((x_iter)+(y_iter+1)*(Nx+1) + (z_iter)*(Nx+1)*(Ny+1));		// 3
							el_nodes[3] = ((x_iter+1)+(y_iter+1)*(Nx+1) + (z_iter+1)*(Nx+1)*(Ny+1));	// 6
							Elem* el1 = add_elem(TET4, el_nodes, _elem.size());
							el_nodes[0] = ((x_iter)+(y_iter)*(Nx+1) + (z_iter)*(Nx+1)*(Ny+1));			// 0
							el_nodes[1] = ((x_iter)+(y_iter+1)*(Nx+1) + (z_iter)*(Nx+1)*(Ny+1));		// 3
							el_nodes[2] = ((x_iter)+(y_iter+1)*(Nx+1) + (z_iter+1)*(Nx+1)*(Ny+1));		// 7
							el_nodes[3] = ((x_iter+1)+(y_iter+1)*(Nx+1) + (z_iter+1)*(Nx+1)*(Ny+1));	// 6
							Elem* el2 = add_elem(TET4, el_nodes, _elem.size());
							el_nodes[0] = ((x_iter)+(y_iter)*(Nx+1) + (z_iter)*(Nx+1)*(Ny+1));			// 0
							el_nodes[1] = ((x_iter)+(y_iter+1)*(Nx+1) + (z_iter+1)*(Nx+1)*(Ny+1));		// 7
							el_nodes[2] = ((x_iter)+(y_iter)*(Nx+1) + (z_iter+1)*(Nx+1)*(Ny+1));		// 4
							el_nodes[3] = ((x_iter+1)+(y_iter+1)*(Nx+1) + (z_iter+1)*(Nx+1)*(Ny+1));	// 6
							Elem* el3 = add_elem(TET4, el_nodes, _elem.size());
							el_nodes[0] = ((x_iter)+(y_iter)*(Nx+1) + (z_iter)*(Nx+1)*(Ny+1));			// 0
							el_nodes[1] = ((x_iter+1)+(y_iter)*(Nx+1) + (z_iter+1)*(Nx+1)*(Ny+1));		// 5
							el_nodes[2] = ((x_iter+1)+(y_iter+1)*(Nx+1) + (z_iter+1)*(Nx+1)*(Ny+1));	// 6
							el_nodes[3] = ((x_iter)+(y_iter)*(Nx+1) + (z_iter+1)*(Nx+1)*(Ny+1));		// 4
							Elem* el4 = add_elem(TET4, el_nodes, _elem.size());
							el_nodes[0] = ((x_iter+1)+(y_iter)*(Nx+1) + (z_iter)*(Nx+1)*(Ny+1));		// 1
							el_nodes[1] = ((x_iter+1)+(y_iter)*(Nx+1) + (z_iter+1)*(Nx+1)*(Ny+1));		// 5
							el_nodes[2] = ((x_iter+1)+(y_iter+1)*(Nx+1) + (z_iter+1)*(Nx+1)*(Ny+1));	// 6
							el_nodes[3] = ((x_iter)+(y_iter)*(Nx+1) + (z_iter)*(Nx+1)*(Ny+1));			// 0
							Elem* el5 = add_elem(TET4, el_nodes, _elem.size());
							el_nodes[0] = ((x_iter+1)+(y_iter)*(Nx+1) + (z_iter)*(Nx+1)*(Ny+1));		// 1
							el_nodes[1] = ((x_iter+1)+(y_iter+1)*(Nx+1) + (z_iter+1)*(Nx+1)*(Ny+1));	// 6
							el_nodes[2] = ((x_iter+1)+(y_iter+1)*(Nx+1) + (z_iter)*(Nx+1)*(Ny+1));		// 2
							el_nodes[3] = ((x_iter)+(y_iter)*(Nx+1) + (z_iter)*(Nx+1)*(Ny+1));			// 0
							Elem* el6 = add_elem(TET4, el_nodes, _elem.size());
							if (x_iter==0)
							{
								_sidesets["LEFT"].insert(std::pair<id_type, id_type>(el2->get_id(), 0));
								_sidesets["LEFT"].insert(std::pair<id_type, id_type>(el3->get_id(), 0));
							}
							if (x_iter==(Nx-1))
							{
								_sidesets["RIGHT"].insert(std::pair<id_type, id_type>(el5->get_id(), 0));
								_sidesets["RIGHT"].insert(std::pair<id_type, id_type>(el6->get_id(), 0));
							}
							if (y_iter==0)
							{
								_sidesets["FRONT"].insert(std::pair<id_type, id_type>(el4->get_id(), 1));
								_sidesets["FRONT"].insert(std::pair<id_type, id_type>(el5->get_id(), 1));
							}
							if (y_iter==(Ny-1))
							{
								_sidesets["BACK"].insert(std::pair<id_type, id_type>(el1->get_id(), 2));
								_sidesets["BACK"].insert(std::pair<id_type, id_type>(el2->get_id(), 2));
							}
							if (z_iter==0)
							{
								_sidesets["BOTTOM"].insert(std::pair<id_type, id_type>(el1->get_id(), 0));
								_sidesets["BOTTOM"].insert(std::pair<id_type, id_type>(el6->get_id(), 3));
							}
							if (z_iter==(Nz-1))
							{
								_sidesets["TOP"].insert(std::pair<id_type, id_type>(el3->get_id(), 2));
								_sidesets["TOP"].insert(std::pair<id_type, id_type>(el4->get_id(), 2));
							}
						}
					}
				}
			}
		}

		// Set all _elem_to_material entries to -1 so that we know which elements have had a material set for it. This will be partitioned and the values will be maintained
		_elem_to_material.resize(_elem.size());
		for(id_type i=0; i<_elem.size(); ++i)
			_elem_to_material[i] = -1;  // Set as the not assigned value

		// If the mesh is periodic, create the periodic nodesets
		if (periodic())
			build_periodic_nodesets(Nx, Ny, Nz);
	}

	// Set boolean
	_was_generated = true;

	PartitionerSerial partitioner(this);
	partitioner.partition();
	_n_original_nodes = n_local_nodes();
	_elem_activity.resize(n_local_elem());
	std::fill(_elem_activity.begin(), _elem_activity.end(), true);
	_active_elem = _elem;
	_n_original_elements = n_local_elem();
	_mesh_topology = 0; // This mesh has not been refined
	generate_node_elem();
}


void Mesh::getDomainLimits(std::vector<double>& mins, std::vector<double>& maxes)
{
	if (!_was_generated)
		err_message("Can only get domain limits of a mesh for a generated mesh!");

	else
	{
		mins = _minimum_bounds;
		maxes = _maximum_bounds;
	}
}



/*
 * Function to build the periodic nodesets associated with a serial generated mesh
 */
void Mesh::build_periodic_nodesets(id_type Nx, id_type Ny, id_type Nz)
{
	if (_dim == 1) // For some reason dim() doesn't work here...
	{
		_n_global_non_periodic_nodes = Nx;
		std::vector<id_type> edges = {0, Nx};
		_periodic_nodesets.push_back(edges);
		_periodic_id_to_set[0] = 0;
		_periodic_id_to_set[Nx] = 0;
	}
	else if (_dim == 2)
	{
		_n_global_non_periodic_nodes = Nx*Ny;
		_periodic_nodesets.resize(1+(Nx-1)+(Ny-1));
		id_type count = 0;

		// Add the corner nodeset
		std::vector<id_type> corners = {0, Nx, (Nx+1)*Ny, (Nx+1)*(Ny+1)-1};
		for (id_type i=0; i<corners.size(); ++i)
			_periodic_id_to_set[corners[i]] = count;
		_periodic_nodesets[count].swap(corners);
		count++;

		// Loop over the bottom and add the link to the top
		for (id_type n=1; n<Nx; ++n)
		{
			std::vector<id_type> vec = {n, (Nx+1)*Ny+n};
			for (id_type i=0; i<vec.size(); ++i)
				_periodic_id_to_set[vec[i]] = count;
			_periodic_nodesets[count].swap(vec);
			count++;
		}

		// Loop over the left and add the link to the right
		for (id_type n=1; n<Ny; ++n)
		{
			std::vector<id_type> vec = {n*(Nx+1), (n+1)*(Nx+1)-1};
			for (id_type i=0; i<vec.size(); ++i)
				_periodic_id_to_set[vec[i]] = count;
			_periodic_nodesets[count].swap(vec);
			count++;
		}
	}
	else if (_dim == 3)
	{
		_n_global_non_periodic_nodes = Nx*Ny*Nz;
		_periodic_nodesets.resize(1 +
								  (Nx-1) + (Ny-1) + (Nz-1) +
								  (Nx-1)*(Ny-1) + (Nx-1)*(Nz-1) + (Ny-1)*(Nz-1));
		id_type count = 0;
		id_type top_start = (Nx+1)*(Ny+1)*Nz;

		// Add the corner nodeset
		std::vector<id_type> corners = {0, Nx, (Nx+1)*Ny, (Nx+1)*(Ny+1)-1,
										top_start+0, top_start+Nx, top_start+(Nx+1)*Ny, top_start+(Nx+1)*(Ny+1)-1};
		for (id_type i=0; i<corners.size(); ++i)
			_periodic_id_to_set[corners[i]] = count;
		_periodic_nodesets[count].swap(corners);
		count++;

		// Add the edges parallel to the x-axis
		for (id_type n=1; n<Nx; ++n)
		{
			std::vector<id_type> vec = {n, (Nx+1)*Ny+n, n+top_start, (Nx+1)*Ny+n+top_start};
			for (id_type i=0; i<vec.size(); ++i)
				_periodic_id_to_set[vec[i]] = count;
			_periodic_nodesets[count].swap(vec);
			count++;
		}

		// Add the edges parallel to the y-axis
		for (id_type n=1; n<Ny; ++n)
		{
			std::vector<id_type> vec = {n*(Nx+1), (n+1)*(Nx+1)-1, n*(Nx+1)+top_start, (n+1)*(Nx+1)-1+top_start};
			for (id_type i=0; i<vec.size(); ++i)
				_periodic_id_to_set[vec[i]] = count;
			_periodic_nodesets[count].swap(vec);
			count++;
		}

		// Add the edges parallel to the z-axis
		for (id_type n=1; n<Nz; ++n)
		{
			id_type level_start = n*(Nx+1)*(Ny+1);
			std::vector<id_type> vec = {level_start+0, level_start+Nx, level_start+(Nx+1)*Ny, level_start+(Nx+1)*(Ny+1)-1};
			for (id_type i=0; i<vec.size(); ++i)
				_periodic_id_to_set[vec[i]] = count;
			_periodic_nodesets[count].swap(vec);
			count++;
		}

		// Loop over the bottom and add the link to the top
		for (id_type n1=1; n1<Nx; ++n1)
		{
			for (id_type n2=1; n2<Ny; ++n2)
			{
				std::vector<id_type> vec = {n2*(Nx+1)+n1, n2*(Nx+1)+n1+top_start};
				for (id_type i=0; i<vec.size(); ++i)
					_periodic_id_to_set[vec[i]] = count;
				_periodic_nodesets[count].swap(vec);
				count++;
			}
		}

		// Loop over the left and add the link to the right
		for (id_type n1=1; n1<Ny; ++n1)
		{
			for (id_type n2=1; n2<Nz; ++n2)
			{
				std::vector<id_type> vec = {n1*(Nx+1)+n2*(Nx+1)*(Ny+1), n1*(Nx+1)+n2*(Nx+1)*(Ny+1)+Nx};
				for (id_type i=0; i<vec.size(); ++i)
					_periodic_id_to_set[vec[i]] = count;
				_periodic_nodesets[count].swap(vec);
				count++;
			}
		}

		// Loop over the front and add the link to the back
		for (id_type n1=1; n1<Nx; ++n1)
		{
			for (id_type n2=1; n2<Nz; ++n2)
			{
				std::vector<id_type> vec = {n1+n2*(Nx+1)*(Ny+1), n1+n2*(Nx+1)*(Ny+1) + Ny*(Nx+1)};
				for (id_type i=0; i<vec.size(); ++i)
					_periodic_id_to_set[vec[i]] = count;
				_periodic_nodesets[count].swap(vec);
				count++;
			}
		}
	}

	_n_local_non_periodic_nodes = _n_global_non_periodic_nodes;
	_n_local_owned_non_periodic_nodes = _n_global_non_periodic_nodes;
}



// Read a mesh from an input file (Abaqus)
//	Taken from libMesh
void Mesh::read_mesh(const std::string& fname)
{
	if (_rank==0)
	{
		AbaqusIO io_object;
		io_object.set_mesh(this);
		io_object.read(fname);

		// Set all _elem_to_material entries to -1 so that we know which elements have had a material set for it. This will be partitioned and the values will be maintained
		_elem_to_material.resize(_elem.size());
		for(id_type i=0; i<_elem.size(); ++i)
			_elem_to_material[i] = -1;  // Set as the not assigned value
	}
	

	PartitionerSerial partitioner(this);
	partitioner.partition();
	_n_original_nodes = n_local_nodes();
	_elem_activity.resize(n_local_elem());
	std::fill(_elem_activity.begin(), _elem_activity.end(), true);
	_active_elem = _elem;
	_n_original_elements = n_local_elem();
	_mesh_topology = 0; // This mesh has not been refined
	_periodic = false; // I don't wanna figure out how to make an Abaqus mesh periodic
	generate_node_elem();
}


// Rotate the mesh about the origin. Angles are about the fixed z-y-x coordinate axes
void Mesh::rotate_mesh(double z_angle, double y_angle, double x_angle)
{
	// Rotate about the z-axis
	if (z_angle != 0.0)
	{
		for(Mesh::node_iterator it=nodes_begin(), end=nodes_end(); it!=end; ++it)
		{
			double x = (*(*it))(0);
			double y = (*(*it))(1);
			double angle = atan2(y, x);
			angle = Utilities::constrain_angle(angle + z_angle);
			double mag = sqrt(x*x + y*y);
			(*it)->set_coord(mag*cos(angle), 0);
			(*it)->set_coord(mag*sin(angle), 1);
		}
	}

	// Rotate about the y-axis
	if (y_angle != 0.0)
	{
		for(Mesh::node_iterator it=nodes_begin(), end=nodes_end(); it!=end; ++it)
		{
			double x = (*(*it))(0);
			double z = (*(*it))(2);
			double angle = atan2(x, z);
			angle = Utilities::constrain_angle(angle + y_angle);
			double mag = sqrt(x*x + z*z);
			(*it)->set_coord(mag*cos(angle), 2);
			(*it)->set_coord(mag*sin(angle), 0);
		}
	}

	// Rotate about the x-axis
	if (x_angle != 0.0)
	{
		for(Mesh::node_iterator it=nodes_begin(), end=nodes_end(); it!=end; ++it)
		{
			double y = (*(*it))(0);
			double z = (*(*it))(2);
			double angle = atan2(z, y);
			angle = Utilities::constrain_angle(angle + y_angle);
			double mag = sqrt(y*y + z*z);
			(*it)->set_coord(mag*cos(angle), 1);
			(*it)->set_coord(mag*sin(angle), 2);
		}
	}
}




// FUNCTIONS THAT HANDLE THE BOUNDARY (NODESETS AND SIDESETS SPECIFICALLY)
// Add a new nodeset. To do this you need a name and a set of the global node ids that will belong to the nodeset
void Mesh::add_nodeset(std::string name, std::set<id_type>& nodes)
{	
	// Check to see if a nodeset with the same name has been created before
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if( nodesetExists(name) )
		err_message("Cannot create multiple nodesets with the same name.");
	else
	{
		// Figure out which nodes are actually in the local mesh and add them
		std::set<id_type> nodes_to_add;
		for(auto it=nodes.begin(), end=nodes.end(); it!=end; ++it)
			if(node_in_local_mesh(*it))
				nodes_to_add.insert(*it);

		_nodesets[name].swap(nodes_to_add);
	}
}

// Add a new nodeset without any checks so its faster
void Mesh::add_nodeset_total(std::string name, std::set<id_type>& nodes)
{
	// Check to see if a nodeset with the same name has been created before
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if ( nodesetExists(name) )
		err_message("Cannot create multiple nodesets with the same name.");
	else
		_nodesets[name].swap(nodes);
}

void Mesh::add_to_nodeset(std::string name, id_type val)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	std::map<std::string, std::set<id_type> >::iterator it = _nodesets.find(name);
	if (it != _nodesets.end())
	{
		if (node_in_local_mesh(val))
			it->second.insert(val);
	}
	else
		err_message("Attempted to add a node id to a node set that does not exist in the mesh.");
}


// Add a new sideset. To do this you need a name and a vector of the global element ids that will belong to the elemset
void Mesh::add_sideset(std::string name, std::set<std::pair<id_type, id_type> >& sides)
{
	// Check to see if a sideset with the same name has been created before
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if( sidesetExists(name) )
		err_message("Cannot create multiple sidesets with the same name.");
	else
	{
		// Figure out which elements are actually in the local mesh and add them
		std::set<std::pair<id_type, id_type> > sides_to_add;
		for(auto it=sides.begin(), end=sides.end(); it!=end; ++it)
			if(elem_in_local_mesh(it->first))
				sides_to_add.insert(*it);

		_sidesets[name].swap(sides_to_add);
	}
}

// Add a new sideset without checks so its faster
void Mesh::add_sideset_total(std::string name, std::set<std::pair<id_type, id_type> >& sides)
{
	// Check to see if a sideset with the same name has been created before
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if ( sidesetExists(name) )
		err_message("Cannot create multiple sidesets with the same name.");
	else
		_sidesets[name].swap(sides);
}

void Mesh::add_to_sideset(std::string name, std::pair<id_type, id_type> val)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	std::map<std::string, std::set<std::pair<id_type, id_type> > >::iterator it = _sidesets.find(name);
	if (it != _sidesets.end())
	{
		if(elem_in_local_mesh(val.first))
			it->second.insert(val);
	}
	else
		err_message("Attempted to add an element id to a side set that does not exist in the mesh.");
}


// Function to return the given nodeset
std::set<id_type> Mesh::get_nodeset(id_type idx)
{
	if(idx>=n_nodesets())
		err_message("Please select a valid nodeset.");
	else
	{
		id_type count=0;
		for(auto it=_nodesets.begin(), end=_nodesets.end(); it!=end; ++it, ++count)
			if(count==idx)
				return it->second;
	}
	err_message("We'll never get here!");
}


// Function to return the given nodeset
std::set<id_type> Mesh::get_nodeset(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if ( nodesetExists(name) )
		return _nodesets[name];
	else
		err_message("Please select a valid nodeset.");
}


// Function to return the given elemset
std::set<std::pair<id_type, id_type> > Mesh::get_sideset(id_type idx)
{
	if (idx>=n_sidesets())
		err_message("Please select a valid sideset.");
	else
	{
		id_type count=0;
		for(auto it=_sidesets.begin(), end=_sidesets.end(); it!=end; ++it, ++count)
			if(count==idx)
				return it->second;
	}
	err_message("We'll never get here!");
}


// Function to return the given elemset
std::set<std::pair<id_type, id_type> > Mesh::get_sideset(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if ( sidesetExists(name) )
		return _sidesets[name];
	else
		err_message("Please select a valid sideset.");
}

bool Mesh::nodesetExists(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if (_nodesets.find(name) != _nodesets.end())
		return true;
	else
		return false;
}
bool Mesh::sidesetExists(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if (_sidesets.find(name) != _sidesets.end())
		return true;
	else
		return false;
}






















void Mesh::delete_elem(id_type id)
{
	if (!elem_in_local_mesh(id))
	{
		id_type idx = _global_to_local_elem[id];
		_elem.erase(_elem.begin()+idx);
		_elem_node.erase(_elem_node.begin()+idx);
		for(id_type i=idx; i<n_local_elem(); ++i)
			_global_to_local_elem[_elem[i]->get_id()]--;

		_n_global_elem--; // Note, this doesn't update on remote processors... However, the only time that I use this is in read_mesh which is serial so I guess its ok
	}
}


Elem* Mesh::add_elem(elem_type etype, std::vector<id_type>& node_ids, id_type global_id)
{
	// Make sure this doesn't exist in the mesh already
	if (elem_in_local_mesh(global_id))
	{
		char buf[100], buf2[100];
		strcpy(buf, "An element with the id %");
		strcat(buf, SPEC);
		strcat(buf, " already exists in the mesh.");
		sprintf(buf2, buf, global_id);
		err_message( buf2 );
	}

	// Build the actual element
	Elem* el = Elem::build(etype);
	
	// Check to make sure that this function call is valid
	if (node_ids.size()!=el->n_nodes())
	{
		delete el;
		err_message("While adding element to the mesh, the number of nodes does not match the number of nodes in the element.");
	}

	// Set the mesh dimension
	_dim = std::max(_dim, el->dim());
	
	// Set the element id and add it to the mesh
	bool do_en_table = false;
	if (_elem.size() == _elem_node.size()) // Actually adding a new element to the mesh, not just one that was partitioned.
		do_en_table = true;
	el->set_id( global_id );
	_elem.push_back(el);
	if (do_en_table)
		_elem_node.push_back( std::vector<id_type>(el->n_nodes()) );
	id_type e_size = _elem.size() - 1;
	
	// Set the entry for the elem_node table and set the elemental node pointers
	for(id_type i=0; i<el->n_nodes(); ++i)
	{
		id_type l_node = global_to_local_node(node_ids[i]);
		Node* node = get_node_local(l_node);
		if (do_en_table)
			_elem_node[e_size][i] = l_node;

		el->set_node(i) = node;
	}
	
	// Update the entry in the global to local element mapping
	_global_to_local_elem.insert(std::pair<id_type, id_type>(global_id, _elem.size()-1));

	// Update the global number of elements
	_n_global_elem++;

	return el;
}



Node* Mesh::add_node(double x, double y, double z, id_type global_id, int owner)
{
	// Make sure this doesn't exist in the mesh already
	if (node_in_local_mesh(global_id))
	{
		char buf[250], buf2[250];
		strcpy(buf, "A node with the id ");
		strcat(buf, SPEC);
		strcat(buf, "already exists in the mesh.");
		sprintf(buf2, buf, global_id);
		err_message( buf2 );
	}
	
	// Add the new Node to the mesh
	_nodes.push_back( new Node(x, y, z, global_id) );
	
	// Update the global to local node mapping
	_global_to_local_node.insert(std::pair<id_type, id_type>(global_id, _nodes.size()-1));
	
	// Update the Node owners array
	_node_owners.push_back(owner);

	// Update the number of global nodes
	_n_global_nodes++;

	return _nodes[_nodes.size()-1];
}





// Function that return the mesh partition number that corresponds to the owner of the given local node number
int Mesh::get_node_owner_local(id_type idx)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling get_node_owner_local.");
	if (!_partitioned) err_message("Must partition the mesh prior to calling get_node_owner_local.");
	
	if (idx < ENRICH_START)
	{
		if (idx<n_local_nodes())
			return _node_owners[idx];
		else
			err_message("Attempted to access a node owner corresponding to a local node not in the local mesh.");
	}
	else if ((idx-ENRICH_START) < n_local_enrich_nodes())
		return _enrich_owners[idx-ENRICH_START];
	else
		err_message("Attempted to access a node owner corresponding to a local node not in the local mesh.");
}

// Function that return the mesh partition number that corresponds to the owner of the given global node id
int Mesh::get_node_owner_global(id_type id)
{
	id_type idx = global_to_local_node(id);
	return get_node_owner_local(idx);
}

// Returns the full node-elem result, which is both the node_elem_local and remote_node_elem. Based on global node number
std::vector<id_type> Mesh::get_node_elem(id_type id)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling get_node_elem.");
	if (!_partitioned) err_message("Must partition the mesh prior to calling get_node_elem.");
	if (!_node_elem_generated) err_message("Must generate the node_elem connectivity table prior to calling get_node_elem.");

	if (id<ENRICH_START) // A normal node
	{
		id_type l_node = global_to_local_node( id );
		std::vector<id_type> ne = get_node_elem_local(l_node);

		std::unordered_map<id_type, std::map<int, std::vector<id_type> > >::iterator it2 = _remote_node_elem.find(id);
		if (it2 != _remote_node_elem.end()) // There are remote neighbors to this node
			for(auto it3=it2->second.begin(), end=it2->second.end(); it3!=end; ++it3)
				ne.insert(ne.end(), it3->second.begin(), it3->second.end()); // Append the retrieved vector to the total list
			
		return ne;
	}
	else if ((id-ENRICH_START)<n_global_enrich_nodes())
		return get_node_elem_global(id-ENRICH_START); // FIXME: This is so simple because currently I dont store the remote node_elem result for enrichment nodes...
	else
		err_message("Attempted to access node_elem table of a global node id not in the global mesh.");
}

// Returns the global element numbers of all local elements neighboring the node with the given local node id
// NOTE: local elements means that this will only return element ids tha exist on the local processor. To get element ids on another processor use get_remote_node_elem.
//		 Should really abstract that away so that there's a function to get all neighboring global element ids regardless of which processor they're on.
std::vector<id_type>& Mesh::get_node_elem_local(id_type idx)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling get_node_elem_local.");
	if (!_partitioned) err_message("Must partition the mesh prior to calling get_node_elem_local.");
	if (!_node_elem_generated) err_message("Must generate the node_elem connectivity table prior to calling get_node_elem_local.");
	
	if (idx < ENRICH_START)
	{
		if (idx<n_local_nodes())
			return _node_elem[idx];
		else
			err_message("Cannot return elements neighboring node, local node number does not exist on the local mesh.");
	}
	else if ((idx-ENRICH_START) < n_local_enrich_nodes())
		return _enrich_node_elem[idx-ENRICH_START];
	else
		err_message("Cannot return elements neighboring node, local node number does not exist on the local mesh.");
}

// Returns the global element numbers of all elements neighboring the node with the given global node id
// NOTE: local elements means that this will only return element ids tha exist on the local processor. To get element ids on another processor use get_remote_node_elem.
//		 Should really abstract that away so that there's a function to get all neighboring global element ids regardless of which processor they're on.
std::vector<id_type>& Mesh::get_node_elem_global(id_type id)
{
	id_type idx = global_to_local_node(id);
	return get_node_elem_local(idx);
}

std::vector<id_type> Mesh::get_remote_node_elem(id_type id, int rank)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling get_remote_node_elem.");
	if (!_partitioned) err_message("Must partition the mesh prior to calling get_remote_node_elem.");
	if (!_node_elem_generated) err_message("Must generate the node_elem connectivity table prior to calling get_remote_node_elem.");
	
	if (id < ENRICH_START) // This is a normal node
	{
		if (id < n_global_nodes()) // This node exists in the global mesh
		{
			std::unordered_map<id_type, std::map<int, std::vector<id_type> > >::iterator it = _remote_node_elem.find(id);
			if (it != _remote_node_elem.end()) // node is in the map
			{
				std::map<int, std::vector<id_type> >::iterator it2 = it->second.find(rank);
				if (it2 != it->second.end())  // node exists on the given rank as well
					return it2->second;
				else
					return std::vector<id_type>();
			}
			else
				return std::vector<id_type>();
		}
		else
			err_message("Attempted to access remote node_elem table of a global node id not in the global mesh.");
	}
	else if ((id-ENRICH_START) < n_global_enrich_nodes())
		return std::vector<id_type>();
	else
		err_message("Attempted to access remote node_elem table of a global node id not in the global mesh.");
}

// Find the global elements (local or remote) that the list of global nodes has in common
std::vector<id_type> Mesh::elem_in_common(const std::vector<id_type>& nodes)
{
	if (!_node_elem_generated)
		err_message("Must call generate_node_elem() before calling elem_in_common()");
	
	if (nodes.size() == 0)
		err_message("Must pass in nodes to find the common partitions between them.");

	std::vector<id_type> intersect = get_node_elem(nodes[0]);
	std::sort(intersect.begin(), intersect.end());
	std::vector<id_type> temp_intersect;
	for(id_type i=1; i<nodes.size(); ++i)
	{
		std::vector<id_type> elem_set = get_node_elem(nodes[i]);
		std::sort(elem_set.begin(), elem_set.end());
		std::set_intersection(intersect.begin(), intersect.end(), 
							  elem_set.begin(), elem_set.end(),
							  std::back_inserter(temp_intersect));
		intersect.swap(temp_intersect);
		temp_intersect.clear();
	}

	return intersect;
}

// Returns global element numbers that share the global nodes in the nodes vector
std::vector<id_type> Mesh::local_elem_in_common(const std::vector<id_type>& nodes)
{
	if (!_node_elem_generated)
		err_message("Must call generate_node_elem() before calling local_elem_in_common()");
	
	if (nodes.size() == 0)
		err_message("Must pass in nodes to find the common partitions between them.");

	std::vector<id_type> intersect = get_node_elem_global(nodes[0]);
	std::sort(intersect.begin(), intersect.end());
	std::vector<id_type> temp_intersect;
	for(id_type i=1; i<nodes.size(); ++i)
	{
		std::vector<id_type> elem_set = get_node_elem_global(nodes[i]);
		std::sort(elem_set.begin(), elem_set.end());
		std::set_intersection(intersect.begin(), intersect.end(), 
							  elem_set.begin(), elem_set.end(),
							  std::back_inserter(temp_intersect));
		intersect.swap(temp_intersect);
		temp_intersect.clear();
	}

	return intersect;
}

std::map<int, std::vector<id_type> > Mesh::remote_elem_in_common(const std::vector<id_type>& nodes)
{
	if (!_node_elem_generated)
		err_message("Must call generate_node_elem() before calling elem_in_common()");
	
	if (nodes.size() == 0)
		err_message("Must pass in nodes to find the common partitions between them.");

	// If any of the nodes aren't on a partition interface then there will never be any remote nodes in common between the nodes
	std::map<int, std::vector<id_type> > common;
	for (id_type n=0; n<nodes.size(); ++n)
		if (!node_on_part_interface(nodes[n]))
			return common;

	// Get the result for the first node
	common = _remote_node_elem[nodes[0]];

	std::set<int> empties;
	std::vector<id_type> temp_intersect;
	for (auto it=common.begin(), end=common.end(); it!=end; ++it) // FIXME: probably best memory access-wise to reverse these loops... But this way makes the most sense to me right now
	{
		int part = it->first;
		std::vector<id_type>& intersect = it->second;
		std::sort(intersect.begin(), intersect.end());
		for (id_type n=1; n<nodes.size(); ++n)
		{
			auto it2 = _remote_node_elem.find(nodes[n]);
			if (it2 != _remote_node_elem.end())
			{
				auto it3 = (*it2).second.find(part);
				if (it3 != (*it2).second.end())
				{
					std::vector<id_type> elem_set = it3->second;
					std::sort(elem_set.begin(), elem_set.end());
					std::set_intersection(intersect.begin(), intersect.end(), 
										  elem_set.begin(), elem_set.end(),
										  std::back_inserter(temp_intersect));
					intersect.swap(temp_intersect);
					temp_intersect.clear();
				}
				else
				{
					intersect.clear();
					break;
				}
			}
			else // This node doesn't even exist on a partition interface so it will never have anything in common
			{
				common.clear();
				return common;
			}
		}

		// Mark the partitions with no entries as emty to delete later
		if (intersect.size() == 0)
			empties.insert(part);
	}

	// Remove entries from the map with no elements in common
	for (auto it=empties.begin(), end=empties.end(); it!=end; ++it)
		common.erase( *it );

	return common;
}

// Function returns the raw pointer to the local node idx. Good for node modification. Use carefully.
Node* Mesh::get_node_local(id_type idx)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling get_node_local.");

	if (idx < ENRICH_START)
	{
		if (idx<n_local_nodes())
			return _nodes[idx];
		else
			err_message("Attempting to access a node not existing on the local mesh.");
	}
	else if ((idx-ENRICH_START) < n_local_enrich_nodes())
		return _enrich_nodes[idx-ENRICH_START];
	else
		err_message("Attempting to access a node not existing on the local mesh.");
}

// Function returns the raw pointer to the local node idx. Good for node modification. Use carefully.
Node* Mesh::get_node_global(id_type id)
{
	id_type idx = global_to_local_node(id);
	return get_node_local(idx);
}

// Nice wrapper function to determine if a global node id exists in the local mesh
bool Mesh::node_in_local_mesh(id_type id)
{
	if (id < ENRICH_START)
	{
		if (_global_to_local_node.find(id)!=_global_to_local_node.end())
			return true;
		else
			return false;
	}
	else if ((id-ENRICH_START) < n_global_enrich_nodes())
	{
		if (_global_to_local_enrich_node.find(id)!=_global_to_local_enrich_node.end())
			return true;
		else
			return false;
	}
	else
		err_message("Attempting to access a node not existing on the global mesh.");
}

bool Mesh::own_node_global(id_type id)
{
	return get_node_owner_global(id)==_rank;
}
bool Mesh::own_node_local(id_type idx)
{
	return get_node_owner_local(idx)==_rank;
}

// Return true if the given node id is on a  periodic boundary
bool Mesh::node_periodic(id_type id)
{
	auto it = _periodic_id_to_set.find( id );
	if (it == _periodic_id_to_set.end())
		return false;
	else
		return true;
}

// Returns the set of node ids that is associated with this periodic node id
std::vector<id_type> Mesh::get_node_periodic_nodeset(id_type id)
{
	auto it = _periodic_id_to_set.find( id );
	if (it == _periodic_id_to_set.end())
		err_message("Attempting to access the periodic nodeset of a node not on a periodic boundary.");
	else
		return _periodic_nodesets[it->second];
}

// Returns whether or not the given node id is the primary node of its periodic nodeset (the lowest id)
bool Mesh::primary_periodic_node(id_type id)
{
	auto it = _periodic_id_to_set.find( id );
	if (it == _periodic_id_to_set.end())
		err_message("Attempting to access the periodic nodeset of a node not on a periodic boundary.");
	else
		return (_periodic_nodesets[it->second][0] == id);
}

// Checks whether or not this mesh partition is responsible for assigning the dofs for this node id
bool Mesh::check_node_responsibility(id_type n_id)
{
	if (own_node_global(n_id))
	{
		if (node_periodic(n_id))
		{
			if (primary_periodic_node(n_id))
				return true;
			else
				return false;
		}
		else
			return true;
	}
	else
		return false;
}



void Mesh::store_shape_functions()
{
	// Clear any old mesh data
	_shape.clear();
	_shape_grad.clear();
	_J.clear();
	_W.clear();
	_Rot.clear();
	_volumetric_qps.clear();

	// Initialize data structures
	id_type nelem = n_local_elem();
	_shape.resize( nelem );
	_shape_grad.resize( nelem );
	_J.resize( nelem );
	_W.resize( nelem );
	_Rot.resize( nelem );
	_volumetric_qps.resize( nelem );

	// Store all shape functions, even for non-active elements (In case of coarsening)
	for(Mesh::element_iterator it=elements_begin(), end=elements_end(); it!=end; ++it)
	{
		Elem* el =(*it);
		id_type local_e = global_to_local_elem(el->get_id());

		// If the element isn't intersected
		if (!el->is_intersected())
		{
			id_type nqp = el->n_q_points();
			_shape[local_e].resize(nqp);
			_shape_grad[local_e].resize(nqp);
			_J[local_e].resize(nqp);
			_W[local_e].resize(nqp);
			_volumetric_qps[local_e] = nqp;
			for (id_type qp=0; qp<nqp; ++qp)
			{
				std::vector<double> xi;
				el->q_point(xi, _W[local_e][qp], qp);
				el->ShapeFunctions(xi, _shape[local_e][qp], _shape_grad[local_e][qp], _J[local_e][qp]);

				// Check for a negative Jacobian (Allow for some noise?)
				if (_J[local_e][qp] < -1e-16)
				{
					std::stringstream ss;
					ss << "Negative element Jacobian (element=" << local_e << " J=" << _J[local_e][qp] << ") detected.";
					std::string out = ss.str();
					err_message( out.data() );
				}
			}
		}

		// The element is intersected
		else
		{
			id_type ngqp = 0;
			for (id_type ie=0; ie<el->n_integration_elem(); ++ie)
				ngqp += el->get_integration_elem(ie)->n_q_points();
			_volumetric_qps[local_e] = ngqp;
			id_type nvqp = ngqp;
			for (id_type ce=0; ce<el->n_cohesive_elem(); ++ce)
				ngqp += el->get_cohesive_elem(ce)->n_q_points();

			_shape[local_e].resize(ngqp);
			_shape_grad[local_e].resize(ngqp);
			_J[local_e].resize(ngqp);
			_W[local_e].resize(ngqp);
			_Rot[local_e].resize(ngqp - _volumetric_qps[local_e]);

			id_type local_qp = 0;;
			for (id_type ie=0; ie<(*it)->n_integration_elem(); ++ie)
			{
				Elem* int_el = el->get_integration_elem(ie);
				id_type nqp = int_el->n_q_points();

				for (id_type qp=0; qp<nqp; ++qp)
				{
					// First get the integration point in child element referenece coordinates 
					std::vector<double> child_rcoords;
					int_el->q_point(child_rcoords, _W[local_e][local_qp], qp);

					// Compute the child element's shape functions, gradients, and mapping Jacobian determinant
					std::vector<double> child_shape(int_el->n_nodes());
					std::vector<std::vector<double> > child_grad_x;
					int_el->ShapeFunctions(child_rcoords, child_shape, child_grad_x, _J[local_e][local_qp]);

					// Check for a negative Jacobian
					if (_J[local_e][local_qp] < 0.0)
						err_message("Negative element Jacobian detected.");

					// Convert these coordinates to the global coordinates now
					std::vector<double> gcoords(el->dim());
					for(id_type d=0; d<dim(); ++d)
						for(id_type n=0; n<int_el->n_nodes(); ++n)
							gcoords[d] += (*int_el)(n)(d)*child_shape[n];

					// Now convert these back to parent element reference coordinates
					bool inside = false;
					std::vector<double> rcoords = el->inverse_map(gcoords, inside);
					if (!inside)
						err_message("Bad inverse map detected!");

					// Use the parent coordinates to compute the parent shape functions and gradients at the new reference coords
					std::vector<double> shape(el->n_nodes());
					std::vector<std::vector<double> > grad_x;
					double parent_J;
					el->ShapeFunctions(rcoords, shape, grad_x, parent_J);

					// Check that we got the correct inverse map
					std::vector<double> gcoords2(el->dim());
					for(id_type d=0; d<dim(); ++d)
						for(id_type n=0; n<el->n_nodes(); ++n)
							gcoords2[d] += (*el)(n)(d)*shape[n];
					double diff_sq = 0.0;
					for(id_type d=0; d<dim(); ++d)
						diff_sq += pow(gcoords2[d]-gcoords[d], 2);
					if (sqrt(diff_sq) > 1e-9)
						err_message("Bad inverse map in assembly detected.");


					// Now, look through the integration element's nodes and see if they match any of the enrichment nodes
					// If they do, then I need to add that shape function and gradient to the appropriate place in the shape and grad_x structures
					shape.resize(el->n_nodes() + el->n_enrich_nodes());
					for(id_type en=0; en<el->n_enrich_nodes(); en++)
						grad_x.push_back(std::vector<double>(dim()));
					
					for(id_type en=0; en<el->n_enrich_nodes(); ++en)
					{
						id_type enrich_id = el->get_enrich_node(en)->get_id();
						for(id_type n=0; n<int_el->n_nodes(); ++n)
						{
							id_type int_node_id = (*int_el)(n).get_id();
							if (int_node_id == enrich_id)
							{
								shape[el->n_nodes()+en] = child_shape[n];
								grad_x[el->n_nodes()+en] = child_grad_x[n];
								break;
							}
						}
					}

					// Set the shape function and gradients
					_shape[local_e][local_qp] = shape;
					_shape_grad[local_e][local_qp] = grad_x;

					// Update the local quadrature point
					local_qp++;
				}
			}
			for (id_type ce=0; ce<(*it)->n_cohesive_elem(); ++ce)
			{
				CohesiveElem* coh_el = (*it)->get_cohesive_elem(ce);
				id_type nqp = coh_el->n_q_points();
				for (id_type qp=0; qp<nqp; ++qp)
				{
					std::vector<double> xi;
					coh_el->q_point(xi, _W[local_e][local_qp], qp);
					coh_el->ShapeFunctions(xi, _shape[local_e][local_qp], _Rot[local_e][local_qp-nvqp], _J[local_e][local_qp]);
					local_qp++;
				}
			}
		}
	}
}

double& Mesh::get_W(id_type local_e, id_type qp)
{
	if (local_e < _W.size())
	{
		if (qp < _W[local_e].size())
			return _W[local_e][qp];
		else
			err_message("Invalid Quadrature point.");
	}
	else
		err_message("Invalid element number.");
}
double& Mesh::get_J(id_type local_e, id_type qp)
{
	if (local_e < _J.size())
	{
		if (qp < _J[local_e].size())
			return _J[local_e][qp];
		else
			err_message("Invalid Quadrature point.");
	}
	else
		err_message("Invalid element number.");
}
std::vector<double>& Mesh::get_shape(id_type local_e, id_type qp)
{
	if (local_e < _shape.size())
	{
		if (qp < _shape[local_e].size())
			return _shape[local_e][qp];
		else
			err_message("Invalid Quadrature point.");
	}
	else
		err_message("Invalid element number.");
}
std::vector<std::vector<double> >& Mesh::get_shape_grad(id_type local_e, id_type qp)
{
	if (local_e < _shape_grad.size())
	{
		if (qp < _shape_grad[local_e].size())
			return _shape_grad[local_e][qp];
		else
			err_message("Invalid Quadrature point.");
	}
	else
		err_message("Invalid element number.");
}
DenseMatrix<double>& Mesh::get_rot_matrix(id_type local_e, id_type qp)
{
	if (local_e < _Rot.size())
	{
		qp -= _volumetric_qps[local_e]; // Didn't store rotation matricies for non-cohesive quadrature points
		if (qp >= 0 && qp < _Rot[local_e].size())
			return _Rot[local_e][qp];
		else
			err_message("Invalid Quadrature point.");
	}
	else
		err_message("Invalid element number.");
}
id_type Mesh::n_volumetric_qps(id_type local_e)
{
	if (local_e < _volumetric_qps.size())
	{
		return _volumetric_qps[local_e];
	}
	else
		err_message("Invalid element number.");
}
std::vector<double> Mesh::get_global_coords_undeformed(id_type l_elem, id_type local_qp)
{
	std::vector<double> ret(dim());
	Elem* el = _elem[l_elem];
	std::vector<double>& N = get_shape(l_elem, local_qp);
	for (id_type n=0; n<el->n_nodes(); ++n)
		for (id_type d=0; d<dim(); ++d)
			ret[d] += N[n]*((*el->get_node(n))(d));

	return ret;
}













































std::vector<id_type>& Mesh::get_elem_node_local(id_type idx)
{
	if (idx < n_local_elem())
		return _elem_node[idx];
	else
		err_message("Attempted to access the Elem-Node connectivity table for an element not in the local mesh.");
}
std::vector<id_type>& Mesh::get_elem_node_global(id_type id)
{
	id_type idx = global_to_local_elem(id);
	return get_elem_node_local(idx);
}


// Function returns the raw pointer to the local element idx. Good for element modification. Use carefully.
Elem* Mesh::get_elem_local(id_type idx)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling get_elem_local.");
	
	if (idx >=_elem.size()) // Element is not in the local mesh
		return NULL;
	else
		return _elem[idx];
}

// Function returns the raw pointer to the global element id. Good for element modification. Use carefully.
Elem* Mesh::get_elem_global(id_type id)
{
	id_type idx = global_to_local_elem(id);
	return get_elem_local(idx);
}






/*
 * A function to add a new material to the mesh
 * To use this you should create a local instance of the material and then pass in the address to this function.
 * You must define all parameters and the name of the material. This function will handle the assignment of material ids
 * Ex:
 *	LinearElasticIsotropicMaterial material;
 *	material.set_parameter("E", 69000000000);
 *	material.set_parameter("nu", 0.334);
 *	material.set_name("Aluminum");
 *	mesh.add_material(&material);
*/
 int Mesh::find_mat_number(std::string name)
 {
	std::map<std::string, int>::iterator it = _mat_name_to_idx.find(name);
	if (it != _mat_name_to_idx.end())
		return it->second;
	else
		return -1;
 }

void Mesh::add_material(Material* mat)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling add_material.");
	
	if (!mat->name_set()) err_message("Please give your material a name before attempting to add it to the mesh."); // Because materials use the name as a unique identifier

	int mat_num = find_mat_number(mat->get_name());
	if (mat_num != -1)
		err_message("Cannot create multiple materials with the same name.");

	id_type max=0;
	for(id_type i=0; i<_materials.size(); ++i)
	{
		id_type id = _materials[i]->get_id();
		if (id > max)
			max = id;
	}
	// FIXME: Do some sort of reduction here to determine the max overall id?
	mat->set_id(max+1);
	_materials.push_back(mat->allocate_and_copy());
	_mat_name_to_idx.insert(std::pair<std::string, int>(mat->get_name(), _materials.size()-1));

	if (mat->is_cohesive())
	{
		_cohesive_mats.push_back(_materials.size()-1);
		_cohesive_to_mat_pair.push_back(std::pair<int, int>(-1,-1)); // Put a wildcard pair as a placeholder
		_is_cohesive = true;
	}
}


// Function that returns the raw pointer to the ith material object. This will allow for modification of current materials in the mesh
Material* Mesh::get_material(id_type idx)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling get_material.");
	
	if (idx < _materials.size())
		return _materials[idx];
	else
		err_message("Material index does not match a material. Please select a valid material.");
}


// Function that returns the raw pointer to the material with the matching name. This will allow for modification of current materials in the mesh
Material* Mesh::get_material(std::string name)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling get_material.");
	
	int mat_num = find_mat_number(name);

	if (mat_num<0) // Didn't find a valid matching material
		err_message("Material name does not match a material. Please select a valid material.");
	else
		return _materials[mat_num];
}


// Function that a user should call if they wish to remain ignorant of mesh partitioning. Simply calls get_element_material_global
Material* Mesh::get_elem_material(id_type id, id_type int_el)
{
	return get_element_material_global(id, int_el);
}


// Function that returns a pointer to the material contained in the elemnt with global id id (FIXME: this function doesn't really work agnostic of partitioning
Material* Mesh::get_element_material_global(id_type id, id_type int_el)
{
	id_type idx = global_to_local_elem( id );
	return get_element_material_local(idx, int_el);
}


// Function that returns a pointer to the material contained in the local element idx
Material* Mesh::get_element_material_local(id_type idx, id_type int_el)
{
	if (idx < _elem_to_material.size()) // Local mesh containes the element
	{
		Elem* el = get_elem_local(idx);
		if (el->is_intersected())
		{
			if (int_el < el->n_integration_elem())
				return el->get_int_elem_mat( int_el );
			else
				err_message("Attempted to access the element material for an integraion element not existing in the given element.");
		}
		else
		{
			if (_elem_to_material[idx]>=0) // Contains a valid material
				return _materials[_elem_to_material[idx]];
			else
				err_message("Attempted to access the element material for an element for which a material has not been assigned.");
		}
	}
	else
		err_message("Attempted to access the material for an element not existing in the local mesh.");
}


// Add a new elemset. To do this you need a name and a vector of the global element ids that will belong to the elemset
void Mesh::add_elemset(std::string name, std::set<id_type>& elems)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling add_elemset.");
	
	// Check to see if a elemset with the same name has been created before
	if (_elemsets.find(name) != _elemsets.end())
		err_message("Cannot create multiple elemsets with the same name.");
	else
	{
		// Figure out which elements are actually in the local mesh and add them
		std::set<id_type> elem_to_add;
		for(auto it=elems.begin(), end=elems.end(); it!=end; ++it)
			if (elem_in_local_mesh(*it))
				elem_to_add.insert(*it);

		_elemsets[name].swap(elem_to_add);
	}
}


void Mesh::add_to_elemset(std::string name, id_type val)
{
	std::map<std::string, std::set<id_type> >::iterator it = _elemsets.find(name);
	if (it != _elemsets.end())
	{
		if (elem_in_local_mesh(val))
			it->second.insert(val);
	}
	else
		err_message("Attempted to add an element id to an element set that does not exist in the mesh.");
}


// Function to return the given elemset
std::set<id_type> Mesh::get_elemset(id_type idx)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling get_elemset.");
	if (idx>=n_elemsets())
		err_message("Please select a valid elemset.");
	else
	{
		id_type count=0;
		for(auto it=_elemsets.begin(), end=_elemsets.end(); it!=end; ++it, ++count)
			if (count==idx)
				return it->second;
	}
	err_message("We'll never get here!");
}


// Function to return the given elemset
std::set<id_type> Mesh::get_elemset(std::string name)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling get_elemset.");

	if (_elemsets.find(name) != _elemsets.end())
		return _elemsets[name];
	else
		err_message("Please select a valid elemset.");
}


/*
 * Functions to assign materials to elements
 *	1. Assigns the material with the given material_id to every element of the mesh
 *	2. Assigns the material with the given material_name to every element of the mesh
 *	3. Assigns the material with the given material_id to the elements with the given list of global element ids
 *	4. Assigns the material with the given material_name to the elements with the given list of global element ids
 *	5. Assigns the material with the given material_id to the elements in the elemset with the given name
 *	6. Assigns the material with the given material_name to the elements in the elemset with the given name
*/

// 1. Assigns the material with the given material_id to every element of the mesh
void Mesh::set_material(id_type id)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling set_material.");

	// Determine which material id belongs to
	int mat = -1;
	for(id_type i=0; i<_materials.size(); ++i)
	{
		if (id==_materials[i]->get_id())
		{
			mat = i;
			break;
		}
	}
	if (mat<0)
		err_message("Given material id does not match any known material.");
	
	// Add that material to all elements
	_elem_to_material.resize(_elem.size());
	for(id_type i=0; i<_elem.size(); ++i)
		_elem_to_material[i] = mat;

	_materials_assigned = true;
}
// 2. Assigns the material with the given material_name to every element of the mesh
void Mesh::set_material(std::string name)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling set_material.");
	
	// Determine which material name belongs to
	int mat = -1;
	std::transform(name.begin(), name.end(), name.begin(), ::toupper);
	for(id_type i=0; i<_materials.size(); ++i)
	{
		std::string name_test = _materials[i]->get_name();
		std::transform(name_test.begin(), name_test.end(), name_test.begin(), ::toupper);
		if (name_test==name)
		{
			mat = i;
			break;
		}
	}
	if (mat<0)
		err_message("Given material name does not match any known material.");
	
	// Add that material to all elements
	_elem_to_material.resize(_elem.size());
	for(id_type i=0; i<_elem.size(); ++i)
		_elem_to_material[i] = mat;
	
	_materials_assigned = true;
}
// 3. Assigns the material with the given material_id to the elements with the given list of global element ids
void Mesh::set_material(id_type id, std::vector<id_type>& elems)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling set_material.");
	
	// Determine which material id belongs to
	int mat = -1;
	for(id_type i=0; i<_materials.size(); ++i)
	{
		if (id==_materials[i]->get_id())
		{
			mat = i;
			break;
		}
	}
	if (mat<0)
		err_message("Given material id does not match any known material.");
	
	_elem_to_material.resize(_elem.size());
	for(id_type i=0; i<elems.size(); ++i)
		if ( elem_in_local_mesh(elems[i]) )
			_elem_to_material[_global_to_local_elem[elems[i]]] = mat;
	
	// Assigns the boolean variable
	_materials_assigned = true;
	for(id_type i=0; i<_elem_to_material.size(); ++i)
		if (_elem_to_material[i]<0)
			_materials_assigned = false;
		
}
// 4. Assigns the material with the given material_name to the elements with the given list of global element ids
void Mesh::set_material(std::string name, std::vector<id_type>& elems)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling set_material.");
	
	// Determine which material name belongs to
	int mat = -1;
	std::transform(name.begin(), name.end(), name.begin(), ::toupper);
	for(id_type i=0; i<_materials.size(); ++i)
	{
		std::string name_test = _materials[i]->get_name();
		std::transform(name_test.begin(), name_test.end(), name_test.begin(), ::toupper);
		if (name_test==name)
		{
			mat = i;
			break;
		}
	}
	if (mat<0)
		err_message("Given material name does not match any known material.");
	
	_elem_to_material.resize(_elem.size());
	for(id_type i=0; i<elems.size(); ++i)
		if ( elem_in_local_mesh(elems[i]) )
			_elem_to_material[_global_to_local_elem[elems[i]]] = mat;
	
	// Assigns the boolean variable
	_materials_assigned = true;
	for(id_type i=0; i<_elem_to_material.size(); ++i)
		if (_elem_to_material[i]<0)
			_materials_assigned = false;
}
// 6. Assigns the material with the given material_name to the elements in the elemset with the given name
void Mesh::set_material_from_elemset(std::string mat_name, std::string elem_name)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling set_material_from_elemset.");
	
	// Determine which material mat_name belongs to
	int mat = -1;
	for(id_type i=0; i<_materials.size(); ++i)
	{
		std::string name_test = _materials[i]->get_name();
		std::transform(name_test.begin(), name_test.end(), name_test.begin(), ::toupper);
		std::transform(mat_name.begin(), mat_name.end(), mat_name.begin(), ::toupper);
		if (name_test==mat_name)
		{
			mat = i;
			break;
		}
	}
	if (mat<0)
		err_message("Given material name does not match any known material.");

	if (_elemsets.find(elem_name) == _elemsets.end())
		err_message("Please select a valis elemset.");
	else
	{
		_elem_to_material.resize(_elem.size());
		for(auto it=_elemsets[elem_name].begin(), end=_elemsets[elem_name].end(); it!=end; ++it)
			if (elem_in_local_mesh(*it))
				_elem_to_material[_global_to_local_elem[*it]] = mat; // No need to check if the elemnt is in the local mesh as it couldn't be in the elemset otherwise
	}

	// Assigns the materials assigned boolean variable
	_materials_assigned = true;
	for(id_type i=0; i<_elem_to_material.size(); ++i)
		if (_elem_to_material[i]<0)
			_materials_assigned = false;
}






// Function that makes it easy to determine if a vector of global node ids are all in the local mesh.
// Returns true if all ids exist in the local mesh. Returns false if any node id does not exist in the local mesh.
bool Mesh::nodes_in_local_mesh(std::vector<id_type>& nodes)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling nodes_in_local_mesh.");
	
	for(id_type i=0; i<nodes.size(); ++i)
	{
		if ( !node_in_local_mesh(nodes[i]) ) // Couldn't find the global id in the local mesh
			return false;
	}
	return true;
}


// Function that makes it easy to determine if a vector of global node ids are all in the local mesh.
// Returns true if all ids exist in the local mesh. Returns false if any node id does not exist in the local mesh.
bool Mesh::nodes_in_local_mesh(std::set<id_type>& nodes)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling nodes_in_local_mesh.");

	for(auto it=nodes.begin(), end=nodes.end(); it!=end; ++it)
	{
		if ( !node_in_local_mesh(*it) ) // Couldn't find the global id in the local mesh
			return false;
	}
	return true;
}


// Function that makes it easy to determine if a vector of global element ids are all in the local mesh.
// Returns true if all ids exist in the local mesh. Returns false if any element id does not exist in the local mesh.
bool Mesh::elems_in_local_mesh(std::vector<id_type>& elems)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling elems_in_local_mesh.");
	
	for(id_type i=0; i<elems.size(); ++i)
	{
		if ( !elem_in_local_mesh(elems[i]) ) // Couldn't find the global id in the local mesh
			return false;
	}
	return true;
}


// Function that makes it easy to determine if a vector of global element ids are all in the local mesh.
// Returns true if all ids exist in the local mesh. Returns false if any element id does not exist in the local mesh.
bool Mesh::elems_in_local_mesh(std::set<id_type>& elems)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling elems_in_local_mesh.");
	
	for(auto it=elems.begin(), end=elems.end(); it!=end; ++it)
	{
		// if (!_global_to_local_elem.count(*i))
		if ( !elem_in_local_mesh(*it) ) // Couldn't find the global id in the local mesh
			return false;
	}
	return true;
}


// Nice wrapper function to determine if a global elem id exists in the local mesh
bool Mesh::elem_in_local_mesh(id_type id)
{
	if (_global_to_local_elem.find(id)==_global_to_local_elem.end()) // Couldn't find the global id in the local mesh
		return false;
	else
		return true;
}





// Finds all of the partitions that the list of global node ids has in common
std::vector<int> Mesh::partitions_in_common(const std::vector<id_type>& nodes)
{
	if (nodes.size() == 0)
		err_message("Must pass in nodes to find the common partitions between them.");

	for (id_type n=0; n<nodes.size(); ++n)
		if (!node_in_local_mesh(nodes[n]))
			err_message("Attempting to find common partitions of nodes not in the local mesh");

	// Find the set intersection of each of the partition lists for this node
	bool all_on_interface = true;
	for (id_type n=0; n<nodes.size(); ++n)
		if (!node_on_part_interface(nodes[n]))
			all_on_interface = false;

	std::vector<int> intersect;
	if (!all_on_interface)
		intersect.push_back( _rank );
	else
	{
		intersect = _node_partition_interface[nodes[0]];
		std::sort(intersect.begin(), intersect.end());
		std::vector<int> temp_intersect;
		for (id_type n=1; n<nodes.size(); ++n)
		{
			std::vector<int> vec = _node_partition_interface[nodes[n]];
			std::sort(vec.begin(), vec.end());
			std::set_intersection(intersect.begin(), intersect.end(), vec.begin(), vec.end(), std::back_inserter(temp_intersect));
			intersect.swap(temp_intersect);
			temp_intersect.clear();
		}
	}

	return intersect;
}


// Function to determine if the given global node id is located on a partition interface
bool Mesh::node_on_part_interface(id_type id)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling node_on_part_interface.");
	if (!_partitioned) err_message("Must partition the mesh prior to calling node_on_part_interface.");
	
	if (id < ENRICH_START)
	{
		if (_node_partition_interface.find(id)==_node_partition_interface.end())
			return false;
		else
			return true;
	}
	else
	{
		if (_enrich_nodes_interface_lists.find(id)==_enrich_nodes_interface_lists.end())
			return false;
		else
			return true;
	}
}


// Function to determine if a group of nodes are all on a partition interface
bool Mesh::nodes_on_part_interface(const std::vector<id_type>& node_group)
{
	bool on_interface = true;
	for (id_type n=0; n<node_group.size(); ++n)
		if (!node_on_part_interface(node_group[n]))
		{
			on_interface = false;
			break;
		}

	return on_interface;
}

// Get the vector of partitions that have a copy of this node
std::vector<int> Mesh::get_node_parts(id_type id)
{
	if (id < ENRICH_START)
	{
		if (node_on_part_interface(id))
			return _node_partition_interface[id];
		else
			return std::vector<int>(1, _rank);
	}
	else
	{
		if (node_on_part_interface(id))
			return _enrich_nodes_interface_lists[id];
		else
			return std::vector<int>(1, _rank);
	}
}


// Function to determine how many patitions have a copy of this node. Return 1 if the node is not n a partition interface.
id_type Mesh::n_parts_with_node(id_type id) // Doesn't include the current partition
{
	if (!_init) err_message("Must initiailize the mesh prior to calling n_parts_with_node.");
	if (!_partitioned) err_message("Must partition the mesh prior to calling n_parts_with_node.");

	
	if (node_on_part_interface(id))
	{
		if (id < ENRICH_START)
			return _node_partition_interface[id].size();
		else
			return _enrich_nodes_interface_lists[id].size();
	}
	else
		return 1;
}


// Used to iterate through which partitions have a copy of this node. Usually would be called after calling n_parts_with_node(id).
id_type Mesh::get_part_int_from_node(id_type id, id_type interface)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling get_part_int_from_node.");
	if (!_partitioned) err_message("Must partition the mesh prior to calling get_part_int_from_node.");
	
	if (node_on_part_interface(id))
	{
		if (id < ENRICH_START)
		{
			if (interface < _node_partition_interface[id].size())
				return _node_partition_interface[id][interface];
			else
				err_message("Invalid interface selection.");
		}
		else
		{
			if (interface < _enrich_nodes_interface_lists[id].size())
				return _enrich_nodes_interface_lists[id][interface];
			else
				err_message("Invalid interface selection.");
		}
	}
	else
		return -1;
}

id_type Mesh::global_to_local_node(id_type id)
{
	if (id < ENRICH_START) // Normal node
	{
		std::unordered_map<id_type, id_type>::iterator it = _global_to_local_node.find(id);
		if (it != _global_to_local_node.end())
			return it->second;
		else
			err_message("Requested local node index of a node that does not exist in the local mesh.");
	}
	else if ((id-ENRICH_START) < n_global_enrich_nodes())  // Enrichment node
	{
		std::unordered_map<id_type, id_type>::iterator it = _global_to_local_enrich_node.find(id);
		if (it != _global_to_local_enrich_node.end())
			return it->second;
		else
			err_message("Requested local node index of a node that does not exist in the local mesh.");
	}
	else
		err_message("Requested local node index of a node that does not exist in the global mesh.");
}
id_type Mesh::global_to_local_elem(id_type id)
{
	std::unordered_map<id_type, id_type>::iterator it = _global_to_local_elem.find(id);
	if (it != _global_to_local_elem.end())
		return it->second;
	else
	{
		std::stringstream ss;
		ss << id;
		std::string out = "Requested global element id " + ss.str() + " which does not exist in the local mesh.";
		err_message( out.data() );
	}
}






/*
 * Create the _node_elem table from the _elem_node table
 * NOTE: Whereas the _elem_node table is local elem numbers to local node numbers for conveniences sake,
 *		 the _node_elem table is local node numbers to GLOBAL elem numbers (stores as the elem id's) because they can span across partitions
*/
void Mesh::generate_node_elem()
{
	if (!_init) err_message("Must initiailize the mesh prior to calling generate_node_elem.");
	if (!_partitioned) err_message("Must partition the mesh prior to calling generate_node_elem.");

	// Clear the date structures. Will be regenerated for a refined mesh
	_node_elem.clear();
	_remote_node_elem.clear();

	// Handle serial case, much simpler
	if (_serial)
	{
		int n_nodes = n_local_nodes();
		for(int i=0; i<n_nodes; ++i)
		{
			_node_elem.push_back(std::vector<id_type>());
		}
		for(element_iterator it=active_elements_begin(), end=active_elements_end(); it!=end; ++it)
		{
			id_type l_elem = global_to_local_elem((*it)->get_id());
			int neln = _elem_node[l_elem].size();
			for(int n=0; n<neln; ++n)
			{
				_node_elem[_elem_node[l_elem][n]].push_back((*it)->get_id());
			}
		}

		_node_elem_generated = true;
		return;
	}
	
	/*
	Steps to actually generate the _node_elem table
	1. Loop over all the current elements in the elem_node table and add an entry for that element to each of its nodes
	2. For any nodes that are on an interface (use _node_partition_interface map), send the current _node_elem entry to all other partitions that have that node
	3. Recieve the corresponding lists from all the other partitions that have that node.
		- Would it be better to store remote elements separetely with some method for determining which partition they're on?
		- For communication purposes, since the number of elements for each node is unknown, some form of
			ptr/ind structure will have to be set up that also includes the (global) node numbers
	*/
	
	// 1. Loop over all the current elements in the elem_node table and add an entry for that element to each of its nodes
	int n_nodes = n_local_nodes();
	for(int i=0; i<n_nodes; ++i)
	{
		_node_elem.push_back(std::vector<id_type>());
	}
	for(element_iterator it=active_elements_begin(), end=active_elements_end(); it!=end; ++it)
	{
		id_type l_elem = global_to_local_elem((*it)->get_id());
		int neln = _elem_node[l_elem].size();
		for(int n=0; n<neln; ++n)
		{
			_node_elem[_elem_node[l_elem][n]].push_back((*it)->get_id());
		}
	}
	
	// 2. For any nodes that are on an interface (use _node_partition_interface map), send the current _node_elem entry to all other partitions that have that node
	int nranks;
	MPI_Comm_size(_communicator, &nranks);
	std::vector<std::vector<id_type> > neptr_send_lists;
	std::vector<std::vector<id_type> > neind_send_lists;
	std::unordered_map<id_type, std::vector<int> >::iterator j = _node_partition_interface.begin();
	std::unordered_map<id_type, std::vector<int> >::iterator end = _node_partition_interface.end();
	for(int i=0; i<nranks; ++i)
	{
		neptr_send_lists.push_back(std::vector<id_type>());
		neind_send_lists.push_back(std::vector<id_type>());
	}
	std::vector<int> curr_inds(nranks,0); // This can also act as my n_expected (if the entry is 0, then I'm not sending any nodes ot that partition)
	for(;j!=end;++j)
	{
		id_type gnode = (*j).first;
		id_type lnode = _global_to_local_node[gnode];
		for(id_type k=0; k<(*j).second.size(); ++k) // Loop over partitions node exists on
		{
			int part = (*j).second[k];
			neptr_send_lists[part].push_back(gnode);
			neptr_send_lists[part].push_back(curr_inds[part]);
			for(id_type l=0; l<_node_elem[lnode].size(); ++l)
			{
				neind_send_lists[part].push_back(_node_elem[lnode][l]);
				curr_inds[part]++;
			}
			neptr_send_lists[part].push_back(curr_inds[part]);
		}
	}
	// Send all of the node_elem info
	int n_sends = 0;
	for(int i=0; i<nranks; ++i)
	{
		if ((i!=_rank) && (curr_inds[i] != 0))
			++n_sends;
	}
	MPI_Request reqs[n_sends*2];
	int n_sent = 0;
	for(int i=0; i<nranks; ++i)
	{
		if ((i!=_rank) && (curr_inds[i] != 0))
		{
			MPI_Isend(&neptr_send_lists[i][0], neptr_send_lists[i].size(), MPI_ID, i, 0, _communicator, &reqs[2*n_sent]);
			MPI_Isend(&neind_send_lists[i][0], neind_send_lists[i].size(), MPI_ID, i, 1, _communicator, &reqs[2*n_sent+1]);
			n_sent++;
		}
	}
	
	// 3. Recieve the corresponding lists from all the other partitions that have that node.
	MPI_Status status[2];
	for(int i=0; i<nranks; ++i)
	{
		if ((i!=_rank) && (curr_inds[i]!=0))
		{
			// Recieve the node_elem lists
			int count0, count1;
			std::vector<id_type> neptr_recv_vec;
			std::vector<id_type> neind_recv_vec;
			MPI_Probe(i, 0, MPI_COMM_WORLD, &status[0]);
			MPI_Probe(i, 1, MPI_COMM_WORLD, &status[1]);
			MPI_Get_count(&status[0], MPI_ID, &count0);
			MPI_Get_count(&status[1], MPI_ID, &count1);
			neptr_recv_vec.resize(count0);
			neind_recv_vec.resize(count1);
			MPI_Recv(&neptr_recv_vec[0], count0, MPI_ID, i, 0, _communicator, &status[0]);
			MPI_Recv(&neind_recv_vec[0], count1, MPI_ID, i, 1, _communicator, &status[1]);
			
			// Put this node_elem info in my data structure
			for(int j=0; j<count0; j=j+3)
			{
				id_type node_id = neptr_recv_vec[j];
				id_type start  = neptr_recv_vec[j+1];
				id_type e  = neptr_recv_vec[j+2];
				for(id_type k=start; k<e; ++k)
				{
					if (_remote_node_elem.find(node_id)==_remote_node_elem.end()) // Don't have any remote info for this node yet
					{
						std::vector<id_type> v(1);
						v[0] = neind_recv_vec[k];
						std::map<int, std::vector<id_type> > tmp_map;
						tmp_map.insert(std::pair<int, std::vector<id_type> >(i, v));
						_remote_node_elem.insert(std::pair<id_type, std::map<int, std::vector<id_type> > >(node_id, tmp_map));
					}
					else
					{
						if (_remote_node_elem[node_id].find(i)==_remote_node_elem[node_id].end()) // Have some info for the node but not for this partition
						{
							std::vector<id_type> v(1);      
							v[0] = neind_recv_vec[k];
							_remote_node_elem[node_id].insert(std::pair<id_type, std::vector<id_type> >(i, v));
						}
						else
							_remote_node_elem[node_id][i].push_back(neind_recv_vec[k]);
					}
				}
			}
		}
	}
	
	_node_elem_generated = true;
	MPI_Waitall(n_sends*2, reqs, MPI_STATUSES_IGNORE);
}











// EVERYTHING FROM HERE ON OUT IS SPECIFIC TO IGFEM
// ============================================================================================================================

// Function which resets the mesh to the Pre-IGFEM state
void Mesh::clearEnrichments()
{
	// IGFEM stuff
	for(id_type n=0; n< _enrich_nodes.size(); ++n)
		delete _enrich_nodes[n];
	for (Mesh::element_iterator it=elements_begin(), end=elements_end(); it!=end; ++it)
		(*it)->clearEnrichments();
	_enrich_nodes.clear();
	_global_to_local_enrich_node.clear();
	_enrich_owners.clear();
	_node_detect.clear();
	_n_global_enrich_nodes = 0;
	_n_local_owned_enrich_nodes = 0;
	_enrich_nodes_interface_lists.clear();
	_enrich_node_elem.clear();
	_nodes_detected = false;
	_elements_detected = false;
}

void Mesh::add_cohesive_material_pair(std::string coh_mat, std::string mat1, std::string mat2)
{
	int coh_mat_num = find_mat_number(coh_mat);
	coh_mat_num = std::find(_cohesive_mats.begin(), _cohesive_mats.end(), coh_mat_num) - _cohesive_mats.begin();
	int mat1_num = find_mat_number(mat1);
	int mat2_num = find_mat_number(mat2);

	if (coh_mat_num==-1 || mat1_num==-1 || mat2_num==-1)
		err_message("Material name does not match a material. Please select a valid material.");
	else
	{
		if (mat1_num < mat2_num)
			_cohesive_to_mat_pair[coh_mat_num] = std::pair<int, int>(mat1_num, mat2_num);
		else
			_cohesive_to_mat_pair[coh_mat_num] = std::pair<int, int>(mat2_num, mat1_num);
	}
}


void Mesh::add_inclusion(Inclusion* inc)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling add_inclusion.");
	if (_nodes_detected) err_message("Cannot add any more inclusions after detecting nodes.");

	_IGFEM = true;

	// Check to see if the material already exists in the mesh
	Material * outside_mat = inc->get_material();
	Material * mesh_mat_ptr = NULL;
	std::string name = outside_mat->get_name();
	int mat = -1;
	std::transform(name.begin(), name.end(), name.begin(), ::toupper);
	for(id_type i=0; i<_materials.size(); ++i)
	{
		std::string name_test = _materials[i]->get_name();
		std::transform(name_test.begin(), name_test.end(), name_test.begin(), ::toupper);
		if (name_test==name)
		{
			mat = i;
			break;
		}
	}
	if (mat<0) // Didn't find a valid matching material so we have to make one and add it to the mesh
	{
		add_material(outside_mat);
		mesh_mat_ptr = _materials[n_materials()-1];
	}
	else		// Found a matching material so I'll use that pointer to in the inclusion
		mesh_mat_ptr = _materials[mat];

	inc->set_id(n_inclusions());
	_inclusions.push_back(inc->allocate_and_copy());
	_inclusions[n_inclusions()-1]->set_material(mesh_mat_ptr);
}

Inclusion* Mesh::get_inclusion(id_type idx)
{
	if (!_init) err_message("Must initiailize the mesh prior to calling get_inclusion.");

	if (idx < n_inclusions())
		return _inclusions[idx];
	else
		err_message("Inclusion index does not match an inclusion. Please select a valid inclusion.");
}

// Takes the local node number of an enrichment node and returns a pair of global node numbers that make the ede this enrichment node lies on
std::pair<id_type, id_type> Mesh::get_enriched_edge_local(id_type idx)
{
	if (idx >= ENRICH_START)
	{
		if ((idx-ENRICH_START) < n_local_enrich_nodes())
			return _enrich_node_to_edge[idx-ENRICH_START];
		else
			err_message("Attempting to access the enriched edge of a node not existing in the local mesh");
	}
	else
		err_message("Can only find enriched edges of enrichment nodes.");
}

// Takes the local node number of an enrichment node and returns the local inclusion number 
id_type Mesh::get_enriched_inclusion_local(id_type idx)
{
	if (idx >= ENRICH_START)
	{
		if ((idx-ENRICH_START) < n_local_enrich_nodes())
			return _enrich_node_to_inclusion[idx-ENRICH_START];
		else
			err_message("Attempting to access the enriched inclusion of a node not existing in the local mesh");
	}
	else
		err_message("Can only find enriched inclusion of enrichment nodes.");
}


int Mesh::get_nodal_detection_local(id_type idx)
{
	if (idx < n_local_nodes())
		return _node_detect[idx];
	else
		err_message("Attempting to access the nodal detection value of a node which does not exist on the local mesh.");
}

// Analyze all nodes and determine which inclusion they lie in
void Mesh::analyze_nodes()
{
	if (!IGFEM())
		return;

	// Serial, easier case
	if (_serial)
	{
		_node_detect.resize( n_local_nodes() );

		// Loop through all nodes and determine which inclusion they lie in, if any
		for(id_type n=0; n<n_local_nodes(); ++n)
			_node_detect[n] = detect_node(_nodes[n]);
	}

	else
	{
		_node_detect.resize( n_local_nodes() );
		
		// This is a container that is filled with nodes that are on a partition interface and also lie on an inclusion surface. Very unlucky
		// But they need to be moved by the process that ownes them.
		std::unordered_map<id_type, Inclusion*> nodes_to_reprocess;
		
		// Loop through all nodes and determine which inclusion they lie in, if any
		for(id_type n=0; n<n_local_nodes(); ++n)
		{
			_node_detect[n] = detect_node(_nodes[n]);

			// If the node lies on an inclusion surface and is also on a partition interface we must move it later
			if (_node_detect[n] >= (int)n_inclusions())
				nodes_to_reprocess[n] = _inclusions[_node_detect[n] - n_inclusions()];
		}
		
		// Function to move nodes that are on an inclusion surface and also on a partition interface
		reprocess_nodes(nodes_to_reprocess);
	}

	_nodes_detected = true;
}

// Moves node on an inclusion surface and also on a partition interface
// local node number and a pointer to the incusion it's on
void Mesh::reprocess_nodes(std::unordered_map<id_type, Inclusion*>& nodes_to_reprocess)
{
	// Now we can handle those unfortunate nodes that are on a partition interface and also on an inclusion surface
	// Node id lists will contain the global node ids as well as the eventual detection id for each node.
	// Coords send vector will contain the coordinates of each node. For every 3 entries in the coords vec, there will be 2 in the ids vec.
	std::map<int, std::vector<int> > node_id_send_lists;
	std::map<int, std::vector<double> > node_coord_send_lists;
	std::map<int, int> n_expected;
	for(auto it=nodes_to_reprocess.begin(), end=nodes_to_reprocess.end(); it!=end; ++it)
	{
		id_type n = (*it).first;
		Node* node = _nodes[n];
		id_type id = node->get_id();
		Inclusion* inc = (*it).second;
		int times_moved = 0;
		// If I own this node, then modify it and add the new location to the send lists of all partitions that have a copy of it
		if (get_node_owner_local(n)== _rank)
		{
			// Try to modify the node location and check to see if was moved onto another partition
			id_type ret;
			while(times_moved < 5)
			{
				modify_node_location(node, inc);
				ret = detect_node(node);
				if ((ret>=n_inclusions()) && (ret!=2*n_inclusions()))
					times_moved++;
				else
					break;
			}
			// If it was moved too many times then just throw an error, this thing won't be fixed
			if (times_moved==5)
			{
				char buf[100], buf2[100];
				strcpy(buf, "Unable to move node %");
				strcat(buf, SPEC);
				strcat(buf, " off of an inclusion surface");
				sprintf(buf2, buf, id);
				err_message( buf2 );
			}
			// Otherwise, if the node is on a partition interface pack the node into the buffers and get ready to send them
			if (_node_partition_interface.find(id) != _node_partition_interface.end())
			{
				for(id_type i=0; i<_node_partition_interface[id].size(); ++i)
				{
					int part = _node_partition_interface[id][i];
					if (part != _rank)
					{
						// Make sure I at least have empty vectors to insert into
						if (node_id_send_lists.find(part)==node_id_send_lists.end())
						{
							node_id_send_lists[part] = std::vector<int>();
							node_coord_send_lists[part] = std::vector<double>();
						}
						node_id_send_lists[part].push_back( id );
						node_id_send_lists[part].push_back( ret );
						node_coord_send_lists[part].push_back( (*node)(0) );
						node_coord_send_lists[part].push_back( (*node)(1) );
						node_coord_send_lists[part].push_back( (*node)(2) );
					}
				}
			}
		}
		// Otherwise, I need to add it to my n_expected array
		else
		{
			if (_node_partition_interface.find(id) != _node_partition_interface.end())
				n_expected[get_node_owner_local(n)]++;
		}
	}
	
	// Now send those lists
	int n_sends = node_id_send_lists.size()*2;
	MPI_Request reqs[n_sends];
	int n_sent = 0;
	for(auto it=node_id_send_lists.begin(), end = node_id_send_lists.end(); it!=end; ++it)
	{
		int part = (*it).first;
		MPI_Isend(&(*it).second[0], (*it).second.size(), MPI_ID, part, 0, _communicator, &reqs[n_sent*2]);
		MPI_Isend(&node_coord_send_lists[part][0], node_coord_send_lists[part].size(), MPI_DOUBLE, part, 1, _communicator, &reqs[n_sent*2 + 1]);
		n_sent++;
	}
	
	// And now recieve them
	MPI_Status status[2];
	for(auto it=n_expected.begin(), end=n_expected.end(); it!=end; ++it)
	{
		int part = (*it).first;
		int count0, count1;
		std::vector<int> id_recv_vec;
		std::vector<double> coord_recv_vec;
		MPI_Probe(part, 0, MPI_COMM_WORLD, &status[0]);
		MPI_Probe(part, 1, MPI_COMM_WORLD, &status[1]);
		MPI_Get_count(&status[0], MPI_ID, &count0);
		MPI_Get_count(&status[1], MPI_DOUBLE, &count1);
		if (count0/2 != (*it).second)
			err_message("Number of inclusion surface/partition interface nodes expected did not match the number recieved.");
		id_recv_vec.resize(count0);
		coord_recv_vec.resize(count1);
		MPI_Recv(&id_recv_vec[0], count0, MPI_ID, part, 0, _communicator, &status[0]);
		MPI_Recv(&coord_recv_vec[0], count1, MPI_DOUBLE, part, 1, _communicator, &status[1]);
		
		// Now actually use the info recieved to modify the nodes
		for(int n=0; n<(count0/2); ++n)
		{
			int idx = n*2;
			int idx2 = n*3;
			id_type id = id_recv_vec[idx];
			id_type node = _global_to_local_node[id];
			int inc = id_recv_vec[idx+1];
			double x = coord_recv_vec[idx2];
			double y = coord_recv_vec[idx2+1];
			double z = coord_recv_vec[idx2+2];
			_nodes[node]->set_coords(x, y, z);
			_node_detect[node] = inc;
			
		}
	}
}

// Helper function for nodal detection
int Mesh::detect_node(Node* node)
{
	// Loop through all of the inclusions in the mesh, checking each one
	// We loop through the inclusions in reverse order because we will allow inclusions to be placed
	// on top of other incusions and thus take precedence over the inclusios that were place earlier.
	// For example: an ellipsoidal inclusion placed within another ellpsoidal inclusion could be used to represent some sort of "cored" structure
	int times_moved = 0;
	for(int j=(n_inclusions()-1); j>=0; --j)
	{
		// Get the limits on the size of the inclusion
		std::vector<double> dl = _inclusions[j]->get_domain_lims();
		// If this node does a least lie within the domain limits
		if (dl[0]<=(*node)(0) && dl[1]<=(*node)(1) && dl[2]<=(*node)(2) && (*node)(0)<=dl[3] && (*node)(1)<=dl[4] && (*node)(2)<=dl[5])
		{
			// Determining if a node is within an inclusion could vary from inclusion
			// to inclusion so this will be left up to the individual inclusion objects
			id_type res = _inclusions[j]->detect_node(node);
			if (res == 0)  // Node does not fall within inclusion
				continue; // point is outside this inclusion, checking point in the next inclusion
				
			// Node is on the inclusion surface
			else if (res == 1)
			{
				// If the node is not on a partition interface then just move it now (Also should do this in serial)
				if (_node_partition_interface.find(node->get_id()) == _node_partition_interface.end())
				{
					// This has happened too many times, we're just gonna throw an error
					if (times_moved >= 5)
					{
						char buf[100], buf2[100];
						strcpy(buf, "Unable to move node %");
						strcat(buf, SPEC);
						strcat(buf, " off of an inclusion surface");
						sprintf(buf2, buf, node->get_id());
						err_message( buf2 );
					}
					modify_node_location(node, _inclusions[j]);
					++times_moved;
					j = n_inclusions(); // Used to restart the loop over the inclusions
				}
				else
					return n_inclusions() + j;
				
			}
			else // Node is contained in inclusion
			{
				return j;
			}
		}
	}
	
	// If the node was not found to lie within any of the inclusions, then it must lie in a matrix material
	return -1;
}


// Helper function for actually moving the Nodal position
void Mesh::modify_node_location(Node* node, Inclusion* inc)
{
	id_type l_node = global_to_local_node(node->get_id());
	
	// First get the outward facing unit normal to the surface at the given node location
	// This is the direction we weill move along.
	std::vector<double> surf_norm = inc->get_surface_normal(node);
	
	// Next, I need a characteristic length scale for the elements that this node belongs to
	// To do this, I will pick one of the element it beongs to and one of that elemnts nodes and find the distance to that node
	// FIXME: I should then probably check to see if all element Jacobians are still positive...
	// NOTS: This method assumes that all elements that have this node are of similar size
	Node* other = _elem[ _global_to_local_elem[ _node_elem[l_node][0] ] ]->get_node(0);
	if (node == other) // If I happened to pick the same node, just pick the next one instead
		other = _elem[ _global_to_local_elem[ _node_elem[l_node][0] ] ]->get_node(1);
	
	// Get the distance from node to other and scale it by 5%
	double move_dist = 0.05*( sqrt(pow(((*node)(0) - (*other)(0)), 2) + pow(((*node)(1) - (*other)(1)), 2) + pow(((*node)(2) - (*other)(2)), 2)) );
	
	// Actually move the node position
	node->set_coords(((*node)(0) + move_dist*surf_norm[0]), ((*node)(1) + move_dist*surf_norm[1]), ((*node)(2) + move_dist*surf_norm[2]));
	
	// Here is where I would check to see if the Jacobians were still positive...
}


// Analyze all elements and determine if/how they are intersected by the inclusions
void Mesh::analyze_elements()
{
	if (!IGFEM())
		return;

	if (!_nodes_detected) err_message("Must detect nodes the mesh prior to calling analyze_elements.");

	// Same in serial and parallel
	// Build the vector of inclusions materials
	std::vector<Material*> mats;
	std::vector<id_type> mat_nums;
	mats.push_back(NULL); // Placeholder for the material already contained in the element
	mat_nums.push_back( 0 );
	for(id_type i=0; i<n_inclusions(); ++i)
	{
		Material* mat = _inclusions[i]->get_material();
		mats.push_back(mat);
		mat_nums.push_back( find_mat_number(mat->get_name()) );
	}

	for(Mesh::element_iterator it=active_elements_begin(), end=active_elements_end(); it!=end; ++it)
	{
		// Get the local element number
		id_type l_elem = global_to_local_elem((*it)->get_id());

		// Add the subdomain material that this element contains (Makes sure integration element materials get assigned appropriately)
		mats[0] = get_element_material_local(l_elem);
		mat_nums[0] = find_mat_number( mats[0]->get_name() );

		// First build up the list of nodal detection values
		std::vector<int> nodal_detection((*it)->n_nodes());
		for(id_type n=0; n<nodal_detection.size(); ++n)
			nodal_detection[n] = _node_detect[_elem_node[l_elem][n]];
		
		// Perform the element detection
		(*it)->detection(nodal_detection, mats, is_cohesive());

		// Update the element's material if it contained entirely inside an inclusions
		if (!(*it)->is_intersected())
		{
			int inc_num = nodal_detection[0]; // All detection values should be the same
			if (inc_num != -1) // Entirely contained in an inclusion
				_elem_to_material[l_elem] = mat_nums[inc_num + 1];
		}

		if ((*it)->is_intersected())
		{
			_intersected_elem.push_back(l_elem);

			// If the mesh is cohesive and we've gotten this far then I need to add the cohesive material to the element
			if (is_cohesive())
			{
				// Get the cohesive element structure
				std::vector<std::vector<short_id_type> >& coh_elems = (*it)->getCohesiveElemStructure();
				id_type n_nodes = (*it)->n_nodes();

				for (id_type ce=0; ce<coh_elems.size(); ++ce)
				{
					// Get the inclusion number associated with this cohesive element
					id_type inclusion = (*it)->get_inclusion_from_enrich(coh_elems[ce][0] - n_nodes); // Any node on cohesive element will do since they're all enrichment nodes
					id_type inc_mat_num = mat_nums[inclusion + 1];
					id_type matrix_mat_num = mat_nums[0];

					// Create a pair with the material associated with the cohesive element and this element's matrix material
					// NOTE: this makes the assumption that inclusions do not intersect within an element (cored inclusion ide no longer works)
					std::pair<int, int> cohesive_pair;
					if (matrix_mat_num < inc_mat_num)
						cohesive_pair = std::pair<int, int>(matrix_mat_num, inc_mat_num);
					else
						cohesive_pair = std::pair<int, int>(inc_mat_num, matrix_mat_num);

					// Find the corrcet cohesive material (If no matching pair is found and there's a default cohesive material (-1,-1), the default material will be used)
					std::pair<int, int> def(-1,-1);
					int coh_mat_num = -1;
					int default_coh_mat = -1;
					Material* coh_mat = NULL;
					for(id_type cm=0; cm<_cohesive_to_mat_pair.size(); ++cm)
					{
						if (cohesive_pair == _cohesive_to_mat_pair[cm])
						{
							coh_mat_num = _cohesive_mats[cm];
							break;
						}
						if (default_coh_mat==-1 && _cohesive_to_mat_pair[cm]==def)
							default_coh_mat = _cohesive_mats[cm];
					}
					if (coh_mat_num!=-1 || default_coh_mat!=-1)
					{
						if (coh_mat_num!=-1) // found a matching cohesive material
							coh_mat = _materials[coh_mat_num];
						else 				// using the first default material in the mesh (-1,-1)
							coh_mat = _materials[default_coh_mat];
					}
					else
						err_message("No matching cohesive material found for an element.");

					// Add the cohesive material to the element
					(*it)->set_cohesive_material(coh_mat, ce);

				} // End cohesive element loop
			} // End is cohesive
		} // End is intersected
	} // End element iteratr

	_elements_detected = true;
}






// Function to add all enrichment nodes to the mesh, keeping consistent across processors
// Loop through all elements and add their cut edges to a set of node-pairs.
// This ensures that only one enrichment node is created for each unique edge.
// Then find all the intersections, create the nodes, and add the enrichment nodes to the elements
void Mesh::add_enrichments()
{
	if (!IGFEM())
		return;

	if (!_elements_detected) err_message("Must detect elements the mesh prior to calling add_enrichments.");

	// Steps that this function takes
	// 1. Compiles a list of all unique node-pairs that constitute edges cut by an inclusion interface
	//	-One list of all edges that are not on a partition interface and another list
	//	 of all edges that are on a partition interface. These will require special handling
	// 2. Find the actual intersection points of the inclusions with the given edges.
	//	-Do this for all edges not on a partition interface as well as those on a
	//	 partition interface where I own the node with the smaller global ID of the pairs
	//	-Create new nodes at these points
	// 3. Assign global enrichment node numbers
	//	-Do a parallel prefix sum to determine number of enrichment node on processors before
	//	 mine and then label similar to dofs
	// 4. Communicate partition interface nodes to all other partitions that have copies of both of the edge's ends.
	// 5. Add the node pointers to the elements
	//	-Make sure order of addition matches the order of the _cutEdge vector for each element
	// 6. Build the integration elements for each element

	// Ok, this function is really dificult to do a separate serial section and parallel section but for the sake of consistency
	// I guess we'll keep it as separate sections and not just add some if statements in the middle. Results in a lot of
	// code duplication but oh well.
	if (_serial)
	{
		std::set<std::pair<std::pair<id_type, id_type>, id_type> > cut_edges;
		std::unordered_map<std::pair<id_type, id_type>, std::pair<id_type, id_type> > parent_to_refined_edge;
	
		// 1. Compiles a list of all unique node-pairs that constitute edges cut by an inclusion interface
		// Loop over all elements and add their cut edges to the sets
		for(id_type e=0; e<_intersected_elem.size(); ++e)
		{
			Elem* el = _elem[_intersected_elem[e]];
			std::vector<short_id_type>& int_edges = el->getIntersectedEdges();
			std::vector<id_type>& inc_nums = el->getInclusionNumbers();
			
			// Loop over all edges, creating pairs of global IDs for each edge and adding it to the appropriate set
			for(id_type i=0; i<int_edges.size(); ++i)
			{
				// Should I figure out some way of extending this to quadratic edges...?
				id_type n0 = (*el)(el->edge_nodes_map(int_edges[i], 0)).get_id();
				id_type n1 = (*el)(el->edge_nodes_map(int_edges[i], 1)).get_id();
				std::pair<id_type, id_type> edge;
				if (n0<n1)
					edge = std::make_pair(n0, n1);
				else
					edge = std::make_pair(n1, n0);

				auto it0 = _edge_hanging_nodes.find(n0);
				auto it1 = _edge_hanging_nodes.find(n1);
				if (it0 != _edge_hanging_nodes.end() || it1 != _edge_hanging_nodes.end())
				{
					// Dont defer if both are hanging nodes this means that the edge belongs only to a refined element?)
					if (it0 == _edge_hanging_nodes.end() || it1 == _edge_hanging_nodes.end())
					{
						// Figure out which is the hanging node and get the corresponding neighbors
						id_type hang, other;
						if (it0 != _edge_hanging_nodes.end()) // Node 0 is the hanging node
							{hang = n0; other = n1;}
						else
							{hang = n1; other = n0;}
						// Get the neighbors of this hanging node
						std::vector<id_type>& neighbors = get_hanging_node_neighbors( hang );

						// If this is an edge hanging node then I need to check if it needs to be substituted with the parent edge
						if (neighbors.size() == 2)
						{
							// If the other node from the refined edge is in the neighbors then this parent edge is the same as the refined edge
							// Otherwise, this is an edge that belongs only to a refined element
							if (std::find(neighbors.begin(), neighbors.end(), other) != neighbors.end()) // std::find here is fine since its just a two elment long vector
							{
								std::pair<id_type, id_type> refined_edge = edge;
								edge = {neighbors[0], neighbors[1]}; // Should already be sorted

								// Add entry that maps parent edges to the refined edges for hanging node cases
								parent_to_refined_edge[edge] = refined_edge;
							}
						}
					}
				}

				cut_edges.insert( std::pair<std::pair<id_type, id_type>, id_type>(edge, inc_nums[i]) );
			}
		}

		// 2. Find the actual intersection points of the partition interfaces with the given edges.
		//	-Create new nodes at these points
		//	-Create a map from the local enrichment node numbers back to the edge it belongs to
		_enrich_nodes.resize(cut_edges.size());
		_enrich_owners.resize(cut_edges.size());
		_enrich_node_to_edge.resize(cut_edges.size());
		_enrich_node_to_inclusion.resize(cut_edges.size());
		id_type count = 0;
		for(auto it=cut_edges.begin(), end=cut_edges.end(); it!=end; ++it)
		{
			std::pair<id_type, id_type> edge = it->first;
			id_type inclusion = it->second;
			// Form a vector of nodes that define the edge that is intersected. This can later be extended to 3+ nodes
			std::vector<Node> edge_nodes(2);
			edge_nodes[0] = *_nodes[ _global_to_local_node[edge.first] ];
			edge_nodes[1] = *_nodes[ _global_to_local_node[edge.second] ];
			
			// Actually find the intersection point and return a new node
			if (inclusion >= n_inclusions() || inclusion < 0)
				std::cout << "Error impending here.\n";
			std::vector<double> coords = _inclusions[inclusion]->find_intersection(edge_nodes);
			_enrich_nodes[count] = new Node(coords[0], coords[1], coords[2], 0); // Don't worry about global ID yet, we'll worry about that next
			_enrich_owners[count] = _rank;
			_enrich_node_to_edge[count] = edge;
			_enrich_node_to_inclusion[count] = inclusion;
			count++;
		}
		// If this mesh is cohesive, then need to add the mirror nodes to the mesh
		id_type n_initial = _enrich_nodes.size();
		if (is_cohesive())
		{
			_enrich_nodes.resize( n_initial * 2 );
			_enrich_owners.resize( n_initial * 2 );
			_enrich_node_to_edge.resize( n_initial * 2 );
			_enrich_node_to_inclusion.resize( n_initial * 2 );
			for(id_type n=0; n<n_initial; ++n)
			{
				Node* copy = _enrich_nodes[n];
				_enrich_nodes[n+n_initial] = new Node((*copy)(0), (*copy)(1), (*copy)(2), 0);
				_enrich_owners[n+n_initial] = _rank;
				_enrich_node_to_edge[n+n_initial] = _enrich_node_to_edge[n];
				_enrich_node_to_inclusion[n+n_initial] = _enrich_node_to_inclusion[n];
			}
		}

		// 3. Assign global enrichment node numbers
		_n_local_owned_enrich_nodes = _enrich_nodes.size();
		_n_global_enrich_nodes = _n_local_owned_enrich_nodes;
		// Now give all of them the actual numbering
		for(id_type local = 0; local<_enrich_nodes.size(); ++local)
		{
			_enrich_nodes[local]->set_id( local + ENRICH_START );
			_global_to_local_enrich_node[local + ENRICH_START] = local + ENRICH_START; // This is really kinda useless in serial but in parallel its useful
		}

		// 5. Add the node pointers to the elements
		//	-Make sure order of addition matches the order of the _cutEdge vector for each element
		_enrich_node_elem.resize(n_local_enrich_nodes());
		for(id_type n=0; n<n_initial; ++n)
		{
			std::pair<id_type, id_type> edge = _enrich_node_to_edge[n];
			id_type inclusion = _enrich_node_to_inclusion[n];
			add_enrich_node_to_elem(edge, inclusion, n, n_initial);

			// If this edge contains a hanging node, then add the enrichment node to the refined edge as well
			auto it = parent_to_refined_edge.find( edge );
			if (it != parent_to_refined_edge.end())
				add_enrich_node_to_elem(it->second, inclusion, n, n_initial);
		}

		// 6. Build the integration elements for each element
		for(element_iterator it=active_elements_begin(), end=active_elements_end(); it!=end; ++it)
			(*it)->generate_integration_elem();

		// Assign the materials added boolean because IGFEM takes care of that by having each element have a pointer to its material
		_materials_assigned = true;
	}
	



	// Now do parallel version, significantly more complicated
	// ---------------------------------------------------------------------------------
	else
	{
		std::set<std::pair<std::pair<id_type, id_type>, id_type> > cut_edges;
		std::set<std::pair<std::pair<id_type, id_type>, id_type> > cut_edges_on_interface;
		std::unordered_map<std::pair<id_type, id_type>, std::vector<int> > edge_int_map; // Maps an edge to what partitions actually have that edge (defined by have at least one element that shares both of the nodes on the edge)
		std::unordered_map<std::pair<id_type, id_type>, std::pair<id_type, id_type> > parent_to_refined_edge;
		
		
		// 1. Compiles a list of all unique node-pairs that constitute edges cut by an inclusion interface
		//	-One list of all edges that are not on a partition interface and another list
		//	 of all edges that are on a partition interface. These will require special handling
		// Loop over all elements and add their cut edges to the sets
		for(id_type e=0; e<_intersected_elem.size(); ++e)
		{
			Elem* el = _elem[_intersected_elem[e]];
			std::vector<short_id_type>& int_edges = el->getIntersectedEdges();
			std::vector<id_type>& inc_nums = el->getInclusionNumbers();
			
			// Loop over all edges, creating pairs of global IDs for each edge and adding it to the appropriate set
			for(id_type i=0; i<int_edges.size(); ++i)
			{
				// Should I figure out some way of extending this to quadratic edges...?
				id_type n0 = (*el)(el->edge_nodes_map(int_edges[i], 0)).get_id();
				id_type n1 = (*el)(el->edge_nodes_map(int_edges[i], 1)).get_id();
				std::pair<id_type, id_type> edge;
				if (n0<n1)
					edge = std::make_pair(n0, n1);
				else
					edge = std::make_pair(n1, n0);

				// Check to see if this edge contains a hanging node
				auto it0 = _edge_hanging_nodes.find(n0);
				auto it1 = _edge_hanging_nodes.find(n1);
				if (it0 != _edge_hanging_nodes.end() || it1 != _edge_hanging_nodes.end())
				{
					// Dont defer if both are hanging nodes this means that the edge belongs only to a refined element?)
					if (it0 == _edge_hanging_nodes.end() || it1 == _edge_hanging_nodes.end())
					{
						// Figure out which is the hanging node and get the corresponding neighbors
						id_type hang, other;
						if (it0 != _edge_hanging_nodes.end()) // Node 0 is the hanging node
							{hang = n0; other = n1;}
						else
							{hang = n1; other = n0;}
						// Get the neighbors of this hanging node
						std::vector<id_type>& neighbors = get_hanging_node_neighbors( hang );

						// If this is an edge hanging node then I need to check if it needs to be substituted with the parent edge
						if (neighbors.size() == 2)
						{
							// If the other node from the refined edge is in the neighbors then this parent edge is the same as the refined edge
							// Otherwise, this is an edge that belongs only to a refined element
							if (std::find(neighbors.begin(), neighbors.end(), other) != neighbors.end()) // std::find here is fine since its just a two elment long vector
							{
								std::pair<id_type, id_type> refined_edge = edge;
								edge = {neighbors[0], neighbors[1]}; // Should already be sorted

								// Add entry that maps parent edges to the refined edges for hanging node cases
								parent_to_refined_edge[edge] = refined_edge;
							}
						}
					}
				}

				// If I've already processed this edge then it will be in one of the sets
				std::pair<std::pair<id_type, id_type>, id_type> insertion(edge, inc_nums[i]);
				if (cut_edges.find(insertion)!=cut_edges.end() || cut_edges_on_interface.find(insertion)!=cut_edges_on_interface.end())
					continue;

				// See if there are other elements on remote processors that also have this edge
				std::vector<id_type> nodes = {edge.first, edge.second};
				std::vector<int> common_parts = partitions_in_common( nodes );

				// If the only partition is myself, then I am inside my partition
				if (common_parts.size() == 1)
					cut_edges.insert( insertion );

				// Otherwise I might be on an interface
				else
				{
					// I have to check if I actually have elements on remote processors that have an edge containing this element
					std::map<int, std::vector<id_type> > r_elem = remote_elem_in_common(nodes);
					if (r_elem.size() != 0)
					{
						cut_edges_on_interface.insert( insertion );
						edge_int_map[edge] = std::vector<int>(1, _rank);
						for (auto it=r_elem.begin(), end=r_elem.end(); it!=end; ++it)
							edge_int_map[edge].push_back(it->first);
					}
					else // Otherwise I'm not actually on an interface
						cut_edges.insert( insertion );
				}
				
			} // End edge loop
		} // End element loop
		
		
		
		
		// 2. Find the actual intersection points of the inclusions with the given edges.
		//	-Do this for all edges not on a partition interface as well as those on a
		//	 partition interface where I own the node with the smaller global ID of the pairs
		//	-Create new nodes at these points
		//	-Create a map from the local enrichment node numbers back to the edge it belongs to
		_enrich_nodes.resize(cut_edges.size() + cut_edges_on_interface.size());
		_enrich_owners.resize(cut_edges.size() + cut_edges_on_interface.size());
		_enrich_node_to_edge.resize(cut_edges.size() + cut_edges_on_interface.size());
		_enrich_node_to_inclusion.resize(cut_edges.size() + cut_edges_on_interface.size());
		id_type count = 0;
		for(auto it=cut_edges.begin(), end=cut_edges.end(); it!=end; ++it)
		{
			std::pair<id_type, id_type> edge = it->first;
			id_type inclusion = it->second;
			// Form a vector of nodes that define the edge that is intersected. This can later be extended to 3+ nodes
			std::vector<Node> edge_nodes(2);
			edge_nodes[0] = *_nodes[ _global_to_local_node[edge.first] ];
			edge_nodes[1] = *_nodes[ _global_to_local_node[edge.second] ];
			
			// Actually find the intersection point and return a new node
			std::vector<double> coords = _inclusions[inclusion]->find_intersection(edge_nodes);
			_enrich_nodes[count] = new Node(coords[0], coords[1], coords[2], 0); // Don't worry about global ID yet, we'll worry about that next
			_enrich_owners[count] = _rank;
			_enrich_node_to_edge[count] = edge;
			_enrich_node_to_inclusion[count] = inclusion;
			count++;
		}
		// Set the number nodes that are not on an interface so I know which ones to send later
		id_type n_local_interior_e_nodes = count;
		
		std::map<int, id_type> n_expected;
		for(auto it=cut_edges_on_interface.begin(), end=cut_edges_on_interface.end(); it!=end; ++it)
		{
			std::pair<id_type, id_type> edge = it->first;
			id_type inclusion = it->second;
			// If the owner of the node with the lowest id also has this edge, then that partition is in charge of making the enrichment node
			// If it does not have this edge, then the lowest partition number that does have this edge is in charge of it
			int owner = get_node_owner_global(edge.first);
			bool lowest_owner = (std::find(edge_int_map[ edge ].begin(), edge_int_map[ edge ].end(), owner) != edge_int_map[ edge ].end()); // The owner of the lowest node does have this edge
			int lowest_part = *(std::min_element(edge_int_map[ edge ].begin(), edge_int_map[ edge ].end()));

			// If I own the node with the lower global node ID
			if ((lowest_owner && owner==_rank) || (!lowest_owner && lowest_part==_rank))
			{
				// Form a vector of nodes that define the edge that is intersected. This can later be extended to 3+ nodes
				std::vector<Node> edge_nodes;
				edge_nodes.push_back( *_nodes[ _global_to_local_node[edge.first] ] );
				edge_nodes.push_back( *_nodes[ _global_to_local_node[edge.second] ] );
				
				// Actually find the intersection point and return a new node
				std::vector<double> coords = _inclusions[inclusion]->find_intersection(edge_nodes);
				_enrich_nodes[count] = new Node(coords[0], coords[1], coords[2], 0); // Don't worry about global ID yet, we'll worry about that next
				_enrich_owners[count] = _rank;
				_enrich_node_to_edge[count] = edge;
				_enrich_node_to_inclusion[count] = inclusion;
				count++;
			}
			else
			{
				if (lowest_owner)
					n_expected[ owner ]++;
				else
					n_expected[ lowest_part ]++;
			}
		}

		// Create the mirror nodes for all of the nodes just created
		_n_local_owned_enrich_nodes = count;
		std::vector<Node*> mirror_nodes; // This is the vector of all the mirror nodes that will exist in this mesh (Will be concatenated w/ _enrich_nodes later)
		if (is_cohesive())
		{
			mirror_nodes.resize(_enrich_nodes.size());
			for(id_type n=0; n<_n_local_owned_enrich_nodes; ++n)
			{
				Node* copy = _enrich_nodes[n];
				mirror_nodes[n] = new Node((*copy)(0), (*copy)(1), (*copy)(2), 0);
			}
		}
		
		
		
		
		// 3. Assign global enrichment node numbers
		//	-Do a parallel prefix sum to determine number of enrichment node on processors before
		//	 mine and then label similar to dofs
		id_type partial_sum;
		MPI_Scan(&_n_local_owned_enrich_nodes, &partial_sum, 1, MPI_ID, MPI_SUM, _communicator);
		MPI_Allreduce(&_n_local_owned_enrich_nodes, &_n_global_enrich_nodes, 1, MPI_ID, MPI_SUM, _communicator);
		partial_sum -= _n_local_owned_enrich_nodes;
		
		// Now give all of them the actual numbering
		for(id_type local = 0; local<_n_local_owned_enrich_nodes; ++local)
		{
			_enrich_nodes[local]->set_id( local + partial_sum + ENRICH_START );
			_global_to_local_enrich_node[local + partial_sum + ENRICH_START] = local + ENRICH_START;
		}

		// Now if there are mirror nodes, given them numbering after all of the regular enrichment nodes
		if (is_cohesive())
		{
			for(id_type local = 0; local<_n_local_owned_enrich_nodes; ++local)
				mirror_nodes[local]->set_id( local + partial_sum + ENRICH_START + _n_global_enrich_nodes );
				// Note: Don't set the _global_to_local map here becuase I don't know exactly how many new nodes I'll be recieving (I guess I do, but I don't feel like figuring it out)
		}
		
		
		
		
		// 4. Communicate partition interface nodes to all other partitions that have copies of both of the edge's ends.
		int ids_per_coords = 1;
		if (is_cohesive()) ids_per_coords = 2; // If I'm sending the mirror nodes, then I need to send 2 ids for every coordinate set
		std::map<int, std::vector<id_type> > id_send_lists;
		std::map<int, std::vector<double> > coord_send_lists;
		std::map<int, std::vector<id_type> > pptr_send_lists;
		std::map<int, std::vector<int> > pind_send_lists;
		std::map<int, id_type> curr_inds;
		for(id_type n=n_local_interior_e_nodes; n<_n_local_owned_enrich_nodes; ++n)
		{
			Node* e_node = _enrich_nodes[n];
			
			// Send the edge information to all partitions that were found to contain the edge previously (edge_int_map)
			std::pair<id_type, id_type> edge = _enrich_node_to_edge[n];
			id_type inclusion = _enrich_node_to_inclusion[n];
			std::vector<int>& parts = edge_int_map[ edge ];
			_enrich_nodes_interface_lists.insert(std::pair<id_type, std::vector<int> >(e_node->get_id(), parts));

			for(id_type p=0; p<parts.size(); ++p)
			{
				int part = parts[p];
				if (part != _rank) // Only need to do this for ranks that are not my own
				{
					// Make sure there are at least empty vectors to insert into
					if (id_send_lists.find(part)==id_send_lists.end())
					{
						id_send_lists.insert(std::pair<int, std::vector<id_type> >(part, std::vector<id_type>()));
						coord_send_lists.insert(std::pair<int, std::vector<double> >(part, std::vector<double>()));
						pptr_send_lists.insert(std::pair<int, std::vector<id_type> >(part, std::vector<id_type>(1, 0)));
						pind_send_lists.insert(std::pair<int, std::vector<int> >(part, std::vector<int>()));
					}
					
					// Add the enrichment node id, as well as the regular edge node ids to the id_send_lists
					id_send_lists[part].push_back( e_node->get_id() );
					if (is_cohesive()) id_send_lists[part].push_back( mirror_nodes[n]->get_id() );
					id_send_lists[part].push_back( edge.first );
					id_send_lists[part].push_back( edge.second );
					id_send_lists[part].push_back( inclusion );
					
					// Add the enrichment node's coordinates to the coord_snend_list
					coord_send_lists[part].push_back( (*e_node)(0) );
					coord_send_lists[part].push_back( (*e_node)(1) );
					coord_send_lists[part].push_back( (*e_node)(2) );

					// Add the interface information
					pind_send_lists[part].insert(pind_send_lists[part].end(), parts.begin(), parts.end());
					curr_inds[part] += parts.size();
					pptr_send_lists[part].push_back(curr_inds[part]);
				}
			}
		}
		
		// Send the information
		MPI_Request reqs[id_send_lists.size()*4];
		id_type n_sent = 0;
		for(auto it=id_send_lists.begin(), end=id_send_lists.end(); it!=end; ++it)
		{
			int part = (*it).first;
			MPI_Isend(&(*it).second[0], (*it).second.size(), MPI_ID, part, 0, _communicator, &reqs[4*n_sent]);
			MPI_Isend(&coord_send_lists[part][0], coord_send_lists[part].size(), MPI_DOUBLE, part, 1, _communicator, &reqs[4*n_sent + 1]);
			MPI_Isend(&pptr_send_lists[part][0], pptr_send_lists[part].size(), MPI_ID, part, 2, _communicator, &reqs[4*n_sent+2]);
			MPI_Isend(&pind_send_lists[part][0], pind_send_lists[part].size(), MPI_INT, part, 3, _communicator, &reqs[4*n_sent+3]);
			n_sent++;
		}
		
		// Recieve the enrichment node information
		int msg_count = 2+ids_per_coords + 1; // 2 for the 2 edge nodes, 1 for the inclusion number
		for(auto it=n_expected.begin(), end=n_expected.end(); it!=end; ++it)
		{
			int part = (*it).first;
			std::vector<id_type> id_recv_vec;		Utilities::RecieveUnknown(id_recv_vec, part, 0, MPI_ID, get_comm());
			std::vector<double> coord_recv_vec;		Utilities::RecieveUnknown(coord_recv_vec, part, 1, MPI_DOUBLE, get_comm());
			std::vector<id_type> pptr_recv_vec;		Utilities::RecieveUnknown(pptr_recv_vec, part, 2, MPI_ID, get_comm());
			std::vector<int> pind_recv_vec;			Utilities::RecieveUnknown(pind_recv_vec, part, 3, MPI_INT, get_comm());
			if ((id_recv_vec.size()/msg_count)!=(*it).second) // This is really just a safety check
				err_message("The number of enrichment nodes recieved did not match the number of enrichment nodes expected");
			
			// Parse the recieved vectors
			for(id_type n=0; n<(id_recv_vec.size()/msg_count); ++n)
			{
				id_type e_id = id_recv_vec[n*msg_count];
				id_type mirror_id = 0;
				id_type n0_id = 0;
				id_type n1_id = 0;
				id_type inclusion = 0;
				if (is_cohesive())
				{
					mirror_id = id_recv_vec[n*msg_count+1];
					n0_id = id_recv_vec[n*msg_count+2];
					n1_id = id_recv_vec[n*msg_count+3];
					inclusion = id_recv_vec[n*msg_count+4];
				}
				else
				{
					n0_id = id_recv_vec[n*msg_count+1];
					n1_id = id_recv_vec[n*msg_count+2];
					inclusion = id_recv_vec[n*msg_count+3];
				}
				std::pair<id_type, id_type> edge(n0_id, n1_id);
				std::pair<std::pair<id_type,id_type>, id_type> edge_inc(edge, inclusion);
				if (cut_edges_on_interface.find(edge_inc) == cut_edges_on_interface.end())
					err_message("Unexpected enriched edge-inclusion combination.");
				double x = coord_recv_vec[n*3];
				double y = coord_recv_vec[n*3+1];
				double z = coord_recv_vec[n*3+2];
				_global_to_local_enrich_node[e_id] = count + ENRICH_START;
				_enrich_nodes[count] = new Node(x, y, z, e_id);
				_enrich_owners[count] = part;
				_enrich_node_to_edge[count] = edge;
				_enrich_node_to_inclusion[count] = inclusion;

				id_type start = pptr_recv_vec[n];
				id_type end = pptr_recv_vec[n+1];
				_enrich_nodes_interface_lists.insert(std::pair<id_type, std::vector<int> >(e_id, std::vector<int>(pind_recv_vec.begin()+start, pind_recv_vec.begin()+end)));

				// Create the mirror nodes
				if (is_cohesive())
					mirror_nodes[count] = new Node(x, y, z, mirror_id);

				// Upate the enrichent node counter
				count++;
			}
		}
		// Add all of the mirror nodes
		id_type n_initial = _enrich_nodes.size();
		if (is_cohesive())
		{
			// Copy the actual nodes
			_enrich_nodes.insert(_enrich_nodes.end(), mirror_nodes.begin(), mirror_nodes.end());

			// Copy the vector objects
			_enrich_node_to_edge.resize(2*n_initial);
			std::copy_n(_enrich_node_to_edge.begin(), n_initial, _enrich_node_to_edge.begin()+n_initial);
			_enrich_node_to_inclusion.resize(2*n_initial);
			std::copy_n(_enrich_node_to_inclusion.begin(), n_initial, _enrich_node_to_inclusion.begin()+n_initial);
			_enrich_owners.resize(2*n_initial);
			std::copy_n(_enrich_owners.begin(), n_initial, _enrich_owners.begin()+n_initial);
		
			// Copy the map objects
			for(id_type n=n_initial; n<_enrich_nodes.size(); ++n)
				_global_to_local_enrich_node[_enrich_nodes[n]->get_id()] = n + ENRICH_START;

			for(id_type n=n_initial; n<_enrich_nodes.size(); ++n)
				_enrich_nodes_interface_lists.insert(std::pair<id_type, std::vector<int> >
						(_enrich_nodes[n]->get_id(), _enrich_nodes_interface_lists[_enrich_nodes[n-n_initial]->get_id()]));

			// Update some number storage
			_n_local_owned_enrich_nodes *= 2;
			_n_global_enrich_nodes *= 2;
		}
		
		
		
		
		// 5. Add the node pointers to the elements
		//	-Make sure order of addition matches the order of the _cutEdge vector for each element
		_enrich_node_elem.resize(n_local_enrich_nodes());
		for(id_type n=0; n<n_initial; ++n)
		{
			std::pair<id_type, id_type> edge = _enrich_node_to_edge[n];
			id_type inclusion = _enrich_node_to_inclusion[n];
			add_enrich_node_to_elem(edge, inclusion, n, n_initial);

			// If this edge contains a hanging node, then add the enrichment node to the refined edge as well
			auto it = parent_to_refined_edge.find( edge );
			if (it != parent_to_refined_edge.end())
				add_enrich_node_to_elem(it->second, inclusion, n, n_initial);
		}

		// 6. Build the integration elements for each element (And the cohesive elements)
		for(Mesh::element_iterator it=active_elements_begin(), end=active_elements_end(); it!=end; ++it)
			(*it)->generate_integration_elem();


		// Assign the materials added boolean because IGFEM takes care of that by having each element have a pointer to its material
		_materials_assigned = true;

		MPI_Waitall(id_send_lists.size()*2, reqs, MPI_STATUSES_IGNORE);
	}
}



// Asigns the element node connectivity for the n'th local enrichment node
void Mesh::add_enrich_node_to_elem(std::pair<id_type, id_type> edge, id_type inclusion, id_type n, id_type n_initial)
{
	std::vector<id_type> nodes = {edge.first, edge.second};
	std::vector<id_type> common = local_elem_in_common(nodes);

	// Add these common elements to the node's elem node entries
	std::vector<id_type> vec_union;
	std::set_union(_enrich_node_elem[n].begin(), _enrich_node_elem[n].end(),
					common.begin(), common.end(), std::back_inserter(vec_union));
	_enrich_node_elem[n].swap(vec_union);
	if (is_cohesive())
		_enrich_node_elem[n+n_initial] = _enrich_node_elem[n];

	for(id_type e=0; e<common.size(); ++e)
	{
		id_type local_e = _global_to_local_elem[common[e]];
		Elem* el = _elem[local_e];

		std::vector<id_type> nodes = {edge.first, edge.second};
		int correct_edge = el->get_edge_from_nodes(nodes);
		if (correct_edge == -1)
			err_message("Matching edge was not found in element when adding enrichment node.");
	
		// Find which element local enrichment node this is
		// Find which entry of the cutEdges vector matches the edge where the correponding entry in the inc_nums array also matches the inclusion
		std::vector<short_id_type>& cutEdges = el->getIntersectedEdges();
		std::vector<id_type>& inc_nums = el->getInclusionNumbers();
		std::vector<short_id_type>::iterator iter = cutEdges.begin();
		while (1)
		{
			iter = std::find(iter, cutEdges.end(), correct_edge);
			if (iter==cutEdges.end())
				err_message("Matching cut edge was not found in the element.");
			else
			{
				int node = iter-cutEdges.begin();

				// Found the correct enrichment node!
				if (inc_nums[node] == inclusion)
				{
					el->set_enrich_node(node, _enrich_nodes[n]);

					// If the mesh is cohesive, then set the upper surface node as the mirror node of the current edge
					if (is_cohesive())
						el->set_enrich_node(node+cutEdges.size(), _enrich_nodes[n+n_initial]);

					// Break out of the while loop
					break;
				}

				// Found the right edge but not the correct inclusion. Keep searching the rest of the cutEdge vector for another matching edge
				else
					iter++;
			}
		}
	}
}




// Computes how far along the enrichment node specified by the id is along its edge
// 0 is defined as on top of the node with the lowest is, 1 is on the other side of the edge
double Mesh::compute_enrich_node_interpolation(id_type id)
{
	id_type l_node = global_to_local_node( id );
	std::pair<id_type, id_type> edge = get_enriched_edge_local(l_node);
	std::vector<double> coords1 = get_node_global( edge.first )->get_coords();
	std::vector<double> coords2 = get_node_global( edge.second )->get_coords();
	std::vector<double> ecoords = get_node_global( id )->get_coords();
	double ratio = 0.0;
	int dim_count = 0;
	for(id_type d=0; d<dim(); ++d)
	{
		double diff = coords2[d] - coords1[d];
		// For some reason the absolute value function isn't working here... So I'll just hardcode it
		double check1 = diff/coords1[d]; if (check1 < 0.0) check1=-1.0*check1;
		double check2 = diff/coords2[d]; if (check2 < 0.0) check2=-1.0*check2;
		if (check1 > 1e-10 ||
		   check2 > 1e-10)
		{
			ratio += (ecoords[d] - coords1[d])/diff;
			dim_count++;
		}
	}
	ratio = ratio/dim_count;

	if (ratio > 1.0 || ratio < 0.0)
		err_message("Invalid result for an enrichment node interpolation factor.");

	return ratio;
}


std::pair<id_type, id_type> Mesh::get_edge_from_enrich_node(id_type enrich_id)
{
	id_type l_node = global_to_local_node( enrich_id );
	if ((l_node-ENRICH_START) < n_local_enrich_nodes())
		return _enrich_node_to_edge[l_node - ENRICH_START];
	else
		err_message("Attempted to access the edge of an enrichment node not existing on the local mesh.");
}














// FUNCTIONS HAVING TO DO WITH A REFINED MESH
// ===================================================================================


void Mesh::set_active_elements()
{
	// Clear the current list of active elements
	_active_elem.clear();

	// Get the number of active elements
	id_type n_active = std::count(_elem_activity.begin(), _elem_activity.end(), true);
	_active_elem.resize(n_active);

	id_type count = 0;
	for (id_type e=0; e<_elem.size(); ++e)
	{
		if (_elem_activity[e])
		{
			_active_elem[count] = _elem[e];
			count++;
		}

	}
}


// Searches the _refine_node_neighbors structure for amatching global node id for the given neighbors
id_type Mesh::get_refine_node_from_neighbors(const std::vector<id_type>& neighbors)
{
	auto it = _refine_node_neighbors.find(neighbors);
	if (it != _refine_node_neighbors.end())
		return it->second;
	else
		err_message("Attempted to find the refinement node for a set of neighbors not existing in the local mesh.");
}


unsigned char Mesh::get_p_level_local(id_type idx)
{
	if (idx < n_local_elem())
		return _p_level[idx];
	else
		err_message("Attempting to access the p-level of an element not in the local mesh.");
}



void Mesh::detect_hanging_nodes()
{
	// Loop over all of the nodes that have been added since the original generation of the mesh (refinement nodes)
	// Look at the node neighbors of each refienment node and find the set intersection of all of the node_elem results
	// of the refined node neighbors. If the set intersection is not empty then this is a hanging node.
	// This is because full refinement of all elements around the given edge or face where a refinement node is located
	// should have fully separated from each other

	// Clear the current hanging node structure
	_edge_hanging_nodes.clear();
	_face_hanging_nodes.clear();
	_hanging_node_neighbors.clear();

	// Create the inverse of the _refine_node_neighbors structure
	std::vector<std::vector<id_type> > refine_neighbors_inv(n_local_nodes() - _n_original_nodes);
	for (auto it=_refine_node_neighbors.begin(), end=_refine_node_neighbors.end(); it!=end; ++it)
		refine_neighbors_inv[global_to_local_node(it->second) - _n_original_nodes] = it->first;

	// Loop over all nodes that have been added since the beginning of refinement
	for (id_type n=_n_original_nodes; n<n_local_nodes(); ++n)
	{
		std::vector<id_type>& neighbors = refine_neighbors_inv[n - _n_original_nodes];
		std::vector<id_type> common_elements = elem_in_common(neighbors);

		// If the size of the intersection isn't 0, then this is a hanging node
		if (common_elements.size() != 0)
		{
			id_type id = _nodes[n]->get_id();
			if ( neighbors.size() == 2)
				_edge_hanging_nodes.insert( id );
			else
				_face_hanging_nodes.insert( id );

			_hanging_node_neighbors.insert( std::pair<id_type, std::vector<id_type> >(id, neighbors) );
		}
	}

	// Clear the _refine_node_neighbors structure of results having to do with non-hanging node?
	// I don't think non-hanging nodes can ever become hanging nodes
	// Unless I decide to implement mesh coarsening... I guess I'll leave it for now
}


// Returns the neighbors of the hanging node specified by the id if it is a hanging node
std::vector<id_type>& Mesh::get_hanging_node_neighbors(id_type id)
{
	auto it = _hanging_node_neighbors.find(id);
	if (it != _hanging_node_neighbors.end())
		return it->second;
	else
		err_message("Attempted to access the neighbors of non-hanging node.");
}

// Finds the nodes and 
void Mesh::set_constraint_values()
{
	// Clear the current structure
	_hanging_node_constraints.clear();

	// Determine the erichment nodes associated with my hanging nodes
	std::unordered_map<id_type, id_type> edge_hanging_to_enrich_nodes;
	std::unordered_map<id_type, std::vector<id_type> > face_hanging_to_enrich_nodes;
	find_hanging_enriched_neighbors(edge_hanging_to_enrich_nodes, face_hanging_to_enrich_nodes);

	// Loop over all of the edge hanging nodes and find out if the edge they are on is enriched as well
	for (auto it=_edge_hanging_nodes.begin(), end=_edge_hanging_nodes.end(); it!=end; ++it)
	{
		id_type id = *it;
		std::vector<id_type>& neighbors = get_hanging_node_neighbors( *it );

		// Add the basic constrint
		std::pair<id_type, double> entry(neighbors[0], 0.5);
		std::vector<std::pair<id_type, double> > insertion(2, entry); // Simply insert the constraint node twice
		insertion[1].first = neighbors[1]; // And change the id of the second one (fewer lines of code this way)

		// If this hanging node is associated with an enrichment node, add that contribution
		auto it2 = edge_hanging_to_enrich_nodes.find( id );
		if (it2 != edge_hanging_to_enrich_nodes.end())
		{
			id_type enrich_id = it2->second;
			double interpolate = compute_enrich_node_interpolation( enrich_id );
			double e_coeff = 0.5 / abs(1.0-interpolate);
			if (interpolate < 0.5)
				e_coeff = 0.5 / (1.0 - interpolate);
			else
				e_coeff = 0.5 / interpolate;

			insertion.push_back( std::pair<id_type, double>(enrich_id, e_coeff) );
		}

		// Add the constraint
		_hanging_node_constraints.insert( std::pair<id_type, std::vector<std::pair<id_type, double> > >(id, insertion) );

	}

	// Something with face hanging nodes, No idea how to deal with this yet
	for (auto it=_face_hanging_nodes.begin(), end=_face_hanging_nodes.end(); it!=end; ++it)
	{
		// Need to create a way to map the face hanging nodes to the individual integration elemnt that it lies within in the parent element.
		// This integration element is then used to find the shape function for all of the enrichment nodes along the boundary of the face to which this hanging node belongs
		// These shape function values are then used as the coeffiecients for the constraints on those enriched dofs

		id_type id = *it;
		std::vector<id_type>& neighbors = get_hanging_node_neighbors( id );

		// Add the basic constrint
		std::vector<std::pair<id_type, double> > insertion;
		double coeff = 1.0 / neighbors.size();
		for (id_type n=0; n<neighbors.size(); ++n)
			insertion.push_back( std::pair<id_type, double>(neighbors[n], coeff) );
		
		// If this hanging node is associated with an enrichment node, add that contribution
		auto it2 = face_hanging_to_enrich_nodes.find( id );
		if (it2 != face_hanging_to_enrich_nodes.end())
		{
			err_message("Don't know how to handle enriched face hanging nodes yet"); // This will only happen with enriched hexehedrals I think so it shouldn't be a problem yet
		}

		// Add the constraint
		_hanging_node_constraints.insert( std::pair<id_type, std::vector<std::pair<id_type, double> > >(id, insertion) );
	}
}


// If any of the enrichment nodes lie on an enriched edge or face, this detects that and lists whichglobal enrichment nodes are associated with each hanging node
// NOTE: This could have been trivially done in the add_enrichments function but I wasn't sure how to extend that to face anging nodes as well
// This also has the added benefit of keeping the enrichment code mostly separate from the refinement code which hopefully means I can change things without breaking everything
void Mesh::find_hanging_enriched_neighbors(std::unordered_map<id_type, id_type>& edges_map, std::unordered_map<id_type, std::vector<id_type> >& faces_map)
{
	// Find the inverse of the _enriched_node_to_edge map
	std::unordered_map<std::pair<id_type, id_type>, id_type> edge_to_enrich_node;
	id_type nn = n_local_enrich_nodes();
	if (is_cohesive())
		nn = nn/2;
	edge_to_enrich_node.rehash(std::ceil((double)nn / edge_to_enrich_node.max_load_factor()));
	for (id_type n=0; n<nn; ++n)
		edge_to_enrich_node.insert( std::pair<std::pair<id_type, id_type>, id_type>(_enrich_node_to_edge[n], n) );

	// Loop over all the hanging nodes and see if their parent edge is in the list of refined edges
	// (Enriched edges will apear there because we default to the parent edge in the add_enrichments function)
	for (auto it=_edge_hanging_nodes.begin(), end=_edge_hanging_nodes.end(); it!=end; ++it)
	{
		// Get the id of the hanging node
		id_type hanging_id = *it;

		std::vector<id_type>& neighbors = get_hanging_node_neighbors( *it );
		if (neighbors.size() == 2)
		{
			std::pair<id_type, id_type> parent_edge(neighbors[0], neighbors[1]); // Should already be sorted
			auto it2 = edge_to_enrich_node.find(parent_edge);
			if (it2 != edge_to_enrich_node.end())
			{
				id_type enrich_id = _enrich_nodes[it2->second]->get_id();

				// If the mesh isn't cohesive, then I"ve found the matching enrichment node
				if (!is_cohesive())
					edges_map.insert( std::pair<id_type, id_type>(*it, enrich_id) );
				// Otherwise, I need to determeine if I'm supposed to use the mirror node or not
				else
				{
					id_type mirror_id = _enrich_nodes[it2->second + nn]->get_id(); // Get the id of the mirror node which might be needed
					double interpolate = compute_enrich_node_interpolation( enrich_id ); // See how far along the edge the enrichment node is
					if (_node_detect[global_to_local_node(neighbors[0])] < 0) // Base material
					{
						if (interpolate > 0.5)
							edges_map.insert( std::pair<id_type, id_type>(hanging_id, enrich_id) );
						else
							edges_map.insert( std::pair<id_type, id_type>(hanging_id, mirror_id) );
					}
					else // inside the inclusion
					{
						if (interpolate > 0.5)
							edges_map.insert( std::pair<id_type, id_type>(hanging_id, mirror_id) );
						else
							edges_map.insert( std::pair<id_type, id_type>(hanging_id, enrich_id) );
					}
				}
			}
		}
		else
			err_message("Invalid number of neighbors for an edge haing node.");
	}


	// Determine which face hanging nodes lie on an enriched face
	// Similar idea but have to loop over every combo of node pairs
	for (auto it=_face_hanging_nodes.begin(), end=_face_hanging_nodes.end(); it!=end; ++it)
	{
		// Get the id of the hanging node
		id_type hanging_id = *it;

		std::vector<id_type>& neighbors = get_hanging_node_neighbors( *it );

		// Loop over every combination of node pair and see if that combination is in the enrichment node map
		for (id_type n1=0; n1<(neighbors.size()-1); ++n1)
		{
			for (id_type n2=(n1+1); n2<neighbors.size(); ++n2)
			{
				std::pair<id_type, id_type> parent_edge(neighbors[n1], neighbors[n2]); // Should already be sorted

				auto it2 = edge_to_enrich_node.find(parent_edge);
				if (it2 != edge_to_enrich_node.end())
				{
					// Make sure there's at least an empty vector to insert into
					if (faces_map.find(hanging_id) == faces_map.end())
						faces_map.insert( std::pair<id_type, std::vector<id_type> >(hanging_id, std::vector<id_type>()) );

					// Get the id of the base enrichment node
					id_type enrich_id = _enrich_nodes[it2->second]->get_id();

					// If the mesh isn't cohesive, then I"ve found the matching enrichment node
					if (!is_cohesive())
						faces_map[*it].push_back( enrich_id );
					// Otherwise, I need to determeine if I'm supposed to use the mirror node or not
					else
					{
						// This is a little more complicated on a face
						err_message("Don't know how to handle cohesive enrichment on a face anging node yet!");
					}
				}
			}
		}
	}
}



// Returns whether or not the global node id is a hanging node
bool Mesh::is_hanging_node(id_type id)
{
	if (_edge_hanging_nodes.find(id) != _edge_hanging_nodes.end())
		return true;
	else if (_face_hanging_nodes.find(id) != _face_hanging_nodes.end())
		return true;
	else
		return false; 
}

std::vector<std::pair<id_type, double> >& Mesh::get_hanging_node_constraint(id_type id)
{
	auto it = _hanging_node_constraints.find(id);
	if (it != _hanging_node_constraints.end())
		return it->second;
	else
		err_message("Attempted to access the hanging node contraints of a node hat is not a hanging node.");
}

