// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// C++ includes
#include <string>
#include <cstdlib> // std::strtol
#include <sstream>
#include <ctype.h> // isspace

// Local includes
#include "abaqus_io.h"
#include "Mesh.h"
#include "node.h"
#include "elem.h"


// Anonymous namespace to hold mapping Data for Abaqus/libMesh element types


/**
 * Data structure used for mapping Abaqus IDs to libMesh IDs, and
 * eventually (possibly) vice-versa.
 */
struct ElementDefinition
{
	// Maps (zero-based!) Abaqus local node numbers to libmesh local node numbers
	std::vector<unsigned> abaqus_zero_based_node_id_to_libmesh_node_id;
	
	// Maps (zero-based!) Abaqus side numbers to libmesh side numbers
	std::vector<unsigned short> abaqus_zero_based_side_id_to_libmesh_side_id;
};

/**
 * Locally-available map containing all element data.
 */
std::map<elem_type, ElementDefinition> eletypes;

/**
 * Helper function to fill up eletypes map
 */
void add_eletype_entry(elem_type libmesh_elem_type,
                       const id_type* node_map,
                       unsigned node_map_size,
                       const unsigned short* side_map,
                       unsigned side_map_size)
{
	// If map entry does not exist, this will create it
	ElementDefinition& map_entry = eletypes[libmesh_elem_type];
	
	
	// Use the "swap trick" from Scott Meyer's "Effective STL" to swap
	// an unnamed temporary vector into the map_entry's vector.  Note:
	// the vector(iter, iter) constructor is used.
	std::vector<unsigned>
		(node_map, node_map+node_map_size).swap
		(map_entry.abaqus_zero_based_node_id_to_libmesh_node_id);
	
	std::vector<unsigned short>
		(side_map, side_map+side_map_size).swap
		(map_entry.abaqus_zero_based_side_id_to_libmesh_side_id);
}


/**
 * Helper function to initialize the eletypes map.
 */
void init_eletypes ()
{
	// This should happen only once.  The first time this method is
	// called the eletypes data struture will be empty, and we will
	// fill it.  Any subsequent calls will find an initialized
	// eletypes map and will do nothing.
	if (eletypes.empty())
	{
		{
			// EDGE2
			const id_type  node_map[] = {0,1}; // identity
			const unsigned short side_map[] = {0,1}; // identity
			add_eletype_entry(EDGE2, node_map, 2, side_map, 2);
		}
		
		{
			// TRI3
			const id_type  node_map[] = {0,1,2}; // identity
			const unsigned short side_map[] = {0,1,2}; // identity
			add_eletype_entry(TRI3, node_map, 3, side_map, 3);
		}
		
		{
			// QUAD4
			const id_type  node_map[] = {0,1,2,3}; // identity
			const unsigned short side_map[] = {0,1,2,3}; // identity
			add_eletype_entry(QUAD4, node_map, 4, side_map, 4);
		}
		
		{
			// TET4
			const id_type  node_map[] = {0,1,2,3}; // identity
			const unsigned short side_map[] = {0,1,2,3}; // identity
			add_eletype_entry(TET4, node_map, 4, side_map, 4);
		}
		
		{
			// TET10
			const id_type  node_map[] = {0,1,2,3,4,5,6,7,8,9}; // identity
			const unsigned short side_map[] = {0,1,2,3};             // identity
			add_eletype_entry(TET10, node_map, 10, side_map, 4);
		}
		
		{
			// HEX8
			const id_type  node_map[] = {0,1,2,3,4,5,6,7}; // identity
			const unsigned short side_map[] = {0,5,1,2,3,4};     // inverse = 0,2,3,4,5,1
			add_eletype_entry(HEX8, node_map, 8, side_map, 6);
		}
		
		{
			// HEX20
			const id_type   node_map[] = // map is its own inverse
			{0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15};
			const unsigned short side_map[] = // inverse = 0,2,3,4,5,1
			{0,5,1,2,3,4};
			add_eletype_entry(HEX20, node_map, 20, side_map, 6);
		}
		
		{
			// HEX27
			const id_type  node_map[] = // inverse = ...,21,23,24,25,26,22,20
			{0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15,26,20,25,21,22,23,24};
			const unsigned short side_map[] = // inverse = 0,2,3,4,5,1
			{0,5,1,2,3,4};
			add_eletype_entry(HEX27, node_map, 27, side_map, 6);
		}
		
		{
			// PRISM6
			const id_type  node_map[] = {0,1,2,3,4,5}; // identity
			const unsigned short side_map[] = {0,4,1,2,3};   // inverse = 0,2,3,4,1
			add_eletype_entry(PRISM6, node_map, 6, side_map, 5);
		}
		
		{
			// PRISM15
			const id_type  node_map[] = // map is its own inverse
			{0,1,2,3,4,5,6,7,8,12,13,14,9,10,11};
			const unsigned short side_map[] = // inverse = 0,2,3,4,1
			{0,4,1,2,3};
			add_eletype_entry(PRISM15, node_map, 15, side_map, 5);
		}
		
		{
			// PRISM18
			const id_type  node_map[] = // map is its own inverse
			{0,1,2,3,4,5,6,7,8,12,13,14,9,10,11,15,16,17};
			const unsigned short side_map[] = // inverse = 0,2,3,4,1
			{0,4,1,2,3};
			add_eletype_entry(PRISM18, node_map, 18, side_map, 5);
		}
	} // if (eletypes.empty())
}





AbaqusIO::AbaqusIO () :
	_already_seen_part(false)
{
}




AbaqusIO::~AbaqusIO ()
{
}




void AbaqusIO::read (const std::string& fname)
{
	// Clear any existing mesh data (Actually just assume that the mesh is empty for now)
	// the_mesh->clear();
	
	// Open stream for reading
	_in.open(fname.c_str());
	if(!_in.good())
		err_message("Bad input file stream.");
	
	// Initialize the elems_of_dimension array.  We will use this in a
	// "1-based" manner so that elems_of_dimension[d]==true means
	// elements of dimension d have been seen.
	elems_of_dimension.resize(4, false);
	
	// Read file line-by-line... this is based on a set of different
	// test input files.  I have not looked at the full input file
	// specs for Abaqus.
	std::string s;
	while (true)
	{
		// Try to read something.  This may set EOF!
		std::getline(_in, s);
		
		if (_in)
		{
			// Process s...
			//
			// There are many sections in Abaqus files, we read some
			// but others are just ignored...  Some sections may occur
			// more than once.  For example for a hybrid grid, you
			// will have multiple *Element sections...
			
			// Some Abaqus files use all upper-case for section names,
			// so we will just convert s to uppercase
			std::string upper(s);
			std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
			
			// 0.) Look for the "*Part" section
			if (upper.find("*PART") == static_cast<std::string::size_type>(0))
			{
				if (_already_seen_part)
					err_message("We currently don't support reading Abaqus files with multiple PART sections");
				
				_already_seen_part = true;
			}
			
			// 1.) Look for the "*Nodes" section
			if (upper.find("*NODE") == static_cast<std::string::size_type>(0))
			{
				// Some sections that begin with *NODE are actually
				// "*NODE OUTPUT" sections which we want to skip.  I
				// have only seen this with a single space, but it would
				// probably be more robust to remove whitespace before
				// making this check.
				if (upper.find("*NODE OUTPUT") == static_cast<std::string::size_type>(0))
					continue;
				
				// Some *Node sections also specify an Nset name on the same line.
				// Look for one here.
				std::string nset_name = this->parse_label(s, "nset");
				
				// Process any lines of comments that may be present
				this->process_and_discard_comments();
				
				// Read a block of nodes
				this->read_nodes(nset_name);
			}
			
			
			
			// 2.) Look for the "*Element" section
			else if (upper.find("*ELEMENT,") == static_cast<std::string::size_type>(0))
			{
				// Some sections that begin with *ELEMENT are actually
				// "*ELEMENT OUTPUT" sections which we want to skip.  I
				// have only seen this with a single space, but it would
				// probably be more robust to remove whitespace before
				// making this check.
				if (upper.find("*ELEMENT OUTPUT") == static_cast<std::string::size_type>(0))
					continue;
				
				// Some *Element sections also specify an Elset name on the same line.
				// Look for one here.
				std::string elset_name = this->parse_label(s, "elset");
				
				// Process any lines of comments that may be present
				this->process_and_discard_comments();
				
				// Read a block of elements
				this->read_elements(upper, elset_name);
			}
			
			
			
			// 3.) Look for a Nodeset section
			else if (upper.find("*NSET") == static_cast<std::string::size_type>(0))
			{
				std::string nset_name = this->parse_label(s, "nset");
				
				// I haven't seen an unnamed nset yet, but let's detect it
				// just in case...
				if (nset_name == "")
					err_message("Unnamed nset encountered!");

				// Check to see if abaqus used the 'generate' syntax for a set
				bool generate = false;
				if(s.find("generate") != std::string::npos)
					generate = true;
				
				// Process any lines of comments that may be present
				this->process_and_discard_comments();
				
				// Read the IDs, storing them in _nodeset_ids
				this->read_ids(nset_name, _nodeset_ids, generate);
			} // *Nodeset
			
			
			
			// 4.) Look for an Elset section
			else if (upper.find("*ELSET") == static_cast<std::string::size_type>(0))
			{
				std::string elset_name = this->parse_label(s, "elset");
				
				// I haven't seen an unnamed elset yet, but let's detect it
				// just in case...
				if (elset_name == "")
					err_message("Unnamed elset encountered!");

				// Check to see if abaqus used the 'generate' syntax for a set
				bool generate = false;
				if(s.find("generate") != std::string::npos)
					generate = true;
				
				// Process any lines of comments that may be present
				this->process_and_discard_comments();
				
				// Read the IDs, storing them in _elemset_ids
				this->read_ids(elset_name, _elemset_ids, generate);
			} // *Elset
			
			
			
			// 5.) Look for a Surface section.  Need to be a little
			// careful, since there are also "surface interaction"
			// sections we don't want to read here.
			else if (upper.find("*SURFACE,") == static_cast<std::string::size_type>(0))
			{
				// Get the name from the Name=Foo label.  This will be the map key.
				std::string sideset_name = this->parse_label(s, "name");
				
				// Process any lines of comments that may be present
				this->process_and_discard_comments();
				
				// Read the sideset IDs
				this->read_sideset(sideset_name, _sideset_ids);
			}
			
			continue;
		} // if (_in)
		
		// If !file, check to see if EOF was set.  If so, break out
		// of while loop.
		if (_in.eof())
			break;
		
		// If !in and !in.eof(), stream is in a bad state!
		err_message("Stream is bad! Perhaps the file does not exist?");
	} // while
	
	// Set the Mesh dimension based on the highest dimension element seen.
	the_mesh->set_mesh_dimension(this->max_elem_dimension_seen());
	
	// Set element IDs based on the element sets.
	//this->assign_subdomain_ids();
	
	// Assign nodeset values to the BoundaryInfo object
	//this->assign_boundary_node_ids();

	// Add nodesets and element sets to the mesh
	assign_nodesets_and_elemsets();
	
	// Assign sideset values in the BoundaryInfo object
	this->assign_sideset_ids();
	
	// Abaqus files only contain nodesets by default.  To be useful in
	// applying most types of BCs in libmesh, we will definitely need
	// sidesets.  So we can call the new BoundaryInfo function which
	// generates sidesets from nodesets.
	build_sidesets_from_nodesets();

	// Delete lower-dimensional elements from the Mesh.  We assume these
	// were only used for setting BCs, and aren't part of the actual
	// Mesh.
	bool delete_elem = false;
	for(id_type dim=max_elem_dimension_seen()-1; dim>0; dim--)
		if(elems_of_dimension[dim])
			delete_elem = true;

	if(delete_elem)
	{
		unsigned char max_dim = this->max_elem_dimension_seen();
		
		for (id_type e=0; e<the_mesh->n_elem(); ++e)
		{
			Elem* elem = the_mesh->get_elem_global(e);
			
			if (elem->dim() < max_dim)
				the_mesh->delete_elem(e);
		}
	}
}







void AbaqusIO::read_nodes(std::string nset_name)
{
	// In the input files I have, Abaqus neither tells what
	// the mesh dimension is nor how many nodes it has...
	//
	// The node line format is:
	// id, x, y, z
	// and you do have to parse out the commas.
	// The z-coordinate will only be present for 3D meshes
	
	// Temporary variables for parsing lines of text
	char c;
	std::string dummy;

	// Defines the sequential node numbering used by libmesh.  Since
	// there can be multiple *NODE sections in an Abaqus file, we always
	// start our numbering with the number of nodes currently in the
	// Mesh.
	id_type node_id = the_mesh->n_nodes();
	
	// We need to duplicate some of the read_ids code if this *NODE
	// section also defines an NSET.  We'll set up the id_storage
	// pointer and push back IDs into this vector in the loop below...
	std::vector<id_type>* id_storage = NULL;
	if (nset_name != "")
		id_storage = &(_nodeset_ids[nset_name]);
	
	// We will read nodes until the next line begins with *, since that will be the
	// next section.
	// TODO: Is Abaqus guaranteed to start the line with '*' or can there be leading white space?
	while (_in.peek() != '*' && _in.peek() != EOF)
	{
		// Re-Initialize variables to be read in from file
		id_type abaqus_node_id=0;
		double x=0, y=0, z=0;
		
		// Read in the coordinates
		_in >> abaqus_node_id >> c >> x;
		
		// Peek at the next character.  If it is a comma, then there is another
		// value to read!
		if (_in.peek() == ',')
			_in >> c >> y;
		if (_in.peek() == ',')
			_in >> c >> z;
		
		// Read (and discard) the rest of the line, including the newline.
		// This is required so that our 'peek()' at the beginning of this
		// loop doesn't read the newline character, for example.
		std::getline(_in, dummy);

		// If this *NODE section defines an NSET, also store the abaqus ID in id_storage
		if (id_storage)
			id_storage->push_back(abaqus_node_id);
		
		// Set up the abaqus -> libmesh node mapping.  This is usually just the
		// "off-by-one" map.
		_abaqus_node_mapping[abaqus_node_id] = node_id;
		_libmesh_node_mapping[node_id] = abaqus_node_id;
		
		// Add the point to the mesh using libmesh's numbering,
		// and post-increment the libmesh node counter. 0 is the mesh partition
		the_mesh->add_node(x, y, z, node_id++, 0);
	} // while
}





void AbaqusIO::read_elements(std::string upper, std::string elset_name)
{
	// initialize the eletypes map (eletypes is a file-global variable)
	init_eletypes();
	
	elem_type el_type = INVALID_ELEM;
	unsigned n_nodes_per_elem = 0;
	
	// Within s, we should have "type=XXXX"
	if (upper.find("T3D2") != std::string::npos)
	{
		el_type = EDGE2;
		n_nodes_per_elem = 2;
		elems_of_dimension[1] = true;
	}
	else if (upper.find("CPE4") != std::string::npos ||
			upper.find("CPS4") != std::string::npos)
	{
		el_type = QUAD4;
		n_nodes_per_elem = 4;
		elems_of_dimension[2] = true;
	}
	else if (upper.find("CPS3") != std::string::npos ||
			upper.find("S3R") != std::string::npos ||
			upper.find("CPE3") != std::string::npos)
	{
		el_type = TRI3;
		n_nodes_per_elem = 3;
		elems_of_dimension[2] = true;
	}
	else if (upper.find("C3D8") != std::string::npos)
	{
		el_type = HEX8;
		n_nodes_per_elem = 8;
		elems_of_dimension[3] = true;
	}
	else if (upper.find("C3D4") != std::string::npos)
	{
		el_type = TET4;
		n_nodes_per_elem = 4;
		elems_of_dimension[3] = true;
		}
	else if (upper.find("C3D20") != std::string::npos)
	{
		el_type = HEX20;
		n_nodes_per_elem = 20;
		elems_of_dimension[3] = true;
	}
	else if (upper.find("C3D6") != std::string::npos)
	{
		el_type = PRISM6;
		n_nodes_per_elem = 6;
		elems_of_dimension[3] = true;
	}
	else if (upper.find("C3D15") != std::string::npos)
	{
		el_type = PRISM15;
		n_nodes_per_elem = 15;
		elems_of_dimension[3] = true;
	}
	else if (upper.find("C3D10") != std::string::npos)
	{
		el_type = TET10;
		n_nodes_per_elem = 10;
		elems_of_dimension[3] = true;
	}
	else
		err_message("Unrecognized element type");

	// Insert the elem type we detected into the set of all elem types for this mesh
	_elem_types.insert(el_type);
	
	// Grab a reference to the element definition for this element type
	const ElementDefinition& eledef = eletypes[el_type];
	
	// If the element definition was not found, the call above would have
	// created one with an uninitialized struct.  Check for that here...
	if (eledef.abaqus_zero_based_node_id_to_libmesh_node_id.size() == 0)
		err_message("No Abaqus->LibMesh mapping information for elem_type!");
	
	// We will read elements until the next line begins with *, since that will be the
	// next section.
	id_type elem_id = the_mesh->n_elem();
	while (_in.peek() != '*' && _in.peek() != EOF)
	{
		// Read the element ID, it is the first number on each line.  It is
		// followed by a comma, so read that also.  We will need this ID later
		// when we try to assign subdomain IDs
		id_type abaqus_elem_id = 0;
		char c;
		_in >> abaqus_elem_id >> c;
		
		// Associate the ID returned from libmesh with the abaqus element ID
		//_libmesh_to_abaqus_elem_mapping[elem->id()] = abaqus_elem_id;
		_abaqus_elem_mapping[abaqus_elem_id] = elem_id;
		
		// The count of the total number of IDs read for the current element.
		unsigned id_count=0;
		std::vector<id_type> node_ids(n_nodes_per_elem);

		
		// Continue reading line-by-line until we have read enough nodes for this element
		while (id_count < n_nodes_per_elem)
		{
			// Read entire line (up to carriage return) of comma-separated values
			std::string csv_line;
			std::getline(_in, csv_line);
			
			// Create a stream object out of the current line
			std::stringstream line_stream(csv_line);
			
			// Process the comma-separated values
			std::string cell;
			while (std::getline(line_stream, cell, ','))
			{
				// FIXME: factor out this strtol stuff into a utility function.
				char* endptr;
				id_type abaqus_global_node_id = (id_type)
				(std::strtol(cell.c_str(), &endptr, /*base=*/10));
				
				if (abaqus_global_node_id!=0 || cell.c_str() != endptr)
				{
					// Use the global node number mapping to determine the corresponding libmesh global node id
					id_type libmesh_global_node_id = _abaqus_node_mapping[abaqus_global_node_id];
					
					// Note: id_count is the zero-based abaqus (elem local) node index.  We therefore map
					// it to a libmesh elem local node index using the element definition map
					unsigned libmesh_elem_local_node_id =
					eledef.abaqus_zero_based_node_id_to_libmesh_node_id[id_count];

					// Add node to the elmental vector of nodeids
					node_ids[libmesh_elem_local_node_id] = libmesh_global_node_id;
					
					// Increment the count of IDs read for this element
					id_count++;
				} // end if strtol success
			} // end while getline(',')
		} // end while (id_count)
		
		// Ensure that we read *exactly* as many nodes as we were expecting to, no more.
		if (id_count != n_nodes_per_elem)
		{
			char buf[100];
			sprintf(buf, "Needed to read %i nodes, but read %i instead!", n_nodes_per_elem, id_count);
			err_message( buf );
		}

		// Actually ad the element to the mesh
		// post-increment the elemental id so its ready for the next element
		the_mesh->add_elem(el_type, node_ids, elem_id++);
		
		// If we are recording Elset IDs, add this element to the correct set for later processing.
		// Make sure to add it with the Abaqus ID, not the libmesh one!
		if (elset_name != "")
			_elemset_ids[elset_name].push_back(abaqus_elem_id);
	} // end while (peek)
}




std::string AbaqusIO::parse_label(std::string line, std::string label_name)
{
  // Handle files which have weird line endings from e.g. windows.
  // You can check what kind of line endings you have with 'cat -vet'.
  // For example, some files may have two kinds of line endings like:
  //
  // 4997,^I496,^I532,^I487,^I948^M$
  //
  // and we don't want to deal with this when extracting a label, so
  // just remove all the space characters, which should include all
  // kinds of remaining newlines.  (I don't think Abaqus allows
  // whitespace in label names.)
  line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());

  // Do all string comparisons in upper-case
  std::string
    upper_line(line),
    upper_label_name(label_name);
  std::transform(upper_line.begin(), upper_line.end(), upper_line.begin(), ::toupper);
  std::transform(upper_label_name.begin(), upper_label_name.end(), upper_label_name.begin(), ::toupper);

  // Get index of start of "label="
  size_t label_index = upper_line.find(upper_label_name + "=");

  if (label_index != std::string::npos)
    {
      // Location of the first comma following "label="
      size_t comma_index = upper_line.find(",", label_index);

      // Construct iterators from which to build the sub-string.
      // Note: The +1 while initializing beg is to skip past the "=" which follows the label name
      std::string::iterator
        beg = line.begin() + label_name.size() + 1 + label_index,
        end = (comma_index == std::string::npos) ? line.end() : line.begin() + comma_index;

      return std::string(beg, end);
    }

  // The label index was not found, return the empty string
  return std::string("");
}




void AbaqusIO::read_ids(std::string set_name, container_t& container, bool generate)
{
	// Grab a reference to a vector that will hold all the IDs
	std::vector<id_type>& id_storage = container[set_name];

	// If we're using the generate syntax
	if(generate)
	{
		// Read entire comma-separated line into a string
		// This generate section should only be on one line
		std::string csv_line;
		std::getline(_in, csv_line);

		// Parse the line by commas into the 3? parts it should be in
		std::vector<id_type> gen_line;
		std::string cell;
		std::stringstream line_stream(csv_line);
		while(std::getline(line_stream, cell, ','))
		{
			// If no conversion can be performed by strtol, 0 is returned.
			//
			// If endptr is not NULL, strtol() stores the address of the
			// first invalid character in *endptr.  If there were no
			// digits at all, however, strtol() stores the original
			// value of str in *endptr.
			char* endptr;

			// FIXME - this needs to be updated for 64-bit inputs
			id_type id = (id_type)
				(std::strtol(cell.c_str(), &endptr, /*base=*/10));

			// Note that lists of comma-separated values in abaqus also
			// *end* with a comma, so the last call to getline on a given
			// line will get an empty string, which we must detect.
			if (id != 0 || cell.c_str() != endptr)
			{
				// 'cell' is now a string with an integer id in it
				gen_line.push_back( id );
			}
		}

		// I only have seen this syntax with 3 entries (first, last, stride)
		// Anything else I'll just throw an error
		if(gen_line.size() != 3)
			err_message("While generating a node set, unknown form of set generation string.");

		// Push back all of the entries in the generation set
		for(id_type id=gen_line[0]; id<=gen_line[1]; id+=gen_line[2])
			id_storage.push_back( id );
	}
	else
	{
		// Read until the start of another section is detected, or EOF is encountered
		while (_in.peek() != '*' && _in.peek() != EOF)
		{
			// Read entire comma-separated line into a string
			std::string csv_line;
			std::getline(_in, csv_line);

			// On that line, use std::getline again to parse each
			// comma-separated entry.
			std::string cell;
			std::stringstream line_stream(csv_line);
			while (std::getline(line_stream, cell, ','))
			{
				// If no conversion can be performed by strtol, 0 is returned.
				//
				// If endptr is not NULL, strtol() stores the address of the
				// first invalid character in *endptr.  If there were no
				// digits at all, however, strtol() stores the original
				// value of str in *endptr.
				char* endptr;

				// FIXME - this needs to be updated for 64-bit inputs
				id_type id = (id_type)
					(std::strtol(cell.c_str(), &endptr, /*base=*/10));

				// Note that lists of comma-separated values in abaqus also
				// *end* with a comma, so the last call to getline on a given
				// line will get an empty string, which we must detect.
				if (id != 0 || cell.c_str() != endptr)
				{
					// 'cell' is now a string with an integer id in it
					id_storage.push_back( id );
				}
			}
		}
	}
}




void AbaqusIO::read_sideset(std::string sideset_name, sideset_container_t& container)
{
  // Grab a reference to a vector that will hold all the IDs
  std::vector<std::pair<id_type, unsigned> >& id_storage = container[sideset_name];

  // Variables for storing values read in from file
  id_type elem_id=0;
  unsigned side_id=0;
  char c;
  std::string dummy;

  // Read until the start of another section is detected, or EOF is encountered
  while (_in.peek() != '*' && _in.peek() != EOF)
    {
      // The strings are of the form: "391, S2"

      // Read the element ID and the leading comma
      _in >> elem_id >> c;

      // Read another character (the 'S') and finally the side ID
      _in >> c >> side_id;

      // Store this pair of data in the vector
      id_storage.push_back( std::make_pair(elem_id, side_id) );

      // Extract remaining characters on line including newline
      std::getline(_in, dummy);
    } // while
}

// I don't think that I actually care about any of this...
// I think I'l just comment out the call to this function for now
void AbaqusIO::assign_subdomain_ids()
{
}
// You know what, I don't think that I care about this function either
// I guess I'll just comment out the call to this function as well
void AbaqusIO::assign_boundary_node_ids()
{
}

/*
// I don't think that I actually care about any of this...
// I think I'l just comment out the call to this function for now
void AbaqusIO::assign_subdomain_ids()
{
  // Get a reference to the mesh we are reading
  Mesh& the_mesh = MeshInput<Mesh>::mesh();

  // The number of elemsets we've found while reading
  std::size_t n_elemsets = _elemset_ids.size();

  // Fill in a temporary map with (elem_type, index) pairs based on the _elem_types set.  This
  // will allow us to easily look up this index in the loop below.
  std::map<elem_type, unsigned> elem_types_map;
  {
    unsigned ctr=0;
    for (std::set<elem_type>::iterator it=_elem_types.begin(); it!=_elem_types.end(); ++it)
      elem_types_map[*it] = ctr++;
  }

  // Loop over each Elemset and assign subdomain IDs to Mesh elements
  {
    // The maximum element dimension seen while reading the Mesh
    unsigned char max_dim = this->max_elem_dimension_seen();

    // The elemset_id counter assigns a logical numbering to the _elemset_ids keys
    container_t::iterator it = _elemset_ids.begin();
    for (unsigned elemset_id=0; it != _elemset_ids.end(); ++it, ++elemset_id)
      {
        // Grab a reference to the vector of IDs
        std::vector<id_type>& id_vector = it->second;

        // Loop over this vector
        for (std::size_t i=0; i<id_vector.size(); ++i)
          {
            // Map the id_vector[i]'th element ID (Abaqus numbering) to LibMesh numbering
            id_type libmesh_elem_id = _abaqus_elem_mapping[ id_vector[i] ];

            // Get pointer to that element
            Elem* elem = the_mesh->elem(libmesh_elem_id);

            if (elem == NULL)
              err_message("Mesh returned NULL pointer for Elem " << libmesh_elem_id);

            // We won't assign subdomain ids to lower-dimensional
            // elements, as they are assumed to represent boundary
            // conditions.  Since lower-dimensional elements can
            // appear in multiple sidesets, it doesn't make sense to
            // assign them a single subdomain id... only the last one
            // assigned would actually "stick".
            if (elem->dim() < max_dim)
              break;

            // Compute the proper subdomain ID, based on the formula in the
            // documentation for this function.
            subdomain_id_type  computed_id = cast_int<subdomain_id_type >
              (elemset_id + (elem_types_map[elem->type()] * n_elemsets));

            // Assign this ID to the element in question
            elem->subdomain_id() = computed_id;

            // We will also assign a unique name to the computed_id,
            // which is created by appending the geometric element
            // name to the elset name provided by the user in the
            // Abaqus file.
            //std::string computed_name = it->first + "_" + Utility::enum_to_string(elem->type());
            //the_mesh->subdomain_name(computed_id) = computed_name;
          }
      }
  }
}



// You know what, I don't think that I care about this function either
// I guess I'll just comment out the call to this function as well
void AbaqusIO::assign_boundary_node_ids()
{
  // Get a reference to the mesh we are reading
  Mesh& the_mesh = MeshInput<Mesh>::mesh();

  // Iterate over the container of nodesets
  container_t::iterator it = _nodeset_ids.begin();
  for (unsigned short current_id=0; it != _nodeset_ids.end(); ++it, ++current_id)
    {
      // Associate current_id with the name we determined earlier
      the_mesh->get_boundary_info().nodeset_name(current_id) = it->first;

      // Get a reference to the current vector of nodeset ID values
      std::vector<id_type>& nodeset_ids = it->second;

      for (std::size_t i=0; i<nodeset_ids.size(); ++i)
        {
          // Map the Abaqus global node ID to the libmesh node ID
          id_type libmesh_global_node_id = _abaqus_node_mapping[nodeset_ids[i]];

          // Get node pointer from the mesh
          Node* node = the_mesh->node_ptr(libmesh_global_node_id);

          if (node == NULL)
            err_message("Error! Mesh returned NULL node pointer!");

          // Add this node with the current_id (which is determined by the
          // alphabetical ordering of the map) to the BoundaryInfo object
          the_mesh->get_boundary_info().add_node(node, current_id);
        }
    }
}

*/


void AbaqusIO::assign_nodesets_and_elemsets()
{
	// Loop over nodesets first
	std::set<id_type> nodeset;
	container_t::iterator it = _nodeset_ids.begin();
	container_t::iterator end = _nodeset_ids.end();
	for(; it!=end; ++it)
	{
		// Get a reference to the current vector of nodeset ID values
		std::vector<id_type>& nodeset_ids = it->second;

		// Convert all nodeset info to libmesh numbering and put it in a set
		for(id_type i=0; i<nodeset_ids.size(); ++i)
			nodeset.insert(_abaqus_node_mapping[nodeset_ids[i]]);

		// Add the nodeset to the mesh
		the_mesh->add_nodeset_total(it->first, nodeset);

		// clear the temporary nodeset so I can fill it again
		// (Note: because I use add_nodeset_total which swpas the set, this set should actually be empty now)
		nodeset.clear();
	}

	// Loop over elemsets now
	std::set<id_type> elemset;
	it = _elemset_ids.begin();
	end = _elemset_ids.end();
	for(; it!=end; ++it)
	{
		// Get a reference to the current vector of elemset ID values
		std::vector<id_type>& elemset_ids = it->second;

		// Convert all elemset info to libmesh numbering and put it in a set
		for(id_type i=0; i<elemset_ids.size(); ++i)
			elemset.insert(_abaqus_elem_mapping[elemset_ids[i]]);

		// Add the elemset to the mesh
		the_mesh->add_elemset(it->first, elemset);

		// clear the temporary elemset so I can fill it again
		elemset.clear();
	}
}




void AbaqusIO::assign_sideset_ids()
{
  // initialize the eletypes map (eletypes is a file-global variable)
  init_eletypes();

  // Iterate over the container of sidesets
  {
  	// Create a temporary set that I'll use to store the sidesets before passing them into the mesh
  	std::set<std::pair<id_type, id_type> > sideset;

    sideset_container_t::iterator it = _sideset_ids.begin();
    for (unsigned short current_id=0; it != _sideset_ids.end(); ++it, ++current_id)
      {
        // Get a reference to the current vector of nodeset ID values
        std::vector<std::pair<id_type,unsigned> >& sideset_ids = it->second;

        for (std::size_t i=0; i<sideset_ids.size(); ++i)
          {
            // sideset_ids is a vector of pairs (elem id, side id).  Pull them out
            // now to make the code below more readable.
            id_type  abaqus_elem_id = sideset_ids[i].first;
            unsigned abaqus_side_number = sideset_ids[i].second;

            // Map the Abaqus element ID to LibMesh numbering
            id_type libmesh_elem_id = _abaqus_elem_mapping[ abaqus_elem_id ];

            // Get pointer to that element
            Elem* elem = the_mesh->get_elem_global(libmesh_elem_id);

            // Check that the pointer returned from the Mesh is non-NULL
            if (elem == NULL)
              err_message("Mesh returned NULL pointer for Elem.");

            // Grab a reference to the element definition for this element type
            const ElementDefinition& eledef = eletypes[elem->get_type()];

            // If the element definition was not found, the call above would have
            // created one with an uninitialized struct.  Check for that here...
            if (eledef.abaqus_zero_based_side_id_to_libmesh_side_id.size() == 0)
              err_message("No Abaqus->LibMesh mapping information for elem_type!");

            // Add this node with the current_id (which is determined by the
            // alphabetical ordering of the map).  Side numbers in Abaqus are 1-based,
            // so we subtract 1 here before passing the abaqus side number to the
            // mapping array
          	sideset.insert(std::pair<id_type, id_type>(libmesh_elem_id,
          		eledef.abaqus_zero_based_side_id_to_libmesh_side_id[abaqus_side_number-1]));
          }

          // Actually add the sideset to the mesh
          the_mesh->add_sideset_total(it->first, sideset);
          sideset.clear();
      }
  }


  // Some elsets (if they contain lower-dimensional elements) also
  // define sidesets.  So loop over them and build a searchable data
  // structure we can use to assign sidesets.
  {
    unsigned char max_dim = this->max_elem_dimension_seen();

    // multimap from "vector-of-lower-dimensional-element-node-ids" to subdomain ID which should be applied.
    // We use a multimap because the lower-dimensional elements can belong to more than 1 sideset.
    typedef std::multimap<std::vector<id_type>, std::string> provide_bcs_t;
    provide_bcs_t provide_bcs;

    // Create storage for temporary sideset information before passing it to the mesh
    std::map<std::string, std::set<std::pair<id_type, id_type> > > sidesets;

    // The elemset_id counter assigns a logical numbering to the _elemset_ids keys
    container_t::iterator it = _elemset_ids.begin();
    for (unsigned short elemset_id=0; it != _elemset_ids.end(); ++it, ++elemset_id)
      {
        // Grab a reference to the vector of IDs
        std::vector<id_type>& id_vector = it->second;

        // Loop over this vector
        for (std::size_t i=0; i<id_vector.size(); ++i)
          {
            // Map the id_vector[i]'th element ID (Abaqus numbering) to LibMesh numbering
            id_type libmesh_elem_id = _abaqus_elem_mapping[ id_vector[i] ];

            // Get pointer to that element
            Elem* elem = the_mesh->get_elem_global(libmesh_elem_id);

            if (elem == NULL)
              err_message("Mesh returned NULL pointer for Elem.");

            // If the element dimension is equal to the maximum
            // dimension seen, we can break out of this for loop --
            // this elset does not contain sideset information.
            if (elem->dim() == max_dim)
              break;

            // We can only handle elements that are *exactly*
            // one dimension lower than the max element
            // dimension.  Not sure if "edge" BCs in 3D
            // actually make sense/are required...
            if (elem->dim()+1 != max_dim)
            	err_message("Expected boundary element of dimension " << max_dim-1 << " but got " << elem->dim());

          	// Create storage in the map for this sideset
          	sidesets[it->first] = sideset_set_t();

            // To be pushed into the provide_bcs data container
            std::vector<id_type> elem_node_ids(elem->n_nodes());

            // Save node IDs in a local vector which will be used as a key for the map.
            for (unsigned n=0; n<elem->n_nodes(); n++)
              elem_node_ids[n] = elem->get_node(n)->get_id();

            // Sort before putting into the map
            std::sort(elem_node_ids.begin(), elem_node_ids.end());

            // Insert the (key, id) pair into the multimap
            provide_bcs.insert(std::make_pair(elem_node_ids, it->first));
          }
      }

    // Loop over elements and try to assign boundary information
    {
      for (id_type e=0; e<the_mesh->n_elem(); ++e)
        {
          Elem* elem = the_mesh->get_elem_global(e);

          if (elem->dim() == max_dim)
            {
              // This is a max-dimension element that may require BCs.
              // For each of its sides, including internal sides, we'll
              // see if a lower-dimensional element provides boundary
              // information for it.  Note that we have not yet called
              // find_neighbors(), so we can't use elem->neighbor(sn) in
              // this algorithm...
              for (unsigned short sn=0; sn<elem->n_sides(); sn++)
                {
                  Elem* side = elem->build_side(sn);

                  // Build up a node_ids vector, which is the key
                  std::vector<id_type> node_ids(side->n_nodes());
                  for (unsigned n=0; n<side->n_nodes(); n++)
                    node_ids[n] = side->get_node(n)->get_id();

                  // Sort the vector before using it as a key
                  std::sort(node_ids.begin(), node_ids.end());

                  // Look for this key in the provide_bcs multimap
                  std::pair<provide_bcs_t::const_iterator, provide_bcs_t::const_iterator>
                    range = provide_bcs.equal_range (node_ids);

                  // Add boundary information for each side in the range.
                  for (provide_bcs_t::const_iterator s_it = range.first;
                       s_it != range.second; ++s_it)
                  	sidesets[s_it->second].insert(std::pair<id_type, id_type>(elem->get_id(), sn));

                  delete side;
                }
            }
        }

        // Add all of these sidesets to the mesh
        for(auto it=sidesets.begin(), end=sidesets.end(); it!=end; ++it)
        	the_mesh->add_sideset_total(it->first, it->second);
    }
  }
}



void AbaqusIO::build_sidesets_from_nodesets()
{
	std::map<std::string, sideset_set_t> sidesets;
	// Loop over all the elements in the mesh
	for(id_type e=0; e<the_mesh->n_elem(); ++e)
	{
		Elem* el = the_mesh->get_elem_global(e);
		// Loop over all the sides of this element
		for(id_type s=0; s<el->n_sides(); ++s)
		{
			Elem* side = el->build_side(s);
			// Loop over all nodesets and see if this side belongs to a nodeset
			container_t::iterator it = _nodeset_ids.begin();
			container_t::iterator end = _nodeset_ids.end();
			for(; it!=end; ++it)
			{
				// Grad a reference to the vector that we're searching
				std::vector<id_type>& id_vec = it->second;

				// Determine if this side is in the current nodeset
				bool not_in_nodeset = false;
				for(id_type n=0; n<side->n_nodes(); ++n)
				{
					if(std::find(id_vec.begin(), id_vec.end(), _libmesh_node_mapping[side->get_node(n)->get_id()])==id_vec.end())
					{
						not_in_nodeset = true;
						break;
					}
				}
				
				// If it isn't in the current nodeset then continue
				if(not_in_nodeset)
					continue;
				
				// Otherwise, add this side to the appropriate sideset
				// Make sure there's storage there first
				if(sidesets.find(it->first)==sidesets.end())
				{
					sidesets.insert(std::pair<std::string, sideset_set_t>
						(it->first, sideset_set_t()));
				}
				sidesets[it->first].insert(std::pair<id_type, id_type>(e, s));
			}
		}
	}

	// Add all of the sidesets to the mesh
	for(auto it=sidesets.begin(), end=sidesets.end(); it!=end; ++it)
		the_mesh->add_sideset_total(it->first, it->second);
}







void AbaqusIO::process_and_discard_comments()
{
  std::string dummy;
  while (true)
    {
      // We assume we are at the beginning of a line that may be
      // comments or may be data.  We need to only discard the line if
      // it begins with **, but we must avoid calling std::getline()
      // since there's no way to put that back.
      if (_in.peek() == '*')
        {
          // The first character was a star, so actually read it from the stream.
          _in.get();

          // Peek at the next character...
          if (_in.peek() == '*')
            {
              // OK, second character was star also, by definition this
              // line must be a comment!  Read the rest of the line and discard!
              std::getline(_in, dummy);
            }
          else
            {
              // The second character was _not_ a star, so put back the first star
              // we pulled out so that the line can be parsed correctly by somebody
              // else!
              _in.unget();

              // Finally, break out of the while loop, we are done parsing comments
              break;
            }
        }
      else
        {
          // First character was not *, so this line must be data! Break out of the
          // while loop!
          break;
        }
    }
}



unsigned char AbaqusIO::max_elem_dimension_seen ()
{
  unsigned char max_dim = 0;

  unsigned char elem_dimensions_size = (unsigned char)
    (elems_of_dimension.size());
  // The elems_of_dimension array is 1-based in the UNV reader
  for (unsigned char i=1; i<elem_dimensions_size; ++i)
    if (elems_of_dimension[i])
      max_dim = i;

  return max_dim;
}
