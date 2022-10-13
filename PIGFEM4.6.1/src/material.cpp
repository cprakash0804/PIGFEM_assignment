/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "material_lei.h"
#include <iostream>


Material::Material()
	: _name_set(false), _curr_dim(0)
{
}

void Material::set_name(std::string name)
{
	if(name.length() > 0) // An actual name and not just an empty string
	{
		_name = name;
		_name_set = true;
	}
}






// Unpacks a serial char buffer into a new material pointer
// This returns a dynamically allocated material pointer so we must be sure to delete this in the mesh destructor
Material* Material::unpack(char* buf)
{
	Material* mat;
	material_type mat_type;
	size_t name_size, mat_size;

	size_t offset = 0;
	// Read the material type from the buffer
	memcpy(&mat_type, buf+offset, sizeof(material_type));
	offset += sizeof(material_type);

	// Read the size of the material name
	memcpy(&name_size, buf+offset, sizeof(size_t));
	offset += sizeof(size_t);

	// Read the size of the material in bytes
	memcpy(&mat_size, buf+offset, sizeof(size_t));
	offset += sizeof(size_t);

	// Read in the actual name
	char* name_tmp = new char[name_size];
	memcpy(name_tmp, buf+offset, name_size);
	std::string name(name_tmp, name_size); // Construct a string object from the null terminated char array
	offset += name_size;

	switch(mat_type)
	{
		case LINEAR_ELASTIC_ISOTROPIC:
			mat = new LinearElasticIsotropicMaterial;
			memcpy(mat, buf+offset, mat_size);
			mat->set_name(name);
			return mat;
		case INVALID_MATERIAL:
		default:
			err_message("Attempted to unpack an invalid material type.");
	}
	delete [] name_tmp;
	return NULL;
}






Material* Material::allocate_and_copy()
{
	Material* new_mat = allocate();   // Hopefully this works to call the appropriate allocate function based on what instance calls this function
	copy(new_mat);
	return new_mat;
}
