#include "STL_Reader.h"


void STL_Reader::read(std::string filename)
{
	// Try to open the file
	std::ifstream myfile;
	myfile.open(filename.c_str(), std::ios::in | std::ios::binary); // Assuming we're just reading in binary data for now
	if (!myfile.good())
		err_message("The given STL file name is not valid!");

	// Clear any previously stored data
	_facets.clear();
	_normals.clear();

	// Read in the header and number of facets
	char header[80];
	char ntri[4];
	myfile.read(header, 80);
	myfile.read(ntri, 4);
	unsigned int* nf = (unsigned int*) ntri;
	unsigned int n_facets = *nf;
	_facets.resize(n_facets);
	_normals.resize(n_facets);

	// Read in information for each facet
	for (id_type f=0; f<n_facets; ++f)
	{
		// Read in the normal and coordinate data
		char facet_data[48];
		myfile.read(facet_data, 48);
		float* f_ptr;
		f_ptr = (float*) &facet_data[0];
		double n1 = *f_ptr;
		f_ptr = (float*) &facet_data[4];
		double n2 = *f_ptr;
		f_ptr = (float*) &facet_data[8];
		double n3 = *f_ptr;
		f_ptr = (float*) &facet_data[12];
		double f11 = *f_ptr;
		f_ptr = (float*) &facet_data[16];
		double f12 = *f_ptr;
		f_ptr = (float*) &facet_data[20];
		double f13 = *f_ptr;
		f_ptr = (float*) &facet_data[24];
		double f21 = *f_ptr;
		f_ptr = (float*) &facet_data[28];
		double f22 = *f_ptr;
		f_ptr = (float*) &facet_data[32];
		double f23 = *f_ptr;
		f_ptr = (float*) &facet_data[36];
		double f31 = *f_ptr;
		f_ptr = (float*) &facet_data[40];
		double f32 = *f_ptr;
		f_ptr = (float*) &facet_data[44];
		double f33 = *f_ptr;
		_normals[f] = {n1, n2, n3};
		_facets[f] = {f11, f12, f13, f21, f22, f23, f31, f32, f33};

		// Read in the 2 byte attrbute byte count
		char att_byte_cnt[2];
		myfile.read(att_byte_cnt, 2); // Don't do anything with this...
	}

	// Close the file
	myfile.close();
}
