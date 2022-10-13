#ifndef _STL_READER_H_
#define _STL_READER_H_

// Local includes
#include "common.h"

// C++ includes
#include <fstream>
#include <vector>

class STL_Reader
{
	private:

		std::vector<std::vector<double> > _facets;
		std::vector<std::vector<double> > _normals;
	public:

		void read(std::string filename);

		std::vector<std::vector<double> >& get_facets() {return _facets;};
		std::vector<std::vector<double> >& get_normals() {return _normals;};
};

#endif
