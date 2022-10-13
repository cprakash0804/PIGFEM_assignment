#include "Utilities.h"
#include <cstdlib> 
#include <iostream>


int main(int argc, char* argv[])
{
	DenseMatrix<double> mat(3);
	Utilities::timer timer;
	std::cout << "Simulating 10,000,000 determinants...\n";
	timer.start();
	for (int i=0; i<10000000; ++i)
	{
		for (int j=0; j<3; ++j)
			for (int k=0; k<3; ++k)
				mat(j,k) = rand();

		double det = Utilities::det(mat);
	}
	double time = timer.time_elapsed();
	std::cout << "Time elapsed: " << time << std::endl;

	return 0;
}