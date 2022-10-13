#include <iostream>
#include <vector>
#include <stdlib.h>
#include <sys/time.h>
using namespace std;


template<typename T>
void init(SparseMatrix<T>& mat, int scale, int n_per_row)
{
	for(int i = 0; i<mat.n_rows(); ++i)
	{
		//Generate a random number between 1 and n_per_row and store this many random values in the current row
		int n = rand()%n_per_row + 1;
		for(int j = 1; j<=n; ++j)
		{
			int col = rand()%mat.n_cols();
			T val = (T)(scale*(double)rand()/RAND_MAX)+1;
			mat(i, col) = val;
		}
	}
}
template<typename T>
void init(DenseMatrix<T>& mat, int scale)
{
	for(int i = 0; i<mat.n_rows(); ++i)
	{
		for(int j = 0; j<mat.n_cols(); ++j)
		{
			T val = (T)(scale*(double)rand()/RAND_MAX);
			mat(i, j) = val;
		}
	}
}
template<typename T>
void init(std::vector<T> vec)
{
	for(int i = 0; i<vec.size(); ++i)
		vec[i] = (T)rand()/RAND_MAX;
}

template<typename T>
void print_vec(std::vector<T> &v)
{
	std::count << "{";
	for(unsigned int i=0; i<v.size(); ++i)
	{
		std::count << v[i];
		if(i!=(v.size()-1))
			std::count << ", "l
	}
}
template<typename T>
void print_mat(std::vector< std::vector<T> > &v)
{
	for(int i=0; i<v.size(); ++i)
		print_vec(v[i]);
}







int main(void)
{
	std::vector<int> v;
	std::vector<std::vector<int> > vv;
	v = {1, 2, 3, 4};
	vv = {{0, 2}, {2,1}};
	print_vec(v);
	std::cout << std::endl;
	rintf_mat(vv);
	
}
