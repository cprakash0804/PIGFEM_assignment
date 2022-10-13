/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _DENSE_MATRIX_H
#define	_DENSE_MATRIX_H

#include <cstdlib>
#include <map>
#include <vector>
#include <iostream>
#include <cstring>
#include "common.h"
#include "gsl_cblas.h"


template <class T>
class DenseMatrix
{
	public:
		/*
		 * Constructors for Square and Rectangular matricies
		*/
		DenseMatrix();
		DenseMatrix(size_t i);
		DenseMatrix(size_t i, size_t j);
		~DenseMatrix();

		/*
		 * Copy Constructor
		*/
		DenseMatrix(const DenseMatrix & other_mat);

		/*
		 * Clears all data from the matrix
		*/
		void clear();
	
		/*
		 * Copies data from other_mat into the current matrix
		*/
		void copy(const DenseMatrix<T>& other_mat);

		/*
		 * Function to access the number of rows and columns
		*/
		size_t n_rows() const {return _m;};
		size_t n_cols() const {return _n;};
		T* data() const {return _mat;};

		/*
		 * Function to change number of rows and columns of your matrix
		 * Adding rows/columns is easy, removing, check if there is data in there?
		*/
		void resize(size_t i, size_t j);
		void resize(size_t i, size_t j, bool copy);

		/*
		 * Functions to compute Matrix-Vector and Matrix-Matrix multiplys without copying a return vector
		*/
		void leftMultiply(std::vector<T>& res, const std::vector<T>& x) const;
		void leftMultiply(DenseMatrix<T>& res, const DenseMatrix<T>& x) const;

		void rightMultiply(std::vector<T>& res, const std::vector<T>& x) const;
		void rightMultiply(DenseMatrix<T>& res, const DenseMatrix<T>& x) const;

		/*
		 * Function to compute the transpose. The transpose will be stored in res
		*/
		void transpose(DenseMatrix<T>& res) const;



		/*
		 * Element access operators
		*/
		inline
		T& operator()(size_t i, size_t j)
		{
		    if(i<_m && j<_n) 
		    	return _mat[i*_n + j];
		    else
		    	err_message("Indicies out of bounds");
		}
		inline
		T operator()(size_t i, size_t j) const
		{
		    if(i<_m && j<_n) 
		    	return _mat[i*_n + j];
		    else
		    	err_message("Indicies out of bounds");
		}





		// MATH OPERATORS
		//=======================================================
		// Function to add alpha*X to the current matrix
		inline
		void AXPY(T alpha, DenseMatrix<T>& X)
		{
			if (_m!=X.n_rows() || _n!=X.n_cols())
				err_message("AXPY operations must be performed on matrices with mathing dimensions!");
			T* Xdat = X.data();
			for (id_type i=0; i<(_m*_n); ++i)
				_mat[i] += alpha * Xdat[i];
		}
	
		inline
		DenseMatrix<T> operator+ (const DenseMatrix<T>& rhs) const
		{
			if(rhs.n_rows()!=n_rows() || rhs.n_cols()!=n_cols())
				err_message("Matrix dimensions must match for addition!");
		
			DenseMatrix<T> ret(n_rows(), n_cols());
			T* ret_data = ret.data();
			T* rhs_data = rhs.data();
			for (id_type idx=0; idx<(_m*_n); ++idx)
				ret_data[idx] = _mat[idx] + rhs_data[idx];
			return ret;
		}
	
		inline
		DenseMatrix<T> operator* (const DenseMatrix<T>& rhs) const
		{
			DenseMatrix<T> ret(n_rows(), rhs.n_cols());
			leftMultiply(ret, rhs);
			return ret;
		}
		inline
		std::vector<T> operator* (const std::vector<T>& rhs) const
		{
			std::vector<T> ret(n_rows());
			leftMultiply(ret, rhs);
			return ret;
		}
		inline
		DenseMatrix<T> operator* (T rhs) const
		{
			DenseMatrix<T> ret(n_rows(), n_cols());
			T* ret_data = ret.data();
			for (id_type idx=0; idx<(_m*_n); ++idx)
				ret_data[idx] = rhs*_mat[idx];
			return ret;
		}

		inline
		void operator += (const DenseMatrix<T>& rhs)
		{
			if (n_rows()!=rhs.n_rows() || n_cols()!=rhs.n_cols())
				err_message("Row or Column mismatch in += operator!");

			T* rhs_data = rhs.data();
			for (id_type idx=0; idx<(_m*_n); ++idx)
				_mat[idx] += rhs_data[idx];
		}
		inline
		void operator *= (const T& rhs)
		{
			for (id_type idx=0; idx<(_m*_n); ++idx)
				_mat[idx] *= rhs;
		}
		inline
		void operator /= (const T& rhs)
		{
			for (id_type idx=0; idx<(_m*_n); ++idx)
				_mat[idx] /= rhs;
		}
		inline
		void operator = (const DenseMatrix<T>& other_mat)
		{
			copy(other_mat);
		}





		inline
		friend std::vector<T> operator* (std::vector<T> lhs, DenseMatrix<T>& rhs)
		{
			std::vector<T> ret(rhs.n_cols());
			rhs.rightMultiply(ret, lhs);
			return ret;
		}
		inline
		friend DenseMatrix<T> operator* (T lhs, DenseMatrix<T>& rhs)
		{
			DenseMatrix<T> ret(rhs.n_rows(), rhs.n_cols());
			T* ret_data = ret.data();
			T* rhs_data = rhs.data();

			id_type m = rhs.n_rows();
			id_type n = rhs.n_cols();
			for (id_type idx=0; idx<(m*n); ++idx)
				ret_data[idx] = lhs*rhs_data[idx];
			return ret;
		}

		/*
		 * Functions to print out matrix
		*/
		void print() const;
		void print(std::ostream& os) const;

	private:
		T* _mat;
		size_t _m;
		size_t _n;
};








// Constructors
template <typename T>
inline
DenseMatrix<T>::DenseMatrix()
{
	_mat = NULL;
	_m = 0;
	_n = 0;
}
template <typename T>
inline
DenseMatrix<T>::DenseMatrix(size_t i)
{
	if (i==0)
	{
		_m = 0;
		_n = 0;
		_mat = NULL;
	}
	else
	{
		_mat = new T[i*i];
		memset(_mat, 0, i*i*sizeof(T));
		_m = i;
		_n = i;
	}
}
template <typename T>
inline
DenseMatrix<T>::DenseMatrix(size_t i, size_t j)
{
	if(i==0 || j==0)
	{
		_m = 0;
		_n = 0;
		_mat = NULL;
	}
	else
	{
		_mat = new T[i*j];
		memset(_mat, 0, i*j*sizeof(T));
		_m = i;
		_n = j;
	}
}


// Clears all data from the matrix
template <typename T>
inline
void DenseMatrix<T>::clear()
{
	if (_mat != NULL)
		delete [] _mat;
	_mat = NULL;
	_m = 0;
	_n = 0;
}


// Copies data from other_mat into the current matrix
template <typename T>
inline
void DenseMatrix<T>::copy(const DenseMatrix<T>& other_mat)
{
	// Copy rown and column numbers
	size_t m = other_mat.n_rows();
	size_t n = other_mat.n_cols();
	resize(m, n, false);
	// Iterate over the other matrix and copy over relevant data
	T* other_data = other_mat.data();
	memcpy(_mat, other_data, m*n*sizeof(T));
}


// Copy Constructor
template <typename T>
inline
DenseMatrix<T>::DenseMatrix(const DenseMatrix<T>& other_mat)
{
	_mat = NULL; // Initialize to a NULL value
	_m = 0;
	_n = 0;
	copy(other_mat);
}


// Destructor
template <typename T>
inline
DenseMatrix<T>::~DenseMatrix()
{
	clear();
}



// Function to change number of rows and columns of your matrix
// Adding rows/columns is easy, removing, check if there is data in there?
template <typename T>
inline
void DenseMatrix<T>::resize(size_t i, size_t j)
{
	resize(i, j, true);
}
template <typename T>
inline
void DenseMatrix<T>::resize(size_t i, size_t j, bool copy)
{
	if(i==0 || j==0)
	{
		clear();
		return;
	}
	if (i==_m && j==_n) // Nothing to do
		return;

	T* new_mat = new T[i*j];
	if (copy)
	{
		if (j > _n) // new columns are bigger
		{
			if (i > _m) // new number of rows is bigger
			{
				for (id_type ii=0; ii<_m; ++ii)
				{
					memcpy(&new_mat[ii*j], &_mat[ii*_m], _n*sizeof(T));
					memset(&new_mat[ii*j + _n], 0, (j-_n)*sizeof(T));
				}
				for (id_type ii=_m; ii<i; ++ii)
					memset(&new_mat[ii*j], 0, j*sizeof(T));
			}
			else // new number of rows is equal or smaller
			{
				for (id_type ii=0; ii<_m; ++ii)
				{
					memcpy(&new_mat[ii*j], &_mat[ii*_m], _n*sizeof(T));
					memset(&new_mat[ii*j + _n], 0, (j-_n)*sizeof(T));
				}
			}
		}
		else // New number of columns are equal or smaller
		{
			if (i > _m) // new number of rows is bigger
			{
				for (id_type ii=0; ii<_m; ++ii)
					memcpy(&new_mat[ii*j], &_mat[ii*_m], j*sizeof(T));
				for (id_type ii=_m; ii<i; ++ii)
					memset(&new_mat[ii*j], 0, j*sizeof(T));
			}
			else // new number of rows is equal or smaller
			{
				for (id_type ii=0; ii<_m; ++ii)
					memcpy(&new_mat[ii*j], &_mat[ii*_m], j*sizeof(T));
			}
		}
	}

	_m = i;
	_n = j;
	if (_mat != NULL)
		delete [] _mat;
	_mat = new_mat;
}






























// BASIC MATH FUNCTIONS FOR SPECIFIC T TYPES
// ====================================================================================================
template <>
inline
void DenseMatrix<double>::AXPY(double alpha, DenseMatrix<double>& X)
{
	if (_m!=X.n_rows() || _n!=X.n_cols())
		err_message("AXPY operations must be performed on matrices with mathing dimensions!");

	double* Xdat = X.data();
	cblas_daxpy(_m*_n, alpha, Xdat, 1, _mat, 1);
}
template <>
inline
void DenseMatrix<float>::AXPY(float alpha, DenseMatrix<float>& X)
{
	if (_m!=X.n_rows() || _n!=X.n_cols())
		err_message("AXPY operations must be performed on matrices with mathing dimensions!");

	float* Xdat = X.data();
	cblas_saxpy(_m*_n, alpha, Xdat, 1, _mat, 1);
}

template<>
inline
DenseMatrix<double> DenseMatrix<double>::operator+ (const DenseMatrix<double>& rhs) const
{
	if(rhs.n_rows()!=_m || rhs.n_cols()!=_n)
		err_message("Matrix dimensions must match for addition!");

	DenseMatrix<double> ret = *this; // Copy constructor
	cblas_daxpy(_m*_n, 1.0, ret.data(), 1, rhs.data(), 1);
	return ret;
}
template<>
inline
DenseMatrix<float> DenseMatrix<float>::operator+ (const DenseMatrix<float>& rhs) const
{
	if(rhs.n_rows()!=_m || rhs.n_cols()!=_n)
		err_message("Matrix dimensions must match for addition!");

	DenseMatrix<float> ret = *this; // Copy constructor
	cblas_saxpy(_m*_n, 1.0, ret.data(), 1, rhs.data(), 1);
	return ret;
}

template<>
inline
DenseMatrix<double> DenseMatrix<double>::operator* (const DenseMatrix<double>& rhs) const
{
	DenseMatrix<double> ret(n_rows(), rhs.n_cols());

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, _m, rhs.n_cols(), _n, 1.0, _mat, _n, rhs.data(), rhs.n_cols(), 0.0, ret.data(), rhs.n_cols());
	return ret;
}
template<>
inline
DenseMatrix<float> DenseMatrix<float>::operator* (const DenseMatrix<float>& rhs) const
{
	DenseMatrix<float> ret(n_rows(), rhs.n_cols());
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, _m, rhs.n_cols(), _n, 1.0, _mat, _n, rhs.data(), rhs.n_cols(), 0.0, ret.data(), rhs.n_cols());
	return ret;
}

template<>
inline
std::vector<double> DenseMatrix<double>::operator* (const std::vector<double>& rhs) const
{
	std::vector<double> ret(n_rows());
	cblas_dgemv(CblasRowMajor, CblasNoTrans, _m, _n, 1.0, _mat, _n, rhs.data(), 1, 0.0, ret.data(), 1);
	return ret;
}
template<>
inline
std::vector<float> DenseMatrix<float>::operator* (const std::vector<float>& rhs) const
{
	std::vector<float> ret(n_rows());
	cblas_sgemv(CblasRowMajor, CblasNoTrans, _m, _n, 1.0, _mat, _n, rhs.data(), 1, 0.0, ret.data(), 1);
	return ret;
}

template<>
inline
DenseMatrix<double> DenseMatrix<double>::operator* (double rhs) const
{
	DenseMatrix<double> ret = *this; // Copy constructor
	cblas_dscal(_m*_n, rhs, ret.data(), 1);
	return ret;
}
template<>
inline
DenseMatrix<float> DenseMatrix<float>::operator* (float rhs) const
{
	DenseMatrix<float> ret = *this; // Copy constructor
	cblas_sscal(_m*_n, rhs, ret.data(), 1);
	return ret;
}

template<>
inline
void DenseMatrix<double>::operator += (const DenseMatrix<double>& rhs)
{
	if (n_rows()!=rhs.n_rows() || n_cols()!=rhs.n_cols())
		err_message("Row or Column mismatch in += operator!");

	cblas_daxpy(_m*_n, 1.0, rhs.data(), 1, _mat, 1);
}
template<>
inline
void DenseMatrix<float>::operator += (const DenseMatrix<float>& rhs)
{
	if (n_rows()!=rhs.n_rows() || n_cols()!=rhs.n_cols())
		err_message("Row or Column mismatch in += operator!");

	cblas_saxpy(_m*_n, 1.0, rhs.data(), 1, _mat, 1);
}

template<>
inline
void DenseMatrix<double>::operator *= (const double& rhs)
{
	cblas_dscal(_m*_n, rhs, _mat, 1);
}
template<>
inline
void DenseMatrix<float>::operator *= (const float& rhs)
{
	cblas_sscal(_m*_n, rhs, _mat, 1);
}

template<>
inline
void DenseMatrix<double>::operator /= (const double& rhs)
{
	cblas_dscal(_m*_n, 1.0/rhs, _mat, 1);
}
template<>
inline
void DenseMatrix<float>::operator /= (const float& rhs)
{
	cblas_sscal(_m*_n, 1.0/rhs, _mat, 1);
}
















































/*
 * leftMultiply functions perform res = A*x where A is the current matrix and x is the second argument
 * x can be a std::vector, DenseVector, or DenseMatrix
 * Vectors are assumed to be column vectors
*/
template <typename T>
inline
void DenseMatrix<T>::leftMultiply(std::vector<T>& res, const std::vector<T>& x) const
{
    if(_n != x.size())
    	err_message("Vector size must match matrix dimensions");

	res.resize(n_rows());

    for(size_t ii=0; ii<_m; ++ii){
    	res[ii] = 0;
        for(size_t jj=0; jj<_n; ++jj)
            res[ii] += _mat[ii*_n + jj] * x[jj];
    }
}
template <typename T>
inline
void DenseMatrix<T>::leftMultiply(DenseMatrix<T>& res, const DenseMatrix<T>& x) const
{
    if(n_cols() != x.n_rows())
    	err_message("Matrix inner dimensions must match.");

	res.resize(n_rows(), x.n_cols(), false);

    id_type ncols = x.n_cols();
    T* res_data = res.data();
    T* x_data = x.data();

    for(size_t ii=0; ii<_m; ++ii)
	{
		for(size_t jj=0; jj<x.n_cols(); ++jj)
		{
		    res_data[ii*ncols + jj] = 0;
		    for(size_t it=0; it<_n; ++it)
		        res_data[ii*ncols + jj] += _mat[ii*_n + it] * x_data[it*ncols + jj];
		}
    }
}


/*
 * rightMultiply functions perform res = x*A where A is the current matrix and x is the second argument
 * x can be a std::vector, or DenseMatrix
 * Vectors are assumed to be row vectors
*/
template <typename T>
inline
void DenseMatrix<T>::rightMultiply(std::vector<T>& res, const std::vector<T>& x) const
{
    if(_m != x.size())
    	err_message("Vector size must match matrix dimensions");

	res.resize(n_cols());

	for(size_t jj=0; jj<_n; ++jj)
	{
		res[jj] = 0;
		for(size_t ii=0; ii<_m; ++ii)
		    res[jj] += _mat[ii*_n + jj] * x[ii];
	}
}
template <typename T>
inline
void DenseMatrix<T>::rightMultiply(DenseMatrix<T>& res, const DenseMatrix<T>& x) const
{
    if(n_rows() != x.n_cols())
    	err_message("Matrix inner dimensions must match.");

	res.resize(x.n_rows(), n_cols(), false);

    id_type ncols = x.n_cols();
    T* res_data = res.data();
    T* x_data = x.data();

    for(size_t ii=0; ii<x.n_rows(); ++ii)
	{
		for(size_t jj=0; jj<n_cols(); ++jj)
		{
		    res_data[ii*ncols + jj] = 0;
		    for(size_t it=0; it<_m; ++it)
		        res_data[ii*ncols + jj] += x_data[ii*ncols + it] * _mat[it*_n + jj];
		}
    }
}


template <typename T>
inline
void DenseMatrix<T>::transpose(DenseMatrix<T>& res) const
{
	res.resize(n_cols(), n_rows(), false);
	T* res_data = res.data();

	for(size_t ii=0; ii<_m; ++ii)
		for(size_t jj=0; jj<_n; ++jj)
			res_data[jj*_m + ii] = _mat[ii*_n + jj];
}



template <typename T>
inline
void DenseMatrix<T>::print() const
{
	print(std::cout);
}
template <typename T>
inline
void DenseMatrix<T>::print(std::ostream& os) const
{
	os.precision(16);
	os << "[";
	for(size_t ii=0; ii<_m; ++ii)
	{
		if(ii!=0)
			os << " ";
		for(size_t jj=0; jj<_n; ++jj)
		{
			os << _mat[ii*_n + jj];
			if(jj != (_n-1))
				os << ", ";
		}
		if(ii != (_m-1))
			os << std::endl;
		else
			os << "]" << std::endl;
	}
}

#endif	/* _DENSE_MATRIX_H */
