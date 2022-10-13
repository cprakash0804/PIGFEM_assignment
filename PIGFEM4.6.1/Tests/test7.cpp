#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include "dense_matrix.h"
using namespace std;

template <class T>
void clean_vec(std::vector<T*>& obj, std::vector<int>& comp)
{
	int new_size = 0;
	int i = 0;
	while(i!=obj.size())
	{
		if( std::find(comp.begin(), comp.end(), obj[i]->get_id())!=comp.end() ) // This one does belong here
		{
			std::swap(obj[new_size], obj[i]);
			new_size++;
		}
		i++;
	}
	for(int j=new_size; j<obj.size(); ++j)
	{
		cout << "Deleting element" << obj[j]->get_id() << endl;
		delete obj[j];
	}
	obj.resize(new_size);
}

void manipulate_a(DenseMatrix<double>& a)
{
	std::cout << "\t Size of a is " << a.n_rows() << "x" << a.n_cols() << std::endl;
	a.resize(4, 5);
	std::cout << "\t After resizing size of a is " << a.n_rows() << "x" << a.n_cols() << std::endl;
	for(unsigned int i=0; i<a.n_rows(); ++i)
		for(unsigned int j=0; j<a.n_cols(); ++j)
			a(i, j) = i%(j+1) + j*15/(i+1);
}

int main(int argc, char**argv)
{
	vector<int> dofs = {5, 6, 7, 8, 9, 10};
	vector<double> vals = {1.3, 5.1, -2.1, 0.02, 0.00003, 1.0};

	std::vector<int>::iterator it = std::find(dofs.begin(), dofs.end(), 11);
	int idx = it - dofs.begin();
	if(it!=dofs.end())
		cout << "Corresponding value is " << vals[idx] << endl;
	else
		cout << "Corresponding value not found\n";
/*
	vector<unsigned int> node = {1,2,3,4};
	vector<unsigned int> ids = {1,2,3,4,5,6,7,8,9,10};
	vector<Elem*> elems;
	for(unsigned int i=0; i<ids.size(); ++i)
		elems.push_back( new Quad4(nodes, ids[i]) );

	for(unsigned int i=0; i<elems.size(); ++i)
		cout << *elems[i];
	clean_vec(elems, own);
	for(unsigned int i=0; i<elems.size(); ++i)
		cout << *elems[i];
	for(unsigned int i=0; i<elems.size(); ++i)
		delete elems[i];
*/


/*
	DenseMatrix<double> a;
	manipulate_a(a);
	a.resize(4, 4);
	for(unsigned int i=0; i<a.n_rows(); ++i)
	{
		for(unsigned int j=0; j<a.n_cols(); ++j)
		{
			cout << "a(" << i << "," << j << ") = " << a(i, j) << endl;
		}
	}
*/
/*
	vector<int> a(0);
	set<int> b = {2, 3, 6};
	a.insert(a.end(), b.begin(), b.end());
	vector<int>::iterator it = a.begin();
	for(; it!=a.end(); ++it)
		cout << *it << " ";
	cout << endl << a.size() << endl;
*/
	/*
	set<int> res(a.size()+b.size());
	set<int>::iterator it, end;
	end = merge(a.begin(), a.end(), b.begin(), b.end(), res.begin());
	a.resize(end-a.begin());
	for(it=res.begin(); it!=end; ++it)
		cout << *it << " ";
	cout << endl;
	cout << res.size() << endl;
	*/
	return 0;
}
