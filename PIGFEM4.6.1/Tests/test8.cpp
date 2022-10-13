#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include "node.h"
#include "elem.h"
#include "quad4.h"
#include "dense_matrix.h"
using namespace std;

template <class T, class iterator_type>
void clean_ElemNodeElemMatvec(std::vector<T*>& obj, iterator_type comp_begin, iterator_type comp_end)
{
	int new_size = 0;
	int i = 0;
	while(i!=obj.size())
	{
		if( std::find(comp_begin, comp_end, obj[i]->get_id())!=comp_end ) // This one does belong here
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

int main(int argc, char**argv)
{
	Node* n= new Node(1, 1, 1, 1);
	vector<Node*> nodes = {n, n, n, n};
	vector<unsigned int> ids = {1,2,3,4,5,6,7,8,9,10};
	vector<unsigned int> own = {2,4,3,7,9};
	vector<Elem*> elems;
	for(unsigned int i=0; i<ids.size(); ++i)
		elems.push_back( new Quad4(nodes, ids[i]) );

	for(unsigned int i=0; i<elems.size(); ++i)
		cout << *elems[i];
	clean_ElemNodeMatvec(elems, own.begin(), own.end());
	for(unsigned int i=0; i<elems.size(); ++i)
		cout << *elems[i];
	for(unsigned int i=0; i<elems.size(); ++i)
		delete elems[i];



/*
	DenseMatrix<double> a;
	a.resize(4, 4);
	for(unsigned int i=0; i<a.n_rows(); ++i)
	{
		for(unsigned int j=0; j<a.n_cols(); ++j)
		{
			a(i, j) = i*j;
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
