#include <iostream>
#include <vector>
#include "elem.h"
#include "quad4.h"
#include "node.h"
using namespace std;

int main(void)
{
	Elem* el;
	std::vector<Node*> nodes;
	for(int i=0; i<4; ++i)
	{
		double x, y;
		if(i==0 or i==3)
			x = 0.0;
		else
			x = 1.0;

		if(i==0 or i==1)
			y = 0.0;
		else
			y = 1.0;
		nodes.push_back(new Node(x, y, 0.0, i));
	}
	el = new Quad4(nodes, 0);

	cout << "All assigned!" << endl;

	// Print all the nodes out
	cout << (*el) << endl;

	// Test node data assignment functions
	Node* node = new Node;
	node->set_coords(5, 1, -2);
	node->set_id(10);

	// Assign new node to element in various ways
	el->set_node(2) = node;
	cout << (*el) << endl;

	node->set_coords(5, 5, 5);
	node->set_id(5);
	el->set_node(1, node);
	cout << (*el) << endl;

	el->set_nodes(nodes);
	cout << (*el) << endl;

	std::vector<Node*> nodes2 = el->get_nodes();
	for(int i=0; i<4; ++i)
		cout << (*nodes2[i]);

	
	delete node;
	for(int i=0; i<4; ++i)
		delete nodes[i];

	delete el;
	
}
