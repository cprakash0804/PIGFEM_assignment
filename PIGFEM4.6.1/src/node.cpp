/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "node.h"
#include <iostream>


Node::Node()
	:_x(0), _y(0), _z(0), _global_id(0)
{
}

Node::Node(double x, double y, double z, id_type id)
	: _x(x), _y(y), _z(z), _global_id(id)
{
}

void Node::clear()
{
	_x = 0;
	_y = 0;
	_z = 0;
	_global_id = 0;
}

void Node::copy(const Node & other_node)
{
	clear();
	_x = other_node(0);
	_y = other_node(1);
	_z = other_node(2);
	_global_id = other_node.get_id();
}

Node::Node(const Node & other_node)
{
	copy(other_node);
}

void Node::set_coords(double x, double y=0.0, double z=0.0)
{
	_x = x;
	_y = y;
	_z = z;
}

void Node::set_coord(double val, int dir)
{
	if(dir==0)
		_x = val;
	else if(dir==1)
		_y = val;
	else if(dir==3)
		_z = val;
	else
		err_message("Invalid coordinate direcion.");
}

std::vector<double> Node::get_coords()
{
	std::vector<double> c;
	c.push_back(_x);
	c.push_back(_y);
	c.push_back(_z);
	return c;
}

void Node::get_coords(double &x, double &y, double &z)
{
	x = _x;
	y = _y;
	z = _z;
}

double Node::operator ()(id_type i) const
{
	if(i > 2)
		err_message("Node index must be less than the dimension of the node.");
	else
	{
		if(i==0)
			return _x;
		else if(i==1)
			return _y;
		else
			return _z;
	}
}

bool Node::operator <(Node& other_node)
{
	if(_global_id < other_node.get_id())
		return true;
	else
		return false;
}

bool Node::operator >(Node& other_node)
{
	if(_global_id < other_node.get_id())
		return false;
	else
		return true;
}

bool Node::operator ==(Node& other_node)
{
	if(_global_id==other_node.get_id() &&
	   _x==other_node(0) && _y==other_node(1) && _z==other_node(2))
		return true;
	else
		return false;
}

bool Node::less(Node& node1, Node& node2)
{
	return (node1 < node2);
}


void Node::print(std::ostream & os) const
{
	os << "id: " << _global_id << " Coordinates: ";
	os << "(" << _x << "," << _y << "," << _z << ")" << std::endl;
}

/*
void Node::MPI_Create_node_type(MPI_Comm comm)
{
	const int nitems = 4;
	int blocklengths[4] = {1,1,1,1};
	MPI_Datatype types[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
	MPI_Datatype MPI_Node;
	MPI_Aint offsets[4], base;
	offsets[0] = offsetof(Node, Node::_x);
	offsets[1] = offsetof(Node, Node::_y);
	offsets[2] = offsetof(Node, Node::_z);
	offsets[3] = offsetof(Node, Node::_global_id);
	MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPI_Node);
	MPI_Type_commit(&MPI_Node);
}
*/
