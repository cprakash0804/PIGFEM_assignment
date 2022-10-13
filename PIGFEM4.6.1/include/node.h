/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _NODE_H_
#define _NODE_H_
#include "common.h"
#include <iostream>
#include <vector>



class Node
{
	private:
		double _x;
		double _y;
		double _z;
		id_type _global_id;
	public:
		Node();
		Node(double x, double y, double z, id_type id);
		Node(const Node & other_node);

		void clear();
		void copy(const Node & other_node);

		void set_coords(double x, double y, double z);
		void set_coord(double val, int dir);
		std::vector<double> get_coords();
		void get_coords(double &x, double &y, double &z);
		
		void set_id(id_type id) {_global_id = id;};
		id_type get_id() const {return _global_id;};

		double operator ()(id_type dir) const;
		bool operator <(Node& other_node);
		bool operator ==(Node& other_node);
		bool operator >(Node& other_node);
		bool less(Node& node1, Node& node2);

		friend std::ostream& operator<<(std::ostream& os, const Node& node)
		{
			node.print(os);
			return os;
		}
		void print(std::ostream & os) const;
};


#endif
