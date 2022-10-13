#include <iostream>
#include <cstddef>
#include <sys/time.h>
#include <vector>
#include "node.h"
#include "mpi.h"
#include "common.h"
using namespace std;

struct Node_data{
	double _x;
	double _y;
	double _z;
	unsigned int _global_id;
};

int main (int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Datatype MPI_Node;
	MPI_Type_contiguous(sizeof(Node), MPI_CHAR, &MPI_Node);
	MPI_Type_commit(&MPI_Node);
	

	vector<Node> send_vec;
	vector<Node> recv_vec;
	if(rank==0)
	{
		for(int i=0; i<10; ++i)
		{
			Node n(i/3., i/5., i/7., i);
			send_vec.push_back(n);
		}
	}
	if(rank==0)
		cout << "vector<Node>.max_size " << send_vec.max_size() << endl;

	MPI_Status status;
	MPI_Request req;

	struct timeval tp;
	gettimeofday(&tp, NULL);
	long int ms_start = tp.tv_sec * 1000 + tp.tv_usec / 1000;
	if(rank==0)
		MPI_Isend(&send_vec[0], send_vec.size(), MPI_Node, 1, 0, MPI_COMM_WORLD, &req);
	if(rank==1)
	{
		int count;
		MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_Node, &count);
		cout << "count: " << count << endl;
		recv_vec.resize(count);
		MPI_Recv(&recv_vec[0], count, MPI_Node, 0, 0, MPI_COMM_WORLD, &status);
	}
	gettimeofday(&tp, NULL);
	long int ms_end = tp.tv_sec * 1000 + tp.tv_usec / 1000;

	if(rank==1)
	{
		for(int i=0; i<recv_vec.size(); ++i)
			cout << rank << ": " << recv_vec[i];
		
		cout << endl;
		std::vector<Node*> _nodes;
		for(int i=0; i<recv_vec.size(); ++i)
			_nodes.push_back(new Node(recv_vec[i]));

		for(int i=0; i<_nodes.size(); ++i)
			cout << rank << ": " << (*_nodes[i]);

		for(int i=0; i<_nodes.size(); ++i)
			delete _nodes[i];
	}
	//if(rank==0)
	//	cout << rank << ": Total Time: " << (ms_end-ms_start) << "ms" << endl;


	MPI_Finalize();
}
