#include <iostream>
#include <cstddef>
#include <sys/time.h>
#include <vector>
#include "mpi.h"
#include "common.h"

using namespace std;


int main (int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
	int size, rank;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

/*
	MPI_Datatype MPI_Node;
	MPI_Type_contiguous(sizeof(Node), MPI_CHAR, &MPI_Node);
	MPI_Type_commit(&MPI_Node);
*/	

	vector<elem_type> send_vec;
	vector<elem_type> recv_vec;
	elem_type e;
	if(rank==0)
	{
		for(int i=0; i<10; ++i)
		{
			e = HEX8;
			send_vec.push_back(e);
		}
	}

	if(rank == 0)
	{
		cout << "Vector<int>::max_size " << send_vec.max_size() << endl;
		cout << "Elem: " << e << endl << endl;
	}
	switch(e)
	{
		case EDGE2:
			cout << "EDGE2" << endl;
			break;
		case QUAD4:
			cout << "QUAD4" << endl;
			break;
		case HEX8:
			cout << "HEX8" << endl;
			break;
		case TET4:
			cout << "TET4" << endl;
			break;
	}

	MPI_Status status;
	MPI_Request req;

	struct timeval tp;
	gettimeofday(&tp, NULL);
	long int ms_start = tp.tv_sec * 1000 + tp.tv_usec / 1000;
	if(rank==0)
		MPI_Isend(&send_vec[0], send_vec.size(), MPI_INT, 1, 0, MPI_COMM_WORLD, &req);
	if(rank==1)
	{
		int count;
		MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_INT, &count);
		cout << "count: " << count << endl;
		recv_vec.resize(count);
		MPI_Recv(&recv_vec[0], count, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	}
	gettimeofday(&tp, NULL);
	long int ms_end = tp.tv_sec * 1000 + tp.tv_usec / 1000;

	if(rank==1)
	{
		for(int i=0; i<recv_vec.size(); ++i)
			cout << rank << ": " << recv_vec[i] << endl;
	}
	//if(rank==0)
	//	cout << rank << ": Total Time: " << (ms_end-ms_start) << "ms" << endl;


	MPI_Finalize();
}
