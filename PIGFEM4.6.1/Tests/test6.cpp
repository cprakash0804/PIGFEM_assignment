#include <iostream>
#include <string.h>
#include <vector>
#include "material.h"
#include "material_lei.h"
#include "mpi.h"
using namespace std;


int main(int argc, char**argv)
{
	MPI_Init(&argc, &argv);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Datatype MPI_MATERIAL;
	MPI_Type_vector(MAX_MATERIAL_SIZE, 1, 1, MPI_CHAR, &MPI_MATERIAL);
	MPI_Type_commit(&MPI_MATERIAL);
	if(rank==0)
	{
		std::vector<char> tot_vec;
		std::vector<char> vec(MAX_MATERIAL_SIZE);
		char* buf;
		Material* mat = new LinearElasticIsotropicMaterial;
		mat->set_parameter("E", 1000056);
		mat->set_parameter("nu", 0.25);
		mat->set_name("Something");
		Material* mat2 = new LinearElasticIsotropicMaterial;
		mat2->set_parameter("E", 12346234);
		mat2->set_parameter("nu", 0.33);
		mat2->set_name("Something2");
	
		// Pack
		buf = &vec[0];
		mat->pack(buf);
		for(int i=0; i<MAX_MATERIAL_SIZE; ++i)
			tot_vec.push_back(vec[i]);
		mat2->pack(buf);
		for(int i=0; i<MAX_MATERIAL_SIZE; ++i)
			tot_vec.push_back(vec[i]);
		buf = &tot_vec[0];
		MPI_Send(buf, 2, MPI_MATERIAL, 1, 0, MPI_COMM_WORLD);
		delete mat;
		delete mat2;
	}
	else
	{
		std::vector<Material*> _materials;
		char* buf;
		Material* mat;
		material_type mat_type;
		size_t mat_size;
		size_t offset = 0;
		int count;
		MPI_Status stat;
		std::vector<char> recv_vec;
		MPI_Probe(0, 0, MPI_COMM_WORLD, &stat);
		MPI_Get_count(&stat, MPI_MATERIAL, &count);
		recv_vec.resize(count*MAX_MATERIAL_SIZE);
		buf = &recv_vec[0];
		MPI_Recv(buf, count, MPI_MATERIAL, 0, 0, MPI_COMM_WORLD, &stat);
		cout << "Recieve Vector size: " << recv_vec.size() << endl;
		for(int i=0; i<count; ++i)
		{
			offset = MAX_MATERIAL_SIZE*i;
			cout << "Offset: " << offset;
			memcpy(&mat_type, buf+offset, sizeof(material_type));
			cout << " Mat Type: " << mat_type;
			offset += sizeof(material_type);
			memcpy(&mat_size, buf+offset, sizeof(size_t));
			cout << " Mat Size: " << mat_size << endl;
			offset += sizeof(size_t);
			if(mat_type==LINEAR_ELASTIC_ISOTROPIC)
			{
				mat = new LinearElasticIsotropicMaterial;
				memcpy(mat, buf+offset, mat_size);
				_materials.push_back(mat);
			}
		}
		for(int i=0; i<count; ++i)
		{
			cout << _materials[i] << endl;
			cout << " Material " << _materials[i]->get_name() << ": mat_type: " << _materials[i]->get_type() << " mat_E: " << _materials[i]->get_parameter("E") << " mat2_K: " << _materials[i]->get_parameter("K") << endl;
			delete _materials[i];
		}
	}


	MPI_Finalize();
}
