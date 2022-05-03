#include <mpi.h>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]){
	int rc, rank, numprocs;
	if (rc = MPI_Init(&argc, &argv)){
		cout<<"Ошибка запуска MPI"<<endl;
		MPI_Abort(MPI_COMM_WORLD, rc);
	}

	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	cout<<"Hello world from "<<rank<<"-th of "<<numprocs<<endl;

	MPI_Finalize();
	return 0;
}
