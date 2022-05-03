#include <mpi.h>
#include <iostream>

using namespace std;

void abortAll(int rc, string str){
	cout<<str<<endl;
	MPI_Abort(MPI_COMM_WORLD, rc);
	exit(0);
}

int main(int argc, char* argv[]){
/*Setup part*/
	int rc, rank, numprocs;
	MPI_Status status;
	if (rc = MPI_Init(&argc, &argv))
		abortAll(-1, "Ошибка запуска MPI");
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double startTime;
	if (rank==0)
		startTime = MPI_Wtime();

	int ball=0;
	if (rank==0){
		MPI_Send(&ball, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
		MPI_Recv(&ball, 1, MPI_INT, numprocs-1, 0, MPI_COMM_WORLD, &status);
		cout<<ball<<endl;
		MPI_Finalize();
		return 0;
	}

	MPI_Recv(&ball, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &status);
	ball+=1;
	MPI_Send(&ball, 1, MPI_INT, rank==numprocs-1?0:rank+1, 0, MPI_COMM_WORLD);

	MPI_Finalize();
	return 0;
}
