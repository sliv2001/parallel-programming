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
	if (rc = MPI_Init(&argc, &argv))
		abortAll(-1, "Ошибка запуска MPI");
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double startTime;
	if (rank==0)
		startTime = MPI_Wtime();

/*Transmission*/
	int N = atoi(argv[1]);
	long double complete;
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

/*Parallel part*/
	long double res=0.0;
	for (int i=rank; i<N; i+=numprocs){
		res+=(1/((double)i+1));
	}
/*Collection*/
	MPI_Reduce(&res, &complete, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (rank==0){
		cout<<complete<<endl;
		cout<<"time elapsed: "<<MPI_Wtime()-startTime<<endl;
	}
	MPI_Finalize();
	return 0;
}
