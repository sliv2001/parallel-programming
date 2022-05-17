#include <mpi.h>
#include <iostream>
#include <unistd.h>
#include "matplotlibcpp.h"

#define MCW MPI_COMM_WORLD

using namespace std;

void abortAll(int rc, string str){
	cout<<str<<endl;
	MPI_Abort(MCW, rc);
	exit(0);
}

double seconds(struct timespec *start, struct timespec *end){
	return end->tv_sec-start->tv_sec+(end->tv_nsec-start->tv_nsec)/
			(double)1000000000;
}

void checkTime(int rank){
	MPI_Status status;
	char buffer=0;
	struct timespec start, end;
	if (rank==0){
		sleep(1);
		clock_gettime(CLOCK_REALTIME, &start);
		MPI_Send(&buffer, 1, MPI_CHAR, 1, 0, MCW);
		MPI_Recv(&buffer, 1, MPI_CHAR, 1, 0, MCW, &status);
		clock_gettime(CLOCK_REALTIME, &end);
		cout<<"Time for transmission "<<seconds(&start, &end)/2<<endl;
	}
	else if (rank==1){
		MPI_Recv(&buffer, 1, MPI_CHAR, 0, 0, MCW, &status);
		MPI_Send(&buffer, 1, MPI_CHAR, 0, 0, MCW);
	}
}

namespace plt = matplotlibcpp;

int main(int argc, char* argv[]){
/*Setup part*/
	int rc, rank, numprocs;
	MPI_Status status;
	if (rc = MPI_Init(&argc, &argv))
		abortAll(-1, "Ошибка запуска MPI");
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	if (numprocs<2)
		abortAll(-2, "Not enough executors");
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double startTime;
	if (rank==0)
		startTime = MPI_Wtime();

	checkTime(rank);
/*Transmission*/
	

/*Parallel part*/
	

/*Collection*/
	if (rank==0){
		cout<<"time elapsed: "<<MPI_Wtime()-startTime<<endl;
	}
	MPI_Finalize();
	return 0;
}
