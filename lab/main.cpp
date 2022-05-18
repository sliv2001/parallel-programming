#include <mpi.h>
#include <iostream>
#include <unistd.h>
#include "matplotlibcpp.h"

#define MCW MPI_COMM_WORLD

using namespace std;
namespace plt = matplotlibcpp;

double t0(double t){
	return 1;
}

double x0(double t){
	return 2;
}

double f(double t, double x){
	return 3;
}

void abortAll(int rc, string str){
	cout<<str<<endl;
	MPI_Abort(MCW, rc);
	exit(0);
}

double seconds(struct timespec *start, struct timespec *end){
	return end->tv_sec-start->tv_sec+(end->tv_nsec-start->tv_nsec)/
			(double)1000000000;
}

void checkTime(int rank, int length){
	MPI_Status status;
	char buffer[length];
	struct timespec start, end;
	if (rank==0){
		sleep(1);
		clock_gettime(CLOCK_REALTIME, &start);
		MPI_Send(buffer, length, MPI_CHAR, 1, 0, MCW);
		MPI_Recv(buffer, length, MPI_CHAR, 1, 0, MCW, &status);
		clock_gettime(CLOCK_REALTIME, &end);
		cout<<"Time for transmission "<<seconds(&start, &end)/2<<endl;
	}
	else if (rank==1){
		MPI_Recv(buffer, length, MPI_CHAR, 0, 0, MCW, &status);
		MPI_Send(buffer, length, MPI_CHAR, 0, 0, MCW);
	}
}

double** allocate(int N, int M){
	double** res = (double**)malloc(sizeof(double*)*N);
	if (res==NULL)
		return res;
	for (int i=0; i<N; i++)
		if ((res[i]=(double*)malloc(sizeof(double)*M))==NULL)
			return NULL;
	return res;
}

void show(double** field, int N, int M, double h, double tau, double t){
	int td=t/tau;
//	cout<<td<<endl;
	if (td>N)
		abortAll(-5, "excessive index");
	vector<double> f(M);
	vector<double> x(M);
	for (int i=0; i<M; i++){
		f[i]=field[td][i];
		x[i]=i*h;
	}
	plt::plot(x, f);
	plt::show();
}

void release(double** field, int N, int M){
	for (int i=0; i<N; i++)
		free(field[i]);
	free(field);
}

void init_table(	double** field,
			int N,
			int M,
			double h,
			double tau,
			double (*x0)(double x),
			double (*t0)(double t)){
	for (int i=0; i<N; i++)
		field[i][0]=t0(tau*i);
	for (int i=0; i<M; i++)
		field[0][i]=x0(h*i);
}

void solveBrick(double** u,
		int N,
		int M,
		double h,
		double tau,
		int rank,
		int numprocs){
	int brickSize, thisSize;
	if ((M-1)%numprocs==0){
		brickSize=(M-1)/numprocs;
		thisSize=brickSize;
	}
	else{
		brickSize=(M-1)/(numprocs-1);
		thisSize=brickSize;
		if (rank+1==numprocs)
			thisSize=(M-1)%(numprocs-1);
	}
	int a=1+rank*brickSize;
	int b=1+rank*brickSize+thisSize;
	cout<<"a="<<a<<endl;
	cout<<"b="<<b<<endl;
	for (int n=1; n<N-1; n++){
//		cout<<"n="<<n<<endl;
		for (int m=a; m<(rank+1==numprocs?b-1:b); m++){
			//cout<<m<<endl;
			u[n+1][m]=f(n*tau, m*h)+u[n][m]-
					tau/2/h*(u[n][m+1]-u[n][m-1])+
					tau*tau/2/h/h*(u[n][m+1]-2*u[n][m]+
					u[n][m-1]);
		}
		if (rank+1==numprocs)
			u[n+1][b-1]=u[n][b-1]+f(n*tau, (b-1)*h)*tau-
					tau/h*(u[n][b-1]-u[n][b-2]);
	}
}

int sequential(	double h,
		double tau,
		double X,
		double T){
	int M=X/h;
	int N=T/tau;
	cout<<"N="<<N<<endl;
	cout<<"M="<<M<<endl;
	double** field = allocate(N, M);
	if (field==NULL){
	 	abortAll(-4, "allocation error");
	}
	init_table(field, N, M, h, tau, x0, t0);
	solveBrick(field, N, M, h, tau, 0, 1);
	show(field, N, M, h, tau, 9.9);
	release(field, N, M);
}

int main(int argc, char* argv[]){
/*Setup part*/
	int rc, rank, numprocs;
	double h, tau, T, X, t;
	MPI_Status status;
	if (rc = MPI_Init(&argc, &argv))
		abortAll(-1, "Ошибка запуска MPI");
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	/*if (numprocs<2)
		abortAll(-2, "Not enough executors");*/
	if (argc<5)
		abortAll(-3, "Not enough args");
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double startTime;
	if (rank==0)
		startTime = MPI_Wtime();
	h=atof(argv[1]);
	tau = atof(argv[2]);
	X = atof(argv[3]);
	T = atof(argv[4]);
//	checkTime(rank, 128);

/*Transmission*/
	

/*Parallel part*/
	

/*Collection*/
	if (rank==0){
		cout<<"time elapsed: "<<MPI_Wtime()-startTime<<endl;
	}

/*Sequential part*/
	if (rank==0)
		sequential(h, tau, X, T);

	MPI_Finalize();
	return 0;
}
