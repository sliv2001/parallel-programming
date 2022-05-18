#include <mpi.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <vector>
#include "matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;

MPI_Status status;
MPI_Request request;

int N, M;
double X, T, tau, h;
int brickSize, lowSize;

double t0(double x){
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
	MPI_Abort(MPI_COMM_WORLD, rc);
	exit(0);
}

vector<vector<double>> solve(int rank, int numprocs){
	cout<<rank<<" "<<numprocs<<endl;
	if ((M-1)%numprocs==0){
		brickSize=(M-1)/numprocs;
		lowSize = brickSize;
	}
	else {
		brickSize = (M-1)/(numprocs-1);
		lowSize = (M-1)%(numprocs-1);
	}
	vector<vector<double>> field(N);
	int sz;
	if (rank==0)
		sz=brickSize+1;
	else if (rank==numprocs-1)
		sz=lowSize;
	else
		sz=brickSize;
	for (int i=0; i<N; i++)
		field[i] = vector<double>(sz);
	for (int i=0; i<(rank==numprocs-1?lowSize:brickSize); i++)
		field[0][i]=t0(i*h);
	if (rank==0)
		for (int i=0; i<N; i++){
			field[i][0]=x0(i*tau);
		}
	for (int n=0; n<N-1; n++){
		double left, right;
		if (rank!=numprocs-1){
			if (rank!=0)
				MPI_Recv(&left, 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
			else{
				left=x0(n*tau);
				field[n][0]=left;
			}
			field[n+1][0]=f(n*tau, (rank*brickSize)*h)+field[n][0]-
					tau/2/h*(field[n][1]-left)+
					tau*tau/2/h/h*(field[n][1]-2*field[n][0]+
					left);
			cout<<field[n+1][0]<<endl;
			if (rank!=0)
				MPI_Isend(&field[n+1][0], 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &request);
			for (int m=1; m<brickSize-1; m++){
				field[n+1][m]=f(n*tau, (rank*brickSize+m)*h)+field[n][m]-
					tau/2/h*(field[n][m+1]-field[n][m-1])+
					tau*tau/2/h/h*(field[n][m+1]-2*field[n][m]+
					field[n][m-1]);
			}
			MPI_Recv(&right, 1, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &status);
			int m=brickSize-1;
			field[n+1][m]=f(n*tau, (rank*brickSize+m)*h)+field[n][m]-
					tau/2/h*(right-field[n][m-1])+
					tau*tau/2/h/h*(right-2*field[n][m]+
					field[n][m-1]);
			MPI_Isend(&field[n+1][m], 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request);
		}
		else {
			if (rank!=0){
				MPI_Recv(&left, 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
				field[n+1][0]=f(n*tau, (rank*brickSize)*h)+field[n][0]-
						tau/2/h*(field[n][1]-left)+
						tau*tau/2/h/h*(field[n][1]-2*field[n][0]+
						left);

			}
			else{
				left=x0(n*tau);
				field[n][0]=left;
			}
			if (rank!=0)
				MPI_Isend(&field[n+1][0], 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &request);
			for (int m=1; m<lowSize-1; m++){
				field[n+1][m]=f(n*tau, (rank*brickSize+m)*h)+field[n][m]-
					tau/2/h*(field[n][m+1]-field[n][m-1])+
					tau*tau/2/h/h*(field[n][m+1]-2*field[n][m]+
					field[n][m-1]);
//				cout<<field[n+1][m]<<endl;
			}
			int m=lowSize-1;
			field[n+1][m]=field[n][m]+f(n*tau, (m)*h)*tau-
					tau/h*(field[n][m]-field[n][m-1]);
		}
	}
	return field;
}

void save(vector<vector<double>> field, int N, int M){
	ofstream fout("output.txt");
	for (int m=0; m<M; m++)
		for (int n=0; n<N; n++)
			fout<<field[n][m]<<" ";
	fout.close();
}

void draw(int N, int M){
	ifstream fin("output.txt");
	std::vector<std::vector<double>> x, y, z;
	for (int m = 0; m < M;  m++) {
		std::vector<double> x_row, y_row, z_row;
		for (int n = 0; n < N; n++) {
			x_row.push_back(n*tau);
			y_row.push_back(m*h);
			double o;
			fin>>o;
			z_row.push_back(o);
		}
		x.push_back(x_row);
		y.push_back(y_row);
		z.push_back(z_row);
	}

    plt::plot_surface(x, y, z);
    plt::show();
}

int main(int argc, char* argv[]){
	int numprocs, rank, rc;
	if (rc = MPI_Init(&argc, &argv))
		abortAll(-1, "Ошибка запуска MPI");
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (argc<5)
		abortAll(-1, "not enough arguments");
	T = atof(argv[1]);
	tau = atof(argv[2]);
	X = atof(argv[3]);
	h = atof(argv[4]);
	N = T/tau;
	M = X/h;
	if (rank==0){
		vector<vector<double>> field = solve(0, 1);
		save(field, N, M);
		draw(N, M);
	}
	MPI_Barrier(MPI_COMM_WORLD);
/*	 solve(rank, numprocs);*/
/*	if (rank==0){
	}
*/	MPI_Finalize();
	return 0;
}
