#include <iostream>
#include <pthread.h>
#include <unistd.h>
#include <cmath>

#ifndef	THREAD_MAX
#define	THREAD_MAX		128
#endif

#define	RECURSION_MAX		1024
#define INIT_ITERATION		0
#define INIT_POINT_NUMBER	128
#define INIT_EPSILON		0.001
using namespace std;
int threads=16;
long double epsilon=INIT_EPSILON;

struct thread_handle {
	pthread_t tid;
	long double a_local, b_local;
	int recursion_counter;
	long double cumulative=0;
};

long double function(long double x){
	return sin(x);
}

long double integrate(long double (*f)(long double x),
		long double a,
		long double b,
		long long int n)
{
	long double ip=a, cumulative=0, step=(b-a)/((long double)n);
	for (long double i=a; i<=b; i+=step){
		cumulative+=(f(i)+f(ip))*step/2;
		ip=i;
	}
	return cumulative;
}

long double integrate_recursive(long double (*f)(long double x),
		struct thread_handle* arg,
		long double previous, long long int n)
{
	arg->recursion_counter++;
	if (arg->recursion_counter>RECURSION_MAX){
		cerr<<"Max recursion reached."<<endl;
		return previous;
	}
	long double intermed=integrate(f,
		arg->a_local, arg->b_local, 2*n);
	if (abs(intermed - previous)>epsilon)
		return integrate_recursive(f, arg, intermed, 2*n);
	return intermed;
}

void* worker(void* arg){
	((struct thread_handle*)arg)->cumulative+=
			integrate_recursive(function,
			(struct thread_handle*)arg,
			INIT_ITERATION, INIT_POINT_NUMBER);
	return NULL;
}

int main(int argc, char* argv[]){
	if (argc!=3)
		return -1;
	threads=atoi(argv[1]);
	epsilon = atof(argv[2]);
	if (threads<0 || threads>THREAD_MAX)
		return -2;
	struct thread_handle a={
		.tid=0,
		.a_local=0.0,
		.b_local=3.0,
		.recursion_counter=0
	};
	worker((void*)&a);
	cout<<a.cumulative<<endl;
	return 0;
}
