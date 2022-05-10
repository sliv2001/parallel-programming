#include <iostream>
#include <pthread.h>
#include <unistd.h>
#include <cmath>

#ifndef	THREAD_MAX
#define	THREAD_MAX		128
#endif

#define	RECURSION_MAX		1024
#define BRICK_NUMBER		256
#define INIT_ITERATION		0
#define INIT_POINT_NUMBER	128
#define INIT_EPSILON		0.001

using namespace std;
int threads=16;
long double epsilon=INIT_EPSILON, ress=0;
int new_brick=0;

struct thread_handle {
	pthread_t tid;
	pthread_mutex_t mutex;
	long double a_local, b_local;
	int recursion_counter;
	int status; /*0=in progress; 1=finished; 2=ready for calculation; 3=done, result added*/
	long double cumulative=0;
};

void prinths(int i, struct thread_handle* hs){
//	cout<<"i="<<i;
//	cout<<"\t tid="<< hs->tid;
	cout<<"a="<< hs->a_local;
	cout<<"\t b="<< hs->b_local;
	cout<<"\t recursion="<< hs->recursion_counter;
	cout<<"\t status="<< hs->status;
	cout<<"\t cumulative="<<hs->cumulative<<endl;
}

long double function(long double x){
	return sin(2*M_PI*1000*x);
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
//	prinths(&arg);
	return intermed;
}

void* worker(void* arg){
	struct thread_handle* a=((struct thread_handle*)arg);
	while (1){
		if (a->status==2){
			pthread_mutex_lock(&a->mutex);
			a->status=0;
			pthread_mutex_unlock(&a->mutex);
			long double res=integrate_recursive(function,
				a, INIT_ITERATION, INIT_POINT_NUMBER);
			pthread_mutex_lock(&a->mutex);
			a->cumulative=res;
			a->status=1;
			pthread_mutex_unlock(&a->mutex);
		}
		else {
			usleep(10);
		}
	}
	return NULL;
}

void manage(int threads, long double epsilon, long double a, long double b){
	struct thread_handle hs[threads];
	long double result=0;
	long double brick=(b-a)/BRICK_NUMBER;
	for (int i=0; i<threads; i++){
		hs[i].a_local=a+new_brick*brick;
		hs[i].b_local=a+(new_brick+1)*brick>b?b:
			a+(new_brick+1)*brick;
		hs[i].recursion_counter=0;
		hs[i].status=2;
		pthread_mutex_init(&hs[i].mutex, NULL);
		pthread_create(&hs[i].tid, NULL, worker, &hs[i]);
		new_brick++;
	}
	while (new_brick<BRICK_NUMBER){
		for (int i=0; i<threads; i++){
			if (pthread_mutex_trylock(&hs[i].mutex)==0 &&
					hs[i].status==1){
				result+=hs[i].cumulative;
				prinths(i, &hs[i]);
				hs[i].a_local=a+new_brick*brick;
				hs[i].b_local=a+(new_brick+1)*brick>b?b:
						a+(new_brick+1)*brick;
				hs[i].recursion_counter=0;
				if (new_brick>=BRICK_NUMBER)
					hs[i].status=3;
				else{
					hs[i].status=2;
					new_brick++;
				}
			}
			pthread_mutex_unlock(&hs[i].mutex);
		}
		usleep(10);
	}

	int done;
	do{
		done=0;
		for (int i=0; i<threads; i++){
			if (hs[i].status==3)
				done+=1;
			if (hs[i].status==1){
				hs[i].status=3;
				result+=hs[i].cumulative;
				done +=1;
			}
		}
		usleep(10);
	} while (done != threads);
	cout<<result<<endl;
}

int main(int argc, char* argv[]){
	long double a, b;
	if (argc!=5)
		return -1;
	threads=atoi(argv[1]);
	epsilon = atof(argv[2]);
	epsilon = epsilon/sqrt(BRICK_NUMBER);
	a= atof(argv[3]);
	b=atof(argv[4]);
	if (threads<0 || threads>THREAD_MAX)
		return -2;
	if (epsilon<0)
		return -3;
	if (a>b)
		return -4;
	manage(threads, epsilon, a, b);
	return 0;
}
