#include <iostream>
#include <pthread.h>
#include <unistd.h>

using namespace std;

#define NTHREADS 32

pthread_t ntid[NTHREADS];

int threads=16;
long long int N;
int finished=0;
pthread_mutex_t mutex;

struct args {
	int id;
	long double result;
};

void* func(void* arg){
	((struct args*)arg)->result=0.0;
	for (int i=((struct args*)arg)->id; i<N; i+=threads){
//		cout<<((struct args*)arg)->result;
		((struct args*)arg)->result+=1/((long double)i+1);
	}
	pthread_mutex_lock(&mutex);
	finished+=1;
	pthread_mutex_unlock(&mutex);
	return NULL;
}

int main(int argc, char* argv[]){

	int err;
	long double complete=0.0;
	threads=atoi(argv[1]);
	N = atol(argv[2]);
	if (pthread_mutex_init(&mutex, NULL)!=0)
		return 0;

	struct args nums[threads];
	for (int i=0; i<threads; i++){
		nums[i].id=i;
		err=pthread_create(&ntid[i], NULL, func, &nums[i]);
	}

	while (1)
		if (finished<threads)
			usleep(1000);
		else
			break;
	for (int i=0; i< threads; i++)
		complete+=nums[i].result;
	cout<<complete<<endl;

	return 0;
}
