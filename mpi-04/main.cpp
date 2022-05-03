#include <iostream>
#include <pthread.h>
#include <unistd.h>

using namespace std;

#define NTHREADS 32

pthread_t ntid[NTHREADS];

int threads=16;

void* func(void* arg){
	pthread_t tid = pthread_self();
	cout<<"hello world from "<<*(int*)arg<<"-th of "<<threads<<endl;
	return NULL;
}

int main(int argc, char* argv[]){

	int err;
	threads=atoi(argv[1]);
	int nums[threads+1];
	for (int i=0; i<threads; i++){
		nums[i]=i+1;
		err=pthread_create(&ntid[i], NULL, func, &nums[i]);
	}
	sleep(1);
	return 0;
}
