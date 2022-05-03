#include <iostream>
#include <pthread.h>
#include <unistd.h>

using namespace std;

int n=16;
int ball=0, current=0;

void* func(void* arg){
	while (*(int*)arg!=current)
		usleep(10);
	ball++;
	current++;
	return NULL;
}

int main(int argc, char* argv[]){

	n=atoi(argv[1]);
	pthread_t threads[n];
	int nums[n];
	for (int i=0; i<n; i++){
		nums[i]=i;
		if (pthread_create(&threads[i], NULL, func, &nums[i]))
			return -2;
	}
	while (current<n)
		usleep(10);
	cout<<ball<<endl;
	return 0;
}
