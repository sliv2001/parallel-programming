mpi: build-mpi run-mpi

build-mpi:
	mpicxx ./main.cpp -O0 -ggdb -I./matplotlib-cpp -I/usr/include/python3.8 -lpython3.8 -o main

run-mpi:
	mpirun -np 2 --oversubscribe ./main 1 0.01 1 0.01

clean:
	rm main
