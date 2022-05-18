mpi: build-mpi run-mpi

build-mpi:
	mpicxx ./main.cpp -O0 -ggdb -I./matplotlib-cpp -I/usr/include/python3.8 -lpython3.8 -o main

run-mpi:
	mpirun -np 16 ./main 0.001 0.001 10 10

clean:
	rm main
