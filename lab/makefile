mpi: build-mpi run-mpi

build-mpi:
	mpicxx ./main.cpp -I./matplotlib-cpp -I/usr/include/python2.7 -lpython2.7 -o main

run-mpi:
	mpirun -np 16 ./main

clean:
	rm main
