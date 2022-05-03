mpi: build-mpi run-mpi

build-mpi:
	mpicxx ./main.cpp -o main

run-mpi:
	mpirun -np 16 ./main

clean:
	rm main
