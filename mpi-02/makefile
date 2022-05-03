all: build run

build:
	mpicxx ./main.cpp -o main

run:
	mpirun -np 16 ./main 500000000000000
	mpirun -np 1 ./main 500000000000000

clean:
	rm main
