all: build run

build:
	mpicxx ./main.cpp -o main

run:
	mpirun -np 16 ./main 1000000000
	mpirun -np 1 ./main 1000000000

clean:
	rm main
