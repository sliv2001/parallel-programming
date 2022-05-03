all: build run

build:
	mpicxx ./main.cpp -o main

run:
	mpirun -np 16 ./main
