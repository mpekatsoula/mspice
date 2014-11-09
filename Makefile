all:
	nvcc error_func.c csparse.c solution_functions.c main.c parser.c hash.c construct_matrixes.c cleanup.c -g -O3 -lgsl -lgslcblas -lm -o mspice
