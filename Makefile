all:
	nvcc cuda.cu error_func.c csparse.c solution_functions.c main.c parser.c hash.c construct_matrixes.c cleanup.c -g -O3 -lgsl -I ../cusplibrary-0.4.0 -lgslcblas -arch=sm_20 -lm -o cudaspice
