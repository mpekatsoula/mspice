#ifndef _SOLUTION_FUNC_H
#define _SOLUTION_FUNC_H

#include <math.h>
#include <complex.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include "csparse.h"

#define EPS 1e-11
#define PI 3.14159265

/* Sparce matrixes */
cs *MNA_sparse;
cs *C_sparse;
double *b_sparse_vector;
cs *A;

/* Non sparce matrixes */
gsl_matrix *MNA_matrix;
gsl_vector *S_vector;
gsl_vector *Ukn_vector;
gsl_matrix *C;


int MNA_matrix_size;
int GNUPLOT;
double **Scan_Values; // temp buffer

void CG ( int DC_scan, char **unkown_vars );
void Bi_CG ( int DC_scan, char **unkown_vars );
void LU_solve(int DC_scan, char **unknown_vars);
void Cholesky_solve(int DC_scan, char **unknown_vars);
void BE( char **unknown_vars );
double exponential ( struct _element * node, double t );
double _sin ( struct _element *node, double t );
double pulse ( struct _element *node, double t, int k );
double pwl ( struct _element *node, double t );



#endif
