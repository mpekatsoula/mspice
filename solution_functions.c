#include <stdio.h>
#include "hms/error_header.h"
#include "hms/parser.h"
#include "hms/hash.h"
#include "hms/solution_functions.h"
#include "hms/csparse.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <math.h>
#include <omp.h>
#include <time.h>
/* File handlers for data, and plot */
FILE *fp, *fp2,*fp3, *pipe;


/* Transient analysis functions */
void BE_sparse ( char **unknown_vars )
{

  double i;
  int j,k=0;
  gsl_matrix *_C;
  gsl_vector *temp;
  gsl_vector *x;
  gsl_vector *e;
  double h = NetOptions->TRAN_step;

}

void BE ( char **unknown_vars ) 
{

  double i;
  int j,k=0;
  gsl_matrix *_C;
  gsl_vector *temp;
  gsl_vector *x;
  gsl_vector *e;
  double h = NetOptions->TRAN_step;

  e = gsl_vector_calloc ( MNA_matrix_size );
  x = gsl_vector_alloc ( MNA_matrix_size );

  /* Temporary C matrix */
  _C = gsl_matrix_alloc ( MNA_matrix_size, MNA_matrix_size );
  temp = gsl_vector_alloc ( MNA_matrix_size );
  gsl_matrix_memcpy (_C, C);

  // C = C * 1/h
  gsl_matrix_scale (_C ,1.0/h);

  // _C = (G+1/h *C)
  // G = G + C
  gsl_matrix_add ( _C, MNA_matrix );

  // For the first iteration x[tk-1] = solution from DC
  gsl_vector_memcpy( x, Ukn_vector );  

  /* Overwrite mna, we don't need him */
  gsl_matrix_memcpy ( MNA_matrix, _C );

  /* Begin transient analysis */
  for ( i = 0; i < NetOptions->TRAN_fin_time; i = i + NetOptions->TRAN_step ) {
    
    // C * x[tk-1]
    gsl_blas_dgemv( CblasNoTrans, 1.0, C, x, 0.0, temp );

    // 1/h ( C * x[tk-1] )
    gsl_vector_scale( temp, 1.0/h );

    // Calculate b_now array
    for ( j = 0; j < MNA_matrix_size; j++ ) {

      if ( b_info[j]->EXP ) {

        gsl_vector_set ( e, j, exponential( b_info[j], i )  );
        
      }
      else if ( b_info[j]->SIN ) {

        gsl_vector_set ( e, j, _sin( b_info[j], i )  );

      }
      else if ( b_info[j]->PULSE ) {

       gsl_vector_set ( e, j, pulse( b_info[j], i, k )  );

      }
      else if ( b_info[j]->PWL ) {

        gsl_vector_set ( e, j, pwl( b_info[j], i )  );

      }
      else {
  
        gsl_vector_set ( e, j, b_info[j]->value );

      }

    }

    //  e + 1/h ( C * x[tk-1] )
    gsl_vector_add ( e, temp );
    
    /* Overwrite b/S_vector */
    gsl_vector_memcpy( S_vector ,e);

    // Solve again..
    NetOptions->DC = 0;
    if ( NetOptions->SPD && !NetOptions->ITER ) { // Solve with Cholesky

      Cholesky_solve(-1, unknown_vars);
      
      printf("\tSaved in Cholesky.txt \n");
    }
    else if ( NetOptions->ITER && !NetOptions->SPD) {

      Bi_CG(-1, unknown_vars);

      printf("\tSaved in Bi_CG_solution.txt \n");
    }
    else if ( NetOptions->ITER && NetOptions->SPD) {

      CG(-1, unknown_vars);
      printf("\tSaved in CG_solution.txt \n");
    }
    else { // Solve with LU

      LU_solve(-1, unknown_vars);
      printf("\tSaved in LU_solution.txt \n");
    } 

    gsl_vector_memcpy ( x, Ukn_vector );

    k++;

  }

  /* Free mem */
  gsl_matrix_free( _C );
  gsl_vector_free( x );

}

void TR ( char **unknown_vars ) 
{

  double i;
  int j,k=0;
  gsl_matrix *_C;
  gsl_matrix *_C1;
  gsl_vector *temp;
  gsl_vector *x;
  gsl_vector *e;
  gsl_vector *prev_e;
  double h = NetOptions->TRAN_step;

  e = gsl_vector_calloc ( MNA_matrix_size );
  x = gsl_vector_alloc ( MNA_matrix_size );
  prev_e = gsl_vector_calloc ( MNA_matrix_size );

  /* Temporary C matrix */
  _C = gsl_matrix_alloc ( MNA_matrix_size, MNA_matrix_size );
  _C1 = gsl_matrix_alloc ( MNA_matrix_size, MNA_matrix_size );
  temp = gsl_vector_alloc ( MNA_matrix_size );
  gsl_matrix_memcpy (_C, C);

  // C = C * 2/h
  gsl_matrix_scale (_C ,2.0/h);

  // _C1 = C * -2/h
  gsl_matrix_scale (_C1 ,-2.0/h);

  // _C = (G+2/h *C)
  gsl_matrix_add ( _C, MNA_matrix );

  // _C1 = (G+2/h *C)
  gsl_matrix_add ( _C1, MNA_matrix );

  // For the first iteration x[tk-1] = solution from DC
  gsl_vector_memcpy( x, Ukn_vector );  

  /* Overwrite mna, we don't need him */
  gsl_matrix_memcpy ( MNA_matrix, _C );

  gsl_vector_memcpy ( prev_e, S_vector );
  /* Begin transient analysis */
  for ( i = 0; i < NetOptions->TRAN_fin_time; i = i + NetOptions->TRAN_step ) {
    
    // _C1 * x[tk-1]
    gsl_blas_dgemv( CblasNoTrans, 1.0, _C1, x, 0.0, temp );

    // Calculate e array
    for ( j = 0; j < MNA_matrix_size; j++ ) {

      if ( b_info[j]->EXP ) {

        gsl_vector_set ( e, j, exponential( b_info[j], i )  );
        
      }
      else if ( b_info[j]->SIN ) {

        gsl_vector_set ( e, j, _sin( b_info[j], i )  );

      }
      else if ( b_info[j]->PULSE ) {

       gsl_vector_set ( e, j, pulse( b_info[j], i, k )  );

      }
      else if ( b_info[j]->PWL ) {

        gsl_vector_set ( e, j, pwl( b_info[j], i )  );

      }
      else {
  
        gsl_vector_set ( e, j, b_info[j]->value );

      }

    }

    // 

    gsl_vector_add( e, prev_e);
    gsl_vector_sub ( e, temp );
    gsl_vector_memcpy(prev_e, e);
    
    /* Overwrite b/S_vector */
    gsl_vector_memcpy( S_vector ,e);

    // Solve again..
    NetOptions->DC = 0;
    if ( NetOptions->SPD && !NetOptions->ITER ) { // Solve with Cholesky

      Cholesky_solve(-1, unknown_vars);
      
      printf("\tSaved in Cholesky.txt \n");
    }
    else if ( NetOptions->ITER && !NetOptions->SPD) {

      Bi_CG(-1, unknown_vars);

      printf("\tSaved in Bi_CG_solution.txt \n");
    }
    else if ( NetOptions->ITER && NetOptions->SPD) {

      CG(-1, unknown_vars);
      printf("\tSaved in CG_solution.txt \n");
    }
    else { // Solve with LU

      LU_solve(-1, unknown_vars);
      printf("\tSaved in LU_solution.txt \n");
    } 

    gsl_vector_memcpy ( x, Ukn_vector );

    k++;

  }

  /* Free mem */
  gsl_matrix_free( _C );
  gsl_vector_free( x );

}

void TR_sparse( char **unknown_vars )
{

}

double exponential ( struct _element *node, double t )
{
  double res;

  if ( (t >= 0) && (t < node->td1) ) {

    res = node->i1;

  }
  else if ( (t >= node->td1) && (t < node->td2) ) {
  
    res = node->i1 + ( node->i2 - node->i1 )*( 1.0 - exp( -( t - node->td1 )/node->tc1 ) );

  }
  else if ( (t >= node->td2) && (t < NetOptions->TRAN_fin_time) ) {
  
    res = node->i1 + ( node->i2 - node->i1 )*( exp( -( t - node->td2 )/node->tc2 - exp( -( t - node->td1 )/node->tc1 ) ) );

  }

  return res;

}

double _sin ( struct _element *node, double t )
{
  double res;

  if ( (t >= 0) && (t < node->td) ) {

    res = node->i1 + node->ia * sin( 2*PI*node->ph/360 );

  }
  else if ( (t >= node->td) && (t < NetOptions->TRAN_fin_time) ) {
  
    res = node->i1 + node->ia * sin( 2*PI*node->fr*( t - node->td ) + 2*PI*node->ph/360 )*exp( -(t - node->td)*node->df );

  }

  return res;

}

double pwl ( struct _element *node, double t )
{

  int i;

  double res;

  if ( t >= 0 && t <= node->pwl_t[0] ) {
  
    res = node->pwl_i[0];

  }
  else if ( t > node->pwl_t[node->n] && t <= NetOptions->TRAN_fin_time ) {
  
    res = node->pwl_i[node->n];

  }
  else {
    for ( i = 0; i < node->n; i++ ) {

      if ( t > node->pwl_t[i] &&  t > node->pwl_t[i + 1] ) {
        res = node->pwl_i[i];
        break;
      }        

    }
  }

  return res;

}

double pulse ( struct _element *node, double t, int k )
{
  double res;

  if ( (t >= 0) && (t < node->td) ) {

    res = node->i1;

  }
  else if ( (t >= node->td + k*node->per) && (t < node->td + node->tr + k*node->per) ) {
  
    res = node->i1 + t*(node->i2 - node->i1)/((node->td + node->tr + k*node->per) - (node->td + k*node->per));
    
    if (res > node->i2)
      res = node->i2;

  }
  else if ( (t >= node->td + node->tr + k*node->per) && (t < node->td + node->tr + node->pw + k*node->per) ) {
  
    res = node->i2;

  }
  else if ( (t >= node->td + node->tr + node->pw + k*node->per) && (t < node->td + node->tr + node->pw + node->tf + k*node->per) ) {
  
    res = node->i2 - t*(node->i2 - node->i1)/((node->td + node->tr + node->pw + node->tf + k*node->per) - (node->td + node->tr + node->pw + k*node->per));
    
    if (res < node->i1)
      res = node->i1;

  }
  else if ( (t >= node->td + node->tr + node->pw + node->tf + k*node->per) && (t < node->td + node->per + k*node->per) ) {
  
    res = node->i1;

  }

  return res;


}

/* Calculate the norm of a vector */
double norm ( const double *v )
{
  
  int i;
  double m_sum = 0.0;
  double max = fabs(v[0]);
  for ( i = 1; i < MNA_matrix_size; i++ )
    if ( fabs(v[i]) > max )
      max = fabs(v[i]);
    //m_sum  += v[i]*v[i];

  //return sqrt(m_sum);
  return max;

}

/* Plot data in gnuplot */
void gnuplot_and_print_data(int iterations, char **unknown_vars)
{

  int k,j;
  float i;

  /* Initial command for plot */
  for ( k = 0; k < NetOptions->num_of_dc_plots &&  !GNUPLOT ; k++ )
  {
    if ( !k && NetOptions->num_of_dc_plots == 1 )
      fprintf(pipe, "plot '-' title 'node %s' lt %d  with linespoints\n", NetOptions->PLOT_scan[k], k + 1);
    else if ( !k )
      fprintf(pipe, "plot '-' title 'node %s' lt %d  with linespoints, ", NetOptions->PLOT_scan[k], k + 1);
    else if ( k == NetOptions->num_of_dc_plots - 1 )
      fprintf(pipe, "'-' title 'node %s' lt %d  with linespoints\n", NetOptions->PLOT_scan[k], k + 1);
    else 
      fprintf(pipe, "'-' title 'node %s' lt %d  with linespoints,", NetOptions->PLOT_scan[k], k + 1 );
  }

  for ( k = 0; k < NetOptions->num_of_dc_plots; k++ )
  {

    int plot_node = ht_get( nodes_hashtable, NetOptions->PLOT_scan[k]) - 1;
    fprintf ( fp2, "--------Value of node %s ( now %s ), for DC Scanning--------\n", NetOptions->PLOT_scan[k],  unknown_vars[plot_node]);
    for ( j = 0 ,i = NetOptions->DC_start; j < iterations; j++, i +=NetOptions->DC_step  )
    {
      fprintf ( fp2, "value: %g \t result: %g\n",i , Scan_Values[k][j]);
      /* Add the plot data */
      if ( !GNUPLOT )
          fprintf(pipe, "%g %g\n", i, Scan_Values[k][j]);
    }
    if ( !GNUPLOT )
      fprintf(pipe, "e\n"); 
  }

  

}

/* Allocate memory for a helper matrix */
void solution_mem_alloc(int iterations)
{

  int i;

  /* Allocate memory for temporary buffer */
  if ( NetOptions->DC && NetOptions->PLOT )
  {
    
    Scan_Values = malloc ( sizeof ( double ) *  NetOptions->num_of_dc_plots );
    for ( i = 0; i < NetOptions->num_of_dc_plots; i++ )
      Scan_Values[i] = calloc ( sizeof(double), sizeof( double ) * (iterations + 1));
  }

}

/* Free allocated memory */
void solution_mem_free()
{

  int i;

  /* Free memory for temporary buffer */
  if ( NetOptions->DC && NetOptions->PLOT )
  {

    for ( i = 0; i < NetOptions->num_of_dc_plots; i++ )
      free ( Scan_Values[i] );
    free ( Scan_Values );
  }

  /* Force write results and close files */
  if ( !GNUPLOT )
  {
    fflush(pipe);
    fclose(pipe); 
  }
  fflush(fp);
  fclose(fp); 
  fflush(fp2);
  fclose(fp2); 

}

/*******************************************/
/************ SPARSE functions *************/
/*******************************************/

/* Find diagonal*/
double * diag(cs *A) {

  int p,j;
  double *d;
  int cnt;
  
  cnt=0;
  d = malloc( MNA_matrix_size * sizeof(double) );
  
  for ( j = 0; j < MNA_matrix_size; j++ )
    d[j] = 1.0;
  
  for ( j = 0; j < MNA_matrix_size; j++ ){
    for( p = A->p[j]; p < A->p[j+1]; p++ ){
      if( A->i[p] == cnt ){
        d[j] = 1.0/A->x[p];
      }
    }
    cnt++;
  }
  
  return(d);
}


/* LU sparse */
void LU_Sparse(char **unknown_vars)
{

  /* File handlers for data, and plot */
  fp = fopen( "LU_solution.txt", "w" );

  /* Temp buffers */
  css *S;
  csn *N;
  double *x;
  S = cs_malloc(MNA_matrix_size,sizeof(css *));
  N = cs_malloc(MNA_matrix_size,sizeof(csn *));
  x = cs_malloc(MNA_matrix_size,sizeof(double));
  
  S = cs_sqr (2, A, 1);
  N = cs_lu ( A, S, 1);
  cs_ipvec ( N->pinv, b_sparse_vector, x, MNA_matrix_size );

  /* Solve LU */
  cs_lsolve ( N->L, x );
  cs_usolve ( N->U, x );
  cs_ipvec ( S->q, x, b_sparse_vector, MNA_matrix_size );



  int i;

  /* Save solution to file */
  for ( i = 0; i < MNA_matrix_size; i++ ) 
    fprintf (fp, "%s: %g\n", unknown_vars[i], b_sparse_vector[i]);


  /* Free mem */
  cs_spfree ( A );
  cs_nfree ( N );
  cs_sfree ( S );
  

}


/* Cholesky sparse */
void Cholesky_Sparse()
{

  /* Temp buffers */
  css *S;
  csn *N;
  double *x = malloc( sizeof(double) * MNA_matrix_size );

  
  S = cs_schol ( 1, A );
  N = cs_chol ( A, S );
  
  /* Solve */ 
  cs_ipvec ( S->pinv, b_sparse_vector, x, MNA_matrix_size );
  cs_lsolve ( N->L, x);
  cs_ltsolve ( N->L, x);
  cs_pvec ( S->pinv, x, b_sparse_vector, MNA_matrix_size );

  int i;
  for ( i = 0; i < MNA_matrix_size; i++ )
    printf("\t%g\n", b_sparse_vector[i] );

  /* Free mem */
  cs_spfree ( A );
  cs_nfree ( N );
  cs_sfree ( S );

}

clock_t start, end;

void start_timer() {

  start = clock();
}

double end_timer() {

  return ((double)clock() - start ) / CLOCKS_PER_SEC;
}

void CG_Sparse( char **unknown_vars)
{

  /* Local vars */
  int i, iter,j,p1;    
  double norm_s, norm_r, rho, beta, rho1, alpha;
  fp = fopen( "CG_solution.txt", "w" );

  double *x = malloc ( sizeof( double ) * MNA_matrix_size );
  double *y = malloc ( sizeof( double ) * MNA_matrix_size );
  double *r = malloc ( sizeof( double ) * MNA_matrix_size );
  double *z = malloc ( sizeof( double ) * MNA_matrix_size );
  double *p = malloc ( sizeof( double ) * MNA_matrix_size );
  double *q = malloc ( sizeof( double ) * MNA_matrix_size );
  double *M_inv;

  // initial guess, x = 0
  #pragma omp parallel for
  for ( i = 0; i < MNA_matrix_size; i++ )
  {
    x[i] = 0;
    y[i] = 0;
    q[i] = 0;
  }
start_timer();

  // y = A*x + y
  cs_gaxpy ( A, x, y );

printf("cs_gaxpy %g\n",end_timer());

  // Get inv diagonal of A
  M_inv = diag(A);

  // r = b - A*x
  #pragma omp parallel for
  for ( i = 0; i < MNA_matrix_size; i++ )
    r[i] = b_sparse_vector[i] - y[i];


  iter = 0;

  // calculate norms, ||r|| and ||b||
  norm_r = norm( r );
  norm_s = norm ( b_sparse_vector );
  if ( !norm_s )
    norm_s = 1;
start_timer();
  /* Main loop */
  while (norm_r/norm_s > atof("1e-6") && iter < MNA_matrix_size )
  {

    iter++;

    // solve M*z = r
    // Multiply M^-1 with r. z[i] = (1/(a[i][i]))*r[i]
    // rho = (r^T)*z
    rho = 0;
    //#pragma omp parallel for reduction(+:rho)
    for ( i = 0; i < MNA_matrix_size; i++ )
    {
      z[i] = M_inv[i]*r[i];
      rho +=  r[i]*z[i];
    }

    if ( iter == 1 ) { 
      for ( i = 0; i < MNA_matrix_size; i++ )
        p[i] = z[i];  
    }
    else {
      beta = rho/rho1;
      // p = z + beta*p
      //#pragma omp parallel for
      for ( i = 0; i < MNA_matrix_size; i++ )
        p[i] = p[i]*beta + z[i];
    }

    rho1 = rho;

    // q = A*p
    //#pragma omp parallel for
    for ( i = 0; i < MNA_matrix_size; i++ )
      q[i] = 0;
    cs_gaxpy ( A, p, q );

    // (p^T)*q
    alpha = 0;
    //#pragma omp parallel for reduction(+:alpha)
    for ( i = 0; i < MNA_matrix_size; i++ )
      alpha +=  p[i]*q[i];

    // alpha = rho/((p^T)*q)
    alpha = rho/alpha;

    // x = x + alpha*p
    //#pragma omp parallel for 
    for ( i = 0; i < MNA_matrix_size; i++ )  
    {  
      x[i] += alpha*p[i];
      r[i] -= alpha*q[i];
    }

    norm_r = norm( r );

    if ( !(iter % 300) ) {
      printf("Iterations done so far: %d\n", iter);
      printf("norm_r: %g, res: %g\n\n\n",norm_r, x[0]);
    }
  }
printf("total_time %gsec.\n",end_timer());
  /* Save solution to file */
  printf("Iterations done: %d\n", iter);
  for ( i = 0; i < MNA_matrix_size; i++ )
    fprintf (fp, "%s: %g\n", unknown_vars[i], x[i]);

  /* Clean Up memory */
  free(x);
  free(r);
  free(p);
  free(q);
  free(z);
  free(y);
  
}


void Bi_CG_Sparse( char **unknown_vars )
{


  double rfo, rho1, rho, k;
  int iter, i, temp1 = 0;
  double norm_s, norm_r, alpha, beta, omega;
  fp = fopen( "Bi_CG_solution.txt", "w" );

  /* Temp buffers */
  double *x = malloc ( sizeof( double ) * MNA_matrix_size );
  double *y = malloc ( sizeof( double ) * MNA_matrix_size );
  double *r = malloc ( sizeof( double ) * MNA_matrix_size );
  double *r_ = malloc ( sizeof( double ) * MNA_matrix_size );
  double *z = malloc ( sizeof( double ) * MNA_matrix_size );
  double *z_ = malloc ( sizeof( double ) * MNA_matrix_size );
  double *p = malloc ( sizeof( double ) * MNA_matrix_size );
  double *p_ = malloc ( sizeof( double ) * MNA_matrix_size );
  double *q = malloc ( sizeof( double ) * MNA_matrix_size );
  double *q_ = malloc ( sizeof( double ) * MNA_matrix_size );
  double *M_inv;

  // initial guess, x = 0
  //#pragma omp parallel for
  for ( i = 0; i < MNA_matrix_size; i++ )
  {
    x[i] = 0;
    y[i] = 0;
    q[i] = 0;
    q_[i] = 0;
  }

  // A*x
  cs_gaxpy ( A, x, y );

  // Get inv diagonal of A
  M_inv = diag(A);

  // r = b - A*x
  // r~ = r
  //#pragma omp parallel for
  for ( i = 0; i < MNA_matrix_size; i++ )
  {
    r[i] = b_sparse_vector[i] - y[i];
    r_[i] = r[i];
  }

  iter = 0;

  // calculate norms, ||r|| and ||b||
  norm_r = norm ( r );
  norm_s = norm ( b_sparse_vector);
  if ( !norm_s )
    norm_s = 1;

  /* main loop */
  while (norm_r/norm_s > NetOptions->ITOL && iter < MNA_matrix_size )
  {

    iter++;

    // solve M*z = r
    // Multiply M^-1 with r. z[i] = (1/(a[i][i]))*r[i]
    // solve (M^T)*z~ = r~
    // z~ = r~*M^-T
    // rho = (z^T)*r_
    rho = 0;
    //#pragma omp parallel for reduction(+:rho)
    for ( i = 0; i < MNA_matrix_size; i++ )
    {
      z[i] = r[i]*M_inv[i];
      z_[i] = r_[i]*M_inv[i];
      rho +=  z[i]*r_[i];
    }

    if ( fabs(rho) < EPS )
      break;

    if ( iter == 1)
    {
      //#pragma omp parallel for
      for ( i = 0; i < MNA_matrix_size; i++ )
      {
        p[i] = z[i];
        p_[i] = z_[i];
      }
    }
    else
    {

      beta = rho/rho1;
      // p = z + beta*p
      // p~ = z~ + beta*p~
      //#pragma omp parallel for
      for ( i = 0; i < MNA_matrix_size; i++ )
      {
        p[i] = p[i]*beta + z[i];
        p_[i] = p_[i]*beta + z_[i];
      }
  
    }
   
    rho1 = rho;

    // q = A*p
    //#pragma omp parallel for
    for ( i = 0; i < MNA_matrix_size; i++ )
      q[i] = 0;
    cs_gaxpy ( A, p, q );
     
    // q~ = A^T*p~
    cs_gaxpy_trans ( A, p_, q_ );

    // omega = (p_^T)*q
    omega = 0;
    //#pragma omp parallel for reduction(+:omega)
    for ( i = 0; i < MNA_matrix_size; i++ )
      omega +=  p_[i]*q[i];

    if ( fabs(omega) < EPS )
      break;

    alpha = rho/omega;
     
    // x = x + alpha*p
    // r = r - alpha*q
    // r~ = r~ - alpha*q~
    //#pragma omp parallel for
    for ( i = 0; i < MNA_matrix_size; i++ )
    {
      x[i] += alpha*p[i];
      r[i] -= alpha*q[i];
      r_[i] -= alpha*q_[i];
    }

    norm_r = norm ( r );
    if ( !(iter % 300) ) {
      printf("Iterations done so far: %d\n", iter);
      printf("norm_r: %g, res: %g\n\n\n",norm_r, x[0]);
    }
  }

  /* Save solution to file */
  for ( i = 0; i < MNA_matrix_size; i++ )
    fprintf (fp, "%s: %g\n", unknown_vars[i], x[i]);

  /* Clean Up memory */
  free(x);
  free(r);
  free(p);
  free(q);
  free(r_);
  free(p_);
  free(q_);
  free(z);
  free(y);

}


/*******************************************/
/********** Non SPARSE functions ***********/
/*******************************************/


/* This function solves the MNA system 
 * using CG. Prints the results to a file. */
void CG ( int DC_scan, char **unknown_vars )
{

  /* File handlers for data, and plot */
  if ( !GNUPLOT )
    pipe = popen("gnuplot -persist", "w");

  fp = fopen( "CG_solution.txt", "a+" );
  fp2 = fopen( "Plot Values.txt", "w" );

  /* Set general options for plot */
  if ( !GNUPLOT )
  fprintf(pipe, "set grid\n \
                 set title 'DC-sweep analysis CG'\n \
                 set xlabel 'Source value'\n \
                 set ylabel 'Element value'\n"  );


  fprintf (fp, "---------------Complete Circuit Solution--------------\n");
  int i,j;
  float _iter;
  int iterations;

  if ( !GNUPLOT)
  {
    iterations = (NetOptions->DC_end - NetOptions->DC_start)/NetOptions->DC_step;
    // Allocate memory for temporary buffers
    solution_mem_alloc(iterations);
  }

  /* Vectors needed for computation of CG */
  gsl_vector *x, *r, *z, *y, *q, *p, *temp, *M_inv;
  x =  gsl_vector_alloc(MNA_matrix_size);
  r =  gsl_vector_alloc(MNA_matrix_size);
  z =  gsl_vector_alloc(MNA_matrix_size);
  y =  gsl_vector_alloc(MNA_matrix_size);
  p =  gsl_vector_alloc(MNA_matrix_size);
  q =  gsl_vector_alloc(MNA_matrix_size);
  temp =  gsl_vector_alloc(MNA_matrix_size);
  M_inv =  gsl_vector_alloc(MNA_matrix_size);

  double rfo, rho1, rho;
  int iter, temp1 = -2;
  float k;
  double norm_s, norm_r, alpha, beta;

  // Initial guess, x=0, and find 1/a[i][i]
  for ( i = 0; i < MNA_matrix_size; i++ )
  {
    gsl_vector_set( x, i, 0);
    double num = gsl_matrix_get(MNA_matrix, i, i);
    if ( !num )
      num = 1;

    gsl_vector_set( M_inv, i, 1.0/num );
  }


  /* DC_Scan loop */
  for ( k = NetOptions->DC_start; k <= NetOptions->DC_end; k += NetOptions->DC_step)
  {
    temp1++;

    // A*x 
    gsl_blas_dgemv( CblasNoTrans, 1.0, MNA_matrix, x, 0.0, y );

    // r = b - A*x
    gsl_vector_memcpy( r , S_vector);
    gsl_vector_sub( r, y );


    iter = 0;

    // calculate norms, ||r|| and ||b||
    norm_r = gsl_blas_dnrm2 ( r );
    norm_s = gsl_blas_dnrm2 ( S_vector);
    if ( !norm_s )
      norm_s = 1;

    /* Main loop */
    while (norm_r/norm_s > NetOptions->ITOL && iter < MNA_matrix_size )
    {
      iter++;

      // solve M*z = r
      // Multiply M^-1 with r. z[i] = (1/(a[i][i]))*r[i]
      gsl_vector_memcpy (z, M_inv);
      gsl_vector_mul (z, r);

      // rho = (r^T)*z
      rho = 0;
      for ( i = 0; i < MNA_matrix_size; i++ )
        rho +=  gsl_vector_get(r,i)*gsl_vector_get(z,i);


      if ( iter == 1 ) // p = z
        gsl_vector_memcpy (p, z);
      else
      {
        beta = rho/rho1;
        // p = z + beta*p
        gsl_vector_scale ( p , beta );
        gsl_vector_add ( p , z );
      }
      rho1 = rho;

      // q = A*p
      gsl_blas_dgemv( CblasNoTrans, 1.0, MNA_matrix, p, 0.0, q );
     
      // (p^T)*q
      alpha = 0;
      for ( i = 0; i < MNA_matrix_size; i++ )
        alpha +=  gsl_vector_get(p,i)*gsl_vector_get(q,i);

      // alpha = rho/((p^T)*q)
      alpha = rho/alpha;

      // x = x + alpha*p
      gsl_vector_memcpy ( temp, p ); 
      gsl_vector_scale ( temp , alpha );
      gsl_vector_add ( x, temp );

      // r = r - alpha*q
      gsl_vector_memcpy (temp, q); 
      gsl_vector_scale ( temp , alpha );
      gsl_vector_sub ( r, temp );

    }

    /* Store solutions */
    for ( j = 0 ; j < NetOptions->num_of_dc_plots && temp1 > -1 ; j++ )
    {
      int plot_node = ht_get( nodes_hashtable, NetOptions->PLOT_scan[j]) - 1;

      if ( NetOptions->DC && NetOptions->PLOT ) 
        Scan_Values[j][temp1] = gsl_vector_get(x, plot_node);
    }

    /* Print Original solution */
    if ( temp1 == -1 ) {

      /* Print solution data */
      for ( j = 0; j < MNA_matrix_size; j++ )
        fprintf (fp, "%s: %g\n", unknown_vars[j], gsl_vector_get(x, j));

    }

    /* Change DC values */
    if ( NetOptions->DC )
      gsl_vector_set( S_vector, DC_scan, k);

  }


  if ( !GNUPLOT && DC_scan != -1 )
  {
    // plot and print data
    gnuplot_and_print_data( iterations, unknown_vars );
  }


  /* Free mem */
  gsl_vector_free(x);
  gsl_vector_free(q);
  gsl_vector_free(p);
  gsl_vector_free(r);
  gsl_vector_free(M_inv);
  gsl_vector_free(y);
  gsl_vector_free(z);
  gsl_vector_free(temp);

  solution_mem_free();

}

/* This function solves the MNA system 
 * using Bi-CG Prints the results to a file. */
void Bi_CG ( int DC_scan, char **unknown_vars )
{

  /* File handlers for data, and plot */
  if ( !GNUPLOT )
    pipe = popen("gnuplot -persist", "w");

  fp = fopen( "Bi_CG_solution.txt", "a+" );
  fp2 = fopen( "Plot Values.txt", "w" );

  /* Set general options for plot */
  if ( !GNUPLOT )
  fprintf(pipe, "set grid\n \
                 set title 'DC-sweep analysis Bi CG'\n \
                 set xlabel 'Source value'\n \
                 set ylabel 'Element value'\n"  );


  fprintf (fp, "---------------Complete Circuit Solution--------------\n");
  int i,j,iterations;

  if ( !GNUPLOT)
  {
    iterations = (NetOptions->DC_end - NetOptions->DC_start)/NetOptions->DC_step;
    // Allocate memory for temporary buffers
    solution_mem_alloc(iterations);
  }

  /* Vectors needed for computation of Bi-CG */
  gsl_vector *x, *r, *r_, *z, *z_, *y, *q, *q_, *p, *p_, *temp, *M_inv;
  x =  gsl_vector_alloc(MNA_matrix_size);
  r =  gsl_vector_alloc(MNA_matrix_size);
  r_ =  gsl_vector_alloc(MNA_matrix_size);
  z =  gsl_vector_alloc(MNA_matrix_size);
  z_ =  gsl_vector_alloc(MNA_matrix_size);
  y =  gsl_vector_alloc(MNA_matrix_size);
  p =  gsl_vector_alloc(MNA_matrix_size);
  p_ =  gsl_vector_alloc(MNA_matrix_size);
  q =  gsl_vector_alloc(MNA_matrix_size);
  q_ =  gsl_vector_alloc(MNA_matrix_size);
  temp =  gsl_vector_alloc(MNA_matrix_size);
  M_inv =  gsl_vector_alloc(MNA_matrix_size);

  double rfo, rho1, rho, k;
  int iter, temp1 = -2;
  double norm_s, norm_r, alpha, beta, omega;


  // Initial guess, x=0, and find 1/a[i][i]
  for ( i = 0; i < MNA_matrix_size; i++ )
  {
    gsl_vector_set( x, i, 0);
    double num = gsl_matrix_get(MNA_matrix, i, i);
    if ( !num )
      num = 1.0;

    gsl_vector_set( M_inv, i, 1.0/num );
  }

  /* DC_Scan loop */
  for ( k = NetOptions->DC_start; k <= NetOptions->DC_end; k += NetOptions->DC_step)
  {
    temp1++;

    // A*x 
    gsl_blas_dgemv( CblasNoTrans, 1.0, MNA_matrix, x, 0.0, y );

    // r = b - A*x
    gsl_vector_memcpy( r , S_vector);
    gsl_vector_sub( r, y );

    // r~ = r
    gsl_vector_memcpy( r_, r );

    iter = 0;

    // calculate norms, ||r|| and ||b||
    norm_r = gsl_blas_dnrm2 ( r );
    norm_s = gsl_blas_dnrm2 ( S_vector);
    if ( !norm_s )
      norm_s = 1;

    /* main loop */
    while (norm_r/norm_s > NetOptions->ITOL && iter < MNA_matrix_size )
    {
      iter++;

      // solve M*z = r
      // Multiply M^-1 with r. z[i] = (1/(a[i][i]))*r[i]
      gsl_vector_memcpy (z, M_inv);
      gsl_vector_mul (z, r);

      // solve (M^T)*z~ = r~
      // z~ = r~*M^-T
      for ( i = 0; i < MNA_matrix_size; i++ )
        gsl_vector_set(z_, i, gsl_vector_get(r_ , i)*gsl_vector_get(M_inv, i));

      // rho = (z^T)*r_
      rho = 0;
      for ( i = 0; i < MNA_matrix_size; i++ )
        rho +=  gsl_vector_get(z,i)*gsl_vector_get(r_,i);

      if ( fabs(rho) < EPS )
        break;

      if ( iter == 1)
      {
        gsl_vector_memcpy (p, z);
        gsl_vector_memcpy (p_, z_);
      }
      else
      {
        beta = rho/rho1;
        // p = z + beta*p
        gsl_vector_scale ( p , beta );
        gsl_vector_add ( p , z );
    
        // p~ = z~ + beta*p~
        gsl_vector_scale ( p_ , beta );
        gsl_vector_add ( p_ , z_ );
      }
      rho1 = rho;

      // q = A*p
      gsl_blas_dgemv( CblasNoTrans, 1.0, MNA_matrix, p, 0.0, q );
     
      // q~ = A^T*p~
      gsl_blas_dgemv( CblasTrans, 1.0, MNA_matrix, p_, 0.0, q_ );

      // omega = (p_^T)*q
      omega = 0;
      for ( i = 0; i < MNA_matrix_size; i++ )
        omega +=  gsl_vector_get(p_,i)*gsl_vector_get(q,i);

      if ( fabs(omega) < EPS )
        break;

      alpha = rho/omega;
     
      // x = x + alpha*p
      gsl_vector_memcpy ( temp, p ); 
      gsl_vector_scale ( temp , alpha );
      gsl_vector_add ( x, temp );
	
      // r = r - alpha*q
      gsl_vector_memcpy (temp, q); 
      gsl_vector_scale ( temp , alpha );
      gsl_vector_sub ( r, temp );

      // r~ = r~ - alpha*q~
      gsl_vector_memcpy (temp, q_); 
      gsl_vector_scale ( temp , alpha );
      gsl_vector_sub ( r_, temp );
    }

    /* Store solutions */
    for ( j = 0 ; j < NetOptions->num_of_dc_plots && !GNUPLOT && temp1 > -1; j++ )
    {
      int plot_node = ht_get( nodes_hashtable, NetOptions->PLOT_scan[j]) - 1;

      if ( NetOptions->DC && NetOptions->PLOT ) 
        Scan_Values[j][temp1] = gsl_vector_get(x, plot_node);
    }
  
    /* Print solution data */
    if ( temp1 == -1 ) {
      
      for ( j = 0; j < MNA_matrix_size; j++ )
        fprintf (fp, "%s: %g\n", unknown_vars[j], gsl_vector_get(x, j));
    }

    /* Change values */
    if ( NetOptions->DC )
      gsl_vector_set( S_vector, DC_scan, k);
  }



  if ( !GNUPLOT && DC_scan != -1  )
  {
    // plot and print data
    gnuplot_and_print_data( iterations, unknown_vars );
  }

  /* Free mem */
  gsl_vector_free(x);
  gsl_vector_free(q);
  gsl_vector_free(p);
  gsl_vector_free(r);
  gsl_vector_free(q_);
  gsl_vector_free(p_);
  gsl_vector_free(r_);
  gsl_vector_free(M_inv);
  gsl_vector_free(y);
  gsl_vector_free(z);
  gsl_vector_free(z_);
  gsl_vector_free(temp);

  solution_mem_free();

}


/* This function solves the MNA system using Cholesky
 * Decomposition. Prints the results to a file.     */
void
Cholesky_solve(int DC_scan, char **unknown_vars)
{
  /* File handlers for data, and plot */
  if ( !GNUPLOT )
    pipe = popen("gnuplot -persist", "w");

  fp = fopen( "Cholesky_solution.txt", "a+" );
  fp2 = fopen( "Plot Values.txt", "w" );

  /* Set general options for plot */
  if ( !GNUPLOT )
  fprintf(pipe, "set grid\n \
                 set title 'DC-sweep analysis Cholesky'\n \
                 set xlabel 'Source value'\n \
                 set ylabel 'Element value'\n"  );

  fprintf (fp, "---------------Complete Circuit Solution--------------\n");
  gsl_linalg_cholesky_decomp(MNA_matrix);
  gsl_linalg_cholesky_solve (MNA_matrix, S_vector, Ukn_vector);

  int k,j, plot_node;

  int iterations = (NetOptions->DC_end - NetOptions->DC_start)/NetOptions->DC_step;

  for ( k = 0; k < MNA_matrix_size; k++ )
    fprintf (fp, "%s: %g\n", unknown_vars[k], gsl_vector_get(Ukn_vector, k));

  if ( !GNUPLOT )
  {
    // Allocate memory for temporary buffers
    solution_mem_alloc(iterations);
  }
  /* Compute DC scanning */
  float i;
  k = -1;
  for ( i = NetOptions->DC_start; i <= NetOptions->DC_end && NetOptions->DC; i += NetOptions->DC_step)
  {
    k++;
    gsl_vector_set( S_vector, DC_scan, i);
    gsl_linalg_cholesky_solve (MNA_matrix, S_vector, Ukn_vector);

    for ( j = 0 ; j < NetOptions->num_of_dc_plots; j++ )
    {
      plot_node = ht_get( nodes_hashtable, NetOptions->PLOT_scan[j]) - 1;

      if ( NetOptions->DC && NetOptions->PLOT ) 
        Scan_Values[j][k] = gsl_vector_get(Ukn_vector, plot_node);

    }
  }

  if ( !GNUPLOT && DC_scan != -1 )
  {
    // plot and print data
    gnuplot_and_print_data( iterations, unknown_vars );

    // Allocate memory for temporary buffers
    solution_mem_free();
  }

}

/* This function solves the MNA system using LU
 * Decomposition. Prints the results to a file. */
void LU_solve(int DC_scan, char **unknown_vars)
{
  /* File handlers for data, and plot */
  fp = fopen( "LU_solution.txt", "a+" );
  fp2 = fopen( "Plot Values.txt", "w" );
  if ( !GNUPLOT )
    pipe = popen("gnuplot -persist", "w");

  /* Set general options for plot */
  if ( !GNUPLOT )
    fprintf(pipe, "set grid\n \
                   set title 'DC-sweep analysis LU'\n \
                   set xlabel 'Source value'\n \
                   set ylabel 'Element value'\n"  );
  int s;
  gsl_permutation * p = gsl_permutation_alloc (MNA_matrix_size);

  /* Solve it once */
  fprintf (fp, "---------------Complete Circuit Solution--------------\n");
  gsl_linalg_LU_decomp (MNA_matrix, p, &s);
  gsl_linalg_LU_solve (MNA_matrix, p, S_vector, Ukn_vector);

  int k,j, plot_node;

  int iterations = (NetOptions->DC_end - NetOptions->DC_start)/NetOptions->DC_step;

  for ( k = 0; k < MNA_matrix_size; k++ )
    fprintf (fp, "%s: %g\n", unknown_vars[k], gsl_vector_get(Ukn_vector, k));

  if ( !GNUPLOT )
  {
    // Allocate memory for temporary buffers
    solution_mem_alloc(iterations);
  }

  /* Compute DC scanning */
  double i;
  k = -1;
  for ( i = NetOptions->DC_start; i <= NetOptions->DC_end && NetOptions->DC; i += NetOptions->DC_step)
  {
    k++;
    gsl_vector_set( S_vector, DC_scan, i);
    gsl_linalg_LU_solve (MNA_matrix, p, S_vector, Ukn_vector);

    for ( j = 0 ; j < NetOptions->num_of_dc_plots; j++ )
    {
      plot_node = ht_get( nodes_hashtable, NetOptions->PLOT_scan[j]) - 1;

      if ( NetOptions->DC && NetOptions->PLOT ) 
        Scan_Values[j][k] = gsl_vector_get(Ukn_vector, plot_node);

    }
  }

  if ( !GNUPLOT && DC_scan != -1  )
  {
    // plot and print data
    gnuplot_and_print_data( iterations, unknown_vars );

    // Allocate memory for temporary buffers
    solution_mem_free();
  }

  gsl_permutation_free(p);
}
