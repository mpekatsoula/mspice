#include <cusp/coo_matrix.h>
#include <cusp/print.h>
#include <cusp/array1d.h>
#include <cusp/krylov/cg.h>
#include <cusp/precond/aggregation/smoothed_aggregation.h>
#include <cusp/krylov/detail/bicg.inl>


extern "C" {
#include "hms/error_header.h"
#include "hms/csparse.h"
#include "hms/solution_functions.h"
}

clock_t start, end;

void start_timer() {

  start = clock();
}

double end_timer() {

  return ((double)clock() - start ) / CLOCKS_PER_SEC;
}
// 44943
// 127565
// 852536
// 954542
// 1618394
// 2506730

int it1;
cusp::coo_matrix<int,double,cusp::host_memory> cuMNA( 2506730,2506730, 7504532); //147315*2

extern "C"
void cuda_init () {

  it1 = 0;
  
}

extern "C"
void cuda_implementation () {
printf("%d\n",it1);
    cusp::array1d<double, cusp::device_memory> b (MNA_matrix_size );
    cusp::array1d<double, cusp::device_memory> x (MNA_matrix_size );

    int i;
    for ( i = 0; i < MNA_matrix_size; i ++ )
      b[i] = b_sparse_vector[i];
  
    cusp::csr_matrix<int,double,cusp::device_memory> A = cuMNA;
    cusp::default_monitor<double> monitor(b, 36520, 1e-6);
    cusp::precond::aggregation::smoothed_aggregation<int, double, cusp::device_memory> M(A);
    start_timer();
    cusp::krylov::bicg(A,A, x,b, monitor, M, M);
   // cusp::krylov::cg(A, x,b, monitor,M);
    printf("time %g\n",end_timer());
   // cusp::print(x);

}

extern "C"
void add_cuda_element ( int idx1, int idx2, double val ) {

 // if ( temp[idx1][idx2] != -1 ) {

   // cuMNA.values[temp[idx1][idx2]] += val;

 // } 
 // else {
//
 //   temp[idx1][idx2] = it1;
    int i = 0;
    for ( i = 0; i < it1; i ++ ) {
      if ( cuMNA.row_indices[i] == idx1 && cuMNA.row_indices[i] == idx2  ) {
   //     cuMNA.row_indices[it1] = idx1;
    //    cuMNA.column_indices[it1] = idx2;
        cuMNA.values[i] += val;
        return;
      }
    }
        cuMNA.row_indices[it1] = idx1;
        cuMNA.column_indices[it1] = idx2;
        cuMNA.values[it1] = val;    
        
    it1++;
    if ( !(it1 %50000) )
      printf("it %d\n", it1);
 // }
  
}

extern "C"
void add_cuda_b ( int idx, double val ) {
  
  //b[idx] = val;

}
