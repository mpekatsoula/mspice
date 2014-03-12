#include <stdio.h>
#include "hms/error_header.h"
#include "hms/parser.h"
#include "hms/hash.h"
#include "hms/csparse.h"
#include "hms/solution_functions.h"
#include <gsl/gsl_matrix.h>
#include <unistd.h>


void
check_gnuplot()
{

  int gnuplot = access("/usr/bin/gnuplot", X_OK);
  gnuplot |= access("/usr/lib/gnuplot/gnuplot_x11", X_OK);
 
  if ( gnuplot != 0 )
    printf("\nWARNING: Gnuplot not found. You need to install gnuplot and gnuplot-x11\n\n");
  
  GNUPLOT = gnuplot;


}

void matrix_operations ( int id1, int id2, double value )
{
  double temp;

  temp = gsl_matrix_get( MNA_matrix, id1 - 1, id2 - 1);
  
  temp += 1.0 / value;
  gsl_matrix_set( MNA_matrix, id1 - 1, id2 - 1, temp); 

}

void sparse_matrix_operations ( int i, int j, double value )
{


  cs_entry( MNA_sparse, i - 1, j - 1, 1.0/value);
  
  
}

void construct_matrixes ( void ) 
{
 
  printf("Constructing MNA matrix\n");

  /* Calculate matrix sizes */
  MNA_matrix_size = nodes_counter + voltage_elements_counter;
  
  double temp;
  unsigned int volt_elements = 0;
  int i,j,k=0;
  int DC_scan, PLOT_scan;
  
  if ( NetOptions->TRAN ) {
    b_info = malloc ( sizeof(struct _element *) * MNA_matrix_size );
    for ( i = 0; i < MNA_matrix_size; i++ )
      b_info[i] = malloc ( sizeof(struct _element) );
  
  }

  /* Sparce matrix calculations */
  if ( NetOptions->SPARSE )
  {
    MNA_sparse = cs_spalloc( MNA_matrix_size, MNA_matrix_size, 10, 1, 1);
    C_sparse = cs_spalloc( MNA_matrix_size, MNA_matrix_size, 10, 1, 1);
    b_sparse_vector = malloc ( sizeof(double) * MNA_matrix_size );
  }
  else
  {
    /* Else non sparse matrix calculations */
    /* Initialize matrix */
    MNA_matrix = gsl_matrix_calloc(MNA_matrix_size, MNA_matrix_size);
    C = gsl_matrix_calloc(MNA_matrix_size, MNA_matrix_size);
    S_vector =  gsl_vector_alloc(MNA_matrix_size);
    Ukn_vector =  gsl_vector_alloc(MNA_matrix_size);
    gsl_vector_set_zero(S_vector);

  }

  char **unknown_vars = malloc ( sizeof(char *) * MNA_matrix_size );
  
  for ( i = 0; i < MNA_matrix_size; i++ )
    unknown_vars[i] = NULL;
 

  /* Traverse the Element list */
  element *next_el = Element_list;
  while ( next_el ) 
  {
    /* Add conductances to the MNA matrix */
    if ( !strncmp( next_el->type, "R", sizeof(char)) )
    {
      if ( next_el->node_1 != 0 )
      {
        if ( !unknown_vars[next_el->node_1 - 1] )
        {
          unknown_vars[next_el->node_1 - 1] = malloc ( sizeof(char) * (floor(log10(abs(next_el->node_1))) + 1 + 4));
          sprintf(unknown_vars[next_el->node_1 - 1], "V(%d)", next_el->node_1);
        }
        if ( !NetOptions->SPARSE )
          matrix_operations ( next_el->node_1, next_el->node_1, next_el->value );
        else
          sparse_matrix_operations ( next_el->node_1, next_el->node_1, next_el->value );
      }

      if ( next_el->node_2 != 0 )
      {
        if ( !unknown_vars[next_el->node_2 - 1] )
        {
          unknown_vars[next_el->node_2 - 1] = malloc ( sizeof(char) * (floor(log10(abs(next_el->node_2))) + 1 + 4));
          sprintf(unknown_vars[next_el->node_2 - 1], "V(%d)", next_el->node_2);
        }
        if ( !NetOptions->SPARSE )
          matrix_operations ( next_el->node_2, next_el->node_2, next_el->value );
        else
          sparse_matrix_operations ( next_el->node_2, next_el->node_2, next_el->value );
      }
      if (( next_el->node_1 != 0 ) && ( next_el->node_2 != 0 ))
      { 
        if ( !NetOptions->SPARSE )
        {
          matrix_operations ( next_el->node_1, next_el->node_2, -(next_el->value) );
          matrix_operations ( next_el->node_2, next_el->node_1, -(next_el->value) );
        }
        else
        {
          sparse_matrix_operations ( next_el->node_1, next_el->node_2, -(next_el->value) );
          sparse_matrix_operations ( next_el->node_2, next_el->node_1, -(next_el->value) );
        }
      }

    }

    /* Construct voltage,L and S matrix */
    if ( !strncmp( next_el->type, "V", sizeof(char)) || !strncmp( next_el->type, "L", sizeof(char)) )
    {
      volt_elements++;
      int index = nodes_counter + volt_elements - 1;

      /* Save info */
      if ( (next_el->EXP || next_el->SIN || next_el->PULSE  || next_el->PWL) && NetOptions->TRAN)
        b_info[index] = next_el;

      /* Placement of the unknown current of V's or L's */
      if ( !unknown_vars[index] )
      {
        unknown_vars[index] = malloc ( sizeof(char) * (floor(log10(abs(index + 1))) + 1 + 4));
        sprintf(unknown_vars[index], "I(%d)", index + 1);
      }

      /* Placement of the correlation between voltages and nodes */
      if ( next_el->node_1 != 0 )
      {
        if ( !NetOptions->SPARSE )
        {
          gsl_matrix_set( MNA_matrix, next_el->node_1 - 1, index, 1);
          gsl_matrix_set( MNA_matrix, index, next_el->node_1 - 1, 1);
        }
        else
        {
          cs_entry( MNA_sparse, next_el->node_1 -1, index, 1);
          cs_entry( MNA_sparse, index, next_el->node_1 - 1, 1);
        }
      }
 
      if ( next_el->node_2 != 0 )  
      {
        if ( !NetOptions->SPARSE )
        {
          gsl_matrix_set( MNA_matrix, next_el->node_2 - 1, index, -1);
          gsl_matrix_set( MNA_matrix, index, next_el->node_2 - 1, -1);
        }
        else
        {
          cs_entry( MNA_sparse, next_el->node_2 - 1, index, -1);
          cs_entry( MNA_sparse, index, next_el->node_2 - 1, -1);
        }
      }
      
      if ( !strncmp( next_el->type, "L", sizeof(char)) )
      {
        if ( !NetOptions->SPARSE ) {
          gsl_vector_set( S_vector, index, 0); 
          gsl_matrix_set ( C, nodes_counter + k , nodes_counter + k, next_el->value  );
          k++;
        }
        else {
          b_sparse_vector[index] = 0;

        }

      }
      else
      {
        if ( !NetOptions->SPARSE )
         gsl_vector_set( S_vector, index, next_el->value);
        else
          b_sparse_vector[index] = next_el->value;
      }

      if ( NetOptions->DC )// Check for .DC
        if ( !strcmp( NetOptions->DC_source_name, next_el->name ) &&  !strcmp( NetOptions->DC_type, "V" ) ) 
          DC_scan = index;
      
    }

    /* C */
    if ( !strncmp( next_el->type, "C", sizeof(char)) ) {

      double temp;

      if ( next_el->node_1 != 0 ) {

        if ( !NetOptions->SPARSE ) {

          temp = gsl_matrix_get( C, next_el->node_1 - 1, next_el->node_1 - 1);
          temp += next_el->value;
          gsl_matrix_set( C, next_el->node_1 - 1, next_el->node_1 - 1, temp); 

        }
        else {

          cs_entry( C_sparse, next_el->node_1 - 1, next_el->node_1 - 1, next_el->value);
        }
   
      }
      if ( next_el->node_2 != 0 ) {

        if ( !NetOptions->SPARSE ) {

          temp = gsl_matrix_get( C, next_el->node_2 - 1, next_el->node_2 - 1);
          temp += next_el->value;
          gsl_matrix_set( C, next_el->node_2 - 1, next_el->node_2 - 1, temp); 

        }
        else {

          cs_entry( C_sparse, next_el->node_2 - 1, next_el->node_2 - 1, next_el->value);
        }
   
      }
      if ( ( next_el->node_2 != 0 ) && ( next_el->node_1 != 0 ) ) {
        
        if ( !NetOptions->SPARSE ) {

          temp = gsl_matrix_get( C, next_el->node_1 - 1, next_el->node_2 - 1);
          temp += -(next_el->value);
          gsl_matrix_set( C, next_el->node_1 - 1, next_el->node_2 - 1, temp); 

          temp = gsl_matrix_get( C, next_el->node_2 - 1, next_el->node_1 - 1);
          temp += -(next_el->value);
          gsl_matrix_set( C, next_el->node_2 - 1, next_el->node_1 - 1, temp); 

        }
        else {

          cs_entry( C_sparse, next_el->node_1 - 1, next_el->node_2 - 1, next_el->value);
          cs_entry( C_sparse, next_el->node_2 - 1, next_el->node_1 - 1, next_el->value);
        }
   

      }

    }

    /* Voltages */
    if ( !strncmp( next_el->type, "I", sizeof(char)) )
    {

      /* Save info */
      if ( (next_el->EXP || next_el->SIN || next_el->PULSE  || next_el->PWL) && NetOptions->TRAN)
        b_info[next_el->node_1 - 1] = next_el;


      if ( next_el->node_1 != 0 ) {

        if ( !NetOptions->SPARSE )
          gsl_vector_set( S_vector, next_el->node_1 - 1, -(next_el->value));
        else
          b_sparse_vector[next_el->node_1 - 1] = -(next_el->value);

        if ( NetOptions->DC ) // Check for .DC
          if ( !strcmp( NetOptions->DC_source_name, next_el->name ) &&  !strcmp( NetOptions->DC_type, "I" )  )
            DC_scan = next_el->node_1 - 1;
      }

      if ( next_el->node_2 != 0 )
      {

        if ( !NetOptions->SPARSE )
          gsl_vector_set( S_vector, next_el->node_2 - 1, next_el->value);
        else
          b_sparse_vector[next_el->node_2 - 1] = next_el->value;

        if ( NetOptions->DC )// Check for .DC
          if ( !strcmp( NetOptions->DC_source_name, next_el->name  ) &&  !strcmp( NetOptions->DC_type, "I" ) ) 
            DC_scan = next_el->node_2 - 1;
      }

    }
 
    next_el = next_el->next;

  }

  for ( i = 0; i < MNA_matrix_size; i++)
  {
    if ( !unknown_vars[i] )
    {
      int pos = i ? i : 1;
      pos = floor(log10(abs(pos)));
      unknown_vars[i] = (char *) malloc ( sizeof(char) * pos + 1 + 4);
      sprintf(unknown_vars[i], "V(%d)", i + 1);
    }
  }

 // printf("\n>C matrix:\n");
 // for ( i =0; i < MNA_matrix_size; i++){
  //  for ( j =0; j < MNA_matrix_size; j++){
     // fprintf ( stdout, "%6.4g\t", gsl_matrix_get( C, i, j));
 //   }
    //fprintf ( stdout, "\n");
 // }

  printf("\n>S vector:\n");
  if ( !NetOptions->SPARSE )
    gsl_vector_fprintf (stdout,  S_vector, "%g");
  //else
   // for ( i = 0; i < MNA_matrix_size; i++ )
   //   fprintf( stdout, "%g\n", b_sparse_vector[i] );
 
  printf("\n>MNA matrix:\n");
  

  if ( !NetOptions->SPARSE )
  {  
    for ( i =0; i < MNA_matrix_size; i++){
      for ( j =0; j < MNA_matrix_size; j++){
        fprintf ( stdout, "%6.4g\t", gsl_matrix_get( MNA_matrix, i, j));
      }
      fprintf ( stdout, "\n");
    }
  }
  else
  {
    int success = cs_print( MNA_sparse, "Sparse MNA matrix", 0 );
    if ( success )
      printf("\tSaved in \"Sparse MNA matrix\".\n");
  }
  
  printf("\n\n");

  /* Solve system */
  printf(">Solution: \n");
    
  /* Check if gnuplot is installed only for DC scan */
  GNUPLOT = 1;
  if ( NetOptions->DC )
    check_gnuplot();

  /***********************************************/
  /*****************DC Analysis*******************/
  /***********************************************/

  if ( NetOptions->SPARSE ) {
    /* Compress MNA */
    A = cs_compress(MNA_sparse);
    
    /* Find duplicates */
    int success = cs_dupl ( A );
    if ( !success )
      printf("Error removing duplicates\n");

    cs_spfree(MNA_sparse);

  }

  if ( NetOptions->SPD && !NetOptions->ITER )
  { // Solve with Cholesky

    if ( NetOptions->SPARSE )
      Cholesky_Sparse( unknown_vars );
    else
      Cholesky_solve(DC_scan, unknown_vars);
    
    printf("\tSaved in Cholesky.txt \n");
  }
  else if ( NetOptions->ITER && !NetOptions->SPD)
  {
    if ( NetOptions->SPARSE )
      Bi_CG_Sparse( unknown_vars );
    else
      Bi_CG(DC_scan, unknown_vars);

    printf("\tSaved in Bi_CG_solution.txt \n");
  }
  else if ( NetOptions->ITER && NetOptions->SPD)
  {
    if ( NetOptions->SPARSE )
      CG_Sparse( unknown_vars  );
    else
      CG(DC_scan, unknown_vars);
    printf("\tSaved in CG_solution.txt \n");
  }
  else
  { // Solve with LU
    if ( NetOptions->SPARSE )
      LU_Sparse(unknown_vars);
    else
      LU_solve(DC_scan, unknown_vars);
    printf("\tSaved in LU_solution.txt \n");
  }

  /* Trasient Analysis */
  if ( NetOptions->TRAN && NetOptions->BE ) {

    if ( !NetOptions->SPARSE )
      BE( unknown_vars );
    else
      BE_sparse( unknown_vars );

  }
  else {

    if ( !NetOptions->SPARSE )
      TR( unknown_vars );
    else
      TR_sparse( unknown_vars );

  }

  /* Cleanup mem */
  if ( !NetOptions->SPARSE )
  {
    gsl_matrix_free (MNA_matrix);
    gsl_vector_free (S_vector);
    gsl_vector_free (Ukn_vector);

    for ( i = 0; i < MNA_matrix_size; i++) 
      free(unknown_vars[i]);
  }


}
