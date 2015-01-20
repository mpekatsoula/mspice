#include "hms/parser.h"
#include "hms/error_header.h"
#include "hms/hash.h"
#include <gsl/gsl_math.h>

/* Add nodes to hash table */
int 
add_node_to_hash_table( char *token ) 
{

  int result;
  
  if ( token[0] == '0' && strlen(token) == 1 ) 
  {
    ground_flag = 1;
    result = ht_set( nodes_hashtable, token, 0 );
    return 0;
  }
  else 
  {
    nodes_counter++;
    result = ht_set( nodes_hashtable, token, nodes_counter );
  }

  /* If the result is not the same, this means that the 
   * node was already imported. So we need to reduce the
   * nodes_counter                                      */
  if ( result != nodes_counter )
    nodes_counter--;

  return result;
}


/* Check if a string is in A-Z,a-z,_. 1 if yes, 0 if no */
int 
check_ascii ( char * string) 
{

  int i;

  for ( i = 0; i < strlen(string); i++ )
    if ( ( string[i] == '\n' ) && (string[i] < 65 || ( string[i] > 90 && string[i] < 97 ) || string[i] > 122 || string[i] != '_' ) )
      return 0;

  return 1;

}

/* Check token */
void 
check_token( char *token, int position, int line_num) 
{

  int i;

  switch ( position ) {

    /* Check circuit element name */
    case 0: 
      if ( token[0] != 'V' && token[0] != 'v' && token[0] != 'I' && token[0] != 'i' &&
           token[0] != 'R' && token[0] != 'r' && token[0] != 'C' && token[0] != 'c' &&
           token[0] != 'L' && token[0] != 'l' && token[0] != 'Q' && token[0] != 'q' &&
           token[0] != 'D' && token[0] != 'd' && token[0] != 'M' && token[0] != 'm' )
        error( LEX_1, line_num );

      if ( token[1] == '\0' )
        error( LEX_2, line_num );

      /* lowercase to uppercase */
      if ( token[0] == 'v' || token[0] == 'V' )
      {
        voltage_elements_counter++;
        token[0] = 'V';
      }
      else if ( token[0] == 'c' )
        token[0] = 'C';
      else if ( token[0] == 'r' || token[0] == 'R' ) {
        token[0] = 'R';
        resistances++;
      }
      else if ( token[0] == 'i' || token[0] == 'I' ) {
        i_elements_counter++;
        token[0] = 'I';
      }
      else if ( token[0] == 'l' || token[0] == 'L' )
      {
	      voltage_elements_counter++;
        token[0] = 'L';
      }
      else if ( token[0] == 'd' )
        token[0] = 'D';
      else if ( token[0] == 'q' )
        token[0] = 'Q';
      else if ( token[0] == 'm' )
        token[0] = 'M';

    break;

    /* Check nodes' validity */
    case 1:
    case 2:
      if ( token[0] == '-'  ) 
        error( LEX_3, line_num );


      for ( i = 0; i < strlen(token); i++ ) 
      {
        if ( ( token[i]  < 48 || ( token[i] > 57 && token[i] < 65 ) || ( token[i] > 90 && token[i] < 97 ) || token[i] > 122 ) && token[i] != '_'  )
          error( LEX_4, line_num );

      }
      

    break;

    case 3:

      if ( token[strlen(token) - 1] == 'e' || token[strlen(token) - 1] == 'E' ) // If 'e' is in the last position, error ( -1 for the \0 )
        error( LEX_5, line_num );  

    break;
  }

}

/* Lexical analysis for the line */
void 
lexical_analysis( char *line ) 
{

  char *token;
  int position = 0;
  element *next_el;

  // Read tokens
  token = strtok( line, "\t\n ");

  /* Assign a special position for the options */
  if ( token[0] == '.' )
    position = OPTIONS_CASE;
  else if ( position == OPTIONS_CASE || position == DC_CASE || position == PLOT_CASE )
    position = position; // do nothing
  else
  {
    /* Allocate memory for element list */
    if ( !Element_list ) 
    {
      Element_list = ( element * ) malloc ( sizeof ( struct _element ) );
      Element_list->next = NULL;
      
      // Init values for error checking
      Element_list->node_3 = -1;
      Element_list->node_4 = -1;
      Element_list->L = -1; 
      Element_list->W = -1;
      Element_list->EXP = 0;
      Element_list->SIN = 0;
      Element_list->PULSE = 0;
      Element_list->PWL = 0;
      next_el = Element_list;
      last_pointer = next_el;
    }
    else 
    {
      // Load last pointer
      next_el = last_pointer;
      next_el->next = ( element * ) malloc ( sizeof ( struct _element ) );
      next_el = next_el->next;

      // Init values for error checking
      next_el->node_3 = -1;
      next_el->node_4 = -1;
      next_el->L = -1;
      next_el->W = -1;
      next_el->next = NULL;
      next_el->EXP = 0;
      next_el->SIN = 0;
      next_el->PULSE = 0;
      next_el->PWL = 0;
      last_pointer = next_el; // set last pointer as new
    }
  }

  while ( token ) 
  {


    /* Fill ELement list with info */
    switch ( position ) {
      case 0:
        check_token( token, position, line_num );
        next_el->type = (char *) malloc ( sizeof (char) );
        strncpy(next_el->type,token,sizeof(char));
        next_el->name = strdup( token + sizeof(char));
      break;
      case 1:
        check_token( token, position, line_num );
        // Add node to hash table and add the info to the element list
        next_el->node_1 = add_node_to_hash_table( token );
      break;
      case 2:
        check_token( token, position, line_num );
        // Add node to hash table and add the info to the element list
        next_el->node_2 = add_node_to_hash_table( token );
      break;
      case 3:
        if ( strncmp(next_el->type,"D",sizeof(char)) && strncmp(next_el->type,"M",sizeof(char)) && strncmp(next_el->type,"Q",sizeof(char)) ) {
          check_token( token, position, line_num );
          next_el->value = atof( token );
          next_el->model_name = NULL; // we have no model name here
        }
        // Copy model name
        if ( !strncmp(next_el->type,"D",sizeof(char)) )
          next_el->model_name = strdup(token);
        
        if ( !strncmp(next_el->type,"Q",sizeof(char)) )
          next_el->node_3 = add_node_to_hash_table( token );
        
      break;
      case 4:        
        /* TRAN spec */
        if ( strncmp(next_el->type,"V",sizeof(char) ) || strncmp(next_el->type,"I",sizeof(char) ) ) {
          /* EXP ( i1 i2 td1 tc1 td2 tc2 ) */
          if ( !strncmp(token, "EXP", sizeof(char) * 3) || !strncmp(token, "exp", sizeof(char) * 3) ) {
            next_el->EXP = 1;

            token = strtok( NULL, "(,\t\n ");
            if ( !token )
              error( SYNTAX_1, line_num );
            
            /* i1 */
            next_el->i1 = atof(token);
            token = strtok( NULL, ",\t\n ");
            if ( !token )
              error( SYNTAX_1, line_num );

            /* i2 */
            next_el->i2 = atof(token);
            token = strtok( NULL, ",\t\n ");
            if ( !token )
              error( SYNTAX_1, line_num );

            /* td1 */
            next_el->td1 = atof(token);
            token = strtok( NULL, ",\t\n ");            
            if ( !token )
              error( SYNTAX_1, line_num );

            /* tc1 */
            next_el->tc1 = atof(token);
            token = strtok( NULL, ",\t\n "); 
            if ( !token )
              error( SYNTAX_1, line_num );

            /* td2 */
            next_el->td2 = atof(token);
            token = strtok( NULL, ",)\t\n "); 
            if ( !token )
              error( SYNTAX_1, line_num );

            /* tc2 */
            next_el->tc2 = atof(token);
            token = strtok( NULL, ",)\t\n "); 
          }
          /* SIN (i1 ia fr td df ph) */
          else if ( !strncmp(token, "SIN", sizeof(char) * 3) || !strncmp(token, "sin", sizeof(char) * 3) ) {
            next_el->SIN = 1;

            token = strtok( NULL, "(,\t\n ");
            if ( !token )
              error( SYNTAX_1, line_num );
            
            /* i1 */
            next_el->i1 = atof(token);
            token = strtok( NULL, ",\t\n ");
            if ( !token )
              error( SYNTAX_1, line_num );

            /* ia */
            next_el->ia = atof(token);
            token = strtok( NULL, ",\t\n ");
            if ( !token )
              error( SYNTAX_1, line_num );

            /* fr */
            next_el->fr = atof(token);
            token = strtok( NULL, ",\t\n ");            
            if ( !token )
              error( SYNTAX_1, line_num );

            /* td */
            next_el->td = atof(token);
            token = strtok( NULL, ",\t\n "); 
            if ( !token )
              error( SYNTAX_1, line_num );

            /* df */
            next_el->df = atof(token);
            token = strtok( NULL, ",)\t\n "); 
            if ( !token )
              error( SYNTAX_1, line_num );

            /* ph */
            next_el->ph = atof(token);
            token = strtok( NULL, ",)\t\n "); 
          }
          /* PULSE i1 i2 td tr tf pw per) */
          else if ( !strncmp(token, "PULSE", sizeof(char) * 5) || !strncmp(token, "pulse", sizeof(char) * 5) ) {
            next_el->PULSE = 1;

            token = strtok( NULL, "(,\t\n ");
            if ( !token )
              error( SYNTAX_1, line_num );
            
            /* i1 */
            next_el->i1 = atof(token);
            token = strtok( NULL, ",\t\n ");
            if ( !token )
              error( SYNTAX_1, line_num );

            /* i2 */
            next_el->i2 = atof(token);
            token = strtok( NULL, ",\t\n ");
            if ( !token )
              error( SYNTAX_1, line_num );

            /* td */
            next_el->td = atof(token);
            token = strtok( NULL, ",\t\n ");            
            if ( !token )
              error( SYNTAX_1, line_num );

            /* tr */
            next_el->tr = atof(token);
            token = strtok( NULL, ",\t\n "); 
            if ( !token )
              error( SYNTAX_1, line_num );

            /* tf */
            next_el->tf = atof(token);
            token = strtok( NULL, ",\t\n "); 
            if ( !token )
              error( SYNTAX_1, line_num );

            /* pw */
            next_el->pw = atof(token);
            token = strtok( NULL, ",\t\n ");
            if ( !token )
              error( SYNTAX_1, line_num );

            /* per */
            next_el->per = atof(token);
            token = strtok( NULL, "),\t\n ");

          }
          /* PWL (t1 i1) ... (tn in) */
          else if ( !strncmp(token, "PWL", sizeof(char) * 3) || !strncmp(token, "pwl", sizeof(char) * 3) ) {
            next_el->PWL = 1;

            token = strtok( NULL, "(,\t\n ");
            if ( !token )
              error( SYNTAX_1, line_num );
            
            /* t1 */
            next_el->pwl_t[0] = atof(token);
            token = strtok( NULL, "),\t\n ");
            if ( !token )
              error( SYNTAX_1, line_num );

            
            /* i1 */
            next_el->pwl_i[0] = atof(token);
            token = strtok( NULL, "(,\t\n ");
            if ( !token )
              error( SYNTAX_1, line_num );

            /* t2 */
            next_el->pwl_t[1] = atof(token);
            token = strtok( NULL, "),\t\n ");
            if ( !token )
              error( SYNTAX_1, line_num );

            /* i2 */
            next_el->pwl_i[1] = atof(token);
            token = strtok( NULL, "(,\t\n ");
            if ( !token )
              break;

            int pairs = 2;
            while ( pairs < MAX_PWL_PAIRS ) {

              /* tn */
              next_el->pwl_t[pairs] = atof(token);
              token = strtok( NULL, "(,\t\n ");
              if ( !token )
                error( SYNTAX_1, line_num );

              /* in */
              next_el->pwl_i[pairs] = atof(token);
              token = strtok( NULL, "(,\t\n ");
              if ( !token )
                break;

              pairs++;

            }
            if ( pairs == MAX_PWL_PAIRS )
              printf("Cannot parse more than %d PWL pairs.\n", MAX_PWL_PAIRS);
           
            next_el->n = pairs - 1;
          }
        }

        if ( !strncmp(next_el->type,"D",sizeof(char)) )
          next_el->area = atoi(token);

        if ( !strncmp(next_el->type,"Q",sizeof(char)) )
          next_el->model_name = strdup(token);

        if ( !strncmp(next_el->type,"M",sizeof(char)) )
          next_el->node_4 = add_node_to_hash_table( token );
      break;
      case 5:
        if ( !strncmp(next_el->type,"Q",sizeof(char)) )
          next_el->area = atoi(token);
        if ( !strncmp(next_el->type,"M",sizeof(char)) )
          next_el->model_name = strdup(token);
      break;
      case 6:
      case 7:
        /* Check for the format L= or W=. Also check for duplicated value, eg: L=2 L=3 ( two L's) */
        if ( token[0] == 'L' && token[1] == '=' && next_el->L == -1 )
          next_el->L = atof(token + 2);
        else if ( token[0] == 'W' && token[1] == '=' && next_el->W == -1 )
          next_el->W = atof(token + 2);
        else
          error( LEX_6, line_num );

      break;
      case OPTIONS_CASE:
        /* Check for .OPTIONS or .DC or PLOT/PRINT */
        if ( !strcmp((token + sizeof(char)), "OPTIONS") )
          NetOptions->OPTIONS = 1;
    
        /* Check for DC */
        if ( !strncmp((token + sizeof(char)), "DC", sizeof(char) * 2) ) {
          NetOptions->DC = 1; 
          position = DC_CASE;
        }

        /* Check for TRAN */
        if ( !strcmp((token + sizeof(char)), "TRAN" )) { 
          NetOptions->TRAN = 1;
          position = TRAN_CASE; 
        }
        /* Check for PLOT */
        if ( !strcmp((token + sizeof(char)), "PLOT") || !strcmp((token + sizeof(char)), "PRINT") ) {
          NetOptions->PLOT = 1;
          position = PLOT_CASE;
        }
        /* Check for SPD */ 
        if ( !strncmp(token, "SPD", strlen(token) ) && NetOptions->OPTIONS )
          NetOptions->SPD = 1; 

        if ( !strncmp(token, "SPARSE", strlen(token) ) && NetOptions->OPTIONS )
          NetOptions->SPARSE = 1; 

        /* Check for ITER */
        if ( !strncmp(token, "ITER", strlen(token) ) && NetOptions->OPTIONS )
          NetOptions->ITER = 1;

        /* Check for ITOL */
        if ( !strncmp(token, "ITOL=", sizeof(char) * 5 ) && NetOptions->OPTIONS )
          NetOptions->ITOL = atof(token + sizeof(char) * 5);      

        /* Check for METHOD */
        if ( !strncmp(token, "METHOD=", sizeof(char) * 7 ) && NetOptions->OPTIONS )
          if ( !strncmp(token + sizeof(char) * 7, "TR", sizeof(char) * 2 ) )
            NetOptions->TR = 1;  
          else if ( !strncmp(token + sizeof(char) * 7, "BE", sizeof(char) * 2 ) )
            NetOptions->BE = 1;  

      break;
      case DC_CASE:
        if ( !token )
          error( SYNTAX_1, line_num );
 
        NetOptions->DC_source_name = strndup(token + sizeof(char), strlen(token) - sizeof(char) );
        NetOptions->DC_type = strndup ( token, sizeof(char) );

        /* Get next token, no need to re-run the loop */
        token = strtok( NULL, "\t\n ");
        if ( !token )
          error( SYNTAX_1, line_num );
        NetOptions->DC_start = atof(token);
        token = strtok( NULL, "\t\n ");
        if ( !token )
          error( SYNTAX_1, line_num );
        NetOptions->DC_end = atof(token);
        token = strtok( NULL, "\t\n ");
        if ( !token )
          error( SYNTAX_1, line_num );
        NetOptions->DC_step = atof(token);

      break;
      case PLOT_CASE:
        if ( flag ) {
          while ( NetOptions->num_of_tran_plots < 256) {

            token = strtok( NULL, "IV()\n\t "); 
            if ( !token )
              break;

            NetOptions->PLOT_scan[NetOptions->num_of_tran_plots] = strndup(token + sizeof(char)*2, strlen(token) - sizeof(char)*3);
            NetOptions->num_of_tran_plots++; // increase counter for dc plots

          }
          if ( NetOptions->num_of_dc_plots == 256 )
            error ( line_num, PLOT_1 );

          flag = 0;
        }
        else {

          while ( NetOptions->num_of_dc_plots < 256 ) {

            token = strtok( NULL, "IV()\n\t "); 
            if ( !token )
              break;

            NetOptions->PLOT_scan[NetOptions->num_of_dc_plots] = strndup(token, sizeof(token));
            NetOptions->num_of_dc_plots++; // increase counter for dc plots
            
            
          }

          if ( NetOptions->num_of_dc_plots == 256 )
            error ( line_num, PLOT_1 );
        }
      break;
      case TRAN_CASE:
        flag = 1;
        if ( !token )
          error( SYNTAX_1, line_num );
        NetOptions->TRAN_step = atof(token);

        token = strtok( NULL, "\t\n ");
        if ( !token )
          error( SYNTAX_1, line_num );
        NetOptions->TRAN_fin_time = atof(token);
      break;     
    }

    token = strtok( NULL, "(\t\n ");
      
    /* Checking if elements per line are missing && if area is pressent */
    if ( !token && position < 3  ) // Less than 3 elements per line, error
      error ( SYNTAX_1, line_num );

    if ( position != OPTIONS_CASE && position != DC_CASE && position != PLOT_CASE && position != TRAN_CASE )
    {
      if ( !strncmp(next_el->type,"D",sizeof(char)) && token && position == 2 && !check_ascii(token) ) // If model name for D is not present or is in not in ascii error
        error( SYNTAX_2, line_num);

      if ( !strncmp(next_el->type,"D",sizeof(char)) && !token && position == 3 ) // If area for D is not present, set to 1
        next_el->area = 1;

      if ( !strncmp(next_el->type,"Q",sizeof(char)) && !token && position < 4 ) // If Q and less than 4 elements error
        error ( SYNTAX_1, line_num );

      if ( !strncmp(next_el->type,"Q",sizeof(char)) && token && position == 3 && !check_ascii(token) ) // If model name for Q is not present or is in not in ascii error
        error( SYNTAX_2, line_num);

      if ( !strncmp(next_el->type,"Q",sizeof(char)) && !token && position == 4 ) // If area for Q is not present, set to 1
        next_el->area = 1;

      if ( !strncmp(next_el->type,"M",sizeof(char)) && !token && position < 7 ) // If M and less than 7 elements error
        error ( SYNTAX_1, line_num );    

      if ( !strncmp(next_el->type,"M",sizeof(char)) && token && position == 4 && !check_ascii(token) ) // If model name for M is not present or is in not in ascii error
        error( SYNTAX_2, line_num);

      /* At last, check the number of elements per line */
      if ( token && ( ( (!strncmp(next_el->type,"V",sizeof(char)) || !strncmp(next_el->type,"R",sizeof(char)) || !strncmp(next_el->type,"I",sizeof(char)) ) && position == 4 ) ||
             (!strncmp(next_el->type,"Q",sizeof(char)) && position == 5) || ( !strncmp(next_el->type,"M",sizeof(char)) && position == 7 ) || ( !strncmp(next_el->type,"D",sizeof(char)) && position == 4 ) ) )
        error ( SYNTAX_3, line_num ); 

      position++;
    }
      

  }

}

void 
parser(char *filename) 
{

  printf("Now parsing file: %s\n", filename);

  int fp;
  char * line = NULL;
  size_t len = 0, file_size;
  ssize_t read;
  line_num = 0;
  struct stat sb;

  ground_flag = 0;

  // Open file
  fp = open( filename, O_RDWR );
  fstat(fp, &sb);

  char * memblock = mmap(NULL, sb.st_size, PROT_WRITE, MAP_PRIVATE, fp, 0);

  if (memblock == MAP_FAILED)
     error( FILE_SIZE_ERROR, -1 );

  if ( fp ) {

    /* Before line parsing create hash table for nodes */
    nodes_hashtable = ht_create( HASH_SIZE );
    nodes_counter = 0;
    voltage_elements_counter = 0;
    resistances = 0;
    i_elements_counter = 0;
    
    /* Allocate memory for NetOptions and init */
    NetOptions = malloc ( sizeof(options) );
    NetOptions->OPTIONS = 0;
    NetOptions->DC = 0;
    NetOptions->SPD = 0;
    NetOptions->DC_end = 1;
    NetOptions->DC_start = 0;
    NetOptions->DC_step = 2;
    NetOptions->PLOT = 0;
    NetOptions->SPARSE = 0;
    NetOptions->ITER = 0;
    NetOptions->TRAN = 0;
    NetOptions->TR = 0;
    NetOptions->BE = 0;
    NetOptions->ITOL = atof("1e-6");
    NetOptions->num_of_dc_plots = 0;
    NetOptions->num_of_tran_plots = 0;
    NetOptions->DC_source_name = NULL;
    NetOptions->DC_type = NULL;
    flag = 0;
    
    // Parse file, get next line
    char *temp;
    line = strtok_r( memblock ,"\n", &temp);

    while ( line ) {

      line_num++;

      // If line starts with * or is an empty line, drop the line
      if ( line[0] == '*' || line[0] == '\n' )  
      {
        line = strtok_r( NULL ,"\n", &temp);
        continue;
      }
      else
      {
        // Lexical analysis for each line
        lexical_analysis( line );
        line = strtok_r( NULL ,"\n", &temp);
      }
    }

    // unmap memory block
    munmap( memblock, sb.st_size );

    printf("Parsing done.\nTotal lines: %d\n\nNodes added (without ground): %d.\n#v: %d\n",line_num,nodes_counter,voltage_elements_counter);
    printf("#i: %d\n\n", i_elements_counter);
    
    if ( !ground_flag )
      error ( MISSING_GROUND, -1 );

  }
  else
    error( FILE_ERROR, -1 );
 

  close( fp );  

}
