#include "hms/error_header.h"

void error( int errcode, int line_num) {

  switch ( errcode ) {
    case ARGUMENTS_ERROR:
      printf("Usage: [executable name] [input_file].\n");
      exit(1);
    break;
    case FILE_ERROR:
      printf("Error opening file.\n");
      exit(1);
    break;    
    case FILE_SIZE_ERROR:
      printf("File is empty. Exiting..\n");
      exit(1);
    break;   
    case LEX_1:
      printf("Lexical error in line %d. Invalid circuit element. Valid elements are: V,I,R,C,L\n", line_num);
      exit(1);
    break;  
    case LEX_2:
      printf("Lexical error in line %d.. Element must have a name.\n", line_num);
      exit(1);
    break;  
    case LEX_3:
      printf("Lexical error in line %d.. Element node cannot be negative.\n", line_num);
      exit(1);
    break; 
    case LEX_4:
      printf("Lexical error in line %d.. Invalid node name.\n", line_num);
      exit(1);
    break; 
    case LEX_5:
      printf("Lexical error in line %d.. Invalid numerical value. No number after e.\n", line_num);
      exit(1);
    break; 
    case LEX_6:
      printf("Lexical error in line %d.. Value must be in the form L=2nm or W=2nm, or you have duplicated value. \n", line_num);
      exit(1);
    break;    
    case MISSING_GROUND:
      printf("Error, ground not present in netlist file.\n");
      exit(1);
    break; 
    case SYNTAX_1:
      printf("Syntax error in line %d. Missing value.\n", line_num);
      exit(1);
    break; 
    case SYNTAX_2:
      printf("Syntax error in line %d. Incorrect/missing model name.\n", line_num);
      exit(1);
    break; 
    case SYNTAX_3:
      printf("Syntax error in line %d. Too many arguments.\n", line_num);
      exit(1);
    break;
    case PLOT_1:
      printf("Error! Cannot plot more than 256 variables.\n");
      exit(1);
    break;
  }



}
