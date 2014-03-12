#include <stdio.h>
#include "hms/error_header.h"
#include "hms/parser.h"
#include "hms/cleanup.h"
#include "hms/hash.h"

int main(int argc, char *argv[]) {

  // Check arguments
  if ( argc < 2 )
    error(ARGUMENTS_ERROR, -1);

  // Call parser
  parser( argv[1] );

  // Construct matrixes
  construct_matrixes ( );

  // Free memory
  CleanUp();

  return 1;

}
