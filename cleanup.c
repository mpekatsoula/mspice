#include <stdio.h>
#include "hms/parser.h"
#include "hms/cleanup.h"
#include "hms/hash.h"

void 
CleanUp( void ) 
{

  int i;

 /* Clean up the Element List */
  element *next_el, *old_el;
  next_el = Element_list;
  while ( next_el ) 
  {
    free(next_el->type);
    free(next_el->name);
    free(next_el->model_name);
    old_el = next_el;
    next_el = next_el->next; 
    free(old_el);
  }

  /* Clean up hash table */
  entry_t *next, *old;
  for ( i = 0; i < HASH_SIZE; i++ )
  {
    next = nodes_hashtable->table[i];
    while ( next ) 
    {
      old = next;
      next = next->next;
      free( old->key );
      free( old );
    }
  
  }
  free ( nodes_hashtable->table );
  free ( nodes_hashtable );

  for ( i = 0; i < NetOptions->num_of_dc_plots; i ++ )
    free(NetOptions->PLOT_scan[i]);

  free(NetOptions->DC_source_name);
  free(NetOptions->DC_type);
  free( NetOptions );
}
