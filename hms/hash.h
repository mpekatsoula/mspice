#ifndef _HASH_H
#define _HASH_H

#define HASH_SIZE 0xC0FFE

struct entry_s {
	char *key;
	unsigned int value;
	struct entry_s *next;
};

typedef struct entry_s entry_t;

struct hashtable_s {
	int size;
	struct entry_s **table;	
};

typedef struct hashtable_s hashtable_t;

/* Nodes hash table */
hashtable_t *nodes_hashtable;

hashtable_t *ht_create( int size );
int ht_hash( hashtable_t *hashtable, char *key );
entry_t *ht_newpair( char *key, unsigned int value );
unsigned int ht_set( hashtable_t *hashtable, char *key, unsigned int value );
int ht_get( hashtable_t *hashtable, char *key );

#endif
