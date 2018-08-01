//                       hash.h
/******************************************************************
 *                                                                *
 * Copyright 2014, Institute for Defense Analyses                 *
 * 4850 Mark Center Drive, Alexandria, VA; 703-845-2500           *
 * This material may be reproduced by or for the US Government    *
 * pursuant to the copyright license under the clauses at DFARS   *
 * 252.227-7013 and 252.227-7014.                                 *
 *                                                                *
 * LARC : Linear Algebra via Recursive Compression                *
 * Authors:                                                       *
 *   - Steve Cuccaro (IDA-CCS)                                    *
 *   - John Daly (LPS)                                            *
 *   - John Gilbert (UCSB, IDA adjunct)                           *
 *   - Jenny Zito (IDA-CCS)                                       *
 *                                                                *
 * Additional contributors are listed in "LARCcontributors".      *
 *                                                                *
 * POC: Jennifer Zito <jszito@super.org>                          *
 * Please contact the POC before disseminating this code.         *
 *                                                                *
 *****************************************************************/

#ifndef LARC_HASH_H
#define LARC_HASH_H

/******************************************************************************
 *     Hash functions for Global Recursive Compression Scheme                 *
 *     for rapid find and insert functions for both the matrix and            *
 *     operations stores                                                      *
 *     last updated:   Fall 2015                                              *
 ******************************************************************************/

#include <stdint.h> // uint64_t
#include "matrix_store.h"


/*
 * allocate a new hash table
 */
hash_table_t *alloc_hash(int exponent);
hash_node_t *alloc_hash_node();
/*
 * lookup/insert/delete an entry
 */
hash_node_t *hash_lookup(hash_table_t *table, record_ptr_t record_ptr, uint64_t hash_val);
hash_node_t *hash_insert_at_tail(hash_table_t *table, record_ptr_t record_ptr, uint64_t hash_val);
hash_node_t *hash_insert_at_head(hash_table_t *table, record_ptr_t record_ptr, uint64_t hash_val);
int hash_node_remove(hash_table_t *table, record_ptr_t record_ptr, uint64_t hash_val);
/*
 * report out
 */
uint64_t hash_report(hash_table_t *table, FILE *fp, int verbose);

/*
 * Utility
 */
static inline uint64_t
hash_max(hash_table_t *table)
{
  return 1<<table->exponent;
}

// creates a mask that grabs the low exponent number of bits
// of the hash_val and then references that entry in the heads array
static inline hash_node_t *
hash_get_chain(hash_table_t *table, uint64_t hash_val)
{
  return table->heads[hash_val&((1<<table->exponent)-1)];
}


/* Returns a 64 bit hash value, where the input is of ScalarType and
 * exponent is the hash table exponent */
uint64_t
hash_from_matrix_scalar(ScalarType scalar, matrix_type_t mat_type, uint64_t exponent);

/* Returns a 64 bit hash value, where the input is an array of 4 mat_add_t 
 * and exponent is the hash table exponent */
uint64_t
hash_from_matrix_panel(mat_add_t sub[4], matrix_type_t mat_type, uint64_t exponent);

// obsolete
// uint64_t hash_from_four_integers 
//   (uint64_t int1, uint64_t int2, uint64_t int3, uint64_t int4, uint64_t exponent);

// obsolete
// uint64_t hash_from_int_list (uint64_t *targets, uint64_t targ_len, uint64_t exponent);

// obsolete
// uint64_t uint64_from_4uint64 (uint64_t int1, uint64_t int2, uint64_t int3, uint64_t int4);

uint64_t 
recursive_hash_from_two_integers
   (uint64_t int1, uint64_t int2, uint64_t exponent);

uint64_t 
recursive_hash_from_three_integers
   (uint64_t int1, uint64_t int2, uint64_t int3, uint64_t exponent);

uint64_t 
recursive_hash_from_four_integers 
   (uint64_t int1, uint64_t int2, uint64_t int3, uint64_t int4, uint64_t exponent);

uint64_t
recursive_hash_from_int_list (uint64_t *tlist, uint64_t tlength, uint64_t exponent);



#ifdef HASHSTATS
// hash table statistics produced when HASHSTATS is defined 
void hashstats_to_files( hash_table_t *hash_table, 
                         char *accesses_file, 
                         char *nodes_file, 
                         char *report_file
			); 

#endif


#endif
