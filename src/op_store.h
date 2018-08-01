//                              op_store.h
//                Matrix Operations Stores for LARC
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


#ifndef LARC_OP_STORE_H
#define LARC_OP_STORE_H

#include <stdint.h>
#include "larc.h"

/* lists out all the op_names and description of them 
   to help users interface with op chain functions  */
void list_op_names();

/* Returns 1 on successfully creating the operations store */
int create_op_store(size_t exponent);

/* Returns 0 on success, or error code if failed  */
int op_set(op_type_t op, mat_add_t in1, mat_add_t in2, mat_add_t out);

/* Returns the matrix ID previously stored, or MATRIX_PTR_INVALID if not found */
mat_add_t op_get(op_type_t op, mat_add_t in1, mat_add_t in2);

/* Prints a summary of operations store usage */
void op_store_report(char *outfilepath);

int64_t op_hashID_by_matrixIDs(int64_t in1_mID, int64_t in2_mID, char *op_name);

// prints a hash chain from the opstore to a file.
int op_hash_chain_info_to_file(uint64_t hash, char *outfilepath, char *comment);

// prints a hash chain from the opstore to screen.
int op_hash_chain_info_to_screen(uint64_t hash, char *comment);

// Cleans a hash chain in the op store for a given hash
int clean_op_hash_chain(uint64_t hashID);


#ifdef HASHSTATS
/* prints hash table statistics for op store, creating files with 
   accesses and nodes for each bucket (hash_val).  It also produces
   as standard hash_report to stdout or file.
   This works only if HASHSTATS is defined.     */
void op_hashstats(char *accesses_file,  char *nodes_file, char *report_file);
#endif

#endif
