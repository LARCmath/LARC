//                  op_store.c
//   Matrix Operations Stores for LARC
/******************************************************************
 *                                                                *
 * Copyright (C) 2014, Institute for Defense Analyses             *
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
 *   - Mark Pleszkoch (IDA-CCS)                                   *
 *   - Jenny Zito (IDA-CCS)                                       *
 *                                                                *
 * Additional contributors are listed in "LARCcontributors".      *
 *                                                                *
 * Questions: larc@super.org                                      *
 *                                                                *
 * All rights reserved.                                           *
 *                                                                *
 * Redistribution and use in source and binary forms, with or     *
 * without modification, are permitted provided that the          *
 * following conditions are met:                                  *
 *   - Redistribution of source code must retain the above        *
 *     copyright notice, this list of conditions and the          *
 *     following disclaimer.                                      *
 *   - Redistribution in binary form must reproduce the above     *
 *     copyright notice, this list of conditions and the          *
 *     following disclaimer in the documentation and/or other     *
 *     materials provided with the distribution.                  *
 *   - Neither the name of the copyright holder nor the names of  *
 *     its contributors may be used to endorse or promote         *
 *     products derived from this software without specific prior *
 *     written permission.                                        *
 *                                                                *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND         *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,    *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF       *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE       *
 * DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER NOR        *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT   *
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;   *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)       *
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN      *
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR   *
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, *
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.             *
 *                                                                *
 *****************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <errno.h>

#include "larc.h"
#include "hash.h"
#include "matrix_store.h"
#include "op_store.h"

/*!
 * \file op_store.c
 * \brief Contains the static data structure for the operations store and the
 * functions which allow user access to the store
 */


/* Define the struct here as it doesn't need to be used elsewhere */
static struct op_store_t {
  uint64_t longest;         // length in characters of longest op name
  uint64_t exponent;        // op_store has hash table of size 2^exponent
  uint64_t size;            // op_store has hash table of size 2^exponent
  uint64_t numsets;         // number of times op_set is called on valid
  uint64_t op_store_hits;   // incremented in op_get when success
  uint64_t op_store_misses; // incremented in op_get when failed
  hash_table_t *hash_table; // pointer to hash table for op_store
  uint64_t oprecords_total_created[INVALID_OP+1];
  uint64_t oprecords_stored_now[INVALID_OP+1];
} op_store = {0};


/* The types of operations to store (OPERATIONS is defined in op_store.h) */
#define X(a, b, c) b,
static char *op_names[] = { OPERATIONS };
#undef X

#define X(a, b, c) c,
static char *op_names_verbose[] = { OPERATIONS };
#undef X

// This function is to help users interface with op chain functions
void list_op_names() {
  printf("The operation names understood by LARC are:\n");
  for (int i=0; i<=INVALID_OP; ++i)
  {
    printf("    %s\n",op_names_verbose[i]);
  }
}

void
op_store_report(char *outfilepath)
{
  FILE *f; 
  if (strcmp(outfilepath,"stdout")) {
    printf("Printing op store report to file %s\n", outfilepath);
    f = fopen(outfilepath, "a"); 
    
  }
  else {
    printf("Printing op store report to screen\n");
    f = stdout;
  }
  
  int64_t size;
  fprintf(f,"\n");
  fprintf(f,"Op store report:\n");
  size = hash_report(op_store.hash_table, f, 0);
  fprintf(f,"Statistics for each operation:\n");
  for (int i=0; i<INVALID_OP; ++i)
  {
    fprintf(f,"%-14s: %12" PRIu64 " records created, %12" PRIu64 " current records\n",
       op_names[i], op_store.oprecords_total_created[i],
       op_store.oprecords_stored_now[i]);
  }
  fprintf(f,"Total Op Store size: %3.2fMB\n", size/(1024.0*1024.0));
  fflush(f);
  if  (strcmp(outfilepath,"stdout")) fclose(f);   // if not stdout then close
}

int
create_op_store(size_t exponent)
{
  if (exponent>=64)
  {
    fprintf(stderr,"You have chosen an op store size of 2^64 or larger\n");
    fprintf(stderr,"Too large... exiting\n");
    exit(0);
  }
  memset(&op_store, 0, sizeof(op_store));
  for (int op = 0; op < INVALID_OP; op++)
  {
    if (strlen(op_names_verbose[op]) > op_store.longest)
    {
      op_store.longest = strlen(op_names_verbose[op]);
    }
    op_store.oprecords_total_created[op] = 0;
    op_store.oprecords_stored_now[op] = 0;
  }
  op_store.exponent = exponent;
  op_store.size = (uint64_t)1 << op_store.exponent;
  op_store.hash_table = alloc_hash(exponent);
  op_store.hash_table->record_size = sizeof(struct larc_op_t);
  return 1;
}

size_t get_op_store_exp(void) {
  return op_store.hash_table->exponent;
}

uint64_t
hash_from_op(uint64_t in1_pID, uint64_t in2_pID, op_type_t op_type) {
  uint64_t hash = recursive_hash_from_three_integers(in1_pID, 
                                          (uint64_t)op_type, 
                                          in2_pID, op_store.exponent);
  return hash;
}

int op_set(op_type_t op_type, uint64_t in1_pID, uint64_t in2_pID,
           uint64_t out_pID, uint64_t hash)
{
// this function should be called whenever 
// a) op_get has returned zero;
// b) the operation has been performed and we have an answer
// It adds the operation to the appropriate hash chain, creating it if necessary
  if (op_type >= INVALID_OP) {
    return EINVAL;
  }

  op_store.numsets++;

  struct larc_op_t *op_rec_ptr = calloc(1,sizeof(struct larc_op_t));
  if (op_rec_ptr == NULL) { ALLOCFAIL(); }
  
  op_rec_ptr->op_type = op_type;
  op_rec_ptr->in1_pID = in1_pID;
  op_rec_ptr->in2_pID = in2_pID;
  op_rec_ptr->out_pID = out_pID;
  op_store.oprecords_total_created[op_type]++;
  op_store.oprecords_stored_now[op_type]++;

#ifndef HASH_CHAIN_GROWS_AT_HEAD
  // INSERT new items at tail of hash collision chain
  hash_insert_at_tail(op_store.hash_table, (record_ptr_t) op_rec_ptr,
                      MATRIX_ID_INVALID, hash);
#endif
#ifdef HASH_CHAIN_GROWS_AT_HEAD
  // INSERT new items at head of hash collision chain
  hash_insert_at_head(op_store.hash_table, (record_ptr_t) op_rec_ptr,
                      MATRIX_ID_INVALID, hash);
#endif

  return 0;
}

int64_t op_get(op_type_t op_type, uint64_t in1_pID,
                 uint64_t in2_pID, uint64_t hash)
{
  // if this operation has been done before,  the
  // matrixID of the result will be returned;
  // otherwise we get MATRIX_ID_INVALID

  int64_t ret_ID = MATRIX_ID_INVALID;

  if (op_type >= INVALID_OP) return ret_ID;

  // this function returns a pointer to the op chain
  // which for this op contains these two inputs

  int verbose = 0;

  op_ptr_t op_rec_ptr = NULL;

  if (verbose) printf("In %s, hash is %" PRIu64 "\n",__func__,hash);

#ifdef HASHSTATS
  (op_store.hash_table->num_accesses[hash])++;
#endif

  // find hash chain for this op containing this hash value,
  // if it exists

  if (verbose) printf("In %s, about to call hash_get_chain with op %d\n",
                      __func__,op_type);

  hash_node_t *node_ptr = hash_get_chain(op_store.hash_table, hash);

  if (verbose) printf("In %s, returned from hash_get_chain\n",__func__);

  //hash_node_t *last_zero_hit = NULL;
  int len = 0;

  if (verbose) printf("In %s, about to walk through chain\n",__func__);

#ifdef HASHSTATS
  (op_store.hash_table->num_accesses[hash])++;
#endif

  // walk through hash chain, checking to see if the inputs
  // are present; if found, return the pointer to the larc_op_t
  // location
  while (node_ptr)
  {
    len++;
    op_rec_ptr = (op_ptr_t) (node_ptr->record_ptr);

    // take opportunity to clean obsolete nodes
    if ( record_is_invalid(get_recordPTR_from_pID(op_rec_ptr->in1_pID,"first",__func__,0))
         || record_is_invalid(get_recordPTR_from_pID(op_rec_ptr->in2_pID,"second",__func__,0))
         || record_is_invalid(get_recordPTR_from_pID(op_rec_ptr->out_pID,"third",__func__,0))
          )
    {
      // if you are going to remove the node, then you need to know what the
      // previous node was  
      hash_node_t *next_ptr = node_ptr->next;
      // remove the node directly, stitching up the chain
      hash_node_remove_node(op_store.hash_table,node_ptr,hash);
      // keep statistics of nodes that have been removed
      if (op_rec_ptr->op_type!=INVALID_OP)
          op_store.oprecords_stored_now[op_rec_ptr->op_type]--;

      //FREE RECORD FOR NODE
      free(op_rec_ptr);
      // change node pointer to the next node and start from top of loop
      node_ptr = next_ptr;
      continue;
    } // if (record contains invalid matrix pointer)

    // node has passed validity checks
    if (in1_pID  == op_rec_ptr->in1_pID
          && in2_pID == op_rec_ptr->in2_pID
          && op_type == op_rec_ptr->op_type)
    {
      op_store.hash_table->hits++;
      if (!node_ptr->hits_maxxed)
      {
        node_ptr->record_hits++;
        if (node_ptr->record_hits == ((uint32_t) (-1)))
              node_ptr->hits_maxxed = 1;
      }
      op_store.op_store_hits++;
      ret_ID = op_rec_ptr->out_pID;
      // node contains correct record, so stop traversing/searching chain
      break;
    } // if (node contains correct record)

    //  This was going to be used so that we could truncate the 
    //  end of the hash chain which would contain all nodes with zero hits
    //if (node_ptr->hits == 0) last_zero_hit = node_ptr;

    // We move to the next record in the chain if the record checked
    // is not the record we wanted.
    node_ptr = node_ptr->next;

  } // END while (node_ptr)

  if (node_ptr == NULL)
  {
    if (verbose) printf("In %s, got through chain\n",__func__);

    // at this point we know the op is not in the op_store
    op_store.hash_table->misses++;
    op_rec_ptr = NULL;
    op_store.op_store_misses++;
  }

  return ret_ID;
}


op_type_t get_op_type_from_string_name(char *op_name)
{
  op_type_t op_type;

  // default value, is we fail (INVALID_OP)
  for (op_type=0; op_type<INVALID_OP; ++op_type) {
    if (!strcmp(op_name,op_names[op_type])) break;
  }
  return (op_type);
}

/****************************************************************************
*                   op_hash_chain_info_to_file                              *
*  This function prints information about op_store hash chain corresponding *
*  to the given the hash value for that chain.                              *
*  To get the appropriate hash value for the argument of this function      *
*     hashID = hash_from_op(in1_pID, in2_pID, op_name);             *
*     op_type can be seen by ...                                            *
*  The output file path and a user comment to be printed in the file are    *
*  also arguments.                                                          * 
****************************************************************************/
int
op_hash_chain_info_to_file(uint64_t hash, char *outfilepath, char *comment)
{
  FILE *f; 
  if (strcmp(outfilepath,"stdout"))  {
    printf("printing op store hash chain with hash value %" PRIu64 " to file %s\n",
	   hash,outfilepath);  
     f = fopen(outfilepath, "w"); 
  }
  else {
    printf("printing op store hash chain with hash value %" PRIu64 " to screen\n",hash);  
    f = stdout;
  }
  
  int ret = 1;
  fprintf(f,"Comment: %s\n",comment);
  fprintf(f,"This is the op_store hash chain info for hash value %" PRIu64 "\n", hash);
  fprintf(f,"\n\n========================(Op Hash Chain)=============================\n");
  fprintf(f,"operation    matrixID In_1 matrixID In_2  matrixID Out\n");

  // get pointer to the head of hash chain from hash table
  hash_table_t *table_ptr = op_store.hash_table;
  hash_node_t *node_ptr = table_ptr->heads[hash];
  
  op_ptr_t record_ptr;
  
  while (node_ptr)  { 
    record_ptr = (op_ptr_t) node_ptr->record_ptr; 
    
    // HEADER: "operation matrixID In_1 matrixID In_2  matrixID Out"
    uint64_t in1_pID = record_ptr->in1_pID;
    uint64_t in2_pID = record_ptr->in2_pID;
    uint64_t out_pID = record_ptr->out_pID;
    op_type_t op_type = record_ptr->op_type;
    fprintf(f,"%-12s %13" PRIu64 " %13" PRIu64 " %13" PRIu64 " \n",
            op_names[op_type], MID_FROM_PID(in1_pID), MID_FROM_PID(in2_pID),
            MID_FROM_PID(out_pID));

    node_ptr = node_ptr->next;
  }
  
  fprintf(f,"\n======================================================================\n");
  
  if  (strcmp(outfilepath,"stdout")) fclose(f);   // if not stdout then close

  return ret;
}


/*************************************************************************** 
 *                   op_hash_chain_info_to_screen                            *
 *  Calls op_hash_chain_info_to_file with stdout.                            *
 ****************************************************************************/
int 
op_hash_chain_info_to_screen(uint64_t hash, char *comment) 
{ 
  char *file_name = "stdout";  //print to screen
  return  op_hash_chain_info_to_file(hash, file_name, comment); 
  
}  


/*************************************************************************** 
 *                   base_op_hash_clean                        *
 *   Cleans a hash chain in the op store for a given hash and op_type      *
 *   by removing nodes with invalid matrices.                               *
 ***************************************************************************/
/*!
 * \ingroup larc
 * \brief Cleans a hash chain in the op store for a given hash and op_type
 * \param hashID Identifies the hash chain
 * \return 1
 */
static int
base_op_hash_clean(uint64_t hashID)
{
	hash_table_t *table_ptr = op_store.hash_table;
	hash_node_t *node_ptr = table_ptr->heads[hashID];
	op_ptr_t current_op_ptr;
	
	// walk through the hash chain, removing nodes along the way if any of the 3 matrices are invalid
	while (node_ptr)  {
		current_op_ptr = (op_ptr_t) node_ptr->record_ptr;
		struct larc_op_t *op_rec_ptr = (struct larc_op_t *) current_op_ptr;
        if (
	    record_is_invalid(get_recordPTR_from_pID(op_rec_ptr->in1_pID,"first",__func__,0))
	    || record_is_invalid(get_recordPTR_from_pID(op_rec_ptr->in2_pID,"second",__func__,0))
	    || record_is_invalid(get_recordPTR_from_pID(op_rec_ptr->out_pID,"third",__func__,0))
	    )
	    {
               hash_node_remove(op_store.hash_table,
                    (record_ptr_t)current_op_ptr, MATRIX_ID_INVALID,hashID);
               // keep statistics of nodes that have been removed 
               if (op_rec_ptr->op_type!=INVALID_OP)
                  op_store.oprecords_stored_now[op_rec_ptr->op_type]--;
                        
               // FREE RECORD FOR NODE;
			free(op_rec_ptr);
            }

            //  TODO: Need to reset last_zero_hit after cleaning 

            node_ptr = node_ptr->next;
        }
	return(1);
}


/*************************************************************************** 
 *                   clean_op_hash_chain                                   *
 *   Cleans a hash chain in the op store for a given hash and op_name      *
 *   Calls base_op_hash_chain with checks.                                 *
 ***************************************************************************/
int
clean_op_hash_chain(uint64_t hashID)
{
	uint64_t max = op_store.size;
	if ((hashID < 0) || (hashID >= max)) {
		fprintf(stderr,"Error: hash value out of range\n");
		return(0);
		}

	base_op_hash_clean(hashID);
	
	return(1);
}



/*************************************************************************** 
 *                   repair_op_store                                        *
 *   cleans/repairs all hash chains in the op store                                *
 ***************************************************************************/
int
clean_op_store()
{
    uint64_t max = op_store.size;
    uint64_t hashID;
    for (hashID = 0; hashID < max; hashID ++){
        base_op_hash_clean(hashID);
        // or 
        //clean_op_hash_chain(hashID);
    }
    return 1;
}


/*************************************************************************** 
 *                   clean_op_store                                        *
 *   Empty all hash chains in the op store                                *
 ***************************************************************************/
int empty_op_store()
{
    uint64_t max = op_store.size;
    uint64_t hashID;
    hash_table_t *table_ptr = op_store.hash_table;
    for (hashID = 0; hashID < max; hashID ++){
	hash_node_t *node_ptr = table_ptr->heads[hashID];
	hash_node_t *next_node_ptr;
	op_ptr_t current_op_ptr;
	
	// walk through the hash chain, removing nodes along the way 
	while (node_ptr) {
            // save the next node - alternatively we could always set this to
            // the head since next get's set to the head during hash_node_remove
            next_node_ptr = node_ptr->next;
            current_op_ptr = (op_ptr_t) node_ptr->record_ptr;
            struct larc_op_t *op_rec_ptr = (struct larc_op_t *) current_op_ptr;
            // remove the node that stores this record in the store
            // but does not free the record itself
            hash_node_remove(op_store.hash_table, (record_ptr_t)current_op_ptr,
                             MATRIX_ID_INVALID, hashID);
            // keep statistics of nodes that have been removed 
            if (op_rec_ptr->op_type!=INVALID_OP)
                 op_store.oprecords_stored_now[op_rec_ptr->op_type]--;


            // FREE RECORD FOR NODE;
            free(op_rec_ptr);

            node_ptr = next_node_ptr;
        }
    }

    // reset counters in op store
    op_store.numsets = 0;
    op_store.op_store_hits = 0;
    op_store.op_store_misses = 0;
    memset(op_store.oprecords_total_created, 0x00, (INVALID_OP + 1)*sizeof(uint64_t));
    memset(op_store.oprecords_stored_now,    0x00, (INVALID_OP + 1)*sizeof(uint64_t));

    return 1;
}

#ifdef HASHSTATS
/************************************************************************
 *  This function calls hashstats_to_file with the op_store             *
 *  hash table creating files with the following                        *
 *  hash table statistics:                                              *
 *    (1) accesses[hash_val] = An array of the number of accesses each  *
 *                             hash bucket (hash_val) has so far. This  *
 *                             is incremented by both hits and misses.  *
 *    (2) nodes[hash_val] = An array of the number of nodes (records)   *
 *                          currently in each hash chain.  These values *
 *                          are incremented by inserts, and decremented *
 *                          by removes.                                 *
 *    (3) And a standard hash_report() is sent to stdout or a file.     *
 *                                                                      *
 ************************************************************************/    
void op_hashstats(char *accesses_file,  char *nodes_file, char *report_file)
{
  hashstats_to_files(op_store.hash_table, accesses_file, nodes_file, report_file);
}
#endif

void free_op_store(void)
{
    // free nodes from op_store hash table, along with the matrix record
    // pointed to by each node. This routine should only be called as
    // part of a LARC shutdown procedure.

    // the op store struct only allocates the hash table. Each hash node
    // allocates a record of type larc_op_t. There are no other allocations.

    fprintf(stderr,"freeing memory for operationsStore structure\n");

    uint64_t hashID;
    hash_table_t *table_ptr = op_store.hash_table;
    for (hashID = 0; hashID < op_store.size; hashID++)
    {
        hash_node_t *node_ptr = table_ptr->heads[hashID];
        while (node_ptr)
        {
           hash_node_t *next_node_ptr = node_ptr->next;
           free((op_ptr_t)(node_ptr->record_ptr));
           free(node_ptr);
           node_ptr = next_node_ptr;
        }
    }
    dealloc_hash(op_store.hash_table);
}
