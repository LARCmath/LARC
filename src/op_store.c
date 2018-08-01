//                  op_store.c
//   Matrix Operations Stores for LARC
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <errno.h>

#include "larc.h"
#include "hash.h"
#include "matrix_store.h"
#include "op_store.h"

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
} store = {0};


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
  size = hash_report(store.hash_table, f, 0);
  fprintf(f,"Statistics for each operation:\n");
  for (int i=0; i<INVALID_OP; ++i)
  {
    fprintf(f,"%-12s: %12ld records created, %12ld current records\n",
       op_names[i], store.oprecords_total_created[i],
       store.oprecords_stored_now[i]);
  }
  fprintf(f,"Total Op Store size: %3.2fMB\n", size/(1024.0*1024.0));
  fflush(f);
  if  (strcmp(outfilepath,"stdout")) fclose(f);   // if not stdout then close
}

int
create_op_store(size_t exponent)
{
  memset(&store, 0, sizeof(store));
  for (int op = 0; op < INVALID_OP; op++)
  {
    if (strlen(op_names_verbose[op]) > store.longest)
    {
      store.longest = strlen(op_names_verbose[op]);
    }
    store.oprecords_total_created[op] = 0;
    store.oprecords_stored_now[op] = 0;
  }
  store.exponent = exponent;
  store.size = 1 << store.exponent;
  store.hash_table = alloc_hash(exponent);
  store.hash_table->record_size = sizeof(struct larc_op_t);
  return 1;
}


static uint64_t 
hash_from_op(mat_add_t in1_mat_ptr, mat_add_t in2_mat_ptr, op_type_t op_type) {
  uint64_t hash = recursive_hash_from_three_integers((uint64_t)in1_mat_ptr, 
                                          (uint64_t)op_type, 
                                          (uint64_t)in2_mat_ptr, 
                                          store.exponent);
  return hash;
}
  

/* Returns a pointer to the struct larc_op_t indicated by the two inputs 
 * If the two inputs don't match at their hash location, walk through
 * the hash chain structure and look until
 * the next available entry is found, or until the entire array has been
 * walked once (looping at the end) */
struct larc_op_t *
op_find(op_type_t op_type, mat_add_t in1_mat_ptr, mat_add_t in2_mat_ptr, int create_flag)
{
  int verbose = 0;

  // test to see that op_type is valid
  if (op_type == INVALID_OP) {
    printf("in %s, attempting to find and INVALID_OP\n",__func__);
    exit(-1);
  }

  uint64_t hash = hash_from_op(in1_mat_ptr, in2_mat_ptr, op_type);
  if (verbose) {
    printf("In %s, hash is %" PRIu64 "\n",__func__,hash);
  }

  // NOTE: misses and hits are incremented in places that call op_find
  //       it would probably be better to 
#ifdef HASHSTATS
  if (create_flag == 0) {
    (store.hash_table->num_accesses[hash])++;
  }
#endif
  
  // find hash chain for this op containing this hash value,
  // if it exists
  
  if (verbose) {
    printf("In %s, about to call hash_get_chain with op %d\n",__func__,op_type);
  }
  
  hash_node_t *node_ptr = hash_get_chain(store.hash_table, hash);
  
  if (verbose) {
    printf("In %s, returned from hash_get_chain\n",__func__);
  }
  
  hash_node_t *last_zero_hit = NULL;
  int len = 0;
  
  if (verbose) {
    printf("In %s, about to walk through chain\n",__func__);
  }
  
#ifdef HASHSTATS
  (store.hash_table->num_accesses[hash])++;
#endif


  // walk through hash chain, checking to see if the inputs
  // are present; if found, return the pointer to the larc_op_t
  // location
  while (node_ptr)
    {
      len++;
      op_add_t current_op_ptr = (op_add_t) (node_ptr->record_ptr);
      struct larc_op_t *op_rec_ptr = (struct larc_op_t *) current_op_ptr;
      if (get_matrixID_from_ptr(in1_mat_ptr) == op_rec_ptr->in1_matID 
	  && get_matrixID_from_ptr(in2_mat_ptr) == op_rec_ptr->in2_matID
          && op_type == op_rec_ptr->op_type)
	{
	  store.hash_table->hits++;
	  node_ptr->hits++;
	  return op_rec_ptr;
	}
      else {
	if (
	    matrix_is_invalid(mat_ptr_from_matrixID(op_rec_ptr->in1_matID,"first",__func__,0))
	    || matrix_is_invalid(mat_ptr_from_matrixID(op_rec_ptr->in2_matID,"second",__func__,0))
	    || matrix_is_invalid(mat_ptr_from_matrixID(op_rec_ptr->out_matID,"third",__func__,0))
	    )
	  {
	    hash_node_remove(store.hash_table, (record_ptr_t)current_op_ptr, hash);
	    //FREE RECORD FOR NODE;
        free(op_rec_ptr);
        }
       }

      //  This was going to be used so that we could truncate the 
      //  end of the hash chain which would contain all nodes with zero hits
      if (node_ptr->hits == 0)
	last_zero_hit = node_ptr;
      
      node_ptr = node_ptr->next;
    }
  
  if (verbose) {
    printf("In %s, got through chain\n",__func__);
  }
  
  // at this point we know the op is not in the op_store
  // we either create a new larc_op_t, insert it in an old
  // hash chain (or create a new hash chain for it), and
  // return it to the calling routine, or if create_flag is
  // zero we return zero.
  store.hash_table->misses++;
  if (create_flag)
    {
      struct larc_op_t *op_rec_ptr = calloc(1, sizeof(struct larc_op_t));
      // setup as invalid
      op_rec_ptr->op_type = INVALID_OP;
      op_rec_ptr->in1_matID=op_rec_ptr->in2_matID=op_rec_ptr->out_matID = MATRIX_ID_INVALID;

#ifndef HASH_CHAIN_GROWS_AT_HEAD
      // INSERT new items at tail of hash collision chain
      hash_insert_at_tail(store.hash_table, (record_ptr_t) op_rec_ptr, hash);
#endif
#ifdef HASH_CHAIN_GROWS_AT_HEAD
      // INSERT new items at head of hash collision chain
      hash_insert_at_head(store.hash_table, (record_ptr_t) op_rec_ptr, hash);
#endif
      
      return op_rec_ptr;
    }
  return NULL;
}


int
op_set(op_type_t op_type, mat_add_t in1_mat_ptr, mat_add_t in2_mat_ptr, mat_add_t out_mat_ptr)
{
// this function should be called whenever 
// a) op_get has returned zero;
// b) the operation has been performed and we have an answer
// It adds the operation to the appropriate hash chain, creating
// it if necessary
	if (op_type >= INVALID_OP) {
		return EINVAL;
	}

	store.numsets++;

        // the "1" instructs op_find to create a new larc_op_t if needed
	struct larc_op_t *op_rec_ptr = op_find(op_type, in1_mat_ptr, in2_mat_ptr, 1);

        // if a new larc_op_t has been created, need to initialize it?
	if (op_rec_ptr->in1_matID != get_matrixID_from_ptr(in1_mat_ptr) 
	    || op_rec_ptr->in2_matID != get_matrixID_from_ptr(in2_mat_ptr)
            || op_rec_ptr->op_type == INVALID_OP) { 
	  op_rec_ptr->in1_matID = get_matrixID_from_ptr(in1_mat_ptr);
	  op_rec_ptr->in2_matID = get_matrixID_from_ptr(in2_mat_ptr);
          op_rec_ptr->op_type = op_type;
          store.oprecords_total_created[op_type]++;
          store.oprecords_stored_now[op_type]++;
	}
	op_rec_ptr->out_matID = get_matrixID_from_ptr(out_mat_ptr);

	return 0;
}

mat_add_t
op_get(op_type_t op_type, mat_add_t in1_mat_ptr, mat_add_t in2_mat_ptr)
{
        // if this operation has been done before,  the
        // matrix_ptr of the result will be returned;
        // otherwise we get MATRIX_PTR_INVALID
	mat_add_t ret_mat_ptr = MATRIX_PTR_INVALID;

	if (op_type >= INVALID_OP) {
		return ret_mat_ptr;
	}

        // this function returns a pointer to the op chain
        // which for this op contains these two inputs
	struct larc_op_t *op_rec_ptr = op_find(op_type, in1_mat_ptr, in2_mat_ptr, 0);

	if (op_rec_ptr != NULL) {
	  store.op_store_hits++;
	  ret_mat_ptr = mat_ptr_from_matrixID(op_rec_ptr->out_matID, "first", __func__,0);
	} else {
	  store.op_store_misses++;
	}

	return ret_mat_ptr;
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


/*****************************************************************
 *               op_hashID_by_matrixIDs                            *
 *  If this function succeeds it will return the hash value      *
 *  associated with hashing the operation                        *
 *         hash_from_op(in1_ptr, in2_ptr, op_type);              *
 *  where the in1_ptr is the matrix pointer associated with      *
 *  the matrix matrixID in1_mID, etc.                       *
 *  This function returns either a -1 if it fails because:       *
 *    - one of the matrixIDs is out of range, or            *
 *    - the matrix associated with a matrixID has been      *
 *      removed from the matrix store                            *
 *  Normally a hash is a uint64_t, but if there is a fail,       *
 *  the hashID can be -1, so the return value is int64_t.        *
 *  The user can call the function list_op_names()               *
 *****************************************************************/
int64_t op_hashID_by_matrixIDs(int64_t in1_mID, int64_t in2_mID, char *op_name)
{
  // look up the matrix pointer corresponding to that matrixID
  mat_add_t in1_ptr  = mat_ptr_from_matrixID(in1_mID, "first", __func__,0);
  mat_add_t in2_ptr  = mat_ptr_from_matrixID(in2_mID, "second", __func__,0);

  // check to see if this matrix has already been removed from the matrix store
  if (in1_ptr == MATRIX_PTR_INVALID || in2_ptr == MATRIX_PTR_INVALID) {
    return(-1);
  }

  // Getting the op_type_t from the op_name
  op_type_t op_type = get_op_type_from_string_name(op_name);
  if (op_type == INVALID_OP)
  {
    printf("In %s: argument op_name\n\t%s\n",__func__,op_name);
    printf(" is not a valid operation\n");
    exit(-1);
  }

  // calculate hash value for the operation of op_type
  uint64_t hash = hash_from_op(in1_ptr, in2_ptr, op_type); 

  return (int64_t)hash;
}


/****************************************************************************
*                   op_hash_chain_info_to_file                              *
*  This function prints information about op_store hash chain corresponding *
*  to the given the hash value for that chain.                              *
*  To get the appropriate hash value for the argument of this function      *
*     hashID = op_hashID_by_matrixIDs(in1_mID, in2_mID, op_name);             *
*     op_type can be seen by ...                                            *
*  The output file path and a user comment to be printed in the file are    *
*  also arguments.                                                          * 
****************************************************************************/
int
op_hash_chain_info_to_file(uint64_t hash, char *outfilepath, char *comment)
{
  FILE *f; 
  if (strcmp(outfilepath,"stdout"))  {
    printf("printing op store hash chain with hash value %ld to file %s\n",
	   hash,outfilepath);  
     f = fopen(outfilepath, "w"); 
  }
  else {
    printf("printing op store hash chain with hash value %ld to screen\n",hash);  
    f = stdout;
  }
  
  int ret = 1;
  fprintf(f,"Comment: %s\n",comment);
  fprintf(f,"This is the op_store hash chain info for hash value %ld\n", hash);
  fprintf(f,"\n\n========================(Op Hash Chain)=============================\n");
  fprintf(f,"operation    matrixID In_1 matrixID In_2  matrixID Out\n");

  // get pointer to the head of hash chain from hash table
  hash_table_t *table_ptr = store.hash_table;
  hash_node_t *node_ptr = table_ptr->heads[hash];
  
  op_add_t record_ptr;
  
  while (node_ptr)  { 
    record_ptr = (op_add_t) node_ptr->record_ptr; 
    
    // HEADER: "operation matrixID In_1 matrixID In_2  matrixID Out"
    uint64_t in1_mID = record_ptr->in1_matID;
    uint64_t in2_mID = record_ptr->in2_matID;
    uint64_t out_mID = record_ptr->out_matID;
    op_type_t op_type = record_ptr->op_type;
    fprintf(f,"%-12s %13lu %13lu %13lu \n",op_names[op_type],in1_mID, in2_mID, out_mID);

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
 ***************************************************************************/
static int
base_op_hash_clean(uint64_t hashID)
{
	hash_table_t *table_ptr = store.hash_table;
	hash_node_t *node_ptr = table_ptr->heads[hashID];
	op_add_t current_op_ptr;
	
	// walk throught the hash chain, removing nodes along the way if any of the 3 matrices are invalid
	while (node_ptr)  {
		current_op_ptr = (op_add_t) node_ptr->record_ptr;
		struct larc_op_t *op_rec_ptr = (struct larc_op_t *) current_op_ptr;
        if (
	    matrix_is_invalid(mat_ptr_from_matrixID(op_rec_ptr->in1_matID,"first",__func__,0))
	    || matrix_is_invalid(mat_ptr_from_matrixID(op_rec_ptr->in2_matID,"second",__func__,0))
	    || matrix_is_invalid(mat_ptr_from_matrixID(op_rec_ptr->out_matID,"third",__func__,0))
	    )
	    {
			hash_node_remove(store.hash_table, (record_ptr_t)current_op_ptr, hashID);
                        // keep statistics of nodes that have been removed 
                        if (op_rec_ptr->op_type!=INVALID_OP)
                          store.oprecords_stored_now[op_rec_ptr->op_type]--;
                        
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
 ***************************************************************************/
int
clean_op_hash_chain(uint64_t hashID)
{
	uint64_t max = store.size;
	if ((hashID < 0) || (hashID >= max)) {
		printf("Error: hash value out of range\n");
		return(0);
		}

	base_op_hash_clean(hashID);
	
	return(1);
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
  hashstats_to_files(store.hash_table, accesses_file, nodes_file, report_file);
}
#endif
