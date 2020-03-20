//                  info_store.c
//   Matrix Information Store for Global Recursive Compression Scheme
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
#include "info_store.h"


/* Define the struct here as it doesn't need to be used elsewhere */
static struct info_store_t {
  uint64_t longest;         // length in characters of longest info name
  uint64_t exponent;        // info_store has hash table of size 2^exponent
  uint64_t size;            // info_store has hash table of size 2^exponent
  uint64_t numsets;         // number of times info_set is called on valid
  uint64_t info_store_hits;   // incremented in info_get when success
  uint64_t info_store_misses; // incremented in info_get when failed
  hash_table_t *hash_table; // pointer to hash table for info_store
  uint64_t inforecords_total_created[INVALID_INFO+1];
  uint64_t inforecords_stored_now[INVALID_INFO+1];
} store = {0};


/* The types of information to store (INFO is defined in info_store.h) */
#define X(a, b, c) b,
static char *info_names[] = { INFO };
#undef X

#define X(a, b, c) c,
static char *info_names_verbose[] = { INFO };
#undef X


// This function is to help users interface with info chain functions
void list_info_names() {
  printf("The info names understood by LARC are:\n");
  for (int i=0; i<=INVALID_INFO; ++i)
  {
    printf("    %s\n",info_names_verbose[i]);
  }
}

// This function returns the name associated with enum info_type value
char *return_info_name(uint64_t info_type) {
  if (info_type > INVALID_INFO) {return("");}
  return info_names[info_type];
}

void
info_store_report(char *outfilepath)
{
  FILE *f; 
  if (strcmp(outfilepath,"stdout")) {
    printf("Printing info store report to file %s\n", outfilepath);
    f = fopen(outfilepath, "a"); 
    
  }
  else {
    printf("Printing info store report to screen\n");
    f = stdout;
  }
  
  int64_t size;
  fprintf(f,"\n");
  fprintf(f,"Info store report:\n");
  size = hash_report(store.hash_table, f, 0);
  fprintf(f,"Statistics for each info type:\n");
  for (int i=0; i<INVALID_INFO; ++i)
  {
    fprintf(f,"%-12s: %12ld records created, %12ld current records\n",
       info_names[i], store.inforecords_total_created[i],
       store.inforecords_stored_now[i]);
  }
  fprintf(f,"Total Info Store size: %3.2fMB\n", size/(1024.0*1024.0));
  fflush(f);
  if  (strcmp(outfilepath,"stdout")) fclose(f);   // if not stdout then close
}

int
create_info_store(size_t exponent)
{
  memset(&store, 0, sizeof(store));
  for (int info = 0; info < INVALID_INFO; info++)
  {
    if (strlen(info_names_verbose[info]) > store.longest)
    {
      store.longest = strlen(info_names_verbose[info]);
    }
    store.inforecords_total_created[info] = 0;
    store.inforecords_stored_now[info] = 0;
  }
  store.exponent = exponent;
  store.size = 1L << store.exponent;
  store.hash_table = alloc_hash(exponent);
  store.hash_table->record_size = sizeof(struct larc_info_t);
  return 1;
}

/*
 * \ingroup larc
 * \brief Generates a hash value based on an input matrixID and an info_type_t value
 * \param matID A matrixID
 * \param info_type An enum value indicating the type of metadata stored
 * \return The hash value
 */
static uint64_t 
hash_from_info(int64_t matID, info_type_t info_type) {
  uint64_t hash = recursive_hash_from_two_integers((uint64_t)matID, 
                                          (uint64_t)info_type, 
                                          store.exponent);
  return hash;
}
  

/* Returns a pointer to the struct larc_info_t indicated by the two inputs 
 * If the two inputs don't match at their hash location, walk through
 * the hash chain structure and look until
 * the next available entry is found, or until the entire array has been
 * walked once (looping at the end) */
/* if create_flag is 1, then a new record will be formed */

struct larc_info_t *
info_find(info_type_t info_type, int64_t matID, int create_flag)
{
  int verbose = 0;

  // test to see that info_type is valid
  if (info_type == INVALID_INFO) {
    fprintf(stderr,"in %s, attempting to find and INVALID_INFO\n",__func__);
    exit(-1);
  }

  uint64_t hash = hash_from_info(matID, info_type);
  if (verbose) {
    printf("In %s, hash is %" PRIu64 "\n",__func__,hash);
  }
  
  // find hash chain for this info containing this hash value,
  // if it exists
  
  if (verbose) {
    printf("In %s, about to call hash_get_chain with info %d\n",__func__,info_type);
  }
  
  hash_node_t *node_ptr = hash_get_chain(store.hash_table, hash);
  
  if (verbose) {
    printf("In %s, returned from hash_get_chain\n",__func__);
  }
  
  //hash_node_t *last_zero_hit = NULL;
  int len = 0;
  
  if (verbose) {
    printf("In %s, about to walk through chain\n",__func__);
  }
  
  // walk through hash chain, checking to see if the inputs
  // are present; if found, return the pointer to the larc_info_t
  // location
  while (node_ptr)
    {
      len++;
      info_ptr_t current_info_ptr = (info_ptr_t) (node_ptr->record_ptr);
      struct larc_info_t *info_rec_ptr = (struct larc_info_t *) current_info_ptr;
      if (matID == info_rec_ptr->my_matID && info_type == info_rec_ptr->info_type)
	{
	  store.hash_table->hits++;
	  node_ptr->hits++;
	  return info_rec_ptr;
	}
      else {
	if (matrix_is_invalid(get_matPTR_from_matID(info_rec_ptr->my_matID,"my",__func__,0)))
	  {
	    hash_node_remove(store.hash_table, (record_ptr_t)current_info_ptr, hash);
	    //FREE RECORD FOR NODE;
	    free(info_rec_ptr);
	  }
       }

      //  This was going to be used so that we could truncate the 
      //  end of the hash chain which would contain all nodes with zero hits
      //if (node_ptr->hits == 0)
	//last_zero_hit = node_ptr;
      node_ptr = node_ptr->next;
    }
  
  if (verbose) {
    printf("In %s, got through chain\n",__func__);
  }
  
  // at this point we know the info is not in the info_store
  // we either create a new larc_info_t, insert it in an old
  // hash chain (or create a new hash chain for it), and
  // return it to the calling routine, or if create_flag is
  // zero we return zero.
  store.hash_table->misses++;
  if (create_flag)
    {
      struct larc_info_t *info_rec_ptr = calloc(1, sizeof(struct larc_info_t));
      // setup as invalid
      info_rec_ptr->info_type = INVALID_INFO;
      info_rec_ptr->my_matID = -1;
      strcpy(info_rec_ptr->info_data,"");

#ifndef HASH_CHAIN_GROWS_AT_HEAD
      // INSERT new items at tail of hash collision chain
      hash_insert_at_tail(store.hash_table, (record_ptr_t) info_rec_ptr, hash);
#endif
#ifdef HASH_CHAIN_GROWS_AT_HEAD
      // INSERT new items at head of hash collision chain
      hash_insert_at_head(store.hash_table, (record_ptr_t) info_rec_ptr, hash);
#endif
      
      return info_rec_ptr;
    }
  return NULL;
}


// This function will create a new larc_info_t record and fills
// it with the variables passed to this function.
// Panic if the record exists and you are trying to save conflicting
// info_type or info_data

int
info_set(info_type_t info_type, int64_t my_matID, const char *info_data)
{
// this function should be called whenever 
// a) info_get has failed to find anything
// b) TODO  .... the info 
//     add the info to the appropriate hash chain, creating
//     the record if necessary
	if (info_type >= INVALID_INFO) {
	  fprintf(stderr,"Error in %s, invalid info_type\n",__func__);
	  return EINVAL;
	}

	store.numsets++;

        // the "1" instructs info_find to create a new larc_info_t if needed
	struct larc_info_t *info_rec_ptr = info_find(info_type, my_matID, 1);

        // if a new larc_info_t has been created, need to initialize it?
	if (info_rec_ptr->my_matID != my_matID 
            || info_rec_ptr->info_type != info_type) { 
	  info_rec_ptr->my_matID = my_matID;
	  if (strlen(info_data) <= 256) {
            strcpy(info_rec_ptr->info_data,info_data);
	  }
          else {
            fprintf(stderr,"string is too long in %s\n",__func__); 
            exit(1);
          }
          info_rec_ptr->info_type = info_type;
          store.inforecords_total_created[info_type]++;
          store.inforecords_stored_now[info_type]++;
	}
        else { // record already exists
	  if (!strcmp(info_rec_ptr->info_data,info_data)) {
            fprintf(stderr,"Panic in %s because info_data %s is changing to %s\n",
         	   __func__, info_rec_ptr->info_data,info_data);
	  }
	}
	return 0;
}

// This function searches for an information record associated with 
// the matrix specified by my_matID and information type info_type
// and returns the string of meta information or the NULL string if
// nothing is found.
char *
info_get(info_type_t info_type, int64_t my_matID)
{

  char *ret_string = calloc(256,sizeof(char));
  
  if (info_type >= INVALID_INFO) {
    return ret_string;
  }
  
  // this function returns a pointer to the info record
  struct larc_info_t *info_rec_ptr = info_find(info_type, my_matID, 0);
  
  if (info_rec_ptr != NULL) {
    store.info_store_hits++;
    strcpy(ret_string,info_rec_ptr->info_data);
  } else {
    store.info_store_misses++;
  }
  
  return ret_string;
}


/* translates a string into the appropriate enum value */
info_type_t get_info_type_from_string_name(const char *info_name)
{
  info_type_t info_type;

  // default value, is we fail (INVALID_INFO)
  for (info_type=0; info_type<INVALID_INFO; ++info_type) {
    if (!strcmp(info_name,info_names[info_type])) break;
  }
  return (info_type);
}



/*****************************************************************
 *               info_hashID_by_matrixIDs                            *
 *  If this function succeeds it will return the hash value      *
 *  associated with hashing the info record
 *         hash_from_info(in1_ptr, in2_ptr, info_type);              *
 *  where the in1_ptr is the matrix pointer associated with      *
 *  the matrix matrixID in1_mID, etc.                       *
 *  This function returns either a -1 if it fails because:       *
 *    - one of the matrixIDs is out of range, or            *
 *    - the matrix associated with a matrixID has been      *
 *      removed from the matrix store                            *
 *  Normally a hash is a uint64_t, but if there is a fail,       *
 *  the hashID can be -1, so the return value is int64_t.        *
 *  The user can call the function list_info_names()               *
 *****************************************************************/
int64_t info_hashID_by_matrixIDs(int64_t my_matID, const char *info_name)
{

  // make sure that this matrixID corresponds to a matrix which still exists
  mat_ptr_t my_matPTR  = get_matPTR_from_matID(my_matID, "my matID", __func__,0);

  // check to see if this matrix has already been removed from the matrix store
  if (my_matPTR == MATRIX_PTR_INVALID) {
    return(-1);
  }

  // Getting the info_type_t from the info_name
  info_type_t info_type = get_info_type_from_string_name(info_name);
  if (info_type == INVALID_INFO)
  {
    fprintf(stderr,"In %s: argument info_name\n\t%s\n",__func__,info_name);
    fprintf(stderr," is not a valid info record\n");
    exit(-1);
  }

  // calculate hash value for the operation of info_type
  uint64_t hash = hash_from_info(my_matID, info_type); 

  return (int64_t)hash;
}



/****************************************************************************
*                   info_hash_chain_info_to_file                              *
*  This function prints information about info_store hash chain corresponding *
*  to the given the hash value for that chain.                              *
*  To get the appropriate hash value for the argument of this function      *
*     hashID = info_hashID_by_matrixIDs(in1_mID, in2_mID, info_name);             *
*     info_type can be seen by ...                                            *
*  The output file path and a user comment to be printed in the file are    *
*  also arguments.                                                          * 
****************************************************************************/
int
info_hash_chain_info_to_file(uint64_t hash, char *outfilepath, char *comment)
{
  FILE *f; 
  if (strcmp(outfilepath,"stdout"))  {
    printf("printing info store hash chain with hash value %ld to file %s\n",
	   hash,outfilepath);  
     f = fopen(outfilepath, "w"); 
  }
  else {
    printf("printing info store hash chain with hash value %ld to screen\n",hash);  
    f = stdout;
  }
  
  int ret = 1;
  fprintf(f,"Comment: %s\n",comment);
  fprintf(f,"This is the info_store hash chain info for hash value %ld\n", hash);
  fprintf(f,"\n\n========================(Info Hash Chain)=============================\n");
  fprintf(f,"info_name   my_matID   info_data\n");

  // get pointer to the head of hash chain from hash table
  hash_table_t *table_ptr = store.hash_table;
  hash_node_t *node_ptr = table_ptr->heads[hash];
  
  info_ptr_t record_ptr;
  
  while (node_ptr)  { 
    record_ptr = (info_ptr_t) node_ptr->record_ptr; 
    
    // HEADER: "info_name  my_matID   info_data\n"
    info_type_t info_type = record_ptr->info_type;
    uint64_t my_matID = record_ptr->my_matID;
    char *info_data = record_ptr->info_data;
    fprintf(f,"%-12s %13lu %s\n",info_names[info_type],my_matID,info_data);

    node_ptr = node_ptr->next;
  }
  
  fprintf(f,"\n======================================================================\n");
  
  if  (strcmp(outfilepath,"stdout")) fclose(f);   // if not stdout then close

  return ret;
}


/*************************************************************************** 
 *                   info_hash_chain_info_to_screen                            *
 *  Calls info_hash_chain_info_to_file with stdout.                            *
 ****************************************************************************/
int 
info_hash_chain_info_to_screen(uint64_t hash, char *comment) 
{ 
  char *file_name = "stdout";  //print to screen
  return  info_hash_chain_info_to_file(hash, file_name, comment); 
  
}  



/*************************************************************************** 
 *                   base_info_hash_clean                        *
 *   Cleans a hash chain in the info store                       *
 ***************************************************************************/
/*
 * \ingroup larc
 * \brief Cleans a hash chain in the info store
 * \param hashID The hash specifing the hash chain to be cleaned
 */
static int
base_info_hash_clean(uint64_t hashID)
{
	hash_table_t *table_ptr = store.hash_table;
	hash_node_t *node_ptr = table_ptr->heads[hashID];
	info_ptr_t current_info_ptr;
	
	// walk throught the hash chain, removing nodes along the way if any of the 3 matrices are invalid
	while (node_ptr)  {
		current_info_ptr = (info_ptr_t) node_ptr->record_ptr;
		struct larc_info_t *info_rec_ptr = (struct larc_info_t *) current_info_ptr;
        // if any of the matrices have been deleted from the store, panic
        if (matrix_is_invalid(get_matPTR_from_matID(info_rec_ptr->my_matID,"my_matID",__func__,0)))
	    {
			hash_node_remove(store.hash_table, (record_ptr_t)current_info_ptr, hashID);
                        // keep statistics of nodes that have been removed 
                        if (info_rec_ptr->info_type!=INVALID_INFO)
                          store.inforecords_stored_now[info_rec_ptr->info_type]--;
                        
			// FREE RECORD FOR NODE;
			free(info_rec_ptr);
			}
		//  TODO: Need to reset last_zero_hit after cleaning 

		node_ptr = node_ptr->next;
		}
	return(1);
}


/*************************************************************************** 
 *                   clean_info_hash_chain                                   *
 *   Cleans a hash chain in the info store for a given hash *
 ***************************************************************************/
int
clean_info_hash_chain(uint64_t hashID)
{
	uint64_t max = store.size;
	if ((hashID < 0) || (hashID >= max)) {
		fprintf(stderr,"Error: hash value out of range\n");
		return(0);
		}

	base_info_hash_clean(hashID);
	
	return(1);
}

