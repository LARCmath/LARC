//                                  hash.c 
//                        Hash functions for LARC
//                  for fast find_or_insert functions 
//     for the Matrix Store and the matrix operation stores                                        
/******************************************************************
 *                                                                *
 * Copyright (C) 2014-2024, Institute for Defense Analyses        *
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

// Standard libraries
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <string.h> //for memcpy

// Our header files structures and functions
#include "larc.h"
#include "scalars.h"
#include "hash.h"
#include "global.h"

/*!
 * \file hash.c
 * \brief This file contains the routines needed for LARC's hash tables
 */

/******************************************************************
 * BACKGROUND:
 *  The MatrixStore has a MatrixRecord for each unique matrix and each
 *  submatrix (down to 1 by 1 matrices).  When a new matrix is stored
 *  it is assigned a new MatrixRecord and a unique MatrixID.  MatrixIDs
 *  are assigned in the order matrices are created and never reused
 *  even if a matrix has been deleted from the store. (LARC also uses
 *  PackedIDs, which are MatrixIDs with additional information 'packed'
 *  into the unused bits; currently the only difference between the PackedID
 *  and the MatrixID is a bit indicating whether the MatrixID is for a
 *  scalar or non-scalar matrix. PackedIDs inherit uniqueness from MatrixIDs.)
 *
 *  The value of a 1 by 1 matrix is just the scalar value and is stored in its
 *  MatrixRecord.  The value of a larger matrix is given by a list of 
 *  the four MatrixIDs of that matrices quadrant submatrices; these
 *  determine the matrix recursively.
 *
 *  LARC has a special locality sensitive hash that we have created
 *  to handle numerical precision issues and to mimic symbolic
 *  computation (or mimic exact computation).  Our locality hash works by
 *  defining a division of space into small regions at compile time
 *  determined by parameters regionbitparam and zeroregionbitparam.
 *  LARC will store the first scalar it sees for any region, and 
 *  that will become the unique representative scalar for that region.
 *  This scalar is given a MatrixRecord and a MatrixID.
 *  When an attempt is made to store any new scalar from a region
 *  which already has a stored representative, LARC instead returns
 *  the MatrixID of the existing representative, with the assumption
 *  that the existing representative is "close enough". This avoids
 *  the situation where numerical precision can spawn many similar
 *  scalars  (which would have be identical if computations were 
 *  exact) and the derivative problem of many nearly identical matrices.
 *
 *  HOW HASHING COMES IN:
 *  LARC is able to quickly find whether there is an existing stored
 *  scalar in a region, because the hash functions below take any
 *  scalars which would be in the same region to the same
 *  hash bucket.  This is accomplished by defining our hash
 *  for scalar s to be H(R(s)) where R calculates the scalar value
 *  in the center of the region containing s, and H is a hash function.
 *  The hash in mult_hash_golden below is the multiplicative
 *  Fibonacci hash.  There are some variations
 *  which are handled for different scalarTypes listed below. For 
 *  example the function hash_longdouble converts all the bits of the
 *  mantissa, sign and exponent, into 64 integer bits, then hashes
 *  that. The basic idea is that if we have more than 64 bits of information
 *  in a scalarType (e.g. Complex has a double for 
 *  hash bucket in the matrix store
 *  The hashes provide a fast front end for storage and retrieval
 *  into the store.  The Links column of the Matrix Table keeps track of
 *  the next row in each hash collision chain, as well as the next
 *  free row in the free list of the Matrix Table.  The hash acts
 *  on the MatrixVals stored in the Matrix Table, which are either
 *  a pair of long doubles corresponding to the real and imaginary part
 *  of a complex scalar (if the matrix is 1 by 1), or which are
 *  a list of four matrix pointers corresponding to the four quadrant
 *  submatrices in the case where the matrix is larger.
 *  Thus the find_or_insert_matrix function takes as input a MatrixVal
 *  and returns the matrixID that has been found, or inserted into
 *  the matrix table.
 *  
 *  Similarly the matrix operations are memoized and the hash for
 *  these takes a pair of matrixIDs for 2 argument operations like
 *  Sum, Product, and Kronecker Product.  And takes one argument,
 *  (and a place saver) for single argument functions like Conjugation
 *  (which has an option to be stored in a table and an option to
 *  be stored in the matrix table), Conjugation, and Adjoint.
 *  Trace is handled differently since it's output is a scalarType
 *  number instead of a matrixID
 *******************************************************************/

/***************************************************************************
 * The function region_hash_from_scalarType() does not actually hash the
 * scalar given to it, but instead hashes the region center for the
 * region the scalar is in. This is then passed to the sca_hash() function,
 * which is compiled to different functions depending on the chosen
 * scalarType. Many of these functions are found in scalars.c, but a brief
 * outline is given here.
 *
 * ** INTEGER: do int_hash -> mult_golden_hash
 * ** REAL: do hash_from_one_longdouble
 * ** COMPLEX: do hash_from_two_longdoubles
 * ** MPINTEGER: do larc_mpz_hash -> calls mult_golden_hash on sign bit,
 *    loops overXORing each limb in then hashing
 * ** MPREAL: do larc_mpfr_hash -> tests for zero with mpfr_regular_p()
 *    if zero, returns mult_golden_hash of zero,
 *    otherwise, behaves as MPINTEGER
 * ** MPCOMPLEX: do larc_mpfr_hash on real and imaginary parts, pass results
 *    to recursive_hash_from_two_integers()
 * ** MPRATIONAL: do larc_mpz_hash on numerator, denominator, pass results to
 *    recursive_hash_from_two_integers()
 * ** MPRATCOMPLEX: do mprational hash on real and imaginary parts, pass
 *    results to recursive_hash_from_two_integers()
 */


// NON SCALAR
uint64_t
hash_from_matrix_subMatList(int64_t sub[4], matrix_type_t mat_type, uint64_t exponent)
{

#ifdef DEBUG_HASH_C
  printf("In routine %s\n",__func__);
  printf("     The hash exponent is %ld\n",exponent);
#endif // #ifdef DEBUG_HASH_C

  // EXCEPTION CHECKING
  if(mat_type==SCALAR) {   
    fprintf(stderr,"In %s mat_type was SCALAR\n",__func__);
    exit(1);
  }   

  // sub is 4 int64_t corresponding to the four submatrices
#ifdef DEBUG_HASH_C
  printf("     sub = [%ld , %ld , %ld , %ld]\n",
	 sub[0], sub[1], sub[2], sub[3]);
#endif // #ifdef DEBUG_HASH_C

    return recursive_hash_from_four_integers(
                (uint64_t)sub[0], (uint64_t)sub[1], 
                (uint64_t)sub[2], (uint64_t)sub[3], exponent);
}


#ifndef MAR  // SPRmode  
// SCALAR
/***************************************************************
 *    region_hash_from_scalarType
 *    Orig Author: Jenny Zito, Modified by Steve, Bill, Jenny, Andy ...
 *        Our hashing function is a locality sensitive hash.
 *        We accomplish this by breaking up space into regions.
 *        This function calls another to return the center point
 *        of the region containing the given input scalar.
 *        Then this central point is hashed.  In this way all
 *        scalars in the same region are hashed to the same
 *        value making it easy for LARC to find a previously
 *        saved scalar, that it can use as the representative.
 ****************************************************************/
uint64_t
region_hash_from_scalarType(scalarType scalar, uint64_t exponent)
{

  int debug = 0;
  
  if (debug) {
    printf("In routine %s\n",__func__);
    printf("     The hash exponent is %" PRIu64 "\n",exponent);
  }
  
  // Test to see if the scalar is NaN or INFINITE
  if (testForNaN(scalar))
  {
    fprintf(stderr,"ERROR: in %s, scalar being hashed is NaN.\n",__func__);
    fprintf(stderr,"LARC does not currently support this. Exiting...\n");
    exit(0);
  }
  if (testForInfinity(scalar))
  {
    fprintf(stderr,"ERROR: in %s, scalar being hashed is infinity.\n",__func__);
    fprintf(stderr,"LARC does not currently support this. Exiting...\n");
    exit(0);
  }

  // marker is the center of the region containing the scalar input
  // the region is used for the locality-sensitive hash to ensure
  // that only one point in any region is stored in the MatrixStore
  if (scratchVars.approx_value_in_use)
    fprintf(stderr,"%s reusing scratchVars.approx_value!\n",__func__);
  scratchVars.approx_value_in_use = 1;

  scalarType *marker = &scratchVars.approx_value;
  return_SPRregion_or_MARtile_label(marker, 0, scalar);
  
  if (debug) {
    char *scalar_string = sca_get_readable_approx_str(scalar);
    char *marker_string = sca_get_readable_approx_str(*marker);
    printf("Before rounding, scalar is %s, afterwards it is %s.\n", scalar_string, marker_string);
    free(marker_string);
    free(scalar_string);
    printf("calling sca_hash\n");
  }    
  
  uint64_t ret_hash = sca_hash(*marker, exponent);
  scratchVars.approx_value_in_use = 0;
  return ret_hash;
}
#endif // SPRmode



//   Authors: Bill Carlson, Jenny Zito, Steve Cuccaro
hash_table_t *
alloc_hash(int exponent)
{
  hash_table_t *t = calloc(1,sizeof(hash_table_t));
  
  if (!t)
    ALLOCFAIL();

  t->exponent = exponent;
  t->nentries = (uint64_t)1 << exponent;
  t->mask = t->nentries - 1;
  t->heads = calloc(t->nentries, sizeof(hash_node_t *));
  if (!t->heads)
    ALLOCFAIL();
  
  t->tails = calloc(t->nentries, sizeof(hash_node_t *));
  if (!t->tails)
    ALLOCFAIL();

#ifdef HASHSTATS
  t->num_accesses = calloc(t->nentries, sizeof(uint64_t));
  t->num_nodes = calloc(t->nentries, sizeof(uint64_t));
#endif

  return t;
}

void dealloc_hash(hash_table_t *table)
{
#ifdef HASHSTATS
  free(table->num_nodes);
  free(table->num_accesses);
#endif
  free(table->tails);
  free(table->heads);
  free(table);
}

//   Authors: Bill, Jenny, Steve
hash_node_t *
alloc_hash_node()
{
        hash_node_t *n = calloc(1,sizeof(hash_node_t));
	if (n != NULL) {
	  //	  n->record_ptr = (void *) MATRIX_PTR_INVALID;
	  n->packedID = 0; // also sets n->record_ptr = NULL;
	  n->next = NULL;
	  n->prev = NULL;
	} else {
		ALLOCFAIL();
	}
	return n;
}

/*
 * lookup/insert/delete an entry
 */

hash_node_t *
hash_lookup(hash_table_t *table, record_ptr_t record_ptr,
      int64_t packedID, uint64_t hash_val)
{
  fprintf(stderr,"WARNING: In %s which has NO CODE yet\n",__func__);
  return NULL;
}

// NOT USED RIGHT NOW
//   Authors: Bill, Jenny, Steve
hash_node_t *
hash_insert_at_head(hash_table_t *table, record_ptr_t record_ptr,
      int64_t packedID, uint64_t hash_val)
{
  if ((record_ptr != RECORD_PTR_INVALID) && (packedID != MATRIX_ID_INVALID))
    {
      fprintf(stderr,"in %s, both record_ptr and packedID not invalid!\n",
              __func__);
      return (hash_node_t *)NULL;
    }
  if ((record_ptr == RECORD_PTR_INVALID) && (packedID == MATRIX_ID_INVALID))
    {
      fprintf(stderr,"in %s, both record_ptr and packedID invalid!\n",
              __func__);
      return (hash_node_t *)NULL;
    }

  hash_node_t *new = alloc_hash_node();

  // might want to put this check in a #ifdef DEBUG  
  if (hash_val >= table->nentries)
    {
      fprintf(stderr, "bad hash value %" PRIu64 "\n", hash_val);
      exit(1);
    }

  hash_node_t *head = table->heads[hash_val];
  
  if (record_ptr != RECORD_PTR_INVALID) new->record_ptr = record_ptr;
  else new->packedID = packedID;
  new->next = head;
  new->prev = NULL;

  if (head)
    head->prev = new;
  else
    {
      table->tails[hash_val] = new;
      table->active_chains++;
    }
  
  table->heads[hash_val] = new;
  table->active_entries++;
  table->inserts++;
#ifdef HASHSTATS
  (table->num_nodes[hash_val])++;
#endif

  return new;
};

//   Authors: Bill, Jenny, Steve
/* record_ptr is used to take a pointer to a matrix record or op_store record
   and turn it into a uint64_t that can be referenced */
hash_node_t *
hash_insert_at_tail(hash_table_t *table, record_ptr_t record_ptr,
                    int64_t packedID, uint64_t hash_val)
{
  if ((record_ptr != RECORD_PTR_INVALID) && (packedID != MATRIX_ID_INVALID))
    {
      fprintf(stderr,"in %s, both record_ptr and packedID not invalid!\n",
              __func__);
      return (hash_node_t *)NULL;
    }
  if ((record_ptr == RECORD_PTR_INVALID) && (packedID == MATRIX_ID_INVALID))
    {
      fprintf(stderr,"in %s, both record_ptr and packedID invalid!\n",
              __func__);
      return (hash_node_t *)NULL;
    }
  hash_node_t *new = alloc_hash_node();

  // might want to put this check in a #ifdef DEBUG  
  if (hash_val >= table->nentries)
    {
      fprintf(stderr, "bad hash value %" PRIu64 "\n", hash_val);
      exit(1);
    }

  hash_node_t *tail = table->tails[hash_val];

  if (record_ptr != RECORD_PTR_INVALID) new->record_ptr = record_ptr;
  else new->packedID = packedID;
  new->next = NULL;

  if (tail) 
    {
       new->prev = tail;
       tail->next = new;
    }
  else // no tail, and no head
    {
       table->heads[hash_val] = new;
       new->prev = NULL;
       table->active_chains++;
    }

  table->tails[hash_val] = new;
  table->active_entries++;
  table->inserts++;
#ifdef HASHSTATS
  (table->num_nodes[hash_val])++;
#endif

  return new;
}


//   Authors: Bill, Jenny, Steve
// record_ptr is a uint64_t * which is either a matns_ptr_t, mats_ptr_t,
// op_ptr_t, or info_ptr_t
int
hash_node_remove(hash_table_t *table, record_ptr_t record_ptr, 
                 int64_t packedID, uint64_t hash_val)
{
  if (hash_val >= table->nentries)
    {
      fprintf(stderr, "in %s, bad hash value %" PRIu64 "\n",__func__,hash_val);
      exit(1);
    }

  if ((record_ptr != RECORD_PTR_INVALID) && (packedID != MATRIX_ID_INVALID))
    {
      fprintf(stderr,"in %s, both record_ptr and packedID not invalid!\n",
              __func__);
      return 0;
    }
  if ((record_ptr == RECORD_PTR_INVALID) && (packedID == MATRIX_ID_INVALID))
    {
      fprintf(stderr,"in %s, both record_ptr and packedID invalid!\n",
              __func__);
      return 0;
    }

  // exactly one of record_ptr and packedID is valid

  // starting at head node of hash chain look for our record
  hash_node_t *n = table->heads[hash_val];
  while (n)
    {
      if ( (n->record_ptr == record_ptr) || (n->packedID == packedID) )
	{
          hash_node_remove_node(table, n, hash_val);
	  return 1;   // successfully removed node
	}
      n = n->next;
    }
  return 0;   // nothing was removed
}


//   Authors: Bill, Jenny, Steve
// n points to the specific node to be removed
hash_node_t *
hash_node_remove_node(hash_table_t *table, hash_node_t *n, uint64_t hash_val)
{
  hash_node_t *next_n = n->next;
  if (n->prev)
    {
      n->prev->next = n->next;
    }
  else
    {
      table->heads[hash_val] = n->next;
      if (!n->next)
        {
          table->active_chains--;
        }
    }

  if (n->next)
    {
      n->next->prev = n->prev;
    }
  else
    { 
      table->tails[hash_val] = n->prev;
    }
	  
  table->active_entries--;

#ifdef HASHSTATS
  (table->num_nodes[hash_val])--;
#endif

  // free node pointer, the calling program must free record pointer!
  free(n); 
  return next_n;
}


/*
 * report out
 */
uint64_t 
hash_report(hash_table_t *table, FILE *fp, int verbose)
{
  // record_size is the size of a matrix record or op_store record
  // buckets: total number of entries in store (empty or not, M=1024*1024)
  // active_chains: number of chains with something in them
  // active_entries: number of hash nodes in the entire hash table 
  // hash chains occupied: 100 * (#active_chains/#buckets) 

  uint64_t size = (table->nentries * sizeof (hash_node_t *)
		   + table->active_entries * (sizeof (hash_node_t) + table->record_size));
  fprintf(fp, "The length of the hash table (num hash buckets) is: %" PRIu64 "M\n",
          table->nentries/(1024*1024));
  fprintf(fp, "The memory used for the hash table is: %3.2fMB\n",
          size/(1024.0*1024.0));
  fprintf(fp, "The current percentage of hash buckets which are occupied is: %7.4f%%\n",
          100.0*table->active_chains/table->nentries);
  fprintf(fp, "The current average depth of hash chains for occupied buckets is: %3.1f\n",
          table->active_chains ? (1.0*table->active_entries/table->active_chains):0);
  fprintf(fp, "Total number of hash nodes (items) ever inserted into the hash table is: %" PRIu64 "\n",
   	  table->inserts);
  fprintf(fp, "Total num of hash nodes (items) ever deleted from hash table is not tracked.\n");
  fprintf(fp, "The ratio (num retrievals from the store)/(num insertions): %g\n",
  	  table->inserts?(1.0*table->hits/table->inserts):0);
// commenting out following line to get cumulative numbers instead of numbers
// since last report was printed
//  table->hits = table->misses = table->inserts = 0;

  // HOW DOES THIS COMPARE TO THEORETICAL EXPECTED NUMBER OF EMPTY BINS, AND VAR
  uint64_t r = table->active_entries;     // number of rocks (balls)
  uint64_t s = table->nentries;           // number of socks (bins)
  // p.5 Random Allocations, Kolchin, Sevast'Yanov and Chistyakov, eqn12: D \mu_0
  long double expected_empty = s * powl((1- (1.0L/s)),r);
  long double expected_mean_falling2 = s * (s-1) * powl((1-2.0L/s),r);
  long double square_expected_mean = expected_empty * expected_empty;
  long double var_empty = expected_mean_falling2 + expected_empty
			       - square_expected_mean;
  uint64_t actual_empty = table->nentries - table->active_chains;
  fprintf(fp, "The actual number of empty bins is %ld.\n", (long int)actual_empty);
  fprintf(fp, "The expected number of empty bins is %ld.\n",(long int)roundl(expected_empty));
  fprintf(fp, "The variance of this is %Lg.\n",var_empty);
  long double stdev_empty = sqrtl(var_empty);
  fprintf(fp, "The standard deviation of this is %Lg.\n",stdev_empty);
  // check whether it is 0 or negative 0
  if (fpclassify(stdev_empty) != FP_ZERO) {
    long double difference = (long double)fabsl(actual_empty - expected_empty);
    long double num_stdevs = difference / stdev_empty;
    fprintf(fp, "The actual number of empty bins is %Lg std deviations from mean.\n",
	    num_stdevs);
  }
#ifdef HASHSTATS  
  uint64_t max_chainlength = 0;
  for (int hash_val = 0; hash_val <= table->nentries; hash_val++) {
    if (table->num_nodes[hash_val] > max_chainlength) {
      max_chainlength = table->num_nodes[hash_val];
    }
  }
  fprintf(fp, "The length of the longest hash chain is %ld\n",max_chainlength);
#endif
  
  return size;

}



#ifdef HASHSTATS


/************************************************************************
 *  hashstats_to_files creates files with the following                 *
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

void hashstats_to_files(
                         hash_table_t *hash_table, 
                         char *accesses_file, 
                         char *nodes_file, 
                         char *report_file
                      ) 
{
  FILE *fa, *fn, *fp; 
  
  fa = fopen(accesses_file, "w"); 
  fn = fopen(nodes_file, "w"); 
  
  if (!strcmp(report_file,"stdout"))
    fp = stdout;
  else
    fp = fopen(report_file, "w"); 

  
  // have option to pass in stdout for file path in hash_report
  hash_report(hash_table, fp, 1);
  
  printf("Printing hash table accesses and nodes stats to files.\n");
  
  for (int hash_val = 0; hash_val <= hash_table->nentries; hash_val++) {
    fprintf(fa,"%ld\n",hash_table->num_accesses[hash_val]);
    fprintf(fn,"%ld\n",hash_table->num_nodes[hash_val]);
  }
  
  fclose(fa);
  fclose(fn);
}

// end #ifdef HASHSTATS
#endif


// Hashes for basic types.

/***************************************************************
 *    mult_golden_hash
 *    Author: Jenny Zito
 * LARC can use mult_golden_hash(key, E) as an index into an 
 * E-bit hash table.
 * 
 * LARC often uses mult_golden_hash(key, 64), a hash function which
 * outputs a unsigned 64-bit integer, to build other hash functions.
 * 
 *****************************************************************/
uint64_t mult_golden_hash(uint64_t key, uint64_t exponent)
{
   // LARC uses hash functions for a variety of reasons including the quick
   // retrieval of previously stored information in possibly huge hash tables
   // such as the ScalarStore and the MatrixStore.  So it is very important
   // that LARC has a very well-behaved hash function.
   // LARC has a simple, but elegant, hash function from which all its other
   // hash functions are built.  This hash function is a multiplicative
   // Fibonacci hash, mult_golden_hash(), which acts on 64-bit integers.
   //
   // A multiplicative hash H(key) uses a constant M (called the multiplier) 
   // which is between 0 and 1 and outputs an integer hash value 
   // between 0 and S-1 (in our case  S = hash table size) as follows:
   //        H(key) = integer part of [S * (fractional part of (M * key))]. 
   // See more in Corman, Leiserson, Rivest and Steine's 
   // Introduction to Algorithms (3rd Edition, pp. 263-264).
   //
   // Fibonacci hashes are particularly well-behaved multiplicative hashes
   // whose multiplicative constant is close to the golden ratio
   //         M = (sqrt(5) - 1)/2 = .6180339886... .  
   // See Knuth's Sorting and Searching (2nd Edition, pp.516-518). 
   // As Knuth mentions, it was shown by Vera Sos that the golden ratio
   // is one of a general class of numbers which give multiplicative hashes
   // with the best possible separations when given sequential values as
   // input.  This property of the golden ratio was known to botanists as
   // early as 1837 and seen in the rotational spacing between plant parts.
   //
   // We use a few tricks to implement the multiplicative Fibonacci hash
   // in LARC mult_golden_hash().  It will act on 64-bit integers and output
   // E-bit integers where E is the exponent of the hash table size S = 2^E.
   // We never expect E to be as large as 64, which means that when
   // implementing a multiplicative hash we only need to know at most
   // the first 64-bits to the right of the decimal place in our computations.
   // In particular, we only need to know the first 64-bits to the right of the
   // decimal of the golden ratio M = (sqrt(5) - 1)/2.  We precompute these
   // 64 bits with multiprecision and store them as an integer
   //               11400714819323198485ULL .    
   // This is equal to the integer part of  [2^{64} * (sqrt(5) - 1)/2].
   //
   // We were able to use this prestored value to implement LARC's
   // multiplicative Fibonacci hash entirely in terms of integers:
   //
	return (11400714819323198485ULL * key) >> (64 - exponent);
   //
   // It takes a minute to convince oneself that this is equivalent to the
   // multiplicative Fibonancci hash H(key) which was described above as
   // H(key) =  the integer part of [(2^E) * (fractional part of (M * key))] 
   // where M =  (sqrt(5) - 1)/2].
   //
   // First, remember that: 
   //    11400714819323198485 =  the integer part of [2^{64} * M], 
   //         where the golden mean M = (sqrt(5) - 1)/2] = .61803... , thus
   //    11400714819323198485 is less than 2^{64}, and we can fit it in
   // an unsigned 64-bit integer.
   // This corresponds to the first 64-bits of the fractional part of M.
   //
   // Because all the arguments are unsigned 64-bit integers, when we
   // compute the product P = (11400714819323198485ULL * key)
   // the result is only the low 64-bits of the integer result.
   // Thus it follows that this product as a uint64_t corresponds to
   // the first 64-bits of the fractional part of (M * key).
   //
   // In LARC we assume  the hash exponent E can not be larger than 64,
   // so the number (64-E) is between 0 and 64.
   // The right shift operation >> (64-E) extracts the top E bits of P
   // which are the first E bits of the fractional part of (M * key).
   // Therefore,  the routine mult_golden_hash returns the intended value,
   //    =   the integer part of [(2^E) * (fractional part of (M * key))] 
   //          where M is the golden mean (sqrt(5) - 1)/2].
   //
   // Without having to use any floating point math, we have calculated
   // the multiplicative Fibonacci hash.  
   //
}


/***************************************************************
 *    recursive_hash_from_two_integers
 *    Author: Jenny Zito 
 *    This function takes as input two uint64_t unsigned integers
 *    a,b, and an hash table size exponent e (2^e = size of table).
 *    We take the first integer and multiply it by the fractional
 *    part of the golden mean, then bit-wise mod 2 add in the next
 *    integer and take the mult_golden_hash, grabbing exponent bits.
 ****************************************************************/
uint64_t recursive_hash_from_two_integers
(uint64_t int1, uint64_t int2, uint64_t exponent)
{
  // alternately take mult_golden_hash with full 64 bits, 
  // and mode 2 add in a new integer, ending with a mult_gold_hash with exponent bits.
  uint64_t ping,pong;
  ping = int1;
  pong = mult_golden_hash(ping,64);
  ping = pong ^ int2;
  pong = mult_golden_hash(ping,exponent);

  return pong;

}


/***************************************************************
 *    recursive_hash_from_three_integers
 *    Author: Jenny Zito
 *    This function takes as input three uint64_t unsigned integers
 *    int1,int2,int3 and an hash table size exponent e (2^e = size of table).
 *    We take the first integer and apply a 64 bit mult_golden_hash
 *    (which multiplies it by the fractional part of the golden mean
 *    and and returns the full 64 bits of the result). 
 *    Then we bit-wise mod 2 add in the next integer, and 
 *    alternatively take the mult_golden_hash with 64 bits of result,
 *    and add in the next integer.  We end with an application of
 *    mult_golden_hash with only returning exponent bits.
 ****************************************************************/
uint64_t recursive_hash_from_three_integers
(uint64_t int1, uint64_t int2, uint64_t int3, uint64_t exponent)
{
  // alternately take mult_golden_hash with full 64 bits, 
  // and mode 2 add in a new integer, ending with a mult_gold_hash with exponent bits.
  uint64_t ping,pong;
  ping = int1;
  pong = mult_golden_hash(ping,64);
  ping = pong ^ int2;
  pong = mult_golden_hash(ping,64);
  ping = pong ^ int3;
  pong = mult_golden_hash(ping,exponent);

  return pong;

}




/***************************************************************
 *    recursive_hash_from_four_integers
 *    Author: Jenny Zito and Steve Kratzer
 *    This function takes as input four uint64_t unsigned integers
 *    a,b,c,d and an hash table size exponent e (2^e = size of table).
 *    We take the first integer and apply a 64 bit mult_golden_hash
 *    (which multiplies it by the fractional part of the golden mean
 *    and and returns the full 64 bits of the result). 
 *    Then we bit-wise mod 2 add in the next integer, and 
 *    alternatively take the mult_golden_hash with 64 bits of result,
 *    and add in the next integer.  We end with an application of
 *    mult_golden_hash with only returning exponent bits.
 ****************************************************************/
uint64_t recursive_hash_from_four_integers
(uint64_t int1, uint64_t int2, uint64_t int3, uint64_t int4, uint64_t exponent)
{
  // alternately take mult_golden_hash with full 64 bits, 
  // and mode 2 add in a new integer, ending with a mult_gold_hash with exponent bits.
  uint64_t ping,pong;
  ping = int1;
  pong = mult_golden_hash(ping,64);
  ping = pong ^ int2;
  pong = mult_golden_hash(ping,64);
  ping = pong ^ int3;
  pong = mult_golden_hash(ping,64);
  ping = pong ^ int4;
  pong = mult_golden_hash(ping,exponent);

  return pong;

}


/***************************************************************
 *    recursive_hash_from_int_list
 *    Author: Jenny Zito
 *    This function takes a tlength-long list called tlist 
 *    of  uint64_t unsigned integers, and a
 *    hash table size exponent e (2^e = size of table).
 *    We take the first integer and apply a 64 bit mult_golden_hash
 *    (which multiplies it by the fractional part of the golden mean
 *    and and returns the full 64 bits of the result). 
 *    Then we bit-wise mod 2 add in the next integer, and 
 *    alternatively take the mult_golden_hash with 64 bits of result,
 *    and add in the next integer.  We end with an application of
 *    mult_golden_hash with only returning exponent bits.
  ****************************************************************/
uint64_t
recursive_hash_from_int_list(uint64_t *tlist, uint64_t tlength, uint64_t exponent)
{
   uint64_t ping, pong;
   int i;
   ping = tlist[0];
   for (i=1;i<tlength;i++) {
     pong = mult_golden_hash(ping,64);
     ping = pong ^ tlist[i];
   }
   pong = mult_golden_hash(ping,exponent);

   return pong;
}

/***************************************************************
 *    hash_from_one_longdouble
 *    Author: Jenny Zito and Laurie Law
 *         This converts entropy of a double into an integer and
 *         pass it to the mult_golden_hash function for hashing
 ****************************************************************/
uint64_t
hash_from_one_longdouble (long double ldoub, uint64_t exponent)
{
  uint64_t int_of_doub;

  // ensure that we aren't dealing with a negative zero
  if (fpclassify(ldoub) == FP_ZERO) ldoub = 0.0;

  // interpret the bits of the double to be an integer
  // int_of_doub = *((uint64_t*)&doub);
  memcpy((void *) &int_of_doub, (void *) &ldoub, sizeof(int_of_doub));

  return mult_golden_hash(int_of_doub, exponent);
}

/***************************************************************
 *    Author: Jenny Zito
 *         This converts entropy of doubles into integers,
 *         then circularly shifts the second integer and mod 2
 *         adds it to the first to get a 64 bit integer to send
 *         to the mult_golden_hash function for hashing
 *         A constant, the "double salt" is added to prevent simple collisions
 *         with real vs. imaginary and with the hash from four integers.
 ****************************************************************/
uint64_t
hash_from_two_longdoubles (long double d1l, long double d2l, uint64_t exponent)
{
  uint64_t int1, int2;
#ifdef DEBUG_HASH_C
  printf("in hash_from_two_longdoubles\n");
#endif // #ifdef DEBUG_HASH_C

  // ensure that we aren't dealing with a negative zero
  if (fpclassify(d1l) == FP_ZERO) d1l = 0.0;
  if (fpclassify(d2l) == FP_ZERO) d2l = 0.0;

  // interpret the bits in each double of the complex number to be an integer
  // int1 = *((uint64_t*)&d1);
  // int2 = *((uint64_t*)&d2);
  memcpy((void *) &int1, (void *) &d1l, sizeof(int1));
  memcpy((void *) &int2, (void *) &d2l, sizeof(int1));

  return recursive_hash_from_two_integers(int1,int2,exponent);
#undef DEBUG_HASH_C
}


#ifdef MAR
uint64_t  hash_tile_index(MAR_tile_index_t *target_tile_index_PTR) {
  // retrieve the exponent
  size_t exponent = get_nonscalar_store_exp();
#ifdef IS_COMPLEX   
    uint64_t real_hash = larc_mpz_hash(target_tile_index_PTR->real_index, 64);
    uint64_t imag_hash = larc_mpz_hash(target_tile_index_PTR->imag_index, 64);
    return recursive_hash_from_two_integers(real_hash, imag_hash, exponent);
#else
    return larc_mpz_hash(target_tile_index_PTR->index, exponent);
#endif
}

// The hash filter consists of hashfilter_bits bits which is a short hash
// function that is independent of the standard hash into MatrixStore. Given
// two different  tile indices A and B that produce the same standard hash from
// the call hash_tile_index, we want their hash filters to not be correlated,
// i.e. hashfilter_tile_index(A) - hashfilter_tile_index(B) should be random.
// In this way, we can use the hash filter as a first round check to see
// whether A and B are the identical tile_indices.  Then if the hash filters
// are the same, we will need to make a secondary check to see if the indices
// are the same.
//
uint64_t hashfilter_of_tileindex(const MAR_tile_index_t *target_tile_index_PTR,
        size_t hashfilter_bits)
{
  
  if (scratchVars.counter_in_use)
      fprintf(stderr,"%s reusing scratchVars.counter!\n",__func__);
  scratchVars.counter_in_use = 1;

  mpz_t *temp_q = &scratchVars.counter;
  //  sca_init(&scalar1);
  mpz_set_ui(*temp_q,1);
  
#ifdef IS_COMPLEX   
    // uint64_t real_hash = larc_mpz_hash(target_tile_index_PTR->real_index, 64);
    mpz_add(*temp_q, *temp_q, target_tile_index_PTR->real_index);
    uint64_t real_hash = larc_mpz_hash(*temp_q, 64);
    scratchVars.counter_in_use = 0;
    uint64_t imag_hash = larc_mpz_hash(target_tile_index_PTR->imag_index, 64);
    return recursive_hash_from_two_integers(real_hash, imag_hash,
	hashfilter_bits);
#else
    mpz_add(*temp_q, *temp_q, target_tile_index_PTR->index);
    uint64_t ret_hash = larc_mpz_hash(*temp_q,hashfilter_bits);
    scratchVars.counter_in_use = 0;
    return ret_hash;
#endif
}




#endif //MAR

