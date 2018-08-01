//                                  hash.c 
//                        Hash functions for LARC
//                  for fast find_or_insert functions 
//     for the Matrix Store and the matrix operation stores                                        
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

// Standard libraries
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <string.h> //for memcpy

// Our header files structures and functions
#include "larc.h"
#include "hash.h"


/******************************************************************
 * BACKGROUND:
 *  The Matrix Store stores matrices in a Matrix Table which uses
 *  for a unique matrixID the integer valued index of the row 
 *  in which the matrix is stored of the Matrix Table (mat_add_t).
 *  Having this as an integer allows check pointing and storage of tables.
 *
 *  The hashes provide a fast front end for storage and retrieval
 *  into the store.  The Links column of the Matrix Table keeps track of
 *  the next row in each hash collison chain, as well as the next
 *  free row in the free list of the Matrix Table.  The hash acts
 *  on the MatrixVals stored in the Matrix Table, which are either
 *  a pair of doubles corresponding to the real and imaginary part
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
 *  Trace is handled differently since it's output is a ScalarType
 *  number instead of a matrixID
 *******************************************************************/


/***************************************************************
 *    mult_golden_hash       
 *    Author: Jenny Zito
 *      Takes in a 64 bit unsigned integer and the exponent for hash table
 *      returns a hash value  between 0 and 2^exponent-1.
 *      The method is based on the multiplicative hash
 *      described in Knuth Vol 3 which uses the golden mean.
 *      The trick is to try to get a full 64 bits to right of decimal.
 *      A straight-forward (but faulty) implementation would only
 *      have 56 bits of entropy from the double golden ratio.
 *      Also see p. 263 multiplicative hash in Corman, Leiserson Rivest and Steine
 *      To get around this we created golden_big_int which contains
 *      the fractional part of the golden mean shifted into a 64-bit
 *      integer. Bill Carlson helped me get the double precision number
 *      for the integer portion of 2^64 * (sqrt(5)-1)/2 using gp calculator.
 *      The magic constant here is ^^^^^^^^^^^^^^^^^^^^ that, equal to
 *      11400714819323198485ULL.
 *****************************************************************/
static inline uint64_t
mult_golden_hash(uint64_t key, uint64_t exponent)
{
	return (11400714819323198485ULL * key) >> (64 - exponent);

}

#ifdef USE_COMPLEX
/***************************************************************
 *    Author: Jenny Zito
 *         This converts entropy of doubles into integers, 
 *         then circularly shifts the second integer and mod 2
 *         adds it to the first to get a 64 bit integer to send
 *         to the mult_golden_hash function for hashing
 *         A constant, the "double salt" is added to prevent simple collisions
 *         with real vs. imaginary and with the hash from four integers. 
 ****************************************************************/
static uint64_t
hash_from_two_doubles (double d1, double d2, uint64_t exponent)
{

  uint64_t int1, int2; 

  // ensure that we aren't dealing with a negative zero
  if (fpclassify(d1) == FP_ZERO) d1 = 0.0;
  if (fpclassify(d2) == FP_ZERO) d2 = 0.0;

  // interpret the bits in each double of the complex number to be an integer
  // int1 = *((uint64_t*)&d1);
  // int2 = *((uint64_t*)&d2);
  memcpy((void *) &int1, (void *) &d1, sizeof(int1));
  memcpy((void *) &int2, (void *) &d2, sizeof(int1));
  

  // circularly shift the second integer and mod 2 add it to the first integer
  //
  // Add in the "double salt", which I selected to be:
  // The integer part of 2^64 * (exp(1)-2)
  // Using the double precision calculator gp:
  //   gp
  //   ? exp(1)
  //   %1  = 2.7182818284590452353602874713526624978
  //   ?  2^64 * (exp(1)-2)
  //   %2 = 13249961062380153450.747823268317876541
  // Thus we set
  //   uint64_t int_salt_from_e = 13249961062380153450ULL;
  // And add this value to the imaginary part of the integer from the 
  // imaginary part of the complex value.  Then  we circularly shift the
  // sum 32 bits and add it to the originals.
  //
  // 1. This keeps simple combinations of imag, real from colliding with each other 
  //    e.g. 0+ 1i would collide with 2^32 + 0i
  // 2. This keeps simple inputs of the four-integers hash from colliding with 
  //    simple inputs of the two-doubles hash 
  //    e.g. 0 + 0i   would collide with (0,0,0,0)


  /**************************************************************
   *   jszito 11/2017 Replacing this code with                  *
   *   my recursive golden multiply hash                        *
   *   recursive_hash_from_two_integers                          *
   **************************************************************/
  // uint64_t shifted_sum = int1 ^ (int2 >> 32) ^ (int2 << 32);   
  // uint64_t int_salt_from_e = 13249961062380153450ULL;          
  // uint64_t result = shifted_sum ^ int_salt_from_e;             
  // return mult_golden_hash(result, exponent);                   

  return recursive_hash_from_two_integers(int1,int2,exponent);  


}
#endif

#ifdef USE_REAL
/***************************************************************
 *    hash_from_one_double         (called by hash_from_matrix_content)
 *    Author: Jenny Zito and Laurie Law
 *         This converts entropy of a double into an integer and
 *         pass it to the mult_golden_hash function for hashing
 ****************************************************************/
static uint64_t
hash_from_one_double (double doub, uint64_t exponent)
{

  uint64_t int_of_doub; 

  // ensure that we aren't dealing with a negative zero
  if (fpclassify(doub) == FP_ZERO) doub = 0.0;

  // interpret the bits of the double to be an integer
  // int_of_doub = *((uint64_t*)&doub);
  memcpy((void *) &int_of_doub, (void *) &doub, sizeof(int_of_doub));

  return mult_golden_hash(int_of_doub, exponent);
}
#endif


// NON SCALAR
/***************************************************************
 *    hash_from_matrix_panel
 *    Author: Jenny Zito
 *         Returns the hash of a set of four matrix pointers
 *    Jenny, Steve rewrote this in Nov 2016 to handle part of 
 *        hash_from_matrix_content
 ****************************************************************/
uint64_t
hash_from_matrix_panel(mat_add_t sub[4], matrix_type_t mat_type, uint64_t exponent)
{

#ifdef DEBUG
  printf("In routine %s\n",__func__);
  printf("     The hash exponent is %d\n",exponent);
#endif

  // EXCEPTION CHECKING
  if(mat_type==SCALAR) {   
    printf("In hash_from_matrix_panel mat_type was SCALAR\n");
    exit(1);
  }   

  // sub is 4  mat_add_t corresponding to the four submatrices
#ifdef DEBUG
  printf("     sub = [%d , %d , %d , %d]\n",
	 (uint64_t)sub[0], (uint64_t)sub[1], (uint64_t)sub[2], (uint64_t)sub[3]);
#endif

#ifdef DEBUG
  // DEBUGGING
    printf("IN HASH FROM MATRIX PANEL: List of addresses = [%d , %d , %d , %d]\n",	(uint64_t)sub[0], (uint64_t)sub[1], (uint64_t)sub[2], (uint64_t)sub[3]);
    printf("  MOD 64  = [%d , %d , %d , %d]\n",	(uint64_t)sub[0]%64, (uint64_t)sub[1]%64, (uint64_t)sub[2]%64, (uint64_t)sub[3]%64);
#endif

    return recursive_hash_from_four_integers(
                (uint64_t)sub[0], 
                (uint64_t)sub[1], 
                (uint64_t)sub[2], 
                (uint64_t)sub[3], 
                exponent);

  // Experiment with using matrixIDs instead of matrixPTRs
    /*
    return recursive_hash_from_four_integers(  
                                 get_matrixID_from_ptr(sub[0]), 
                                 get_matrixID_from_ptr(sub[1]), 
                                 get_matrixID_from_ptr(sub[2]), 
                                 get_matrixID_from_ptr(sub[3]),
                                 exponent);
    */

}


// SCALAR
/***************************************************************
 *    hash_from_matrix_scalar
 *    Author: Jenny Zito
 *         Returns the hash a scalar valued matrix
 *    Jenny, Steve, and Bill modified in March of 2016
 *        to always collapse things that are within SIGHASH 
 *        significant bits to the same value before hashing
 *        and also before hashing to collapse values smaller
 *        than zerorealthresh to zero.
 *    Steve and Jenny modified in Nov 2016 to handle scalar functionality
 *        of old function hash_from_matrix_content
 ****************************************************************/
uint64_t
hash_from_matrix_scalar(ScalarType scalar, matrix_type_t mat_type, uint64_t exponent)
{
  
#ifdef DEBUG
  printf("In routine %s\n",__func__);
  printf("     The hash exponent is %d\n",exponent);
#endif
  
  // EXCEPTION CHECKING
  if(mat_type != SCALAR) {   // begin SCALAR TYPE
    printf("In hash_from_matrix_scalar mat_type was not SCALAR\n");
    exit(1);
  }   
  
#ifdef DEBUG
#ifdef USE_COMPLEX  // mat_val_ptr is a complex number
  printf("    scalar = [%g , %g]\n",creal(scalar),cimag(scalar));
#else
  printf("     mat_val_ptr = [%g]\n",scalar);
#endif
#endif
  
  // hash any items with same locality_approximation to the same hash location
  ScalarType marker = locality_approx(scalar);
  
#ifdef DEBUG
  // The %la prints a double in hex including the 0x designation
  if (scalar != marker) {
    printf("Before rounding creal is %la (%g), afterwards it is %la (%g)\n",
	   creal(scalar),creal(scalar),creal(marker),creal(marker));
    printf("Before rounding cimag is %la (%g), afterwards it is %la (%g)\n",
	   cimag(scalar),cimag(scalar),cimag(marker),cimag(marker));
  }
#endif    
  
  //    return hash_from_two_doubles(creal(marker), cimag(marker), exponent);
#ifdef USE_COMPLEX
  return hash_from_two_doubles(creal(marker), cimag(marker), exponent);
#endif
#ifdef USE_REAL
  return hash_from_one_double(marker, exponent);
#endif
#ifdef USE_INTEGER
  return  mult_golden_hash(marker, exponent);
#endif
  
}



//   Authors: Bill Carlson, Jenny Zito, Steve Cuccaro
hash_table_t *
alloc_hash(int exponent)
{
  hash_table_t *t = calloc(1,sizeof(hash_table_t));
  
  if (!t)
    ALLOCFAIL();

  t->exponent = exponent;
  t->nentries = 1L<<exponent;
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


//   Authors: Bill, Jenny, Steve
hash_node_t *
alloc_hash_node()
{
        hash_node_t *n = calloc(1,sizeof(hash_node_t));
	if (n != NULL) {
	  //	  n->record_ptr = (void *) MATRIX_PTR_INVALID;
	  n->record_ptr = NULL;
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
hash_lookup(hash_table_t *table, record_ptr_t record_ptr, uint64_t hash_val)
{
  printf("WARNING: In %s which has NO CODE yet\n",__func__);
  return NULL;
}

// NOT USED RIGHT NOW
//   Authors: Bill, Jenny, Steve
hash_node_t *
hash_insert_at_head(hash_table_t *table, record_ptr_t record_ptr, uint64_t hash_val)
{
  hash_node_t *new = alloc_hash_node();

  // might want to put this check in a #ifdef DEBUG  
  if (hash_val >= table->nentries)
    {
      fprintf(stderr, "bad hash value %ld\n", hash_val);
      exit(1);
    }

  hash_node_t *head = table->heads[hash_val];
  
  new->record_ptr = record_ptr;
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
hash_insert_at_tail(hash_table_t *table, record_ptr_t record_ptr, uint64_t hash_val)
{
  hash_node_t *new = alloc_hash_node();

  // might want to put this check in a #ifdef DEBUG  
  if (hash_val >= table->nentries)
    {
      fprintf(stderr, "bad hash value %ld\n", hash_val);
      exit(1);
    }

  hash_node_t *tail = table->tails[hash_val];

  new->record_ptr = record_ptr;
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
};


//   Authors: Bill, Jenny, Steve
// record_ptr is a uint64_t * which is either a mat_add_t, op_add_t, or info_add_t
int
hash_node_remove(hash_table_t *table, record_ptr_t record_ptr, uint64_t hash_val)
{
  if (hash_val >= table->nentries)
    {
      fprintf(stderr, "bad hash value %ld\n", hash_val);
      exit(1);
    }

  // starting at head node of hash chain look for our record
  hash_node_t *n = table->heads[hash_val];
  while (n)
    {
      if (n->record_ptr == record_ptr)
	{
	  if (n->prev)
	    n->prev->next = n->next;
	  else
	    {
	      table->heads[hash_val] = n->next;
	      if (!n->next)
		table->active_chains--;
	    }

	  if (n->next)
	    n->next->prev = n->prev;
          else 
            table->tails[hash_val] = n->prev;
	  
	  table->active_entries--;

#ifdef HASHSTATS
          (table->num_nodes[hash_val])--;
#endif

          // free node pointer, the calling program must free record pointer!
	  free(n); 

	  return 1;   // successfully removed node
	}
      n = n->next;
    }
  return 0;   // nothing was removed
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
  fprintf(fp, "The length of the hash table (num hash buckets) is: %ldM\n",
          table->nentries/(1024*1024));
  fprintf(fp, "The memory used for the hash table is: %3.2fMB\n",
          size/(1024.0*1024.0));
  fprintf(fp, "The current percentage of hash buckets which are occupied is: %7.4f%%\n",
          100.0*table->active_chains/table->nentries);
  fprintf(fp, "The current average depth of hash chains for occupied buckets is: %3.1f\n",
          table->active_chains ? (1.0*table->active_entries/table->active_chains):0);
  fprintf(fp, "Total number of hash nodes (items) ever inserted into the hash table is: %ld\n",
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
  double var_empty = (double) (expected_mean_falling2 + expected_empty
			       - square_expected_mean);
  uint64_t actual_empty = table->nentries - table->active_chains;
  fprintf(fp, "The actual number of empty bins is %d.\n", (int)actual_empty);
  fprintf(fp, "The expected number of empty bins is %d.\n",(int)roundl(expected_empty));
  fprintf(fp, "The variance of this is %g.\n",var_empty);
  double stdev_empty = sqrt(var_empty);
  fprintf(fp, "The standard deviation of this is %g.\n",stdev_empty);
  // check whether it is 0 or negative 0
  if (fpclassify(stdev_empty) != FP_ZERO) {
    double difference = abs(actual_empty - expected_empty);
    double num_stdevs = difference / stdev_empty;
    fprintf(fp, "The actual number of empty bins is %g std deviations from mean.\n",
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



