//                       matrix_store.c 
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
#include <inttypes.h>
#include <math.h>
#include <string.h>
#include <stdint.h>   // do not know if we use
#include <float.h>    // do not know if we use

#include "matrix_store.h"
#include "matmath.h"
#include "hash.h"
#include "op_store.h"
#include "json.h"
#include "global.h"

// we will have a table of tables of this size, mat_ptr_table_list
// each entry points to a mat_ptr_table
// WARNING LOG_SIZE_MPT CAN NOT BE BIGGER THAN 32
#define LOG_SIZE_MPT 25
#define SIZE_MPT (1<<LOG_SIZE_MPT)
#define MASK_MPT (SIZE_MPT - 1)

static struct matrix_store_t {
  uint64_t matrixID_next;    //  is the unique identification number (the true MatrixID)
                           //  because MatrixPTRs can be freed and reused
  size_t nmatrices;    // the total number matrices in the store currently
  size_t nscalars;     // the total number of scalars in the store currently
  size_t max_nmatrices;  // the maximum number of matrices that were ever 
                         // in the store at the same time
  size_t max_nscalars;   // the maximum number of scalars that were ever 
                         // in the store at the same time
  size_t **hist;         // histogram of matrices by row_level and col_level
  hash_table_t *hash_table;  
// SERIAL TODO:  might want to change store.hash_table name to store.MatrixRecord_HashTable

  // PRELOADED VALUES IN MATRIX STORE THAT WE NEVER WANT TO ERASE AND USE IN MATH FUNCTIONS
  mat_level_t largest_level;    // largest row_level and col_level allowed in store, 
                                 // also largest level of preloaded  matrices
  int sighash;           // locality-approximation parameter: equate two scalars 
                         // when rounded to sighash significant bits
  double zerorealthresh; // locality-approximation  parameter: equate to zero 
                         //when smaller than zerorealthresh
  mat_add_t **zero;     // doubly indexed matrix IDs for the preloaded zero matrices
  mat_add_t *identity;  // array of matrix IDs for the preloaded identity matrices (square)
  mat_add_t *iHadamard; // array of matrix IDs for the preloaded integer Hadamard matrices

  // POSSIBLY CAN BE REMOVED                      
  uint32_t deepest;    // longest successful traversal down a hash chain to retrieve a matrix 

  // this is a list of tables each entry of each table contains a pointer to a matrix record
  // the index into the table is the matrixID, broken into lower and higher bits
  mat_add_t **mat_ptr_table_list;

} store = {0};   // this creates the variable "store" which is a struct of type matrix_store_t


mat_level_t maximum_level(void) {
  return store.largest_level;
}


int get_sighash(void) {
  return store.sighash;
}


double get_zerorealthresh(void) {
  return store.zerorealthresh;
}

uint64_t num_matrices_in_store ( void )
{
  uint64_t num_matrices = 0;

  for ( size_t i = 0 ; i < store.largest_level ; i++ )
    for ( size_t j = 0 ; j < store.largest_level ; j++ )
      num_matrices += store.hist[i][j];

  return num_matrices;
}
 

void
matrix_store_report(char *outfilepath)
{
  FILE *f; 

  if (strcmp(outfilepath,"stdout")) {
     printf("Printing matrix store report to file %s\n", outfilepath);
     f = fopen(outfilepath, "a"); 

  }
  else {
    printf("Printing matrix store report to screen\n");
    f = stdout;
  }
	fprintf(f,"\n");
	fprintf(f,"Matrix Store Report:\n");
        hash_report(store.hash_table, f, 0);
	fprintf(f,"Longest successful hash chain traversal:  %u\n", store.deepest);   
	fprintf(f,"Largest allowable level for matrices in Matrix Store: %u \n", 
                   store.largest_level);
	fprintf(f,"Largest number of scalar matrices ever in Matrix Store at once: %zd\n", 
                   store.max_nscalars);
	fprintf(f,"Largest number of nonscalar matrices ever in Matrix Store at once: %zd\n", 
                   store.max_nmatrices);
        size_t total_matrices = 0;
	fprintf(f,"Number of matrices per level (printed only when count greater than 1):\n");
        // int ijsum = 0;
	for (int i = 0; i <= store.largest_level; i++) {
	  for (int j = 0; j <= store.largest_level; j++) {
            total_matrices += store.hist[i][j];
            if (store.hist[i][j] > 1) {
  	      // printf("(%d,%d): %10zd%c", i, j, store.hist[i][j],
              //       (ijsum % 4 == 3) ? '\n' : '\t');
  	      fprintf(f,"  (%2d,%2d): %13zd\n", i, j, store.hist[i][j]);
	    }
            // ++ijsum;
	  }
	}
	fprintf(f,"The total number of matrices currently in the Matrix Store is %zd.\n", total_matrices);
	fprintf(f,"\n");
	
	// Also print memory used by matrix pointer table
        // The high order bits of matrixID_next give the table index
	uint64_t bytes_ptr_list = SIZE_MPT * sizeof(mat_add_t);   	   			// size of array of pointers to tables
	uint64_t bytes_per_table = SIZE_MPT * sizeof(mat_add_t);  	   			// size of each table
        uint64_t num_tables = ((store.matrixID_next - 1) >> LOG_SIZE_MPT) + 1;    //since indexing starts at 0
	uint64_t bytes_ptr_total = bytes_ptr_list + (num_tables * bytes_per_table);
	fprintf(f, "Memory used by matrix pointer table: \n");
	fprintf(f, "Per table size:  %3.2fMB, Total size:  %3.2fMB \n", bytes_per_table/(1024.0*1024.0), bytes_ptr_total/(1024.0*1024.0));	
	fflush(f);
	if  (strcmp(outfilepath,"stdout")) fclose(f);   // if not stdout then close
}


int
create_matrix_store(size_t exponent, mat_level_t max_level, int sighash, int zerobitthresh)
{
  int verbose = 0;
  if (verbose) {
    printf("inside routine %s with:\n", __func__);
    printf("   matrix store hash exponent %zd\n", exponent);
    printf("   maximum matrix level %d\n",max_level);
    printf("   locality-approximation rounds scalars to %d significant bits\n",sighash);
    printf("                and collapses values within %d bits of zero to zero\n",zerobitthresh);
  }

  store.matrixID_next = 0;
  
  // two dimensional array of histogram for levels
  store.hist = (uint64_t **)malloc((max_level+1)*(sizeof(uint64_t *)));
  if (store.hist == NULL) {
    ALLOCFAIL();
    return 0;
  }
  int i;
  for (i=0; i<=max_level; i++) 
    {
      store.hist[i] = calloc(max_level+1, sizeof(uint64_t));
      if (store.hist[i] == NULL) {
	ALLOCFAIL();
	return 0;
      }
    }
 
  // two dimensional array of zeros and one dimensional arrays of identities and iHadamards
  store.zero = (mat_add_t **)malloc((1 + max_level) * sizeof(mat_add_t *));
  store.identity = (mat_add_t *)malloc((1 + max_level) * sizeof(mat_add_t));
  if ((store.zero == NULL) || (store.identity == NULL)) {
    ALLOCFAIL();
    exit(1);    
  }
  for (int i= 0; i <= max_level; ++i) {
    store.zero[i] = (mat_add_t *)malloc((1 + max_level) * sizeof(mat_add_t ));
    if (store.zero[i] == NULL) {
      ALLOCFAIL();
      exit(1);
    }
  }
  store.iHadamard = (mat_add_t *)malloc((1 + max_level) * sizeof(mat_add_t));
  if (store.iHadamard == NULL) {
    ALLOCFAIL();
    exit(1);
  }
 
  store.largest_level = max_level;
  if (verbose) {
    printf("inside routine %s with:\n", __func__);
    printf("  set store.largest_level %d\n",store.largest_level);
  }

  store.sighash = sighash;
  double zerorealthresh = pow(2.0,-zerobitthresh);
  store.zerorealthresh = zerorealthresh;

  store.hash_table = alloc_hash(exponent);
  store.hash_table->record_size = sizeof(larc_matrix_t);
  if (store.hash_table == NULL) {
    ALLOCFAIL();
    return 0;
  }

  store.nmatrices = 0;
  store.max_nmatrices = 0;
  store.nscalars = 0;
  store.max_nscalars = 0;

  // create list of tables called mat_ptr_table_list
  // and initialize the first table in this list 
  // with zeros corresponding to mat_ptr MATRIX_PTR_INVALID
  // ToDo: may want to have 2 sizes: one for the number of tables, and one for the table size
  store.mat_ptr_table_list = (mat_add_t **)calloc(SIZE_MPT,sizeof(mat_add_t *));  //list of ptrs to tables
  store.mat_ptr_table_list[0] = (mat_add_t *)calloc(SIZE_MPT,sizeof(mat_add_t));  //the first table

  return 1;
}


/************************************************************
 *                   lock_matrix                            *
 *  This routine is used to prevent important matrices      *
 *  (such as identity matrices, and zero matrices required  *
 *  for the math routines to work) from EVER being removed  *
 *  from the matrix store.                                  *
 *  A lock should never be removed.                         *
 *  NOTE: If you want to temporarily hold a matrix use      *
 *        the function set_hold_matrix, and use function    *
 *        release_hold_matrix to remove the hold.           *
 *                                                          *
 ***********************************************************/
int lock_matrix(mat_add_t m_ptr) {
  if (m_ptr == MATRIX_PTR_INVALID) {
    printf ("ERROR: in lock_matrix, matrix does not exist!\n");
    return (0);
  }

  if (m_ptr->lock == 1) {
    printf ("WARN: in lock_matrix, attempt lock_matrix for matrix already locked\n");
    return (1);
  }

  m_ptr->lock = 1;
  
  return (1);
}   

/************************************************************
 *                   set_hold_matrix_from_matrixID	    *
 *    		Python Interface version:accepts a	    *
 *        	matrixID instead of a pointer          *
 *  This routine is used to temporarily hold a matrix       *
 *  so it is not removed from the matrix store.             *
 *  Use the function release_hold_matrix to undo the hold.  *
 *  NOTE: lock_matrix is a permanent lock which is          *
 *  used to prevent important matrices                      *
 *  (such as identity matrices, and zero matrices required  *
 *  for the math routines to work) from EVER being removed  *
 *  from the matrix store.                                  *
 *  A lock should never be removed, a hold can be.          *
 *                                                          *
 ***********************************************************/
int set_hold_matrix_from_matrixID(int64_t m_mID) {
  // get the matrix pointers from the matrixID, and see if still in store 
  mat_add_t m_ptr = mat_ptr_from_matrixID(m_mID, "first", __func__,0);
  return set_hold_matrix(m_ptr);
}

/************************************************************
 *                   set_hold_matrix                        *
 *  This routine is used to temporarily hold a matrix       *
 *  so it is not removed from the matrix store.             *
 *  Use the function release_hold_matrix to undo the hold.  *
 *  NOTE: lock_matrix is a permanent lock which is          *
 *  used to prevent important matrices                      *
 *  (such as identity matrices, and zero matrices required  *
 *  for the math routines to work) from EVER being removed  *
 *  from the matrix store.                                  *
 *  A lock should never be removed, a hold can be.          *
 *                                                          *
 ***********************************************************/
int set_hold_matrix(mat_add_t m_ptr) {
  if (m_ptr == MATRIX_PTR_INVALID) {
    printf ("ERROR: in set_hold_matrix, matrix does not exist!\n");
    return (0);
  }

  ++m_ptr->hold;
  
  return (1);
}   


/************************************************************
 *          release_hold_matrix_from_matrixID	            *
 *    	Python Interface version of release_hold_matrix	    *
 *      using matrixID instead of a pointer            *
 ************************************************************/
int release_hold_matrix_from_matrixID(int64_t m_mID) {
  // get the matrix pointers from the matrixID, and see if still in store 
  mat_add_t m_ptr = mat_ptr_from_matrixID(m_mID, "first", __func__,0);
  return release_hold_matrix(m_ptr);
}




/****************************************************************
 *                  release_hold_matrix                         *
 *  See explanation above for set_hold_matrix and lock_matrix.  *
 *                                                              *
 ****************************************************************/
int release_hold_matrix(mat_add_t m_ptr) {
  if (m_ptr == MATRIX_PTR_INVALID) {
    printf ("ERROR: in release_hold_matrix, matrix does not exist!\n");
    return (0);
  }

  if (m_ptr->hold == 0) {
    printf ("WARN: attempt release_hold_matrix for matrix with no holds\n");
    return (1);
  }

  --m_ptr->hold;
  
  return (1);
}   


/************************************************************************
 *               preload_matrix_store                                   *
 *   This should not be done until the op_store has been initialized.   *
 *   This function loads all the basic 1 and 2 bit matrices.            *
 *   It then loads all the identity matrices and rectangular zero       *
 *   matrices up to the maximum level size store.largest_level. It      *
 *   also loads all the integer Hadamard matrices up to largest_level.  *
 ***********************************************************************/   
int preload_matrix_store(void) {

  int verbose = 0;
  if (verbose) {
    printf("Inside routine %s with:\n", __func__);
    printf("   about to run init_globals\n");
  }
 
  int top_level = store.largest_level;
  if (verbose) {
    printf("Inside routine preload_matrix_store with:\n");
    printf("   retrieved largest level %d\n", top_level);
  }
  
  // panels of sub matrices  
  mat_add_t sm[4];
  
  // 0, 1, -1 are constants needed for preload of zero, identity, integer Hadamard
  ScalarType scalar0;
  ScalarType scalar1;
  ScalarType scalarM1;

#ifndef USE_COMPLEX //REAL OR INTEGER
  scalar0 = 0;
  scalar1 = 1;
  scalarM1 = -1;  
#else   // COMPLEX
  scalar0 = 0 + 0*I;
  scalar1 = 1 + 0*I;
  scalarM1 = -1 + 0*I;
#endif 

  // LOAD ALL THE ZERO MATRICES
  store.zero[0][0] = matrix_get_ptr_scalar(scalar0);
  if (!lock_matrix(store.zero[0][0])) {
      printf("FAIL: scalar zero preload tried to lock a matrix which does not exist\n");
      exit(1);
    }

  for (int k = 1; k <= top_level; k++) {

    // Case: ROW_VECTOR  (0,2**k)
    sm[0] = store.zero[0][k-1];
    sm[1] = store.zero[0][k-1];
    sm[2] = MATRIX_PTR_INVALID;
    sm[3] = MATRIX_PTR_INVALID;
    store.zero[0][k] = matrix_get_ptr_panel(sm, 0, k);
    if (!lock_matrix(store.zero[0][k])) {
      printf("FAIL: row Z preload tried to lock a matrix which does not exist\n");
      exit(1);
    }

    // Case: COL_VECTOR (2**k,0)
    sm[0] = store.zero[k-1][0];
    sm[1] = MATRIX_PTR_INVALID;
    sm[2] = store.zero[k-1][0];
    sm[3] = MATRIX_PTR_INVALID;
    store.zero[k][0] = matrix_get_ptr_panel(sm, k, 0);
    if (!lock_matrix(store.zero[k][0])) {
      printf("FAIL: col Z preload tried to lock a matrix which does not exist\n");
      exit(1);
    }
  }

  // Case: MATRIX  (2**i,2**j) with (0 < i,j <= k) 
  for (int i = 1; i <= top_level; i++) {
    for (int j = 1; j <= top_level; j++) {
      sm[0] = store.zero[i-1][j-1];
      sm[1] = store.zero[i-1][j-1];
      sm[2] = store.zero[i-1][j-1];
      sm[3] = store.zero[i-1][j-1];
      store.zero[i][j] = matrix_get_ptr_panel(sm, i, j);
      if (!lock_matrix(store.zero[i][j])) {
	printf("FAIL: matrix Z preload tried to lock a matrix which does not exist\n");
	exit(1);
      }
    }
  }

  // Prestore all the IDENTITY matrices (square)
  store.identity[0] = matrix_get_ptr_scalar(scalar1);
  if (!lock_matrix(store.identity[0])) {
      printf("FAIL: scalar 0 preload tried to lock a matrix which does not exist\n");
      exit(1);
    }

  for (int i = 1; i <= top_level; i++) {
    sm[1] = sm[2] = store.zero[i-1][i-1];
    sm[0] = sm[3] = store.identity[i-1];
    store.identity[i] = matrix_get_ptr_panel(sm, i, i);
    if (!lock_matrix(store.identity[i])) {
      printf("FAIL: I preload tried to lock a matrix which does not exist\n");
      exit(1);
    }
  }

  // Prestore all the integer HADAMARD matrices and scalar negative one.

  // Load scalar negative one into the otherwise un-used store.iHadamard[0]
  store.iHadamard[0] = matrix_get_ptr_scalar(scalarM1);
  if (!lock_matrix(store.iHadamard[0])) {
      printf("FAIL: scalar negative one preload tried to lock a nonexistent matrix\n");
      exit(1);
    }

  sm[0] = sm[1] = sm[2] = store.identity[0];
  sm[3] = store.iHadamard[0];
  store.iHadamard[1] = matrix_get_ptr_panel(sm,1,1);
  if (!lock_matrix(store.iHadamard[1])) {
      printf("FAIL: iHadamard[1] preload tried to lock a nonexistent matrix\n");
      exit(1);
    }

  for (int i = 2; i <= top_level; i++) {
    sm[0] = sm[1] = sm[2] = store.iHadamard[i-1];
    sm[3] = scalar_mult(store.iHadamard[0],store.iHadamard[i-1]);
    store.iHadamard[i] = matrix_get_ptr_panel(sm, i, i);
    if (!lock_matrix(store.iHadamard[i])) {
      printf("FAIL: HH preload tried to lock a matrix which does not exist\n");
      exit(1);
    }
  }

  //  Preload the other basic 1 and 2 bit matrices and set global names
  init_globals();
 
  return 1;
}


// Python interface version 
int64_t get_zero_matrixID(mat_level_t row_level, mat_level_t col_level) {
  mat_add_t m_ptr = get_zero_matrix_ptr(row_level, col_level);
  return get_matrixID_from_ptr(m_ptr);
 }  


mat_add_t
get_zero_matrix_ptr(mat_level_t row_level, mat_level_t col_level)  
{
  if ((MAX(row_level,col_level)) > store.largest_level) 
  {
    printf("Error: requested zero matrix larger than max_level input to create store.\n");
    exit(1);
  }
  return(store.zero[row_level][col_level]);
}



// Python interface version 
int64_t get_identity_matrixID(mat_level_t level) {
  mat_add_t m_ptr = get_identity_matrix_ptr(level);
  return get_matrixID_from_ptr(m_ptr);
}



mat_add_t get_identity_matrix_ptr(mat_level_t level) {
  if (level >  store.largest_level) 
    {
      printf("Error: requested identity matrix has size %d\n",level);
      printf("       which is larger than maximum matrix size %d.\n",store.largest_level);
      exit(1);
    }
  
  return(store.identity[level]);
}



// Python interface version 
int64_t get_iHadamard_matrixID(mat_level_t level) {
  mat_add_t  m_ptr = get_iHadamard_matrix_ptr(level);
  return get_matrixID_from_ptr(m_ptr);
}




mat_add_t get_iHadamard_matrix_ptr(mat_level_t level) {

  if (level == 0)
  {
    printf("Error in get_iHadamard_matrix_ptr: there is no sensible definition ");
    printf("for a level 0 Hadamard matrix\n");
    exit(1);
  }

  if (level > store.largest_level) 
  {
    printf("Error: requested integer Hadamard matrix larger than maximum size.\n");
    exit(1);
  }

  return(store.iHadamard[level]);
}


int
matrix_appears_as_sub_count_increment(mat_add_t id)
{
  return ++(id->appears_as_sub_count);
}


int
matrix_appears_as_sub_count_decrement(mat_add_t id)
{
  return --(id->appears_as_sub_count);
}

/// SCALAR
/* Returns 0 if matrix differs from provided content, non-zero if equal */
/* this function uses the round_sig_fig_real function to collapse nearly 
   equivalent scalars to the same value, to sighash bits, given in larc.h */
/* this function uses the collapse_near_zero function to collapse tiny numbers to zero 
   with given threshold in larc.h collapse_near_zero(complex input, zerorealthresh) */
//  WARNING: locality-approximation is used here, and could have dire consequences if
//           you don't preload zeros and identities.
/* The second argument of this routine scalar is a ScalarType (complex, double, int) */
static int
matrix_compare_scalar(mat_add_t m_ptr, ScalarType scalar)
{
  // EASY CASES for failed comparison
   if ((m_ptr == MATRIX_PTR_INVALID) || 
       (matrix_type(m_ptr) != SCALAR) )
  { 
     return 0; 
  }

#ifdef USE_INTEGER
  // compare the two integers (no rounding should be present)
  int iseq = (matrix_trace(m_ptr) == scalar);
#endif
#ifndef USE_INTEGER   // COMPLEX or REAL
  // Two scalars are treated as equivalant, if their locality-approximations are the same
  // we are comparing one element from store and one new element
  ScalarType marker_of_stored_element = locality_approx(matrix_trace(m_ptr));
  ScalarType marker_of_new_element = locality_approx(scalar);

  // TODO: change this check into part of unit testing.
  // This should not happen if locality-approximation and storing are done correctly
  if (marker_of_stored_element == 0.0 && !(matrix_is_zero(m_ptr)) ) {
    printf("ERROR in %s, should not have stored value close to zero!\n",
      __func__);
    exit(1);
  }
  int iseq = (marker_of_stored_element == marker_of_new_element);
#endif
    
  // return 1 if the elements locality-approximate to the same value
  return(iseq);

}

/// NON-SCALAR
/* Returns 0 if matrix differs from provided content, non-zero if equal */
/* The second argument of this routine mat_val_ptr is a
   panel[4] which is a list giving the addresses of four matrix records.
*/
static int
matrix_compare_panel(mat_add_t m_ptr, mat_add_t panel[4], matrix_type_t mat_type)
{
  // EASY CASES for failed comparison
   if ((m_ptr == MATRIX_PTR_INVALID) || (matrix_type(m_ptr) != mat_type) )
  { 
     return 0; 
  }

  return (panel[0] == matrix_sub(m_ptr,0) &&
          panel[1] == matrix_sub(m_ptr,1) &&
          panel[2] == matrix_sub(m_ptr,2) &&
          panel[3] == matrix_sub(m_ptr,3) ); 
}

/* Returns the matrix PTR which matches scalar, starting from the hash */
static mat_add_t
matrix_find_scalar(uint64_t hash, ScalarType scalar)
{
#ifdef DEBUG 
printf("In routine %s\n", __func__);
#endif

	hash_node_t *n = hash_get_chain(store.hash_table, hash);
	uint32_t depth = 0;

#ifdef HASHSTATS
        (store.hash_table->num_accesses[hash])++;
#endif

	while (n) {
	  if (matrix_compare_scalar((mat_add_t)(n->record_ptr), scalar)) {
                        store.hash_table->hits++;
			n->hits++;

			if (depth > store.deepest) {
				store.deepest = depth;
			}

			//THIS was commmented out before
#if 0
			/* TODO: optimize hash chain order */
			if (n != store.hash[hash]) {
				struct node *head = store.hash[hash];
				if (n->prev) n->prev->next = n->next;
				if (n->next) n->next->prev = n->prev;
				n->prev = NULL;
				n->next = head;
				if (head)
				  head->prev = n;
				store.hash[hash] = n;
			}
#endif

			return (mat_add_t) (n->record_ptr);
		}
		n = n->next;
		depth++;
	} // end while

	store.hash_table->misses++;

	return MATRIX_PTR_INVALID;
}

/* Returns the matrix PTR which matches panel, starting from the hash */
static mat_add_t
matrix_find_panel(uint64_t hash, mat_add_t panel[4], matrix_type_t mat_type)
{
#ifdef DEBUG 
printf("In routine %s\n", __func__);
#endif

	hash_node_t *n = hash_get_chain(store.hash_table, hash);
	uint32_t depth = 0;

#ifdef HASHSTATS
        (store.hash_table->num_accesses[hash])++;
#endif


	while (n) {
	  if (matrix_compare_panel((mat_add_t)(n->record_ptr), panel, mat_type)) {
                        store.hash_table->hits++;
			n->hits++;

			if (depth > store.deepest) {
				store.deepest = depth;
			}

			//THIS was commmented out before
#if 0
			/* TODO: optimize hash chain order */
			if (n != store.hash[hash]) {
				struct node *head = store.hash[hash];
				if (n->prev) n->prev->next = n->next;
				if (n->next) n->next->prev = n->prev;
				n->prev = NULL;
				n->next = head;
				if (head)
				  head->prev = n;
				store.hash[hash] = n;
			}
#endif

			return (mat_add_t) (n->record_ptr);
		}
		n = n->next;
		depth++;
	}

	store.hash_table->misses++;

	return MATRIX_PTR_INVALID;
}

static int64_t record_mat_ptr_by_matrixID( mat_add_t m_ptr)
{
  // This function adds an entry to the table indexing matrix pointer values
  // by matrixID values. In the future, when we are more actively purging
  // matrices from the matrix store, this may be replaced with a hash.
  int64_t list_index;
  int64_t table_index;
  int64_t index = get_matrixID_from_ptr(m_ptr);

  // high order bits give index in the list of tables
  list_index = index >> LOG_SIZE_MPT;

  // low order bits give index in one of the tables
  table_index = index & MASK_MPT;

  // TODO: this test may be redundant if we have the test when creating matrixIDs
  if (list_index >= SIZE_MPT) 
    {
      printf("ERROR: LOG_SIZE_MPT needs to be larger.\n");    // or, can modify code to allocate more tables as needed
      exit(1);
    }   

  // if necessary, allocate memory for the next table in the list
  if (store.mat_ptr_table_list[list_index] == 0)    
    {
      store.mat_ptr_table_list[list_index] = (mat_add_t *)calloc(SIZE_MPT,sizeof(mat_add_t));
    }
  
  // with setup complete, add pointer to the appropriate table
  store.mat_ptr_table_list[list_index][table_index] = m_ptr;
  return(index);
}


/* Inserts a new scalar matrix into the store and returns the new matrix PTR*/
static mat_add_t
matrix_insert_scalar(uint64_t hash, ScalarType scalar, mat_add_t adjoint)
{
  
  int verbose = 0;
  if (verbose) {   // GOT HERE
    printf("Inside %s: levels (0,0), hash %" PRIu64 "\n",
        __func__, hash);
  }
  
  // allocate space for the new LARC matrix
  mat_add_t m_ptr = calloc(1, sizeof(larc_matrix_t));
  if (m_ptr == NULL) {
    ALLOCFAIL();
    return MATRIX_PTR_INVALID;
  }
  
  // set the basic parameters for the matrix
  m_ptr->matrixID = store.matrixID_next++;
  m_ptr->row_level = 0;
  m_ptr->col_level = 0;
  m_ptr->adjoint = adjoint;
  
  if (verbose) {   //  GOT HERE
    printf("  allocated space, set row level, col level, matrixID, and adjoint\n");
  }
  
#ifndef HASH_CHAIN_GROWS_AT_HEAD
  // INSERT new items at tail of hash collision chain
  hash_insert_at_tail(store.hash_table, (record_ptr_t)m_ptr, hash);
#endif
#ifdef HASH_CHAIN_GROWS_AT_HEAD
  // INSERT new items at head of hash collision chain
  hash_insert_at_head(store.hash_table, (record_ptr_t)m_ptr, hash);
#endif

  if (verbose) {   //  GOT HERE
    printf("  appended hash chain with 0, 0levels\n");
  }
  
  // update matrix store parameters
  store.hist[0][0]++;
  m_ptr->lock = 0;
  m_ptr->hold = 0;
  
  if (verbose) {   //  GOT HERE
    printf("  set row level, col level, matrixID, and adjoint\n");
  }
  
  if (verbose) {
#ifdef USE_COMPLEX
  printf("ADDING TO MATRIX SCALAR STORE: %ld -> %.25g + i * %.25g\n", m_ptr->matrixID, creal(scalar), cimag(scalar));
#else
#ifdef USE_REAL
  printf("ADDING TO MATRIX SCALAR STORE: %ld -> %.25g\n", m_ptr->matrixID, scalar);
#else // USE_INTEGER
  printf("ADDING TO MATRIX SCALAR STORE: %ld -> %ld\n", m_ptr->matrixID, scalar);
#endif
#endif
  }
    
  // the trace of a 1x1 matrix is the element inside the matrix
  m_ptr->trace_element = scalar;
    
  // test to see if scalars are zero or one, and/or self-adjoint
  // if a scalar is real, it is self-adjoint, since transpose of a 1x1 matrix
  // and complex conjugation of a real number both leave the argument unchanged
#ifdef USE_COMPLEX  // COMPLEX
  if (fpclassify(cimag(scalar)) == FP_ZERO) {
    m_ptr->iszero = (fpclassify(creal(scalar)) == FP_ZERO);
    m_ptr->isid = (fpclassify(creal(scalar) - 1.0) == FP_ZERO);
    m_ptr->adjoint = m_ptr;
  }
#else // not USE_COMPLEX
#ifdef USE_REAL
  m_ptr->iszero = (fpclassify(scalar) == FP_ZERO);
  m_ptr->isid = (fpclassify(scalar - 1.0) == FP_ZERO);
#endif
#ifdef USE_INTEGER
  m_ptr->iszero = (scalar==0);
  m_ptr->isid = (scalar==1);
#endif
  m_ptr->adjoint = m_ptr;
#endif // not USE_COMPLEX

  store.nscalars++;
  if( store.nscalars > store.max_nscalars ) {
    store.max_nscalars = store.nscalars;
  }

  // Record the matrix pointer in table indexed by matrixID
  record_mat_ptr_by_matrixID(m_ptr);
  
  return m_ptr;
}

/* Inserts a new matrix into the store and returns the new matrix PTR */
static mat_add_t
matrix_insert_panel(uint64_t hash, mat_add_t panel[4], mat_level_t row_level, mat_level_t col_level, mat_add_t adjoint)
{
  // EXCEPTION TESTING
  if (MAX(row_level,col_level) > store.largest_level){
    printf("ERROR: In %s levels %d %d of incoming matrix too large\n",
        __func__,row_level,col_level);
    exit(1);
  }
  
  matrix_type_t mat_type = matrix_type_from_levels(row_level,col_level);
  int verbose = 0;
  if (verbose) {   // GOT HERE
    printf("Inside %s: levels (%d,%d), hash %" PRIu64 "\n",
        __func__, row_level, col_level, hash);
  }
  
  // allocate space for the new LARC matrix
  mat_add_t m_ptr = calloc(1, sizeof(larc_matrix_t));
  if (m_ptr == NULL) {
    ALLOCFAIL();
    return MATRIX_PTR_INVALID;
  }
  
  // set the basic parameters for the matrix
  m_ptr->matrixID = store.matrixID_next++;
  m_ptr->row_level = row_level;
  m_ptr->col_level = col_level;
  m_ptr->adjoint = adjoint;
  
  if (verbose) {   //  GOT HERE
    printf("  allocated space, set row level, col level, matrixID, and adjoint\n");
  }
  
#ifndef HASH_CHAIN_GROWS_AT_HEAD
  // INSERT new items at tail of hash collision chain
  hash_insert_at_tail(store.hash_table, (record_ptr_t)m_ptr, hash);
#endif
#ifdef HASH_CHAIN_GROWS_AT_HEAD
  // INSERT new items at head of hash collision chain
  hash_insert_at_head(store.hash_table, (record_ptr_t)m_ptr, hash);
#endif

  if (verbose) {   //  GOT HERE
    printf("  appended hash chain with %d, %d levels\n",row_level,col_level);
  }
  
  // update matrix store parameters
  store.hist[row_level][col_level]++;
  m_ptr->lock = 0;
  m_ptr->hold = 0;
  
  if (verbose) {   //  GOT HERE
    printf("  set row level, col level, matrixID, and adjoint\n");
  }
  
  if (mat_type == SCALAR) {
    // this shouldn't happen
    printf("in matrix_insert_panel: somehow got SCALAR for mat_type\n");
    exit(1);
  }

  else if (mat_type == MATRIX) {   

    // Exception checking
    if( (matrix_row_level(panel[0]) != matrix_row_level(panel[1])) || 
	(matrix_row_level(panel[0]) != matrix_row_level(panel[2])) ||
	(matrix_row_level(panel[0]) != matrix_row_level(panel[3])) )
      {
	printf( "ERROR: In %s trying to insert matrix with non-equal row sub-levels %d %d %d %d\n", __func__,
		matrix_row_level(panel[0]), matrix_row_level(panel[1]),
		matrix_row_level(panel[2]), matrix_row_level(panel[3]) );
	exit(1);
      }

    // define new matrix in terms of its four panels
    for (int i = 0; i < 4; i++) {
      m_ptr->submatrix[i] = panel[i];
      if (panel[i] == NULL) {
	printf("Error looking up submatrices\n");
	exit(1);
	// return MATRIX_PTR_INVALID;
      }
      // each panel is a submatrix of the new matrix
      matrix_appears_as_sub_count_increment(panel[i]);
    }

    // special SQUARE MATRIX stuff
    if (row_level == col_level) {
      // trace of matrix is sum of traces of the two diagonal blocks
      m_ptr->trace_element = panel[0]->trace_element + panel[3]->trace_element;
    
      // testing self-adjointness (fails if adjoint field is MATRIX_PTR_INVALID)
      // this should be sufficient to cause the initially loaded square matrices
      // (zeros, identities and integer Hadamards) to be seen as self-adjoint
      if ((panel[1]->adjoint == panel[2]) && (panel[0]->adjoint == panel[0])
           && (panel[3]->adjoint == panel[3])) { m_ptr->adjoint = m_ptr; }
    } // end square matrix tests

    // testing for zero or identity matrices
    if (matrix_is_zero(panel[1]) && matrix_is_zero(panel[2])) {
      if (matrix_is_zero(panel[0]) && matrix_is_zero(panel[3])) {
	  // [0, 0; 0, 0] 
          m_ptr->iszero = 1;
      } else if (matrix_is_id(panel[0]) && matrix_is_id(panel[3])) {
          // [I, 0; 0, I]; only possible for square matrices
          m_ptr->isid = 1;
      }
    }
  }  // end MATRIX type

  else if (mat_type == COL_VECTOR) {
    
    m_ptr->submatrix[0] = panel[0];
    m_ptr->submatrix[2] = panel[2];
    if ((panel[0] == NULL) || (panel[2] == NULL)) {
      printf("Error looking up submatrices\n");
      exit(1);
      return MATRIX_PTR_INVALID;
    }
    matrix_appears_as_sub_count_increment(panel[0]);
    matrix_appears_as_sub_count_increment(panel[2]);
    
    if (matrix_is_zero(panel[0]) && matrix_is_zero(panel[2])) {
      m_ptr->iszero = 1;
    }

  }  // end of COL_VECTOR type
 
  else if (mat_type == ROW_VECTOR) {


    m_ptr->submatrix[0] = panel[0];
    m_ptr->submatrix[1] = panel[1];
    if ((panel[0] == NULL) || (panel[1] == NULL)) {
      printf("Error looking up submatrices\n");
      exit(1);
      return MATRIX_PTR_INVALID;
    }
    matrix_appears_as_sub_count_increment(panel[0]);
    matrix_appears_as_sub_count_increment(panel[1]);

    if (matrix_is_zero(panel[0]) && matrix_is_zero(panel[1])) {
      m_ptr->iszero = 1;
    }

  }  // end of ROW_VECTOR type

  store.nmatrices++;
  if( store.nmatrices > store.max_nmatrices ) {
      store.max_nmatrices = store.nmatrices;
  }
  
  // Record the matrix pointer in table indexed by matrixID
  record_mat_ptr_by_matrixID(m_ptr);
  
  return m_ptr;
}


/* Calculates the adjoint of a matrix and finds or inserts it into the store */
static mat_add_t
matrix_find_or_insert_adjoint(mat_add_t m_ptr)
{
  // EXCEPTION CHECKING
  if matrix_is_invalid(m_ptr) {
    return MATRIX_PTR_INVALID;
  }

  // SHORTCUT: CHECK TO SEE IF WE ALREADY KNOW THIS MATRIX'S ADJOINT
  // (note that self-adjoint cases not in operations store!)
  if (m_ptr->adjoint != MATRIX_PTR_INVALID) return m_ptr->adjoint;
  
  // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE
  mat_add_t adj_ptr = op_get(ADJOINT,m_ptr,m_ptr);

  if (adj_ptr != MATRIX_PTR_INVALID) {
    // we can get rid of this test once we're confident it all works
    printf("in %s: m_ptr->adjoint not set, yet",__func__);
    printf(" we found the adj_ptr in the operations store\n");
  }

  if (adj_ptr == MATRIX_PTR_INVALID) {
    matrix_type_t mat_type = matrix_type(m_ptr);
    matrix_type_t adjoint_type;
    if (mat_type == ROW_VECTOR) {adjoint_type = COL_VECTOR;}
    else if (mat_type == COL_VECTOR) {adjoint_type = ROW_VECTOR;}
    else  {adjoint_type = mat_type;}
  
    uint64_t hash;
    if (matrix_type(m_ptr) == SCALAR) {
#ifdef USE_COMPLEX
      ScalarType scalar;
      // adjoint of complex number might not be in store
      scalar = conj(matrix_trace(m_ptr));
      hash = hash_from_matrix_scalar(scalar, adjoint_type, store.hash_table->exponent);
      adj_ptr = matrix_find_scalar(hash, scalar);
      if matrix_is_invalid(adj_ptr) {
        adj_ptr = matrix_insert_scalar(hash, scalar, m_ptr);
      }
#else
      // adjoint of real number is that number, guaranteed to be in store
      // in fact, m_ptr->adjoint should have already been set...
      // the following code should never be accessed
      printf("in %s: adj_ptr not set for a real scalar value\n",__func__);
      adj_ptr = m_ptr;
#endif
    } else {
      mat_add_t panel[4];
      panel[0] = matrix_adjoint(matrix_sub(m_ptr,0));
      panel[1] = matrix_adjoint(matrix_sub(m_ptr,2));
      panel[2] = matrix_adjoint(matrix_sub(m_ptr,1));
      panel[3] = matrix_adjoint(matrix_sub(m_ptr,3));
      hash = hash_from_matrix_panel(panel, adjoint_type, store.hash_table->exponent);
      adj_ptr = matrix_find_panel(hash, panel, adjoint_type);
      if matrix_is_invalid(adj_ptr) {
        adj_ptr = matrix_insert_panel(hash, panel, 
                  matrix_col_level(m_ptr), matrix_row_level(m_ptr), m_ptr);
      }
    }

    // STORE RESULT IN OPERATIONS STORE
    op_set(ADJOINT, m_ptr, m_ptr, adj_ptr);
  }  // end matrix was not found in operations store

  // tell the matrices they are adjoints of each other
  m_ptr->adjoint = adj_ptr;
  adj_ptr->adjoint = m_ptr;

  return adj_ptr;
}

// These next two routines are used in the python interface
// because we only let python users see matrixIDs
int64_t
matrix_get_matrixID_from_panel(int64_t A_mID, int64_t B_mID, int64_t C_mID, int64_t D_mID, 
                        mat_level_t row_level, mat_level_t col_level) 
{

  // grab the matrix pointers corresponding to each matrixID, and
  // check that these matrices are still in the matrix store
  mat_add_t  panel[4];
  panel[0] = mat_ptr_from_matrixID(A_mID, "first", __func__,0);
  panel[1] = mat_ptr_from_matrixID(B_mID, "second", __func__,0);
  panel[2] = mat_ptr_from_matrixID(C_mID, "third", __func__,0);
  panel[3] = mat_ptr_from_matrixID(D_mID, "fourth", __func__,0);

  // check to see if these matrix pointers are still in the store 
  for (int i=0;i<4;++i) {
    if (panel[i] == MATRIX_PTR_INVALID) { exit(1); }
  }

  mat_add_t m_ptr = matrix_get_ptr_panel(panel,row_level,col_level); 
  return m_ptr->matrixID;
}

// python interface for scalar type
int64_t 
matrix_get_matrixID_from_scalar(ScalarType val)
{
  mat_add_t m_ptr = matrix_get_ptr_scalar(val);
  return m_ptr->matrixID;
}

/* The first argument of this routine is a ScalarType (complex, double, int) */
mat_add_t matrix_get_ptr_scalar(ScalarType scalar)
{
  int verbose = 0;
  if (verbose) { 
    printf("Inside %s:\n", __func__);
    printf("  (0,0) level matrix \n");
  }

  if (verbose) { 
    printf("  About to calculate the matrix_type, hash, and see if matrix is stored\n");
  }

  matrix_type_t mat_type = SCALAR;
  uint64_t hash = hash_from_matrix_scalar(scalar, mat_type, store.hash_table->exponent);
  mat_add_t m_ptr = matrix_find_scalar(hash, scalar);
  if (verbose) { 
    printf("  About to check if the matrix was not found, m_ptr is %p\n",m_ptr);
  }
  
  // if unable to find attempt to insert
  if (m_ptr == MATRIX_PTR_INVALID) {
    if (verbose) { 
          printf("Inside %s, matrix index not found by matrix_find\n",__func__);
          printf("   CALLING  matrix_insert\n");
    }
    m_ptr = matrix_insert_scalar(hash, scalar, MATRIX_PTR_INVALID);
  }
  if (verbose) {   
    printf("    RETURNED from matrix_insert with matrix index %p\n",m_ptr);
  }
  
  // if couldn't insert, then bail out
  if (m_ptr == MATRIX_PTR_INVALID) {
    abort();
  }
  if (verbose) {
    printf("Inside %s\n", __func__);
  }
  
  return m_ptr;
}

/* The first argument of this routine mat_val_ptr is an address (void *) and
   the element at that address is either a ScalarType (complex, double, int)
   or panel[4] which is a list giving the addresses of four matrix records.
*/
mat_add_t
matrix_get_ptr_panel(mat_add_t panel[4], mat_level_t row_level, mat_level_t col_level)
{
  int verbose = 0;
  if (verbose) { 
    printf("Inside %s:\n", __func__);
    printf("  (%d,%d) level matrix \n",row_level,col_level);
  }

  if (verbose) { 
    printf("  About to calculate the matrix_type, hash, and see if matrix is stored\n");
  }

  matrix_type_t mat_type = matrix_type_from_levels(row_level,col_level);
  uint64_t hash = hash_from_matrix_panel(panel, mat_type, store.hash_table->exponent);
  mat_add_t m_ptr = matrix_find_panel(hash, panel, mat_type);
  if (verbose) { 
    printf("  About to check if the matrix was not found, m_ptr is %p\n",m_ptr);
  }
  
  // if unable to find attempt to insert
  if (m_ptr == MATRIX_PTR_INVALID) {
    if (verbose) { 
          printf("Inside %s, matrix index not found by matrix_find\n",__func__);
          printf("   CALLING  matrix_insert\n");
    }
    m_ptr = matrix_insert_panel(hash, panel, row_level, col_level, MATRIX_PTR_INVALID);
  }
  if (verbose) {   
    printf("    RETURNED from matrix_insert with matrix index %p\n",m_ptr);
  }
  
  // if couldn't insert, then bail out
  if (m_ptr == MATRIX_PTR_INVALID) {
    abort();
  }
  if (verbose) {
    printf("Inside %s\n", __func__);
  }
  
  return m_ptr;
}

size_t
matrix_store_scalarCount(void)
{
	return store.nscalars;
}


size_t
matrix_store_matrixCount(void)
{
	return store.nmatrices;
}

int64_t
matrix_adjoint_matrixID(int64_t m_mID)
{
  // get the matrix pointers from the matrixIDs, and see if still in store 
  mat_add_t m_ptr = mat_ptr_from_matrixID(m_mID, "first", __func__,0);
  if (m_ptr == MATRIX_PTR_INVALID) { exit(1); }

  // call matrix pointer version of function
  mat_add_t adjoint_ptr = matrix_adjoint(m_ptr);
  int64_t adjoint_mID = get_matrixID_from_ptr(adjoint_ptr);
  return adjoint_mID;
}

mat_add_t
matrix_adjoint(mat_add_t m_ptr)
{
	if (m_ptr && m_ptr->adjoint != MATRIX_PTR_INVALID) {
		return m_ptr->adjoint;
	}
	return matrix_find_or_insert_adjoint(m_ptr);
}


/* The total number of matrices created and stored in  */
/* the matrix store (the current number of stored     */
/* may be less if there was matrix store cleaning).     */
uint64_t  
num_matrices_created(void)
{
  uint64_t ret = store.matrixID_next;
  return(ret);
}


// This routine takes a matrixID and either returns the address of the
// matrix labeled by that number or MATRIX_PTR_INVALID. The routine is also
// passed two strings: calling_routine, which lets the user know which routine
// passed the bad matrixID; and arg_no (eg "second") to inform the user
// which argument to that routine was bad. 
// print_flag = 1 to print a warning that there is no matrix with this matrixID
// otherwise, print_flag = 0 to suppress warning
// The routine is largely called from C routines which act as wrappers for
// python calls, which don't always behave well when passing pointers around.
// It can be used to confirm that a particular matrixID is valid.
//

mat_add_t
mat_ptr_from_matrixID(int64_t matrixID, const char* arg_no, 
        const char* calling_routine, int print_flag) {

  // grab highest matrixID
  uint64_t highest_matrixID = num_matrices_created()-1;

  // check bounds
  if ((matrixID > highest_matrixID) || (matrixID < 0) ) {
    printf("ERROR: In %s, the %s argument, matrixID %" PRId64 ", is out of range.\n",
        calling_routine, arg_no, matrixID);
    return(MATRIX_PTR_INVALID);
  }

  // get matrix address
  int64_t list_index = matrixID >> LOG_SIZE_MPT;
  int64_t table_index = matrixID & MASK_MPT;
  mat_add_t mat_val_ptr = store.mat_ptr_table_list[list_index][table_index];

  // warn the user if the matrixID doesn't refer to a valid matrix
  if ((mat_val_ptr == MATRIX_PTR_INVALID) && (print_flag==1)) {
    printf("WARNING: In %s, problem with the %s argument:",
        calling_routine, arg_no);
    printf(" matrix with matrixID %" PRId64 " has been removed from the Matrix Store.\n",
         matrixID);
  }

  return mat_val_ptr;
}

static int  
matrix_info_to_file(uint64_t i, FILE *f) {
  // EXCEPTIONS 
  if (f==NULL) { 
    printf("WARN: %s called with bad file pointer\n", __func__); 
    return(0); 
  } 

  // get the matrix pointers from the matrixID, and see if still in store 
  mat_add_t m_ptr = mat_ptr_from_matrixID(i, "first", __func__,0);

  if (matrix_is_invalid(m_ptr)) { return(1); }  // skip printing line

  // HEADER: "matrixID   Levels    Value    Appears_As_Sub Count   Lock   Hold    Adjoint matrixID"
  //   matrixID 
  fprintf(f,"      %zd",get_matrixID_from_ptr(m_ptr));
  //   Levels
  fprintf(f,"    \t(%2d,%2d)",matrix_row_level(m_ptr),matrix_col_level(m_ptr)); 
  //   Value 
  if ((matrix_row_level(m_ptr) == 0) && (matrix_col_level(m_ptr) == 0)) {  // SCALAR

#ifdef USE_INTEGER
    fprintf(f," \t%ld                   ",
              matrix_trace(m_ptr));
#endif
#ifdef USE_REAL  
    fprintf(f," \t%g                   ",
              matrix_trace(m_ptr));
#endif
#ifdef USE_COMPLEX
    fprintf(f," \t%g + %g I            ",
             creal(matrix_trace(m_ptr)),cimag(matrix_trace(m_ptr)));
#endif
  } 
  else {   // MATRIX
    mat_add_t panel[4];
    fprintf(f,"   [ ");
    for (int i = 0; i < 4; i++) {
      panel[i] = matrix_sub(m_ptr,i); 
      if (panel[i] == MATRIX_PTR_INVALID) {
	fprintf(f,"  X ");
      }
      else {
        fprintf(f,"%" PRId64 " ", get_matrixID_from_ptr(panel[i])); 
      }
    } 
    fprintf(f,"]        ");
  }
  //   Appears_As_Sub Count = Reference Count 
  fprintf(f," \t%8d     %4d    %4d     ", matrix_appears_as_sub_count(m_ptr), matrix_has_lock(m_ptr), matrix_has_hold(m_ptr));
  //   matrixID of Adjoint
  mat_add_t adjoint_field = matrix_adjoint_field(m_ptr);
  if (adjoint_field == MATRIX_PTR_INVALID) {
    fprintf(f," \tX ");
  }
  else {
    fprintf(f," \t%" PRId64 " ", get_matrixID_from_ptr(adjoint_field)); 
  }
  //   End of this line 
  fprintf(f,"\n");

  return(1);

}


/***************************************************************************
*                      matrix_store_info_to_file                           *
*  This function prints information about matrices from matrixID      *
*  "start" to "end"  (0 <= start <= end < SIZE_MPT^2
*  The output file path and a user comment to be printed in the file are   *
*  also arguments.                                                         *
***************************************************************************/
int
matrix_store_info_to_file(uint64_t start, uint64_t end,
                          char *outfilepath, char *comment)
{
  if (start > end) {
    printf("ERROR in matrix_store_info_to_file: impossible range %" PRIu64 " > %" PRIu64 "\n", start, end);
    return(0);
  }

  printf("Printing matrices with matrixIDs from %" PRIu64 " to %" PRIu64 " to file %s\n",
          start,end,outfilepath);
  
  int ret = 1;
  FILE *f = fopen(outfilepath, "w"); 
  fprintf(f,"Comment: %s\n",comment); 
  fprintf(f,"Total number of matrices that have been created is:           %zd\n", 
  	  num_matrices_created()); 

  int64_t max_num_matrices = ((int64_t)1)<<(2*LOG_SIZE_MPT);

  /* // fprintf(f,"Largest number of matrices ever stored simultaneously:        %zd\n",  */
  /* //	  store.max_nmatrices); */
  /* // fprintf(f,"Largest number of scalars ever stored simultaneously:         %zd\n",  */
  /* //	  store.max_nscalars); */
  if ((end > max_num_matrices) || (end >= num_matrices_created()) ) {
  /*   // can only print as many matrices as have been created and indexed to print */
     end = MIN(num_matrices_created()-1,max_num_matrices); 
     fprintf(f,"SORRY: currently end matrixID can't be any bigger than %" PRIu64 "\n",end); 
  } 
  fprintf(f,"Printing matrix store entries with matrixIDs from %" PRIu64 " to %" PRIu64 "\n",start,end);
  fprintf(f,"\n\n==========================(Matrix Store)=============================\n");
  fprintf(f,"matrixIDs       Levels           Value   ");
  fprintf(f,"      Appears_As_SubCount    Lock    Hold     AdjointSN\n"); 
  uint64_t i;
     for (i = start; i <= end; ++i) { 
       // fprintf(f,"   --- call to matrix_info_to_file(%d,f) ---\n",i);
       ret = matrix_info_to_file(i,f);
     }
  fprintf(f,"\n======================================================================\n");

  fclose(f); 

  return ret;
} 


/***************************************************************************
*                   matrix_hash_chain_info_to_file                         *
*  This function prints information about matrices to a file,			   *
*  given the hash value for that chain.                                    *
*  To get the appropriate hash value for the argument of this function     *
*     hashID = matrix_hashID_from_matrixID(matrixID);				   *
*  The output file path and a user comment to be printed in the file are   *
*  also arguments.                                                         *
***************************************************************************/
int
matrix_hash_chain_info_to_file(uint64_t hash, char *outfilepath, char *comment)
{
  FILE *f; 

  if (strcmp(outfilepath,"stdout")) {
     printf("Printing matrix hash chain with hash value %ld to file %s\n",
          hash,outfilepath);
     f = fopen(outfilepath, "w"); 

  }
  else {
    printf("Printing matrix hash chain with hash value %ld to screen\n",hash);
    f = stdout;
  }
  
  int ret = 1;
  fprintf(f,"Comment: %s\n",comment); 
  fprintf(f,"This is the matrix hash chain info for hash value %ld\n", hash); 
  fprintf(f,"\n\n===================(Matrix Hash Chain)=============================\n");
  fprintf(f,"MatrixID        Levels           Value   ");
  fprintf(f,"        Appears_As_SubCount  AdjointSN\n"); 

  hash_table_t *table_ptr = store.hash_table;
  hash_node_t *node_ptr = table_ptr->heads[hash];
  mat_add_t record_ptr;

  while (node_ptr)  {
    record_ptr = (mat_add_t) node_ptr->record_ptr;
    uint64_t matrixID = record_ptr->matrixID;
    ret = matrix_info_to_file(matrixID,f);
    node_ptr = node_ptr->next;
  }
  fprintf(f,"\n======================================================================\n");

  if  (strcmp(outfilepath,"stdout")) fclose(f);   // if not stdout then close

  return ret;
} 


/***************************************************************************
*                   matrix_hash_chain_info_to_screen                       *
*  This function prints information about matrices given the hash value    *
*  for that chain.                                                         *
*  To get the appropriate hash value for the argument of this function     *
*     hashID = matrix_hashID_from_matrixID(matrixID);				   *
*  A user comment to be printed is also an argument   					   *
***************************************************************************/
int
matrix_hash_chain_info_to_screen(uint64_t hash, char *comment)
{
  char *filename = "stdout";
  return matrix_hash_chain_info_to_file(hash, filename, comment);
} 


// python interface version
int remove_matrix_from_mat_store_by_matrixID (int64_t  m_mID) {

  // get the matrix pointers from the matrixIDs, and see if still in store 
  mat_add_t m_ptr = mat_ptr_from_matrixID(m_mID, "first", __func__,0);
  
  // Python users could try to remove the same thing twice, but
  // we chose to not kill the program if they do this, just warn them.
  if (m_ptr == MATRIX_PTR_INVALID) {
    printf("WARNING: In %s, matrix %" PRId64 " was formerly removed from the Matrix Store.\n", __func__, m_mID);
    return(0);
  }

  return remove_matrix_from_mat_store (m_ptr);
}



/******************************************************************
*          remove_matrix_from_mat_store (mat_add_t  m_ptr)
*  This function attempts to remove a larc_matrix_t structure 
*  from the matrix store.
*     1 is returned if the remove succeeded
*     0 is returned if the matrix could not be removed
*       because it has larc_appears_as_subs or has a lock or hold.
*******************************************************************/
int remove_matrix_from_mat_store (mat_add_t  m_ptr)
{
  // Exception checking
  if (matrix_is_invalid(m_ptr)) {
    printf("WARNING: remove_matrix_from_mat_store called for invalid matrix pointer\n");
    return 0;
  }
 
  // Check if matrix has appears_as_subs or has a lock or hold
  if (matrix_appears_as_sub_count(m_ptr) != 0) return 0;
  if (matrix_has_lock(m_ptr) == 1) return 0;
  if (matrix_has_hold(m_ptr) == 1) return 0;
  
  // If the matrix had an adjoint, then need to remove the adjoint's adjoint
  mat_add_t adj_ptr = matrix_adjoint_field(m_ptr);
  if (adj_ptr != MATRIX_PTR_INVALID) {
     matrix_adjoint_field(adj_ptr) = MATRIX_PTR_INVALID;
  }

  // If the matrix is not a SCALAR, then recursively remove it
  //   for each child of matrix, decrement appears_as_sub_count 
  //   and if child is now at zero counts  attempt to remove it
  if (matrix_type(m_ptr) != SCALAR) {
    for (int i=0;i<4;++i) {
      mat_add_t child_ptr = matrix_sub(m_ptr, i);
      if (child_ptr != MATRIX_PTR_INVALID) {
	int new_appears_as_sub_count = matrix_appears_as_sub_count_decrement(child_ptr);
	if (new_appears_as_sub_count == 0) remove_matrix_from_mat_store (child_ptr);   
      }
    }
  }

  // Remove the matrix from the hash table for the matrix store

  // Find the hash chain where the matrix is stored
  uint64_t hash;
  matrix_type_t mat_type = matrix_type(m_ptr);
  if (mat_type == SCALAR) {
    // ScalarType *mat_val_ptr;
    // mat_val_ptr = &(matrix_trace(m_ptr));
    // hash = hash_from_matrix_content(mat_val_ptr, SCALAR, store.hash_table->exponent);
    hash = hash_from_matrix_scalar(matrix_trace(m_ptr), SCALAR, store.hash_table->exponent);
  } else {
    // mat_add_t *mat_val_ptr;  // array of four pointers to matrix records
    // mat_val_ptr = m_ptr->submatrix;
    // hash = hash_from_matrix_content(mat_val_ptr, mat_type, store.hash_table->exponent);
    hash = hash_from_matrix_panel(m_ptr->submatrix, mat_type, store.hash_table->exponent);
  }
  hash_node_remove(store.hash_table, (record_ptr_t)m_ptr, hash);
  
  
  // Decrement the histogram and counts of matrices and scalars
  store.hist[matrix_row_level(m_ptr)][matrix_col_level(m_ptr)]--;
  if (mat_type == SCALAR) {
    store.nscalars--;
#if 0
    if( store.nscalars > store.max_nscalars ) {
      printf("WARNING scalars: We do not expect this to happen in %s\n",__func__);
      store.max_nscalars = store.nscalars;
    }
#endif
  } else {
    store.nmatrices--;
#if 0
    if( store.nmatrices > store.max_nmatrices ) {
      printf("WARNING nonscalars: We do not expect this to happen in %s\n",__func__);
      store.max_nmatrices = store.nmatrices;
    }
#endif
  }

  // Remove the matrix pointer from the table indexed by matrixIDs
  table_mat_ptr_by_matrixID_Remove_entry(m_ptr);

  // free the structure, note that no mallocs occur in the matrix
  free(m_ptr);

  return (1);
}


// Remove the matrix pointer from the table indexed by matrixIDs
// by replacing the value with MATRIX_PTR_INVALID
int
table_mat_ptr_by_matrixID_Remove_entry(mat_add_t m_ptr){

  // EXCEPTION CHECKING
  if (m_ptr == MATRIX_PTR_INVALID) { 
    printf("ERROR: in %s\n", __func__); 
    printf("   asked to remove m_idx =MATRIX_PTR_INVALID\n");
    exit(1);
  }

  int64_t index = get_matrixID_from_ptr(m_ptr);
     
  if (index != MATRIX_ID_INVALID) {
    int64_t list_index = index >> LOG_SIZE_MPT;
    int64_t table_index = index & MASK_MPT;
    store.mat_ptr_table_list[list_index][table_index] =  MATRIX_PTR_INVALID;
  }

  return(1);

}


// Locality Approximation Functions

#ifndef USE_INTEGER
// This relies on IEEE 754
// Rounds a single real number to the specified number 
// of significant figures
/*
IEEE 754 double representation
b = base = 2
q = exponent
s = sign bit, 0 for positive and 1 for negative
c = significand = coefficient
number = (-1)^s * c *  b^(q)  
For double there is "bias" = 1023 = 2^10-1
this allows the exponent to be represented
by a number between 0 and 2048 = 2^11 "exp_rep"
so that actually exponent is exp_rep - bias.
DBL_MANT_DIG is the number of bits in the 
mantissa in IEEE 754, Laurie and Jenny think
this is 53 including the implicit bit.
There is a implicit leading bit in the 
mantissa that is always 1
but is invisible), i.e. mantissa = 1.foo   
where foo is the fractional part
of the mantissa.
UINT64_MAX is the all 1's 64-bit quantity
We are guessing the order of the bits
First 1 bit is the sign bit
Next 11 bits are the bias+exp, 
and the last 52 bits are the fractional part
of the mantissa.
 */
static double
round_sig_fig_real(double input_real, int sig)
{
  //  ERROR CHECKING
  if (sig <= 0) {
    printf("Unwise use of %s keeping no bits (or neg bits) of mantissa\n",
        __func__);
    exit(1);
  }
  // since the default sig == DBL_MANT_DIG, the warning has been moved to 
  // initialize_larc, and only prints for sig > DBL_MANT_DIG
  if (sig >= DBL_MANT_DIG) {
    //printf("Questionable use of %s keeping all bits of mantissa\n", __func__);
    return(input_real);
  }
  if (isnan(input_real) ) {
    printf("Inside %s a NaN showed up\n", __func__);
    exit (1);
  }

  // a string 64 1's (1 sign bit, 11 exp+bias, 52 fractional mantissa)
  // DBL_MANT_DIG = 53
  // e.g. sig 5,  64 1's, shift left by (53-5)  to get 16 1's followed by 48 0's
  //       the "implicit bit" is 1 bit of significance and we want 4 more 
  //       fractional bits, to make sig = 5.
  uint64_t mask = UINT64_MAX << (DBL_MANT_DIG - sig);  
  //       mask_next is 1 bit past the ones we are keeping, to determine rounding
  uint64_t mask_next = 1L << (DBL_MANT_DIG - (sig+1));

  uint64_t our_int;

  // if the input is negative, then make it positive before approximating
  uint64_t negative = (input_real < 0.0);
  double abs_input = (negative) ? -input_real : input_real;

  // interpret the floating point as an integer and grab all the bits
  uint64_t * ptr_to_fp = (uint64_t *) &abs_input;
  our_int = *ptr_to_fp;

  // check to see if the bit following our last bit of accuracy is a 1
  uint64_t roundedBits = our_int & mask;
  if (our_int & mask_next)
    roundedBits += mask_next << 1;
  double * ptr_to_double = (double *) &roundedBits;
  double input_real_return = *ptr_to_double;
  if (input_real_return < 0.0) {
    printf("In %s there was a ripple carry all the way to sign bit\n",__func__);
    exit(1);   // can't represent a number this large in floating point
  }

#if DEBUG
  // If you are using this debug section you should consider printing abs_input
  if (roundedBits != our_int)
    printf("fp: %+17.8la (%g) > %+17.8la %lx (%g) > %lx (%lx) (%lx)\n", input_real, input_real, input_real_return, input_real_return,
	   our_int,roundedBits, mask, mask_next);
#endif

  // restore sign of input after approximation completed
  input_real_return =  (negative)? -input_real_return : input_real_return;
  return input_real_return;
}

// Rounds a single real number to 0 if within threshold of zero
static double
collapse_near_zero_real(double input_real, double threshold)
{
  if (fabs(input_real) < threshold) {
    input_real = 0.0;
  }
  return input_real;
}
#endif


#ifdef USE_COMPLEX
// Rounds each part (re/im) of complex input to no. of sig figs specified in larc.h
static complex
round_sig_fig_complex(complex scalar, int sig)
{
  return (round_sig_fig_real(creal(scalar),sig)
        + I*round_sig_fig_real(cimag(scalar),sig));
}

// Rounds each part (re/im) of complex input to exactly zero
// whenever that part is within threshold of zero
static complex
collapse_near_zero_complex(complex input, double threshold)
{
  return collapse_near_zero_real(creal(input),threshold) +
        I*collapse_near_zero_real(cimag(input),threshold);
}
#endif

ScalarType locality_approx(ScalarType scalar)
{
  // We do three things:  set marker = scalar
  //    1. if fabs(marker) < zerorealthresh then set marker = 0
  //    2. marker = truncate marker to sighash significant bits
  //    3. if fabs(marker) < zerorealthresh then set marker = 0

#ifdef USE_INTEGER
  return scalar;
#else
  double zerorealthresh = get_zerorealthresh();
  int sighash = get_sighash();

  // collapse value to zero if it is small enough   
  // round to sighash significant bits
  // collapse value to zero if it is small enough   
  
#ifdef USE_COMPLEX
  complex marker = collapse_near_zero_complex(scalar,zerorealthresh);
  marker = round_sig_fig_complex(marker,sighash);
  marker = collapse_near_zero_complex(marker,zerorealthresh);
#endif
  
#ifdef USE_REAL
  double marker = collapse_near_zero_real(scalar,zerorealthresh);
  marker = round_sig_fig_real(marker,sighash);
  marker = collapse_near_zero_real(marker,zerorealthresh);
#endif
  
  return (marker);
#endif
}


/*****************************************************************
 *                matrix_hashID_from_matrixID                      *
 * If this function succeeds it will return the hash value       *
 *  associated with hashing the matrix pointer of the            *
 *  the given matrixID
 *  This function returns either a -1 if it fails because:       *
 *    - the matrixID is out of range, or                    *
 *    - the matrix associated with the matrixID has been    *
 *      removed from the matrix store                            *
 *  Normally a hash is a uint64_t, but if there is a fail,       *
 *  the hashID can be -1, so the return value is int64_t.        *
 *****************************************************************/
int64_t matrix_hashID_from_matrixID(int64_t m_mID)
{

  // get the matrix pointers from the matrixID, and see if still in store 
  mat_add_t m_ptr = mat_ptr_from_matrixID(m_mID, "first", __func__,0);
  
  if (m_ptr == MATRIX_PTR_INVALID) { return(-1); }

  // go get the hash value that corresponds to the matrix with that matrixID
  matrix_type_t mat_type = matrix_type(m_ptr);
  uint64_t hash;
  if (mat_type == SCALAR) {
    hash = hash_from_matrix_scalar(matrix_trace(m_ptr), SCALAR, store.hash_table->exponent);
  } else {
    hash = hash_from_matrix_panel(m_ptr->submatrix, mat_type, store.hash_table->exponent);
  }
   
  return hash;
}

/************************************************************************
 * Clean matrix store	                                                *
 *                                                                      *
 *  Removes all eligible matrices from the matrix store.                *
 *      The called function, remove_matrix_from_mat_store, will         *
 *      recursively remove eligible children of any deleted matrix      *      	
 ***********************************************************************/
 int clean_matrix_store()
 {
	 uint64_t max = 1<<(store.hash_table->exponent);
	 uint64_t i;
	 
	 // loop though hash chains for each hashID without calling
	 // clean_op_hash_chain so we don't have to repeat checks for invalid hashID
     for (i = 0; i < max; ++i) {
		 hash_table_t *table_ptr = store.hash_table;
		 hash_node_t *node_ptr = table_ptr->heads[i];
		 mat_add_t record_ptr;
		 
		 while (node_ptr)  {
			 record_ptr = (mat_add_t) node_ptr->record_ptr;
			 remove_matrix_from_mat_store(record_ptr);
			 node_ptr = node_ptr->next;
			 }
		 }
	return (1);
}

 
 
 /***********************************************************************
 *  Clean hash chain                                                    *
 *  Input:  a hash ID for the hash chain to clean                       *
 *      removes all eligible matrices from a hash chain.                *
 *      The called function, remove_matrix_from_mat_store, will         *
 *      recursively remove eligible children of any deleted matrix      *
 * **********************************************************************/
 int clean_matrix_hash_chain(uint64_t hashID)
 {
	 //Is the size of hash table stored somewhere?
	 uint64_t max = 1<<(store.hash_table->exponent);
	 
	 if ((hashID < 0) || (hashID >= max)) {
		 printf("Error: hash value out of range\n");
		 return(0);
	}
	hash_table_t *table_ptr = store.hash_table;
	hash_node_t *node_ptr = table_ptr->heads[hashID];
	mat_add_t record_ptr; 	
	 
	 while (node_ptr)  {
		 record_ptr = (mat_add_t) node_ptr->record_ptr;
		 remove_matrix_from_mat_store(record_ptr);
		 node_ptr = node_ptr->next;
	 }
	 return (1);
 
}




int64_t get_diag_hist(uint64_t i){
  int64_t diag_hist;

  if ((i < store.largest_level+1) && (i >= 0))
    {
      // allocate space
      // diag_hist = (int64_t *) calloc(store.largest_level+1,sizeof(int64_t));
      // // warning this may cause a memory leak
      diag_hist = (int64_t)store.hist[i][i];
      return(diag_hist);
    }
  else return(-1);

}


int64_t
matrix_hash_chain_length(uint64_t i)
{
  
  int counter = 0;
  
  if ((i < store.hash_table->nentries) && (i >= 0))
    {
      hash_table_t *table_ptr = store.hash_table;
      hash_node_t *node_ptr = table_ptr->heads[i];
      // mat_add_t record_ptr;
      
      while (node_ptr)  {
        ++counter;
	// record_ptr = (mat_add_t) node_ptr->record_ptr;
	// uint64_t matrixID =  record_ptr->matrixID;
	// ret = matrix_info_to_file(matrixID,f);
	node_ptr = node_ptr->next;
      }

      return counter;
    }  
  else return(-1);
  
  
} 


#ifdef HASHSTATS
/************************************************************************
 *  This function calls hashstats_to_file with the matrix_store         *
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
void matrix_hashstats(char *accesses_file,  char *nodes_file, char *report_file)
{
  hashstats_to_files(store.hash_table, accesses_file, nodes_file, report_file);
}
#endif
