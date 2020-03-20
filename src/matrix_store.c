//                       matrix_store.c 
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
#include "larc.h"
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
  mat_ptr_t **zero;     // doubly indexed matrix IDs for the preloaded zero matrices
  mat_ptr_t *identity;  // array of matrix IDs for the preloaded identity matrices (square)
  mat_ptr_t *iHadamard; // array of matrix IDs for the preloaded integer Hadamard matrices

  // POSSIBLY CAN BE REMOVED                      
  uint32_t deepest;    // longest successful traversal down a hash chain to retrieve a matrix 

  // this is a list of tables each entry of each table contains a pointer to a matrix record
  // the index into the table is the matrixID, broken into lower and higher bits
  mat_ptr_t **mat_ptr_table_list;

} store = {0};   // this creates the variable "store" which is a struct of type matrix_store_t


// former inline functions - were in matrix_store.h

int matrix_is_zero_matrixID(int64_t mat_mID) {
  mat_ptr_t mat_mPTR = get_matPTR_from_matID(mat_mID, "", __func__, 0);
  return(matrix_is_zero(mat_mPTR));
}

int matrix_is_invalid_matrixID(int64_t mat_mID) {
  mat_ptr_t mat_mPTR = get_matPTR_from_matID(mat_mID, "", __func__, 0);
  return(matrix_is_invalid(mat_mPTR));
}

// matrixID versions of MACRO FUNCTIONS
mat_level_t matrix_row_level_matrixID(int64_t mat_mID)
{
    mat_ptr_t mat_ptr = get_matPTR_from_matID(mat_mID, "", __func__, 0);
    return matrix_row_level(mat_ptr);
}
mat_level_t matrix_col_level_matrixID(int64_t mat_mID)
{
    mat_ptr_t mat_ptr = get_matPTR_from_matID(mat_mID, "", __func__, 0);
    return matrix_col_level(mat_ptr);
}
int64_t matrix_sub_matrixID(int64_t mat_mID, int s)
{
    mat_ptr_t mat_ptr = get_matPTR_from_matID(mat_mID, "", __func__, 0);
    mat_ptr_t sub_mat_ptr = matrix_sub(mat_ptr, s);
    return get_matID_from_matPTR(sub_mat_ptr);
}

char *matrix_trace_matrixID(int64_t m_ID)
{
  mat_ptr_t m_ptr = get_matPTR_from_matID(m_ID, "", __func__, 0);
  return sca_get_str(matrix_trace(m_ptr));
}

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
 
mat_level_t matrix_pair_max_level(mat_ptr_t A_ptr, mat_ptr_t B_ptr)
{
    mat_level_t max_level = 0;
    if (A_ptr != NULL){
        max_level = (matrix_row_level(A_ptr) > max_level) ? matrix_col_level(A_ptr) : max_level;
        max_level = (matrix_col_level(A_ptr) > max_level) ? matrix_col_level(A_ptr) : max_level;
    }
    if (B_ptr != NULL){
        max_level = (matrix_row_level(B_ptr) > max_level) ? matrix_col_level(B_ptr) : max_level;
        max_level = (matrix_col_level(B_ptr) > max_level) ? matrix_col_level(B_ptr) : max_level;
    }
    return max_level;
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
	uint64_t bytes_ptr_list = SIZE_MPT * sizeof(mat_ptr_t);   	   			// size of array of pointers to tables
	uint64_t bytes_per_table = SIZE_MPT * sizeof(mat_ptr_t);  	   			// size of each table
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
  store.zero = (mat_ptr_t **)malloc((1 + max_level) * sizeof(mat_ptr_t *));
  store.identity = (mat_ptr_t *)malloc((1 + max_level) * sizeof(mat_ptr_t));
  if ((store.zero == NULL) || (store.identity == NULL)) {
    ALLOCFAIL();
    exit(1);    
  }
  for (int i= 0; i <= max_level; ++i) {
    store.zero[i] = (mat_ptr_t *)malloc((1 + max_level) * sizeof(mat_ptr_t ));
    if (store.zero[i] == NULL) {
      ALLOCFAIL();
      exit(1);
    }
  }
  store.iHadamard = (mat_ptr_t *)malloc((1 + max_level) * sizeof(mat_ptr_t));
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
  store.mat_ptr_table_list = (mat_ptr_t **)calloc(SIZE_MPT,sizeof(mat_ptr_t *));  //list of ptrs to tables
  store.mat_ptr_table_list[0] = (mat_ptr_t *)calloc(SIZE_MPT,sizeof(mat_ptr_t));  //the first table

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
int lock_matrix(mat_ptr_t m_ptr) {
  if (m_ptr == MATRIX_PTR_INVALID) {
    fprintf (stderr,"ERROR: in lock_matrix, matrix does not exist!\n");
    return (0);
  }

  if (m_ptr->lock == 1) {
    if (VERBOSE>DEBUG) {
      fprintf(stderr,"NOTE: in lock_matrix, attempt to lock already-locked matrix\n");
      fprintf(stderr,"assuming that any submatrices are also locked...\n");
    }
    return(1);
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
  mat_ptr_t m_ptr = get_matPTR_from_matID(m_mID, "", __func__,0);
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
int set_hold_matrix(mat_ptr_t m_ptr) {
  if (m_ptr == MATRIX_PTR_INVALID) {
    fprintf (stderr,"ERROR: in set_hold_matrix, matrix does not exist!\n");
    return (0);
  }
  if (m_ptr->lock) {
    if (VERBOSE>BASIC)
      fprintf (stderr,"NOTE: no need to hold a locked matrix.\n");
    return(1);
  }

  ++(m_ptr->hold);
  return (1);
}   


/************************************************************
 *          release_hold_matrix_from_matrixID	            *
 *    	Python Interface version of release_hold_matrix	    *
 *      using matrixID instead of a pointer            *
 ************************************************************/
int release_hold_matrix_from_matrixID(int64_t m_mID) {
  // get the matrix pointers from the matrixID, and see if still in store 
  mat_ptr_t m_ptr = get_matPTR_from_matID(m_mID, "", __func__,0);
  return release_hold_matrix(m_ptr);
}

/****************************************************************
 *                  release_hold_matrix                         *
 *  See explanation above for set_hold_matrix and lock_matrix.  *
 *                                                              *
 ****************************************************************/
int release_hold_matrix(mat_ptr_t m_ptr) {

  if (m_ptr == MATRIX_PTR_INVALID) {
    fprintf (stderr,"ERROR: in release_hold_matrix, matrix does not exist!\n");
    return (0);
  }

  if (m_ptr->lock) {
    if (VERBOSE>BASIC)
      fprintf (stderr,"NOTE: no holds needed on locked matrices\n");
    return (1);
  }

  if (m_ptr->hold == 0) {
    if (VERBOSE>SILENT) 
      fprintf (stderr,"WARN: attempt release_hold_matrix for matrix with no holds\n");
    return (1);
  }

  --(m_ptr->hold);

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
 
  int top_level = store.largest_level;
  if (verbose) {
    printf("Inside routine preload_matrix_store with:\n");
    printf("   retrieved largest level %d\n", top_level);
  }
  
  // panels of sub matrices  
  mat_ptr_t sm[4];
  
  // 0, 1, -1 are constants needed for preload of zero, identity, integer Hadamard
  // NOTE: scalar ops must be instantiated before preloading matrix store!
  scalarType scalar0;
  scalarType scalar1;
  scalarType scalarM1;
  sca_init(&scalar0);
  sca_init(&scalar1);
  sca_init(&scalarM1);

  sca_set_str(&scalar0, "0");
  sca_set_str(&scalar1, "1");
  sca_set_str(&scalarM1, "-1");
  if (verbose) {
    printf("Inside routine preload_matrix_store with:\n");
    char *scalar0_string = sca_get_str(scalar0);
    char *scalar1_string = sca_get_str(scalar1);
    char *scalarM1_string = sca_get_str(scalarM1);
    printf("loaded %s, %s, %s.\n", scalar0_string, scalar1_string, scalarM1_string);
    free(scalarM1_string);
    free(scalar1_string);
    free(scalar0_string);
  }
  

  // LOAD ALL THE ZERO MATRICES
  store.zero[0][0] = get_valMatPTR_from_val(scalar0);
  if (!lock_matrix(store.zero[0][0])) {
      fprintf(stderr,"FAIL: scalar zero preload tried to lock a matrix which does not exist\n");
      exit(1);
    }
  if (verbose) {
    printf("Inside routine preload_matrix_store:\n");
    printf("   preloaded store.zero[0][0]\n");
  }
  

  for (int k = 1; k <= top_level; k++) {

    // Case: ROW_VECTOR  (0,2**k)
    sm[0] = store.zero[0][k-1];
    sm[1] = store.zero[0][k-1];
    sm[2] = MATRIX_PTR_INVALID;
    sm[3] = MATRIX_PTR_INVALID;
    store.zero[0][k] = get_matPTR_from_array_of_four_subMatPTRs(sm, 0, k);
    if (!lock_matrix(store.zero[0][k])) {
      fprintf(stderr,"FAIL: row Z preload tried to lock a matrix which does not exist\n");
      exit(1);
    }
    // Case: COL_VECTOR (2**k,0)
    sm[0] = store.zero[k-1][0];
    sm[1] = MATRIX_PTR_INVALID;
    sm[2] = store.zero[k-1][0];
    sm[3] = MATRIX_PTR_INVALID;
    store.zero[k][0] = get_matPTR_from_array_of_four_subMatPTRs(sm, k, 0);
    if (!lock_matrix(store.zero[k][0])) {
      fprintf(stderr,"FAIL: col Z preload tried to lock a matrix which does not exist\n");
      exit(1);
    }
  }
  if (verbose) {
    printf("Inside routine preload_matrix_store:\n");
    printf("   preloaded row and column zero vectors into store.zero\n");
  }


  // Case: MATRIX  (2**i,2**j) with (0 < i,j <= k) 
  for (int i = 1; i <= top_level; i++) {
    for (int j = 1; j <= top_level; j++) {
      sm[0] = store.zero[i-1][j-1];
      sm[1] = store.zero[i-1][j-1];
      sm[2] = store.zero[i-1][j-1];
      sm[3] = store.zero[i-1][j-1];
      store.zero[i][j] = get_matPTR_from_array_of_four_subMatPTRs(sm, i, j);
      if (!lock_matrix(store.zero[i][j])) {
	fprintf(stderr,"FAIL: matrix Z preload tried to lock a matrix which does not exist\n");
	exit(1);
      }
    }
  }

  // Prestore all the IDENTITY matrices (square)
  store.identity[0] = get_valMatPTR_from_val(scalar1);
  if (!lock_matrix(store.identity[0])) {
      fprintf(stderr,"FAIL: scalar 0 preload tried to lock a matrix which does not exist\n");
      exit(1);
    }


  for (int i = 1; i <= top_level; i++) {
    sm[1] = sm[2] = store.zero[i-1][i-1];
    sm[0] = sm[3] = store.identity[i-1];
    store.identity[i] = get_matPTR_from_array_of_four_subMatPTRs(sm, i, i);
    if (!lock_matrix(store.identity[i])) {
      fprintf(stderr,"FAIL: I preload tried to lock a matrix which does not exist\n");
      exit(1);
    }
  }

  if (verbose) {
    printf("Inside routine preload_matrix_store:\n");
    printf("   preloaded all identity matrices into store.identity\n");
  }


  // Prestore all the integer HADAMARD matrices and scalar negative one.

  // Load scalar negative one into the otherwise un-used store.iHadamard[0]
  store.iHadamard[0] = get_valMatPTR_from_val(scalarM1);
  // if -1 == 1 (eg. in boolean), avoid warning about locking a locked matrix. 
  if (sca_eq(scalarM1, scalar1) == 0){
      if (!lock_matrix(store.iHadamard[0])) {
          fprintf(stderr,"FAIL: scalar negative one preload tried to lock a nonexistent matrix\n");
          exit(1);
      }
  }

  sm[0] = sm[1] = sm[2] = store.identity[0];
  sm[3] = store.iHadamard[0];
  store.iHadamard[1] = get_matPTR_from_array_of_four_subMatPTRs(sm,1,1);
  if (!lock_matrix(store.iHadamard[1])) {
      fprintf(stderr,"FAIL: iHadamard[1] preload tried to lock a nonexistent matrix\n");
      exit(1);
    }

  if (verbose) {
    printf("Inside routine preload_matrix_store:\n");
    printf("   preloaded iHadamard[1]\n");
  }


  for (int i = 2; i <= top_level; i++) {
    sm[0] = sm[1] = sm[2] = store.iHadamard[i-1];
    sm[3] = scalar_mult(store.iHadamard[0],store.iHadamard[i-1]);
    store.iHadamard[i] = get_matPTR_from_array_of_four_subMatPTRs(sm, i, i);
    if (!lock_matrix(store.iHadamard[i])) {
      fprintf(stderr,"FAIL: HH preload tried to lock a matrix which does not exist\n");
      exit(1);
    }
  }

  //  Preload the other basic 1 and 2 bit matrices and set global names
  if (verbose){printf("Initializing globals.\n");}
  init_globals();
 
  sca_clear(&scalar0);
  sca_clear(&scalar1);
  sca_clear(&scalarM1);

  return 1;
}


void exit_if_matrix_ptrs_invalid(const char *function, int mat_ptr_num, ...)
{
  // uses va_start, va_list and va_arg from stdarg.h to take a variable 
  // number of inputs to routine.
  va_list mat_ptrs;

  int at_least_one_invalid = 0;
  
  va_start(mat_ptrs, mat_ptr_num);

  for (int i = 0; i < mat_ptr_num; i++){
    mat_ptr_t mat_ptr = va_arg(mat_ptrs, mat_ptr_t);
    if (mat_ptr == MATRIX_PTR_INVALID){
        at_least_one_invalid = 1;
        fprintf(stderr,"ERROR: In %s the %dth (0-up) matrix checked was not a valid matrix.\n", function, i);
    }
  }
  va_end(mat_ptrs);

  if (at_least_one_invalid)
    exit(EXIT_FAILURE);
}

// Python interface version 
int64_t get_zero_matrixID(mat_level_t row_level, mat_level_t col_level)
{
  mat_ptr_t m_ptr = get_zero_matrix_ptr(row_level, col_level);
  return get_matID_from_matPTR(m_ptr);
}  

mat_ptr_t get_zero_matrix_ptr(mat_level_t row_level, mat_level_t col_level)  
{
  mat_level_t highestInputLevel = MAX(row_level, col_level);

  if (highestInputLevel > store.largest_level) 
  {
    fprintf(stderr,"Error: requested zero matrix (levels %dx%d) larger than max_level input (%d) to create store.\n", row_level, col_level, store.largest_level);
    exit(1);
  }
  return store.zero[row_level][col_level];
}

// Python interface version 
int64_t get_identity_matrixID(mat_level_t level)
{
  mat_ptr_t m_ptr = get_identity_matrix_ptr(level);
  return get_matID_from_matPTR(m_ptr);
}

mat_ptr_t get_identity_matrix_ptr(mat_level_t level) 
{
  if (level > store.largest_level) 
    {
      fprintf(stderr,"Error: requested identity matrix has size %d\n", level);
      fprintf(stderr,"       which is larger than maximum matrix size %d.\n", store.largest_level);
      exit(1);
    }
  
  return store.identity[level];
}

// Python interface version 
int64_t get_iHadamard_matrixID(mat_level_t level) 
{
  mat_ptr_t m_ptr = get_iHadamard_matrix_ptr(level);
  return get_matID_from_matPTR(m_ptr);
}

mat_ptr_t get_iHadamard_matrix_ptr(mat_level_t level) 
{
  if (level == 0)
  {
    fprintf(stderr,"Error in get_iHadamard_matrix_ptr: there is no sensible definition ");
    fprintf(stderr,"for a level 0 Hadamard matrix\n");
    exit(1);
  }

  if (level > store.largest_level) 
  {
    fprintf(stderr,"Error: requested integer Hadamard matrix larger than maximum size.\n");
    exit(1);
  }

  return store.iHadamard[level];
}


uint32_t matrix_appears_as_sub_count_increment(mat_ptr_t id)
{
  if ( ++(id->appears_as_sub_count) == 0) {
    printf("ERROR: increment of sub_count caused overflow\n");
    exit(1);
  }
  return id->appears_as_sub_count;
}


uint32_t matrix_appears_as_sub_count_decrement(mat_ptr_t id)
{
  if (id->appears_as_sub_count==0) {
    printf("ERROR: trying to decrement sub_count of zero\n");
    exit(1);
  }
  return --(id->appears_as_sub_count);
}


// This function compares the nbhd approx of a scalar obtained from
// a matrix record with the nbhd approx of a (possibly unstored) new scalar
/*!
 * \ingroup larc
 * \brief A utility function which gets the neighborhood approximation to a scalar and compares it with a value already in the matrix store
 * \param m_ptr The pointer to a scalar value
 * \param scalar A scalar value to compare with it
 * \return 1 if the values are the same, 0 if they are different or m_ptr does not point to a scalar
 */
static int compare_nbhd_approx_stored_vs_new(mat_ptr_t m_ptr, scalarType scalar)
{
#ifdef DEBUG
  printf("in %s\n",__func__);
#endif
  // EASY CASES for failed comparison
  if ((m_ptr == MATRIX_PTR_INVALID) || (matrix_type(m_ptr) != SCALAR))
  {
#ifdef DEBUG
     printf(">>invalid pointer or not a scalar!\n");
#endif
     return 0;
  }

#ifdef DEBUG
  printf("%s compares the nbhd approx of a scalar obtained from a matrix record\n", __func__);
  printf("with the nbhd approx of a (possibly unstored) new scalar\n");
  char *matrix_trace_string= sca_get_str(matrix_trace(m_ptr));
  char *scalar_string = sca_get_str(scalar);
  printf("nbhd approx of scalar stored in a record is %s, nbhd approx of the new scalar is %s.\n", matrix_trace_string,scalar_string); 
  free(matrix_trace_string);
  free(scalar_string);
  printf("calling sca_eq_approx:\n");
#endif    
  
  int approxs_are_eq = sca_eq_approx(matrix_trace(m_ptr), scalar);
#ifdef DEBUG
  printf("approxs_are_eq = %d\n",approxs_are_eq);
#endif
  return approxs_are_eq;
#undef DEBUG
}



/// NON-SCALAR
/* Returns 0 if matrix differs from provided content, non-zero if equal */
/*!
 * \ingroup larc
 * \brief
 * \param m_ptr The pointer to a stored matrix
 * \param panel A set of four pointers to submatrices (which define a larger matrix)
 * \param mat_type The type of matrix the panel describes, one of MATRIX, ROW_VECTOR, COL_VECTOR 
 * \return 1 if the matrices are the same, 0 if not
 */
static int matrix_compare_panel(mat_ptr_t m_ptr, mat_ptr_t panel[4], matrix_type_t mat_type)
{
  // EASY CASES for failed comparison
  if ((m_ptr == MATRIX_PTR_INVALID) || (matrix_type(m_ptr) != mat_type))
  { 
     return 0; 
  }

  return (panel[0] == matrix_sub(m_ptr,0) &&
          panel[1] == matrix_sub(m_ptr,1) &&
          panel[2] == matrix_sub(m_ptr,2) &&
          panel[3] == matrix_sub(m_ptr,3) ); 
}

/* Returns the matrix PTR which matches scalar, starting from the hash */
/*!
 * \ingroup larc
 * \brief Finds the matrix pointer for a given scalar value (does not add it to the matrix store if it is not found)
 * \param hash A value indicating which hash chain would contain the scalar value
 * \param scalar A scalarType value
 * \return The matrix pointer to the stored scalar, or MATRIX_PTR_INVALID if not found
 */
static mat_ptr_t matrix_find_scalar(uint64_t hash, scalarType scalar)
{

#ifdef DEBUG 
  printf("In routine %s\n", __func__);
  printf(">>>calling hash_get_chain\n");
#endif

  hash_node_t *n = hash_get_chain(store.hash_table, hash);
  uint32_t depth = 0;

#ifdef HASHSTATS
  (store.hash_table->num_accesses[hash])++;
#endif

  while (n) {

#ifdef DEBUG
    printf(">>>in loop, comparing old and new approx values\n");
#endif
    if (compare_nbhd_approx_stored_vs_new((mat_ptr_t)(n->record_ptr), scalar)){
      store.hash_table->hits++;    
      n->hits++;

      if (depth > store.deepest) {
        store.deepest = depth;
      }

#if 0
    // This code was a start at working on optimizing the hash chain order
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

#ifdef DEBUG 
  printf("Leaving %s with valid matrix pointer\n", __func__);
#endif

      return (mat_ptr_t) (n->record_ptr);
    }
    n = n->next;
    depth++;
  } // end while

  store.hash_table->misses++;

#ifdef DEBUG 
  printf("Leaving %s with invalid matrix pointer\n", __func__);
#endif

  return MATRIX_PTR_INVALID;
#undef DEBUG  
}

/* Returns the matrix PTR which matches panel, starting from the hash */
/*!
 * \ingroup larc
 * \brief Finds the matrix pointer for a given matrix as describe by a panel of submatrices (does not add it to the matrix store if it is not found)
 * \param hash A value indicating which hash chain would contain the panel
 * \param panel An array of four matrix pointers
 * \param mat_type The type of matrix the panel describes, one of MATRIX, ROW_VECTOR, COL_VECTOR 
 * \return The matrix pointer to the stored matrix, or MATRIX_PTR_INVALID if not found
 */
static mat_ptr_t
matrix_find_panel(uint64_t hash, mat_ptr_t panel[4], matrix_type_t mat_type)
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
	  if (matrix_compare_panel((mat_ptr_t)(n->record_ptr), panel, mat_type)) {
                        store.hash_table->hits++;
			n->hits++;

			if (depth > store.deepest) {
				store.deepest = depth;
			}


#if 0
    // This code was a start at working on optimizing the hash chain order
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

			return (mat_ptr_t) (n->record_ptr);
		}
		n = n->next;
		depth++;
	}

	store.hash_table->misses++;

	return MATRIX_PTR_INVALID;
}

/*!
 * \ingroup larc
 * \brief Add an entry to the table indexing matrix pointer values by matrixID
 * \param m_ptr A pointer to a matrix
 * \return The matrixID for that pointer
 */
static int64_t record_mat_ptr_by_matrixID( mat_ptr_t m_ptr)
{
  // This function adds an entry to the table indexing matrix pointer values
  // by matrixID values. In the future, when we are more actively purging
  // matrices from the matrix store, this may be replaced with a hash.
  int64_t list_index;
  int64_t table_index;
  int64_t index = get_matID_from_matPTR(m_ptr);

  // high order bits give index in the list of tables
  list_index = index >> LOG_SIZE_MPT;

  // low order bits give index in one of the tables
  table_index = index & MASK_MPT;

  // TODO: this test may be redundant if we have the test when creating matrixIDs
  if (list_index >= SIZE_MPT) 
  {
    fprintf(stderr,"ERROR: LOG_SIZE_MPT needs to be larger.\n");    // or, can modify code to allocate more tables as needed
    exit(1);
  }   

  // if necessary, allocate memory for the next table in the list
  if (store.mat_ptr_table_list[list_index] == 0)    
  {
    store.mat_ptr_table_list[list_index] = (mat_ptr_t *)calloc(SIZE_MPT,sizeof(mat_ptr_t));
  }
  
  // with setup complete, add pointer to the appropriate table
  store.mat_ptr_table_list[list_index][table_index] = m_ptr;
  return(index);
}


/* Inserts a new scalar matrix into the store and returns the new matrix PTR*/
// NOTE: a new copy of scalar is made for storage matrices (in trace_element)
// so routines above this are responsible for releasing (sca_clear) scalar when
// finished (namely the ones that initialize scalar). Similarly those routines
// can then change scalar without changing the entry that was put into the
// matrix store. 
/*!
 * \ingroup larc
 * \brief Inserts a scalar value into the matrix store
 * \param hash The hash for the scalar value
 * \param scalar The scalarType value to be stored
 * \return The matrix pointer for the stored scalar
 */
static mat_ptr_t matrix_insert_scalar(uint64_t hash, scalarType scalar) {
  
  int verbose = 0;
  if (verbose) {  
    printf("Inside %s: levels (0,0), hash %" PRIu64 "\n", __func__, hash);
  }
  
  // allocate space for the new LARC matrix
  mat_ptr_t m_ptr = calloc(1, sizeof(larc_matrix_t));
  if (m_ptr == NULL) {
    ALLOCFAIL();
    return MATRIX_PTR_INVALID;
  }
  
  // set the basic parameters for the matrix
  m_ptr->matrixID = store.matrixID_next++;
  m_ptr->row_level = 0;
  m_ptr->col_level = 0;

  if (verbose) { 
    printf("  allocated space, set row level, col level, and matrixID\n");
  }
  
#ifndef HASH_CHAIN_GROWS_AT_HEAD
  // INSERT new items at tail of hash collision chain
  hash_insert_at_tail(store.hash_table, (record_ptr_t)m_ptr, hash);
#endif
#ifdef HASH_CHAIN_GROWS_AT_HEAD
  // INSERT new items at head of hash collision chain
  hash_insert_at_head(store.hash_table, (record_ptr_t)m_ptr, hash);
#endif

  if (verbose) {
    printf("  appended hash chain with 0, 0levels\n");
  }
  
  // update matrix store parameters
  store.hist[0][0]++;
  m_ptr->info = 0;
  m_ptr->hold = 0;
  m_ptr->lock = 0;
  m_ptr->appears_as_sub_count = 0;
  
  if (verbose) { 
    printf("  set row level, col level, and matrixID\n");
  }
  
  if (verbose) {
    char *scalar_string = sca_get_str(scalar);
    printf("ADDING TO MATRIX SCALAR STORE: %ld -> %s\n", m_ptr->matrixID, scalar_string);
    free(scalar_string);
  }
    
  // the trace of a 1x1 matrix is the element inside the matrix
  // trace_element is a separate copy of scalar, not a reference to it. 
  sca_init(&(m_ptr->trace_element));
  sca_set(&(m_ptr->trace_element), scalar);
    
  // test to see if scalars are zero or one, and/or self-adjoint
  scalarType *sca_id = &scratchVars.quick_use;
  sca_set_str(sca_id, "0");
  m_ptr->iszero = (0 != sca_eq(m_ptr->trace_element, *sca_id));
  sca_set_str(sca_id, "1");
  m_ptr->isid = (0 != sca_eq(m_ptr->trace_element, *sca_id));
  
  store.nscalars++;
  if( store.nscalars > store.max_nscalars ) {
    store.max_nscalars = store.nscalars;
  }

  // Record the matrix pointer in table indexed by matrixID
  record_mat_ptr_by_matrixID(m_ptr);
  
  return m_ptr;
}

/* Inserts a new matrix into the store and returns the new matrix PTR */
/*!
 * \ingroup larc
 * \brief Inserts a matrix into the matrix store
 * \param hash The hash for the panel describing the matrix
 * \param panel An array of four matrix pointers; the quadrant submatrices of the matrix to be stored
 * \param row_level The row level of the matrix to be stored
 * \param col_level The column level of the matrix to be stored
 * \return The matrix pointer for the stored matrix
 */
static mat_ptr_t
matrix_insert_panel(uint64_t hash, mat_ptr_t panel[4], mat_level_t row_level, mat_level_t col_level)
{
  // EXCEPTION TESTING
  if (MAX(row_level,col_level) > store.largest_level){
    fprintf(stderr,"ERROR: In %s levels %d %d of incoming matrix too large\n",
        __func__,row_level,col_level);
    exit(1);
  }
  
  matrix_type_t mat_type = matrix_type_from_levels(row_level, col_level);
  int verbose = 0;
  if (verbose) { 
    printf("Inside %s: levels (%d,%d), hash %" PRIu64 "\n",
        __func__, row_level, col_level, hash);
  }
  
  // allocate space for the new LARC matrix
  mat_ptr_t m_ptr = calloc(1, sizeof(larc_matrix_t));
  if (m_ptr == NULL) {
    ALLOCFAIL();
    return MATRIX_PTR_INVALID;
  }
  
  // set the basic parameters for the matrix
  m_ptr->matrixID = store.matrixID_next++;
  m_ptr->row_level = row_level;
  m_ptr->col_level = col_level;
  
  if (verbose) { 
    printf("  allocated space, set row level, col level, and matrixID \n");
  }
  
#ifndef HASH_CHAIN_GROWS_AT_HEAD
  // INSERT new items at tail of hash collision chain
  hash_insert_at_tail(store.hash_table, (record_ptr_t)m_ptr, hash);
#endif
#ifdef HASH_CHAIN_GROWS_AT_HEAD
  // INSERT new items at head of hash collision chain
  hash_insert_at_head(store.hash_table, (record_ptr_t)m_ptr, hash);
#endif

  if (verbose) { 
    printf("  appended hash chain with %d, %d levels\n",row_level,col_level);
  }
  
  // update matrix store parameters
  store.hist[row_level][col_level]++;
  m_ptr->info = 0;
  m_ptr->hold = 0;
  m_ptr->lock = 0;
  m_ptr->appears_as_sub_count = 0;
  
  if (verbose) { 
    printf("  set row level, col level, matrixID, and adjoint\n");
  }
  
  if (mat_type == SCALAR) {
    // this shouldn't happen
    fprintf(stderr,"in matrix_insert_panel: somehow got SCALAR for mat_type\n");
    exit(1);
  }

  else if (mat_type == MATRIX) {

    // Exception checking
    if( (matrix_row_level(panel[0]) != matrix_row_level(panel[1])) || 
	(matrix_row_level(panel[0]) != matrix_row_level(panel[2])) ||
	(matrix_row_level(panel[0]) != matrix_row_level(panel[3])) )
      {
	fprintf(stderr,"ERROR: In %s trying to insert matrix with non-equal row sub-levels %d %d %d %d\n", __func__,
		matrix_row_level(panel[0]), matrix_row_level(panel[1]),
		matrix_row_level(panel[2]), matrix_row_level(panel[3]) );
	exit(1);
      }

    if( (matrix_col_level(panel[0]) != matrix_col_level(panel[1])) || 
	(matrix_col_level(panel[0]) != matrix_col_level(panel[2])) ||
	(matrix_col_level(panel[0]) != matrix_col_level(panel[3])) )
      {
	fprintf(stderr,"ERROR: In %s trying to insert matrix with non-equal col sub-levels %d %d %d %d\n", __func__,
		matrix_col_level(panel[0]), matrix_col_level(panel[1]),
		matrix_col_level(panel[2]), matrix_col_level(panel[3]) );
	exit(1);
      }

    // define new matrix in terms of its four panels
    for (int i = 0; i < 4; i++) {
      if (panel[i] == NULL) {
	fprintf(stderr,"Error looking up submatrices\n");
	exit(1);
	// return MATRIX_PTR_INVALID;
      }
      // each panel is a submatrix of the new matrix
      m_ptr->submatrix[i] = panel[i];
      matrix_appears_as_sub_count_increment(panel[i]);
    }

    // special SQUARE MATRIX stuff
    if (row_level == col_level) {
      // trace of matrix is sum of traces of the two diagonal blocks
      //m_ptr->trace_element = panel[0]->trace_element + panel[3]->trace_element;
      sca_init(&(m_ptr->trace_element));
      sca_add(&(m_ptr->trace_element), panel[0]->trace_element, panel[3]->trace_element);
    } // end square matrix tests

    // testing for zero or identity matrices
    if (matrix_is_zero(panel[1]) && matrix_is_zero(panel[2])) {
      if (matrix_is_zero(panel[0]) && matrix_is_zero(panel[3])) {
	  // [0, 0; 0, 0] 
          m_ptr->iszero = 1;
      }
      else if (matrix_is_id(panel[0]) && matrix_is_id(panel[3])) {
          // [I, 0; 0, I]; only possible for square matrices
          m_ptr->isid = 1;
      }
    }
  }  // end MATRIX type

  else if (mat_type == COL_VECTOR) {
    
    // Exception checking
    if( matrix_col_level(panel[0]) != matrix_col_level(panel[2]) )
    {
	fprintf(stderr,"ERROR: In %s trying to insert column vector with non-equal col sub-levels %d %d\n", __func__,
		matrix_col_level(panel[0]), matrix_col_level(panel[2]) );
	exit(1);
    }
    if ((panel[0] == NULL) || (panel[2] == NULL)) {
      fprintf(stderr,"Error looking up submatrices\n");
      exit(1);
      // return MATRIX_PTR_INVALID;
    }

    m_ptr->submatrix[0] = panel[0];
    m_ptr->submatrix[2] = panel[2];
    matrix_appears_as_sub_count_increment(panel[0]);
    matrix_appears_as_sub_count_increment(panel[2]);
    
    if (matrix_is_zero(panel[0]) && matrix_is_zero(panel[2])) {
      m_ptr->iszero = 1;
    }

  }  // end of COL_VECTOR type
 
  else if (mat_type == ROW_VECTOR) {

    // exception checking
    if( matrix_row_level(panel[0]) != matrix_row_level(panel[1]) )
    {
	fprintf(stderr,"ERROR: In %s trying to insert row vector with non-equal row sub-levels %d %d\n", __func__,
		matrix_row_level(panel[0]), matrix_row_level(panel[1]) );
	exit(1);
    }
    if ((panel[0] == NULL) || (panel[1] == NULL)) {
      fprintf(stderr,"Error looking up submatrices\n");
      exit(1);
      // return MATRIX_PTR_INVALID;
    }

    m_ptr->submatrix[0] = panel[0];
    m_ptr->submatrix[1] = panel[1];
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
mat_ptr_t matrix_adjoint(mat_ptr_t m_ptr)
{
  // EXCEPTION CHECKING
  if matrix_is_invalid(m_ptr) {
    return MATRIX_PTR_INVALID;
  }

  // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE
  mat_ptr_t adj_ptr = op_get(ADJOINT,m_ptr,m_ptr);

  if (adj_ptr == MATRIX_PTR_INVALID) {
    matrix_type_t mat_type = matrix_type(m_ptr);
    matrix_type_t adjoint_type;
    if (mat_type == ROW_VECTOR) {adjoint_type = COL_VECTOR;}
    else if (mat_type == COL_VECTOR) {adjoint_type = ROW_VECTOR;}
    else  {adjoint_type = mat_type;}
  
    uint64_t hash;
    if (matrix_type(m_ptr) == SCALAR) {
      if (sca_is_real(matrix_trace(m_ptr))) {
        // adjoint of real number is that number, guaranteed to be in store
        adj_ptr = m_ptr;
      } else {
        scalarType scalar;
        sca_init(&scalar); 
      // adjoint of complex number might not be in store
        sca_conj(&scalar, matrix_trace(m_ptr));
        hash = hash_from_matrix_scalar(scalar, adjoint_type, store.hash_table->exponent);
        adj_ptr = matrix_find_scalar(hash, scalar);
        if matrix_is_invalid(adj_ptr) {
          adj_ptr = matrix_insert_scalar(hash, scalar);
        }
        sca_clear(&scalar);
      }
    } else {
      mat_ptr_t panel[4];
      panel[0] = matrix_adjoint(matrix_sub(m_ptr,0));
      panel[1] = matrix_adjoint(matrix_sub(m_ptr,2));
      panel[2] = matrix_adjoint(matrix_sub(m_ptr,1));
      panel[3] = matrix_adjoint(matrix_sub(m_ptr,3));
      hash = hash_from_matrix_panel(panel, adjoint_type, store.hash_table->exponent);
      adj_ptr = matrix_find_panel(hash, panel, adjoint_type);
      if matrix_is_invalid(adj_ptr) {
        adj_ptr = matrix_insert_panel(hash, panel, 
                  matrix_col_level(m_ptr), matrix_row_level(m_ptr));
      }
    }

    // STORE RESULT IN OPERATIONS STORE
    op_set(ADJOINT, m_ptr, m_ptr, adj_ptr);
  }  // end matrix was not found in operations store

  return adj_ptr;
}

// These next two routines are used in the python interface
// because we only let python users see matrixIDs
int64_t
get_matID_from_four_subMatIDs(int64_t A_mID, int64_t B_mID, int64_t C_mID, int64_t D_mID, 
                        mat_level_t row_level, mat_level_t col_level) 
{

  // grab the matrix pointers corresponding to each matrixID, and
  // checks on validity occur in get_matPTR_from_array_of_four_subMatPTRs
  mat_ptr_t  panel[4];
  panel[0] = get_matPTR_from_matID(A_mID, "first", __func__,0);
  panel[1] = get_matPTR_from_matID(B_mID, "second", __func__,0);
  panel[2] = get_matPTR_from_matID(C_mID, "third", __func__,0);
  panel[3] = get_matPTR_from_matID(D_mID, "fourth", __func__,0);

  mat_ptr_t m_ptr = get_matPTR_from_array_of_four_subMatPTRs(panel,row_level,col_level); 
  return m_ptr->matrixID;
}

// python interface for scalar type
int64_t get_valID_from_valString(char *val)
{
  scalarType sca;
  sca_init(&sca);
  sca_set_str(&sca, val);
  mat_ptr_t m_ptr = get_valMatPTR_from_val(sca);
  sca_clear(&sca);
  return m_ptr->matrixID;
}

/* The first argument of this routine is a scalarType (complex, double, int) */
mat_ptr_t get_valMatPTR_from_val(scalarType scalar)
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
  if (verbose) {
    printf("The hash of the nbhd approximation returned from hash_from_matrix_scalar is %zd\n",hash);
  }
  mat_ptr_t m_ptr = matrix_find_scalar(hash, scalar);
  if (verbose) {
    printf("  called matrix_find_scalar, result (m_ptr) is %p\n", m_ptr);
  }
  
  // if unable to find, attempt to insert
  if (m_ptr == MATRIX_PTR_INVALID) {
    if (verbose) { 
          printf("Inside %s, matrix ptr not found by matrix_find\n", __func__);
          printf("   CALLING  matrix_insert\n");
    }
    m_ptr = matrix_insert_scalar(hash, scalar);
  }
  if (verbose) {   
    printf("    RETURNED from matrix_insert with matrix ptr %p\n", m_ptr);
  }
  
  // if couldn't insert, then bail out
  if (m_ptr == MATRIX_PTR_INVALID) {
    abort();
  }
  if (verbose) {
    printf("Leaving %s\n", __func__);
  }
  
  return m_ptr;
}

/* The first argument of this routine mat_val_ptr is an address (void *) and
   the element at that address is either a scalarType (complex, double, int)
   or panel[4] which is a list giving the addresses of four matrix records.
*/
mat_ptr_t
get_matPTR_from_array_of_four_subMatPTRs(mat_ptr_t panel[4], mat_level_t row_level, mat_level_t col_level)
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

  // check to see if validity of matrix pointers match matrix type
  if (mat_type == MATRIX){
    for (int i = 0; i < 4; i++){
      if (panel[i] == MATRIX_PTR_INVALID){
        fprintf(stderr,"ERROR: attempting to form matrix from panels with invalid ptr in panel %d\n", i);
        exit(EXIT_FAILURE);
      }
      else {
        if (matrix_row_level(panel[i]) + 1 != row_level) 
          fprintf(stderr,"WARNING on panel %d: row levels inconsistent\n",i);
        if (matrix_col_level(panel[i]) + 1 != col_level) 
          fprintf(stderr,"WARNING on panel %d: col levels inconsistent\n",i);
      }
    }
  }
  else if (mat_type == ROW_VECTOR){
    for (int i = 0; i < 4; i++){
      if (i < 2){
        if (panel[i] == MATRIX_PTR_INVALID){
          fprintf(stderr,"ERROR: attempting to form row vector from panels with invalid ptr in panel %d\n", i);
          exit(EXIT_FAILURE);
        }
        else {
          if (matrix_row_level(panel[i]) != row_level) 
            fprintf(stderr,"WARNING on panel %d: row levels inconsistent\n",i);
          if (matrix_col_level(panel[i]) + 1 != col_level) 
            fprintf(stderr,"WARNING on panel %d: col levels inconsistent\n",i);
        }
      }
      else {
        if (panel[i] != MATRIX_PTR_INVALID){
          fprintf(stderr,"ERROR: attempting to form row vector from panels with valid ptr in panel %d\n", i);
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  else if (mat_type == COL_VECTOR){
    for (int i = 0; i < 4; i++){
      if (i % 2 == 0){
        if (panel[i] == MATRIX_PTR_INVALID){
          fprintf(stderr,"ERROR: attempting to form column vector from panels with invalid ptr in panel %d\n", i);
          exit(EXIT_FAILURE);
        }
        else {
          if (matrix_row_level(panel[i]) + 1 != row_level) 
            fprintf(stderr,"WARNING on panel %d: row levels inconsistent\n",i);
          if (matrix_col_level(panel[i]) != col_level) 
            fprintf(stderr,"WARNING on panel %d: col levels inconsistent\n",i);
        }
      }
      else{
        if (panel[i] != MATRIX_PTR_INVALID){
          fprintf(stderr,"ERROR: attempting to form column vector from panels with valid ptr in panel %d\n", i);
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  else {
    fprintf(stderr,"ERROR: attempting to form matrix of type %d from panels\n", mat_type);
    exit(EXIT_FAILURE);
  }

  uint64_t hash = hash_from_matrix_panel(panel, mat_type, store.hash_table->exponent);
  mat_ptr_t m_ptr = matrix_find_panel(hash, panel, mat_type);
  if (verbose) { 
    printf("  About to check if the matrix was not found, m_ptr is %p\n",m_ptr);
  }
  
  // if unable to find attempt to insert
  if (m_ptr == MATRIX_PTR_INVALID) {
    if (verbose) { 
          printf("Inside %s, matrix index not found by matrix_find\n",__func__);
          printf("   CALLING  matrix_insert\n");
    }
    m_ptr = matrix_insert_panel(hash, panel, row_level, col_level);
  }
  if (verbose) {   
    printf("    RETURNED from matrix_insert with matrix index %p\n",m_ptr);
  }
  
  // if couldn't insert, then bail out
  if (m_ptr == MATRIX_PTR_INVALID) {
    fprintf(stderr,"ERROR: inserting matrix from panels into store failed.\n");
    exit(EXIT_FAILURE);
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
  mat_ptr_t m_ptr = get_matPTR_from_matID(m_mID, "", __func__,0);
  if (m_ptr == MATRIX_PTR_INVALID) { exit(1); }

  // call matrix pointer version of function
  mat_ptr_t adjoint_ptr = matrix_adjoint(m_ptr);
  int64_t adjoint_mID = get_matID_from_matPTR(adjoint_ptr);
  return adjoint_mID;
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

mat_ptr_t
get_matPTR_from_matID(int64_t matrixID, const char* arg_no, 
        const char* calling_routine, int print_flag) {

  // grab highest matrixID
  uint64_t highest_matrixID = num_matrices_created()-1;

  // check bounds
  if ((matrixID > highest_matrixID) || (matrixID < 0) ) {
    if (print_flag)
      fprintf(stderr,"ERROR: In %s, the %s argument, matrixID %" PRId64 ", is out of range.\n",
          calling_routine, arg_no, matrixID);
    return(MATRIX_PTR_INVALID);
  }

  // get matrix address
  int64_t list_index = matrixID >> LOG_SIZE_MPT;
  int64_t table_index = matrixID & MASK_MPT;
  mat_ptr_t mat_val_ptr = store.mat_ptr_table_list[list_index][table_index];

  // warn the user if the matrixID doesn't refer to a valid matrix
  if ((mat_val_ptr == MATRIX_PTR_INVALID) && (print_flag==1)) {
    fprintf(stderr,"WARNING: In %s, problem with the %s argument:",
        calling_routine, arg_no);
    fprintf(stderr," matrix with matrixID %" PRId64 " has been removed from the Matrix Store.\n",
         matrixID);
  }

  return mat_val_ptr;
}

/*!
 * \ingroup larc
 * \brief Writes statistic information about a matrix in the matrix_store to a file
 * \param mat_ID The matrixID of the matrix 
 * \param f A file pointer
 * \return 1 
 */
static int  
matrix_info_to_file(uint64_t mat_ID, FILE *f) {
  // EXCEPTIONS 
  if (f==NULL) { 
    fprintf(stderr,"WARN: %s called with bad file pointer\n", __func__); 
    return(0); 
  } 

  // get the matrix pointers from the matrixID, and see if still in store 
  mat_ptr_t m_ptr = get_matPTR_from_matID(mat_ID, "", __func__,0);

  if (matrix_is_invalid(m_ptr)) { return(1); }  // skip printing line

  // HEADER: "matrixID   Levels    Value    Appears_As_Sub Count   Lock   Hold
  //   matrixID 
  fprintf(f,"      %zd",get_matID_from_matPTR(m_ptr));
  //   Levels
  fprintf(f,"    \t(%2d,%2d)",matrix_row_level(m_ptr),matrix_col_level(m_ptr)); 
  //   Value 
  if ((matrix_row_level(m_ptr) == 0) && (matrix_col_level(m_ptr) == 0)) {  // SCALAR
    char *trace_string = sca_get_str(matrix_trace(m_ptr));
    fprintf(f," \t%s                   ", trace_string);
    free(trace_string);
  } 
  else {   // MATRIX
    mat_ptr_t panel[4];
    fprintf(f,"   [ ");
    for (int i = 0; i < 4; i++) {
      panel[i] = matrix_sub(m_ptr,i); 
      if (panel[i] == MATRIX_PTR_INVALID) {
	fprintf(f,"  X ");
      }
      else {
        fprintf(f,"%" PRId64 " ", get_matID_from_matPTR(panel[i])); 
      }
    } 
    fprintf(f,"]        ");
  }
  //   Appears_As_Sub Count = Reference Count 
  fprintf(f," \t%8u     %4d    %4d     ", matrix_appears_as_sub_count(m_ptr), matrix_has_lock(m_ptr), matrix_has_hold(m_ptr));
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
    fprintf(stderr,"ERROR in matrix_store_info_to_file: impossible range %" PRIu64 " > %" PRIu64 "\n", start, end);
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
  fprintf(f,"      Appears_As_SubCount    Lock    Hold\n"); 
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
  fprintf(f,"        Appears_As_SubCount\n"); 

  hash_table_t *table_ptr = store.hash_table;
  hash_node_t *node_ptr = table_ptr->heads[hash];
  mat_ptr_t record_ptr;

  while (node_ptr)  {
    record_ptr = (mat_ptr_t) node_ptr->record_ptr;
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
  mat_ptr_t m_ptr = get_matPTR_from_matID(m_mID, "", __func__,0);
  
  // Python users could try to remove the same thing twice, but
  // we chose to not kill the program if they do this, just warn them.
  if (m_ptr == MATRIX_PTR_INVALID) {
    fprintf(stderr,"WARNING: In %s, matrix %" PRId64 " was formerly removed from the Matrix Store.\n", __func__, m_mID);
    return(0);
  }

  return remove_matrix_from_mat_store_by_matrix_andor_node(m_ptr, NULL, 0);
}



/******************************************************************
*          remove_matrix_from_mat_store_by_matrix_andor_node (mat_ptr_t  m_ptr, hash_node_t *n, uint64_t supplied_hash)
*  This function attempts to remove a larc_matrix_t structure 
*  from the matrix store.
*  If the node pointer n is not null, then n and supplied_hash will be used to
*  delete the hash node from the hash table.
*  If the node pointer n is null, the function will calculate the hash
*  of the matrix from scratch and then search the hash table for
*  the matrix.
*     1 is returned if the remove succeeded
*     0 is returned if the matrix could not be removed
*       because it has larc_appears_as_subs or has a lock or hold.
*******************************************************************/
int remove_matrix_from_mat_store_by_matrix_andor_node (mat_ptr_t  m_ptr, hash_node_t *n, uint64_t supplied_hash)
{
  // Exception checking
  if (matrix_is_invalid(m_ptr)) {
    fprintf(stderr,"WARNING: %s called for invalid matrix pointer\n",__func__);
    return 0;
  }
 
  // Check if matrix has appears_as_subs or has a lock or hold
  if (matrix_appears_as_sub_count(m_ptr) != 0) return 0;
  if (matrix_has_lock(m_ptr)) return 0;
  if (matrix_has_hold(m_ptr)) return 0;
  
  // If the matrix is not a SCALAR, then recursively remove it
  //   for each child of matrix, decrement appears_as_sub_count 
  //   and if child is now at zero counts  attempt to remove it
  matrix_type_t mat_type = matrix_type(m_ptr);
  if (mat_type != SCALAR) {
    //printf("DEBUG: record type not scalar in %s.\n", __func__); //debug
    for (int i=0;i<4;++i) {
      mat_ptr_t child_ptr = matrix_sub(m_ptr, i);
      if (child_ptr != MATRIX_PTR_INVALID) {
	uint32_t new_appears_as_sub_count =
                matrix_appears_as_sub_count_decrement(child_ptr);
	if (new_appears_as_sub_count == 0)
          remove_matrix_from_mat_store_by_matrix_andor_node(child_ptr, NULL, 0);
      }
    }
  }

  // Remove the matrix from the hash table for the matrix store
  if (n) {
    // If n is not null, use it to directly delete the hash node
    hash_node_remove_node(store.hash_table, n, supplied_hash);
  } else {
    // If n is null, find the hash chain where the matrix is stored
    uint64_t hash;
    if (mat_type == SCALAR) {
      //printf("DEBUG: record type scalar in %s.\n", __func__); //debug
      hash = hash_from_matrix_scalar(matrix_trace(m_ptr), SCALAR, store.hash_table->exponent);
    } else {
      hash = hash_from_matrix_panel(m_ptr->submatrix, mat_type, store.hash_table->exponent);
    }
    hash_node_remove(store.hash_table, (record_ptr_t)m_ptr, hash);
  }

  // Get row level and col level for matrix.
  mat_level_t row_level;
  mat_level_t col_level;
  row_level = matrix_row_level(m_ptr);
  col_level = matrix_col_level(m_ptr);

  // Decrement the histogram and counts of matrices and scalars
  store.hist[row_level][col_level]--;
  if (mat_type == SCALAR) {
    store.nscalars--;
  } else {
    store.nmatrices--;
  }

  // special SQUARE MATRIX stuff
  if (row_level == col_level) {
    // Free/clear the scalar
    sca_clear(&(matrix_trace(m_ptr)));
  }

  // Food for thought: if the COMPLEX adjoint is added to the store just
  // for this scalar matrix, then should it be removed as well? 

  // Remove the matrix pointer from the table indexed by matrixIDs
  table_mat_ptr_by_matrixID_Remove_entry(m_ptr);

  // free the structure, note that no mallocs other than scalar occur in the matrix
  free(m_ptr);

  return (1);
}


/******************************************************************
*          remove_matrix_from_mat_store_by_matrix (mat_ptr_t  m_ptr)
*  This function attempts to remove a larc_matrix_t structure 
*  from the matrix store.
*     1 is returned if the remove succeeded
*     0 is returned if the matrix could not be removed
*       because it has larc_appears_as_subs or has a lock or hold.
*******************************************************************/
int remove_matrix_from_mat_store_by_matrix (mat_ptr_t  m_ptr)
{
  return remove_matrix_from_mat_store_by_matrix_andor_node(m_ptr, NULL, 0);
}


/******************************************************************
*          remove_matrix_from_mat_store_by_node (hash_node_t *n, uint64_t supplied_hash)
*  This function attempts to remove a larc_matrix_t structure 
*  from the matrix store.
*     1 is returned if the remove succeeded
*     0 is returned if the matrix could not be removed
*       because it has larc_appears_as_subs or has a lock or hold.
*******************************************************************/
int remove_matrix_from_mat_store_by_node (hash_node_t *n, uint64_t supplied_hash)
{
  mat_ptr_t  m_ptr = (mat_ptr_t) n->record_ptr;
  return remove_matrix_from_mat_store_by_matrix_andor_node(m_ptr, n, supplied_hash);
}


// Remove the matrix pointer from the table indexed by matrixIDs
// by replacing the value with MATRIX_PTR_INVALID
int
table_mat_ptr_by_matrixID_Remove_entry(mat_ptr_t m_ptr){

  // EXCEPTION CHECKING
  if (m_ptr == MATRIX_PTR_INVALID) { 
    fprintf(stderr,"ERROR: in %s\n", __func__); 
    fprintf(stderr,"   asked to remove m_idx =MATRIX_PTR_INVALID\n");
    exit(1);
  }

  int64_t index = get_matID_from_matPTR(m_ptr);
     
  if (index != MATRIX_ID_INVALID) {
    int64_t list_index = index >> LOG_SIZE_MPT;
    int64_t table_index = index & MASK_MPT;
    store.mat_ptr_table_list[list_index][table_index] = MATRIX_PTR_INVALID;
  }
  return(1);

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
  mat_ptr_t m_ptr = get_matPTR_from_matID(m_mID, "", __func__,0);
  
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
    uint64_t max = 1L<<(store.hash_table->exponent);
    uint64_t i;

    // loop though hash chains for each hashID without calling
    // clean_op_hash_chain so we don't have to repeat checks for invalid hashID
    hash_table_t *table_ptr = store.hash_table;
    for (i = 0; i < max; ++i) {
        hash_node_t *node_ptr = table_ptr->heads[i];
        hash_node_t *next_ptr;
        
        while (node_ptr) {
            next_ptr = node_ptr->next;
            remove_matrix_from_mat_store_by_node(node_ptr, i);
            node_ptr = next_ptr;
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
		 fprintf(stderr,"Error: hash value out of range\n");
		 return(0);
	}
	hash_table_t *table_ptr = store.hash_table;
	hash_node_t *node_ptr = table_ptr->heads[hashID];
        hash_node_t *next_ptr;
	 
	 while (node_ptr)  {
		 next_ptr = node_ptr->next;
		 remove_matrix_from_mat_store_by_node(node_ptr, hashID);
		 node_ptr = next_ptr;
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
      // mat_ptr_t record_ptr;
      
      while (node_ptr)  {
        ++counter;
	// record_ptr = (mat_ptr_t) node_ptr->record_ptr;
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


/****************************************************************** 
   The routine get_valString_from_matID_and_coords returns a string version
   of the scalar located in (row,col) of the matrix with matrix mID
   ALERT: It has a memory leak, because it is allocating space for
   the string and not deleting it.
   TODO: perhaps have matrixID in the name.
********************************************************************/
char *get_valString_from_matID_and_coords(int64_t mID, int64_t row, int64_t col)
{
    // get the matrix pointer from the matrixID
    mat_ptr_t mPTR = get_matPTR_from_matID(mID, "", __func__, 0);

    // call matrix pointer version of function with scalarType
    // scalarType val;
    // sca_init(&val);
    scalarType *valPTR = &scratchVars.misc;
    get_valPTR_from_matPTR_and_coords(valPTR, mPTR, row, col);
    char *return_String = sca_get_str(*valPTR);
    // sca_clear(&val);
    return return_String;
}

  
/****************************************************************** 
   The routine get_valPTR_from_matPTR_and_coords is provided a pointer to scalar retScalar.
   It will assign the value to retScalar of the scalar at the matrix position 
   specified by indices (row,col) of the matrix with matrixID mID. 
********************************************************************/
void get_valPTR_from_matPTR_and_coords(scalarType *retScalar, mat_ptr_t mPTR, int64_t row, int64_t col)
{

  mat_ptr_t sPTR = get_valMatPTR_from_matPTR_and_coords(row, col, mPTR);
  sca_set(retScalar, matrix_trace(sPTR));

}



/****************************************************************** 
   The routine get_valID_from_matID_and_coords returns the matrixID of the
   Scalar at the matrix position specified by indices (row,col) 
   of the matrix with matrixID mID. 
   This is the user safe wrapper for get_valMatPTR_from_matPTR_and_coords
********************************************************************/
int64_t get_valID_from_matID_and_coords(int64_t mID, int64_t row, int64_t col)
{
    // get the matrix pointer from the matrixID
    mat_ptr_t mPTR = get_matPTR_from_matID(mID, "", __func__, 0);

    mat_ptr_t scalarPTR;
    
    scalarPTR = get_valMatPTR_from_matPTR_and_coords(row, col, mPTR);

    int64_t scalarID = get_matID_from_matPTR(scalarPTR);

    return(scalarID);
}    
    

// MatrixID wrapper for get_matPTR_from_oldMatPTR_newVal_and_coords. 
int64_t get_matID_from_oldMatID_newValString_and_coords(int64_t mID, int64_t row, int64_t col, char *val)
{
    // get the matrix pointer from the matrixID
    mat_ptr_t m = get_matPTR_from_matID(mID, "", __func__, 0);

    // call matrix pointer version of function with scalarType
    scalarType sca_val;
    sca_init(&sca_val);
    sca_set_str(&sca_val, val);
    mat_ptr_t ret_ptr = get_matPTR_from_oldMatPTR_newVal_and_coords(m, row, col, sca_val);
    sca_clear(&sca_val);

    return get_matID_from_matPTR(ret_ptr);
}


// Create a new matrix that is the same as m but with the entry at
// (row,col) set to v. This requires saving the panel values at each level
// of recursion and modifying the appropriate panel entry for the submatrix
// which has changed due to changing the scalar at (row,col)
mat_ptr_t get_matPTR_from_oldMatPTR_newVal_and_coords(mat_ptr_t m, int64_t row, int64_t col, scalarType val)
{
  // return matrix if it already has val set at (row, col)
  scalarType cur_val;
  sca_init(&cur_val);
  get_valPTR_from_matPTR_and_coords(&cur_val, m, row, col);
  int val_already_at_entry = sca_eq(cur_val, val);
  sca_clear(&cur_val);
  if (val_already_at_entry)
    return m;

  // calculate the number of rows and columns in the matrix
  int64_t num_rows = 1L << matrix_row_level(m);
  int64_t num_cols = 1L << matrix_col_level(m);

  // check to see if the indices for row, col are in span for the matrix
  if (row >= num_rows){
    fprintf(stderr,"ERROR: in get_matPTR_from_oldMatPTR_newVal_and_coords the row %ld is too large for matrix %p\n", row, m);
    exit(1);
  }

  if (col >= num_cols){
    fprintf(stderr,"ERROR: in get_matPTR_from_oldMatPTR_newVal_and_coords the col %ld is too large for matrix %p\n", col, m);
    exit(1);
  }

  mat_ptr_t panels[4];

  // Handle the scalar case
  if (1 == num_rows && 1 == num_cols){
    return get_valMatPTR_from_val(val);
  }

  // Handle the one-row case
  if (1 == num_rows){

    panels[0] = matrix_sub(m, 0);
    panels[1] = matrix_sub(m, 1);
    
    if (2*col < num_cols) 
      panels[0] = get_matPTR_from_oldMatPTR_newVal_and_coords(panels[0], row, col         , val);
    else              
      panels[1] = get_matPTR_from_oldMatPTR_newVal_and_coords(panels[1], row, col - num_cols/2, val);

    return join(panels[0], panels[1]);
  }

  // Handle the one-column case
  if (1 == num_cols){

    panels[0] = matrix_sub(m, 0);
    panels[1] = matrix_sub(m, 2);

    if (2*row < num_rows) panels[0] = get_matPTR_from_oldMatPTR_newVal_and_coords(panels[0], row         , col, val);
    else              panels[1] = get_matPTR_from_oldMatPTR_newVal_and_coords(panels[1], row - num_rows/2, col, val);

    return stack(panels[0], panels[1]);
  }

  // At this point matrices have at least two columns and at least two
  // rows, and the data structure at the next level down looks like four
  // copies of this one.
  for (int quad = 0; quad < 4; quad++) 
    panels[quad] = matrix_sub(m, quad);

  if (2*row < num_rows){
    // element to be changed is in top half of matrix
    if (2*col < num_cols) panels[0] = get_matPTR_from_oldMatPTR_newVal_and_coords(matrix_sub(m, 0), row         , col         , val);
    else              panels[1] = get_matPTR_from_oldMatPTR_newVal_and_coords(matrix_sub(m, 1), row         , col - num_cols/2, val);
  } 
  else {
    // element to be changed is in bottom half of matrix
    if (2*col < num_cols) panels[2] = get_matPTR_from_oldMatPTR_newVal_and_coords(matrix_sub(m, 2), row - num_rows/2, col         , val);
    else panels[3] = get_matPTR_from_oldMatPTR_newVal_and_coords(matrix_sub(m, 3), row - num_rows/2, col - num_cols/2, val);
  }

  return get_matPTR_from_array_of_four_subMatPTRs(panels, matrix_row_level(m), matrix_col_level(m));
}



// ROUTINE TO RETURN A SCALAR FROM A SPECIFIED MATRIX LOCATION
// This is in io.c because it is called from print_naive_by_matPTR,
// write_naive_by_matPTR, and write_matrix_nonzeros_by_matPTR.
// It is not particularly efficient, but for the small matrices
// which are feasible to print, efficiency is not needed.
// TODO: We now also use this to allow python users to find ijth entry
// thus we need reorganize and rename things e.g. getPTR getID ...
/*********************************************************************
 *                  get_valMatPTR_from_matPTR_and_coords                            *
 *                                                                   *
 *   This algorithm finds the ith jth entry of a matrix, by          *
 *   looking at the binary representation of i and j and             *
 *   recursively calling the function on the appropriate sub matrix  *
 *   depending on the (high bit of i, high bit of j):                *
 *     (0,0) -> sub_matrix[0],  (0,1) -> sub_matrix[1],              *
 *     (1,0) -> sub_matrix[1],  (1,1)-> sub_matrix[3]                *
 *                                                                   *
 *   NOTE:  This function returns a mat_ptr_t return_ptr             *
 *          To get the actual value of the i_jth entry  use          *
 *          matrix_trace(return_ptr)                                 *
 *********************************************************************/
mat_ptr_t get_valMatPTR_from_matPTR_and_coords(long int row_i, long int col_j, mat_ptr_t m_ptr) 
{

  // calculate the total number of rows and columns in the matrix m_ptr
  int64_t num_rows = 1L << matrix_row_level(m_ptr);
  int64_t num_cols = 1L << matrix_col_level(m_ptr);

  // check to see if the indices for row, col are in span for the matrix
  if (row_i >= num_rows){
    fprintf(stderr,"ERROR: in %s the row %ld is too large for matrix %p\n",
	    __func__,row_i, m_ptr);
    exit(1);
  }

  if (col_j >= num_cols){
    fprintf(stderr,"ERROR: in %s the col %ld is too large for matrix %p\n",
	    __func__,col_j, m_ptr);
    exit(1);
  }

  
  matrix_type_t mat_type = matrix_type(m_ptr);

  // Case: SCALAR  return the matrix_ptr of that scalar
  if (mat_type == SCALAR) {return m_ptr;
  }

  /*************************************************
   * We want to determine which of the quadrant submatrices
   * contains the row_i,col_j element that we are seeking.
   * If we refer to these four quadrant submatrices as A00, A01, A10, A11
   * then the first bit after A is 0 when the indexed row 
   * is in the first half of the rows of the matrix, and 1
   * when this indexed row is in the second half of rows.
   * Similarly the second bit after A is determined by the columns.
   * Thus A10 for example would be used if the row_i was 
   * at least as big as half the total number of rows,
   * and col_j was less than half the total number of columns.
   * 
   * The new indices used in the recursion are obtained
   * by removing the top bit.
  **************************************************/

  // In all other cases the calculation will proceed recursively
  // using a subpanel of the matrix and new_i and new_j
  long int half_i_dim;  // half size i dimension
  int i_adjust;    // 0 or half_i_dim depending on size of i
  int i_top_bit = 0;
  long int new_i = 0;  // mod off half dim
  long int half_j_dim;  // half size j dimension
  int j_adjust; // (0 or half_j_dim) depending on size of j
  int j_top_bit = 0;
  long int new_j = 0; // mod off half dim

  // Case: mat_type  MATRIX / COL_VECTOR, calculate new_i, and i_top_bit
  int row_level = matrix_row_level(m_ptr);
  if (row_level > 0) {
    half_i_dim = 1 << (row_level-1);  // half size i dimension
    i_adjust = half_i_dim & row_i;   // 0 or half_i_dim depending on size of i
    i_top_bit = i_adjust >> (row_level-1);
    new_i = row_i - i_adjust; // mod off half dim
  } 

  // Case: mat_type  MATRIX / ROW_VECTOR, calculate new_j, and j_top_bit
  int col_level = matrix_col_level(m_ptr);
  if  (col_level > 0) {
    half_j_dim = 1 << (col_level-1);  // half size j dimension
    j_adjust = half_j_dim & col_j;   // (0 or half_j_dim) depending on size of j
    j_top_bit = j_adjust >> (col_level-1);
    new_j = col_j - j_adjust; // mod off half dim
  }
  
  // Calculate the correct submatrix to use in recursive call
  int panel_index = i_top_bit * 2 + j_top_bit;
  return get_valMatPTR_from_matPTR_and_coords(new_i,new_j,matrix_sub(m_ptr,panel_index));

}
 

/************************************************************************
*  The routine get_quad_in_submatrix (by Jenny and Steve)
*  assumes the existence of a big matrix B.  B can be subdivided
*  into a grid of level small_level submatrices.  Only one of 
*  these submatrices of B contains the location (big_row,big_col).  
*  Call this submatrix A. This routine mods off the high bits of 
*  big_row and big_col to get the coordinates (small_row, small_col) 
*  of the location  within A.  Then the routine returns the quad_index 
*  which labels the quadrant submatrix of A that contains the location
*  (small_row,small_col).
*  These quadrant submatrices, as usual, are numbered 0,1,2,3.
************************************************************************/
int get_quad_in_submatrix(int64_t big_row, int64_t big_col,
			  mat_level_t small_level) {
  int64_t small_row = big_row % (1L<<small_level);  /* small_level low bits */
  int64_t small_col = big_col % (1L<<small_level); /* small_level low bits */
  int quad_index = 2*(small_row>>(small_level-1))+(small_col>>(small_level-1));
  return quad_index;
}



/************************************************************************
*  This routine will return a square matrix of the given level with a
*  single nonzero entry specified by the matrix PTR one_scalar, located
*  at position (row,col), where row and col are  zero-indexed.
*  Note: when the user makes this call, i and j need to be smaller than 
*        the dimension of the matrix.  But, when this calls itself, 
*        i and j will not change even as the level is reduced.  This 
*        is all handled gracefully in get_quad_in_submatrix.
************************************************************************/
mat_ptr_t get_matPTR_single_nonzero_using_valPTR_at_coords(mat_level_t level,
	  int64_t row_i, int64_t col_j, scalarType *one_scalar) {
  mat_ptr_t result_ptr;

  // If one_scalar is scalar0 then do not need to do the recursion
  // and can return zero matrix
  if (sca_eq(*one_scalar,scalar0))  return (get_zero_matrix_ptr(level,level));

  // if level> zero, create a matrix from three zero panels, and one panel
  // with the scalar one_scalar in the appropriate position
  mat_ptr_t panel[4];

  // Get the matrix ptr for the zero matrix of one level down
  mat_ptr_t zero_matrix = get_zero_matrix_ptr(level-1,level-1);

  int i;

  // Loop through and intialize all the panels to zero matrices
  for (i=0; i < 4 ; i++) {
    panel[i] = zero_matrix;
  }

  // Find out which panel should contain one_scalar and set by a recursive call
  int panel_number = get_quad_in_submatrix(row_i, col_j, level-1);

  panel[panel_number] = get_matPTR_single_nonzero_using_valPTR_at_coords(level-1, row_i, col_j, one_scalar);
  
  // build a matrix from the four panels
  result_ptr = get_matPTR_from_array_of_four_subMatPTRs(panel, level, level);

  return (result_ptr); 
}

