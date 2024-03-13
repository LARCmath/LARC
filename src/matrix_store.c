//                       matrix_store.c
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


#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <string.h>
#include <stdint.h>   // do not know if we use
#include <float.h>    // do not know if we use

#include "matrix_store.h"
#include "mar.h"
#include "matmath.h"
#include "hash.h"
#include "op_store.h"
#include "json.h"
#include "larc.h"
#include "global.h"
#include "scalars.h"
#include "spr.h"


/*!
 * \file matrix_store.c
 * \brief Contains the static data structure for the matrix store and the
 * functions which allow user access to the store
 */

// we will have a table of tables of this size, ptrs_indexed_by_matrixID
// each entry points to a matns_ptr_table
// WARNING LOG_SIZE_MPT CANNOT BE BIGGER THAN 32
// MPT stands for Matrix PTR Table
#define LOG_SIZE_MPT 25
#define SIZE_MPT ((uint32_t)1<<LOG_SIZE_MPT)
#define MASK_MPT (SIZE_MPT - 1)

//static struct matrixID_pointer_validity
static struct matrixID_tracker
{
   // contains the ID for the next scalarRecord or matrixRecord
   // to be stored - while the Records are different, their
   // matrixIDs come from a common pool and are assigned in order
   uint64_t matrixID_next;

   // if a pointer is freed, the value indexed is set
   // to MATRIX_PTR_INVALID
   record_ptr_t ** ptrs_indexed_by_matrixID;
   
} matrixID_tracker = {0};


static struct matrix_store_t matrix_store = {0};   // this creates the variable "matrix_store" which is a struct of type matrix_store_t

/*!
 * \ingroup larc
 * \brief A utility function which determins the matrix_type enum value given the dimensions of the matrix
 * \param row_level The log base 2 of the number of rows
 * \param col_level The log base 2 of the number of columns
 * \return SCALAR, ROW_VECTOR, COL_VECTOR or MATRIX
 */
static matrix_type_t matrix_type_from_levels(mat_level_t row_level,
                                             mat_level_t col_level)
{
  if ((row_level == 0) && (col_level == 0)) { return SCALAR;}
  if (row_level == 0) { return ROW_VECTOR;}
  if (col_level == 0) { return COL_VECTOR;}
  return MATRIX;
}

matrix_type_t matrix_type(matns_ptr_t m_ptr)
{
  mat_level_t row_level = m_ptr->row_level;
  mat_level_t col_level = m_ptr->col_level;
  return matrix_type_from_levels(row_level,col_level);
}


/*!
 * \ingroup larc
 * \brief finds the indices needed for the ptrs_indexed_by_matrixID table
 *
 * The table relates matrixIDs to their associated pointers. If the matrices
 * are cleared from the matrixStore or scalarStore, the pointers are set to
 * NULL to indicate that this matrixID is no longer valid.
 *
 * \param packedID Contains the matrixID information + some other info
 * \param block_index The index into each allocated block of the table
 * \param index The index within each allocated block of the table
 */
static inline void get_matrixID_table_indices(const int64_t packedID,
	uint64_t *block_index, uint64_t *index)
{
  // the table matrixID_tracker.ptrs_indexed_by_matrixID is doubly dimensioned,
  // with the first dimension having an allocated size of SIZE_MPT and the
  // second dimension allocated only as needed, up to SIZE_MPT (if we try to
  // allocate more than this we get an error). We divide our matrixID value
  // into upper and lower parts to index into this array.
  
  if (packedID == MATRIX_ID_INVALID)
  {
     fprintf(stderr,"in %s, invalid matrixID passed\n",__func__);
     return;
  }

  // first we remove the bit that indicates whether the matrix is a SCALAR
  // or non-SCALAR matrix (set to 1 if SCALAR)
  int64_t matrixID = MID_FROM_PID(packedID);

  // first index into table
  *block_index = matrixID >> LOG_SIZE_MPT;
  // second index into table
  *index = matrixID & MASK_MPT;
  return;
}

int matrix_is_zero(int64_t mat_pID)
{
  if (mat_pID==MATRIX_ID_INVALID) return 0;
  record_ptr_t r_ptr = get_recordPTR_from_pID(mat_pID,"",__func__,0);
  if (IS_SCALAR(mat_pID)) return(((mats_ptr_t)r_ptr)->iszero);
  return(((matns_ptr_t)r_ptr)->iszero);
}

int matrix_is_identity(int64_t mat_pID)
{
  if (mat_pID==MATRIX_ID_INVALID) return 0;
  record_ptr_t r_ptr = get_recordPTR_from_pID(mat_pID,"",__func__,0);
  if (IS_SCALAR(mat_pID)) return(((mats_ptr_t)r_ptr)->isid);
  return(((matns_ptr_t)r_ptr)->isid);
}

int matrix_is_invalid(int64_t mat_pID)
{
  if (mat_pID==MATRIX_ID_INVALID) return 1;
  record_ptr_t r_ptr = get_recordPTR_from_pID(mat_pID,"",__func__,0);
  return (r_ptr == RECORD_PTR_INVALID);
}

mat_level_t matrix_row_level(int64_t mat_pID)
{
    if (IS_SCALAR(mat_pID)) return 0;
    matns_ptr_t mat_ptr = (matns_ptr_t)get_recordPTR_from_pID(mat_pID,
        "", __func__, 0);
    return mat_ptr->row_level;
}

mat_level_t matrix_col_level(int64_t mat_pID)
{
    if (IS_SCALAR(mat_pID)) return 0;
    matns_ptr_t mat_ptr = (matns_ptr_t)get_recordPTR_from_pID(mat_pID,
        "", __func__, 0);
    return mat_ptr->col_level;
}

int64_t get_pID_of_indexed_submatrix(int64_t mat_pID, int s)
{
    if (IS_SCALAR(mat_pID)) return MATRIX_ID_INVALID;
    matns_ptr_t mat_ptr = (matns_ptr_t)get_recordPTR_from_pID(mat_pID,
        "", __func__, 0);
    return mat_ptr->subMatList[s];
}

size_t get_nonscalar_store_exp(void) {
  return matrix_store.nonscalar_hash_table->exponent;
}

size_t get_scalar_store_exp(void) {
  return matrix_store.scalar_hash_table->exponent;
}

mat_level_t max_level_allowed_matrixStore(void) {
  if (VERBOSE==DEBUG)
  {
    printf("in the ONE TRUE max_level_allowed_matrixStore routine\n");
    printf ("matrix_store.largest_level = %d, store is at %p\n", matrix_store.largest_level, &matrix_store);
  }
  return matrix_store.largest_level;
}


int get_regionbitparam(void) {
  return matrix_store.regionbitparam;
}

int get_zeroregionbitparam(void) {
  return matrix_store.zeroregionbitparam;
}

#ifndef MAR // SPRmode

int is_special_zeroregion(void) {
  return matrix_store.special_zeroregion;
}

#ifdef IS_RATIONAL
const mpq_t* const get_zerorealthresh(void) {
  return &(matrix_store.zerorealthresh);
}
#elif defined(IS_MP)
const mpfr_t* const get_zerorealthresh(void) {
  return &(matrix_store.zerorealthresh);
}
#else
double get_zerorealthresh(void) {
  return matrix_store.zerorealthresh;
}
#endif  // IS_RATIONAL

#endif // SPRmode


uint64_t num_matrices_in_store ( void )
{
  uint64_t num_matrices = 0;

  for ( size_t i = 0 ; i <= matrix_store.largest_level ; i++ )
    for ( size_t j = 0 ; j <= matrix_store.largest_level ; j++ )
      num_matrices += matrix_store.hist[i][j];

  return num_matrices;
}

void matrix_store_report(char *outfilepath)
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
	fprintf(f,"Common Store Report:\n");
	fprintf(f,"Longest successful hash chain traversal:  %u\n", matrix_store.deepest);   
	fprintf(f,"Largest allowable level for matrices in Matrix Store: %u \n", 
                   matrix_store.largest_level);
	fprintf(f,"Largest number of scalar matrices ever in Matrix Store at once: %zd\n", 
                   matrix_store.max_num_scalars);
	fprintf(f,"Largest number of nonscalar matrices ever in Matrix Store at once: %zd\n", 
                   matrix_store.max_num_nonscalars);
        size_t total_matrices = 0;
	fprintf(f,"Number of matrices per level (printed only when count greater than 1):\n");
        // int ijsum = 0;
	for (int i = 0; i <= matrix_store.largest_level; i++) {
	  for (int j = 0; j <= matrix_store.largest_level; j++) {
            total_matrices += matrix_store.hist[i][j];
            if (matrix_store.hist[i][j] > 1) {
  	      // printf("(%d,%d): %10zd%c", i, j, matrix_store.hist[i][j],
              //       (ijsum % 4 == 3) ? '\n' : '\t');
  	      fprintf(f,"  (%2d,%2d): %13zd\n", i, j, matrix_store.hist[i][j]);
	    }
            // ++ijsum;
	  }
	}
	fprintf(f,"The total number of matrices currently in the store is %zd.\n", total_matrices);
	fprintf(f,"\n");
	
// TODO: figure out what needs to be updated here to get accurate sizes
	// Also print memory used by matrix pointer table
        // The high order bits of matrixID_next give the table index
	uint64_t bytes_ptr_list = SIZE_MPT * sizeof(matns_ptr_t);   	   			// size of array of pointers to tables
	uint64_t bytes_per_table = SIZE_MPT * sizeof(matns_ptr_t);  	   			// size of each table
        uint64_t num_tables = ((matrixID_tracker.matrixID_next - 1) >> LOG_SIZE_MPT) + 1;    //since indexing starts at 0
	uint64_t bytes_ptr_total = bytes_ptr_list + (num_tables * bytes_per_table);
	fprintf(f, "Memory used by matrix pointer table: \n");
	fprintf(f, "Per table size:  %3.2fMB, Total size:  %3.2fMB \n", bytes_per_table/(1024.0*1024.0), bytes_ptr_total/(1024.0*1024.0));	
	fflush(f);
	fprintf(f,"Matrix Hash Table Report:\n");
        hash_report(matrix_store.nonscalar_hash_table, f, 0);
	fprintf(f,"Scalar Hash Table Report:\n");
        hash_report(matrix_store.scalar_hash_table, f, 0);
	if  (strcmp(outfilepath,"stdout")) fclose(f);   // if not stdout then close
}

int create_matrix_store(size_t exponent, mat_level_t max_level, int regionbitparam, int zeroregionbitparam)
{
  int verbose = 0;
  if (verbose) {
    printf("inside routine %s with:\n", __func__);
    printf("   matrix store hash exponent %zd\n", exponent);
    printf("   maximum matrix level %d\n",max_level);
    printf("   the locality sensitive hash uses regions that are:\n");
    printf("   s = 2^{-%d} wide\n",regionbitparam);
    if (zeroregionbitparam < regionbitparam) {
      printf("   except around zero where they are z= 2^{-%d}-s wide.\n",zeroregionbitparam);
    }
    printf("MATRIX STORE LOCATION = %p\n", &matrix_store);
  }

  if (exponent>=64)
  {
    fprintf(stderr,"You have chosen a matrix store size of 2^64 or larger\n");
    fprintf(stderr,"Too large... exiting\n");
    exit(0);
  }

  matrixID_tracker.matrixID_next = 0;
  
  // two dimensional array of histogram for levels
  matrix_store.hist = (size_t **)malloc((max_level+1)*(sizeof(size_t *)));
  if (matrix_store.hist == NULL) {
    ALLOCFAIL();
    return 0;
  }
  int i;
  for (i=0; i<=max_level; i++) 
    {
      matrix_store.hist[i] = (size_t *) calloc(max_level+1, sizeof(size_t));
      if (matrix_store.hist[i] == NULL) {
	ALLOCFAIL();
	return 0;
      }
    }
 
  // two dimensional array of zeros and one dimensional arrays of identities and iHadamards
  matrix_store.zero = (int64_t **)malloc((1 + max_level) * sizeof(int64_t *));
  matrix_store.identity = (int64_t *)malloc((1 + max_level) * sizeof(int64_t));
  if ((matrix_store.zero == NULL) || (matrix_store.identity == NULL)) {
    ALLOCFAIL();
    exit(1);    
  }
  for (int i= 0; i <= max_level; ++i) {
    matrix_store.zero[i] = (int64_t *)malloc((1 + max_level) * sizeof(int64_t ));
    if (matrix_store.zero[i] == NULL) {
      ALLOCFAIL();
      exit(1);
    }
  }
  matrix_store.iHadamard = (int64_t *)malloc((1 + max_level) * sizeof(int64_t));
  if (matrix_store.iHadamard == NULL) {
    ALLOCFAIL();
    exit(1);
  }
 
  matrix_store.largest_level = max_level;
  if (verbose) {
    printf("inside routine %s with:\n", __func__);
    printf("  set matrix_store.largest_level %d\n",matrix_store.largest_level);
    printf("check routine: %d\n",max_level_allowed_matrixStore());
  }


  matrix_store.regionbitparam = regionbitparam;
  matrix_store.zeroregionbitparam = zeroregionbitparam;

#ifndef MAR // SPRmode
  // store the threshold which is half the width of the zero centered regions
  // this threshold has no imaginary part even with complex scalarTypes
  double zerothreshold;
  if (zeroregionbitparam<regionbitparam-1)
      zerothreshold = pow(2.0,-zeroregionbitparam-1)
			- pow(2.0,-regionbitparam-1);
  else
      zerothreshold = pow(2.0,-regionbitparam-1);

#ifdef IS_RATIONAL
  mpq_init(matrix_store.zerorealthresh);
  mpq_set_d(matrix_store.zerorealthresh, zerothreshold);
#elif defined(IS_MP)
  mpfr_init2(matrix_store.zerorealthresh, mpreal_precision);
  mpfr_set_d(matrix_store.zerorealthresh, zerothreshold, MPFR_RNDN);
#else
  matrix_store.zerorealthresh = zerothreshold;
#endif

  // make special SPR regions in vicinity of zero?
  matrix_store.special_zeroregion = (zeroregionbitparam < regionbitparam);
#endif // SPRmode
  
  matrix_store.nonscalar_hash_table = alloc_hash(exponent);
  matrix_store.scalar_hash_table = alloc_hash(exponent);

  matrix_store.nonscalar_hash_table->record_size = sizeof(larc_nsmatrix_t);
  if (matrix_store.nonscalar_hash_table == NULL) {
    ALLOCFAIL();
    return 0;
  }

  matrix_store.scalar_hash_table->record_size = sizeof(larc_smatrix_t);
  if (matrix_store.scalar_hash_table == NULL) {
    ALLOCFAIL();
    return 0;
  }

  matrix_store.num_nonscalars = 0;
  matrix_store.max_num_nonscalars = 0;
  matrix_store.num_scalars = 0;
  matrix_store.max_num_scalars = 0;

  // create list of tables called ptrs_indexed_by_uniqeID
  // and initialize the first table in this list 
  // with zeros corresponding to mat_ptr MATRIX_PTR_INVALID
  // ToDo: may want to have 2 sizes: one for the number of tables, and one for the table size
  matrixID_tracker.ptrs_indexed_by_matrixID = (record_ptr_t **)calloc(SIZE_MPT,sizeof(record_ptr_t *));  //list of ptrs to tables
  matrixID_tracker.ptrs_indexed_by_matrixID[0] = (record_ptr_t *)calloc(SIZE_MPT,sizeof(record_ptr_t));  //the first table
  if ( (matrixID_tracker.ptrs_indexed_by_matrixID == NULL) ||
       (matrixID_tracker.ptrs_indexed_by_matrixID[0] == NULL) )
  {
    ALLOCFAIL();
    exit(1);
  }

  return 1;
}

/************************************************************
 *                   lock_matrix                            *
 *  These routines are used to prevent important matrices   *
 *  (such as identity matrices, and zero matrices required  *
 *  for the math routines to work) from EVER being removed  *
 *  from the matrix store.                                  *
 *  A lock should never be removed.                         *
 *  NOTE: If you want to temporarily hold a matrix use      *
 *        set_hold_matrix, and use            *
 *        release_hold_matrix to remove       *
 *        the hold.                                         *
 *                                                          *
 ***********************************************************/
int lock_matrix(int64_t m_pID) {
  record_ptr_t rec_ptr;
  check_validity_one_input(m_pID, __func__, &rec_ptr);

  if (IS_SCALAR(m_pID))
  {
    mats_ptr_t s_ptr = (mats_ptr_t)rec_ptr;
    if (s_ptr->lock == 1) 
    {
      if (VERBOSE>DEBUG) {
        fprintf(stderr,"NOTE: in %s, attempt to lock already-locked matrix\n",
                __func__);
      }
    }
    else s_ptr->lock = 1;
    return (1);
  }

  matns_ptr_t m_ptr = (matns_ptr_t)rec_ptr;
  if (m_ptr->lock == 1)
  {
    if (VERBOSE>DEBUG) {
      fprintf(stderr,"NOTE: in %s, attempt to lock already-locked matrix\n",
              __func__);
      fprintf(stderr,"assuming that any submatrices are also locked...\n");
    }
  }
  else m_ptr->lock = 1;
  
  return (1);
}   


/************************************************************
 *                   set_hold_matrix	    *
 *  Python Interface version:accepts a	    *
 *  packedID instead of a pointer          *
 *  This routine is used to temporarily hold a matrix       *
 *  so it is not removed from the matrix store.             *
 *  Use the function release_hold_matrix to undo the hold.  *
 *  NOTE: lock_smatrixPTR, lock_nsmatrixPTRs set permanent  *
 *  locks which are used to prevent important matrices      *
 *  (such as identity matrices, and zero matrices required  *
 *  for the math routines to work) from EVER being removed  *
 *  from the matrix store.                                  *
 *  A lock should never be removed, a hold can be.          *
 *                                                          *
 ***********************************************************/
int set_hold_matrix(int64_t m_pID) {

  record_ptr_t r_ptr;
  check_validity_one_input(m_pID,__func__,&r_ptr);

  if (IS_SCALAR(m_pID))
  {
    mats_ptr_t s_ptr = (mats_ptr_t)r_ptr;
    if (s_ptr->lock == 1)
    {
      if (VERBOSE>DEBUG) fprintf(stderr,"no need to hold a locked matrix\n");
    }
    // check and warning added on 22Mar22 - at this time, hold is a 13-bit
    // subfield of an unsigned int within the scalar/nonscalar matrix structs
    else if (s_ptr->hold == 0x1fff) // custom for 13-bit field
    {
      // for now, we convert the hold into a lock and print a warning
      s_ptr->hold = 0;
      s_ptr->lock = 1;
      fprintf(stderr,"in %s: hold = 0x1fff would be incremented,", __func__);
      fprintf(stderr," converting to lock\n");
      fprintf(stderr,"\tmatrixID = %ld, row_level = %d, col_level = %d\n",
        MID_FROM_PID(m_pID), matrix_row_level(m_pID), matrix_col_level(m_pID));
    }
    else ++(s_ptr->hold);
    return (1);
  }
  // NONSCALAR case
  matns_ptr_t m_ptr = (matns_ptr_t)r_ptr;
  if (m_ptr->lock == 1) 
  {
    if (VERBOSE>DEBUG) fprintf(stderr,"no need to hold a locked matrix\n");
  }
  // check and warning added on 22Mar22 - at this time, hold is a 13-bit
  // subfield of an unsigned int within the scalar/nonscalar matrix structs
  else if (m_ptr->hold == 0x1fff) // custom for 13-bit field
  {
    // for now, we convert the hold into a lock and print a warning
    m_ptr->hold = 0;
    m_ptr->lock = 1;
    fprintf(stderr,"in %s: hold = 0x1fff would be incremented,", __func__);
    fprintf(stderr," converting to lock\n");
    fprintf(stderr,"\tmatrixID = %ld, row_level = %d, col_level = %d\n",
      MID_FROM_PID(m_pID), matrix_row_level(m_pID), matrix_col_level(m_pID));
  }
  else ++(m_ptr->hold);
  return (1);
}

int release_hold_matrix(int64_t m_pID) {

  record_ptr_t r_ptr;
  check_validity_one_input(m_pID,__func__,&r_ptr);

  if (IS_SCALAR(m_pID))
  {
    mats_ptr_t s_ptr = (mats_ptr_t)r_ptr;
    if (s_ptr->lock)
    {
      if (VERBOSE>DEBUG)
        fprintf (stderr,"NOTE: no holds needed on locked matrices\n");
    }
    else if (s_ptr->hold == 0)
    {
      if (VERBOSE>SILENT)
        fprintf (stderr,"WARN: attempt %s for matrix with no holds\n",__func__);
    }
    else --(s_ptr->hold);
    return (1);
  }
  // NONSCALAR case  
  matns_ptr_t m_ptr = (matns_ptr_t)r_ptr;
  if (m_ptr->lock)
  {
    if (VERBOSE>DEBUG)
      fprintf (stderr,"NOTE: no holds needed on locked matrices\n");
  }
  else if (m_ptr->hold == 0)
  {
    if (VERBOSE>SILENT) 
      fprintf (stderr,"WARN: attempt %s for matrix with no holds\n",__func__);
  }
  else --(m_ptr->hold);
  return (1);
}

/*!
 * \ingroup larc
 * \brief Removes a matrix pointer from the matrixID table
 *
 * This is done when cleaning has removed a matrix from the matrix store;
 * the pointer is set to NULL in the table
 *
 * \param packedID the ID for which the pointer is to be NULLed
 * \return 1
 */
int invalidate_recordPTR_in_indexTable(int64_t packedID)
{
  // check that this is not MATRIX_ID_INVALID
  if (packedID == MATRIX_ID_INVALID) return(1);
  uint64_t block_index, index;
  get_matrixID_table_indices(packedID, &block_index, &index);
  if (matrixID_tracker.ptrs_indexed_by_matrixID[block_index][index] ==
           RECORD_PTR_INVALID)
    fprintf(stderr,"in %s, invalidating previously invalid entry\n",__func__);
// else
//    fprintf(stderr,"in %s, invalidating packedID %ld \n",__func__,packedID);
  matrixID_tracker.ptrs_indexed_by_matrixID[block_index][index] =
           RECORD_PTR_INVALID;
  return(1);
}

/// NON-SCALAR

/*!
 * \ingroup larc
 * \brief A utility function which increments the number of times a matrix is a submatrix of a larger matrix (if greater than zero, the matrix will not be cleaned)
 * \param id A pointer to a matrix
 * \returns The number of times that the given matrix is a quadrant submatrix
 */
static uint32_t matrix_appears_as_sub_count_increment_ptr(matns_ptr_t id)
{
  if (id->lock == 1)
  {
     // if the matrix is locked, there is no need to count how many times it
     // is used, since it will not be cleaned
     // could return(id->appears_as_sub_count) but complier will simplify
     // the returned constant more easily
     return(0);
  }
  if ( ++(id->appears_as_sub_count) == 0) {
    printf("ERROR: increment of sub_count caused overflow\n");
    exit(1);
  }
  return id->appears_as_sub_count;
}

/*!
 * \ingroup larc
 * \brief A utility function which increments the number of times a matrix is a submatrix of a larger matrix (if greater than zero, the matrix will not be cleaned)
 * \param id A pointer to a matrix
 * \returns The number of times that the given matrix is a quadrant submatrix
 */
static uint32_t scalar_appears_as_sub_count_increment_ptr(mats_ptr_t id)
{
  if (id->lock == 1)
  {
     // if the matrix is locked, there is no need to count how many times it
     // is used, since it will not be cleaned
     // could return(id->appears_as_sub_count) but complier will simplify
     // the returned constant more easily
     return(0);
  }
  if ( ++(id->appears_as_sub_count) == 0) {
    printf("ERROR: increment of sub_count caused overflow\n");
    exit(1);
  }
  return id->appears_as_sub_count;
}

/*!
 * \ingroup larc
 * \brief A utility function which decrements the number of times a matrix is a submatrix of a larger matrix (if greater than zero, the matrix will not be cleaned)
 * \param id A pointer to a matrix
 * \returns The number of times that the given matrix is a quadrant submatrix
 */
static uint32_t matrix_appears_as_sub_count_decrement_ptr(matns_ptr_t id)
{
  if (id->lock == 1)
  {
     // if the matrix is locked, there is no need to count how many times it
     // is used, since it will not be cleaned
     // could return(id->appears_as_sub_count) but complier will simplify
     // the returned constant more easily
     return(0);
  }
  if (id->appears_as_sub_count==0) {
    printf("ERROR in %s: trying to decrement sub_count of zero\n",__func__);
    printf("id = %p, packedID = %ld\n",id,id->packedID);
    exit(1);
  }
  return --(id->appears_as_sub_count);
}

/*!
 * \ingroup larc
 * \brief A utility function which decrements the number of times a matrix is a submatrix of a larger matrix (if greater than zero, the matrix will not be cleaned)
 * \param id A pointer to a matrix
 * \returns The number of times that the given matrix is a quadrant submatrix
 */
static uint32_t scalar_appears_as_sub_count_decrement_ptr(mats_ptr_t id)
{
  if (id->lock == 1)
  {
     // if the matrix is locked, there is no need to count how many times it
     // is used, since it will not be cleaned
     // could return(id->appears_as_sub_count) but complier will simplify
     // the returned constant more easily
     return(0);
  }
  if (id->appears_as_sub_count==0) {
    printf("ERROR in %s: trying to decrement sub_count of zero\n",__func__);
    printf("id = %p, packedID = %ld\n",id,id->packedID);
    exit(1);
  }
  return --(id->appears_as_sub_count);
}

/* Returns 0 if matrix differs from provided content, non-zero if equal */
/*!
 * \ingroup larc
 * \brief
 * \param m_ptr The pointer to a stored matrix
 * \param subMatList A set of four matrixIDs (which define a larger matrix)
 * \return 1 if the matrices are the same, 0 if not
 */
static int matrix_compare_subMatList(matns_ptr_t m_ptr, int64_t subMatList[4])
{
  // EASY CASE for failed comparison
  if (m_ptr == MATRIX_PTR_INVALID) return 0; 

  return (subMatList[0] == m_ptr->subMatList[0] &&
          subMatList[1] == m_ptr->subMatList[1] &&
          subMatList[2] == m_ptr->subMatList[2] &&
          subMatList[3] == m_ptr->subMatList[3] ); 
}

/*!
 * \ingroup larc
 * \brief Finds the packedID for a given matrix as described by a list of matrixIDs for submatrices (does not add it to the matrix store if not found)
 * \param hash A value indicating which hash chain would contain the record
 * \param subMatList An array of four matrixIDs
 * \return The packedID of the stored matrix, or MATRIX_ID_INVALID if not found
 */
static int64_t matrix_find_subMatList(uint64_t hash, int64_t subMatList[4])
{
    int verbose = 0;  

    if (verbose) printf("In routine %s\n", __func__);

    hash_node_t *n = hash_get_chain(matrix_store.nonscalar_hash_table, hash);
    uint32_t depth = 0;

#ifdef HASHSTATS
    (matrix_store.nonscalar_hash_table->num_accesses[hash])++;
#endif

    while (n) {
       int64_t m_pID = n->packedID;
       matns_ptr_t m_ptr =
             (matns_ptr_t)get_recordPTR_from_pID(m_pID,"",__func__,0);
       if (matrix_compare_subMatList(m_ptr, subMatList))
       {
           matrix_store.nonscalar_hash_table->hits++;
           if (!n->hits_maxxed) {
               n->record_hits++;
               if (n->record_hits == ((uint32_t) (-1))) n->hits_maxxed = 1;
           }

           if (depth > matrix_store.deepest) matrix_store.deepest = depth;

    // This code was a start at working on optimizing the hash chain order
    // It is probably obsolete at this point, but we are leaving it in place
    // in case it is a good start point for a future optimization
//			if (n != matrix_store.hash[hash]) {
//				struct node *head = matrix_store.hash[hash];
//				if (n->prev) n->prev->next = n->next;
//				if (n->next) n->next->prev = n->prev;
//				n->prev = NULL;
//				n->next = head;
//				if (head)
//				  head->prev = n;
//				matrix_store.hash[hash] = n;
//			}

           return m_pID;
        }
	n = n->next;
	depth++;
    }

    matrix_store.nonscalar_hash_table->misses++;

    return MATRIX_ID_INVALID;
}

/*!
 * \ingroup larc
 * \brief Add an entry to the table indexing matrix pointer values by matrixID
 * \param r_ptr A pointer to a matrix or scalar
 * \param packedID The packedID assigned to that matrix or scalar
 * \return The packedID for that pointer
 */
static int64_t add_recordPTR_to_indexTable(record_ptr_t r_ptr, int64_t packedID)
{
  // This function adds an entry to the table indexing matrix pointer values
  // by matrixID values. In the future, when we are more actively purging
  // matrices from the matrix store, this may be replaced with a hash.

  // confirm that the pointer is not MATRIX_PTR_INVALID
  if (r_ptr==RECORD_PTR_INVALID)
  {
    fprintf(stderr,"ERROR: pointer passed to %s is invalid!\n",
            __func__);
    exit(0);
  }

  // confirm that packedID assigned to pointer is valid
  if (packedID == MATRIX_ID_INVALID)
  {
    fprintf(stderr,"ERROR: pointer passed to %s has invalid matrixID!\n",
            __func__);
    exit(0);
  }

  // find indices into table of pointers indexed by matrixID
  uint64_t block_index, index;
  get_matrixID_table_indices(packedID, &block_index, &index);

  // TODO: this test may be redundant if we have the test when creating matrixIDs
  if (block_index >= SIZE_MPT)
  {
    fprintf(stderr,"ERROR: LOG_SIZE_MPT needs to be larger.\n");    // or, can modify code to allocate more tables as needed
    exit(1);
  }   

  // if necessary, allocate memory for the next table in the list
  if (matrixID_tracker.ptrs_indexed_by_matrixID[block_index] == 0)
  {
    matrixID_tracker.ptrs_indexed_by_matrixID[block_index] =
	(record_ptr_t *)calloc(SIZE_MPT,sizeof(record_ptr_t));
    if (matrixID_tracker.ptrs_indexed_by_matrixID[block_index]==NULL)
	{ ALLOCFAIL(); }
  }

  // with setup complete, add pointer to the appropriate table
  matrixID_tracker.ptrs_indexed_by_matrixID[block_index][index] = r_ptr;
  return(packedID);
}


/* Inserts a new matrix into the store and returns the new packedID */
/*!
 * \ingroup larc
 * \brief Inserts a matrix into the matrix store
 * \param hash The hash for the panel describing the matrix
 * \param rpanel An array of four record pointers; the quadrant submatrices of the matrix to be stored
 * \param row_level The row level of the matrix to be stored
 * \param col_level The column level of the matrix to be stored
 * \return The packedID for the stored matrix
 */
static int64_t
matrix_insert_panel(uint64_t hash, record_ptr_t rpanel[4], mat_level_t row_level, mat_level_t col_level)
{
  // EXCEPTION TESTING
  if (MAX(row_level,col_level) > matrix_store.largest_level){
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
  matns_ptr_t m_ptr = calloc(1, sizeof(larc_nsmatrix_t));
  if (m_ptr == NULL) {
    ALLOCFAIL();
    return MATRIX_ID_INVALID;
  }
  
  if (matrixID_tracker.matrixID_next+1==0)
  {
    fprintf(stderr,"ERROR: 64-bit matrixID_next about to overflow!\n");
    exit(0);
  }

  // set the basic parameters for the matrix
  // since this is a nonscalar matrix, the matrixID will not
  // have the scalar bit set
  int64_t nextID = matrixID_tracker.matrixID_next++;
  int64_t packedID = PID_FROM_NONSCALAR_MID(nextID);
  m_ptr->packedID = packedID;
  m_ptr->row_level = row_level;
  m_ptr->col_level = col_level;
  if (VERBOSE>=DEBUG) {
    fprintf(stderr,"adding nonscalar with matrixID %ld", nextID);
    fprintf(stderr," and packedID %ld", m_ptr->packedID);
    fprintf(stderr," to chain with hash %lu\n",hash);
    fprintf(stderr,"\trecord pointer is %p\n",m_ptr);
  }

  if (verbose) { 
    printf("  allocated space, set row level, col level, and matrixID \n");
  }
  
#ifndef HASH_CHAIN_GROWS_AT_HEAD
  // INSERT new items at tail of hash collision chain
  hash_insert_at_tail(matrix_store.nonscalar_hash_table, RECORD_PTR_INVALID,
                      packedID, hash);
#endif
#ifdef HASH_CHAIN_GROWS_AT_HEAD
  // INSERT new items at head of hash collision chain
  hash_insert_at_head(matrix_store.nonscalar_hash_table, RECORD_PTR_INVALID,
                      packedID, hash);
#endif

  if (verbose) { 
    printf("  appended hash chain with %d, %d levels\n",row_level,col_level);
  }
  
  // update matrix store parameters
  matrix_store.hist[row_level][col_level]++;
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

  else if (mat_type == MATRIX)
  {
    if ((row_level==1) && (col_level==1))
    {
        // rpanel has scalar pointers
        mats_ptr_t *panel = (mats_ptr_t *)rpanel;
        // define new matrix in terms of its four panels
        for (int i = 0; i < 4; i++) {
            if (panel[i] == NULL)
            {
                fprintf(stderr,"Error looking up submatrices\n");
                exit(1);
            }
            // each panel is a submatrix of the new matrix
            m_ptr->subMatList[i] = panel[i]->packedID;
            scalar_appears_as_sub_count_increment_ptr(panel[i]);
        }
        // testing for zero or identity matrices
        if (panel[1]->iszero && panel[2]->iszero) {
          if (panel[0]->iszero && panel[3]->iszero) {
	      // [0, 0; 0, 0] 
              m_ptr->iszero = 1;
          }
          else if (panel[0]->isid && panel[3]->isid) {
              // [I, 0; 0, I]; only possible for square matrices
              m_ptr->isid = 1;
          }
        }
        // trace of matrix is sum of traces of the two diagonal blocks
        sca_init(&(m_ptr->trace_element));
        sca_add(&(m_ptr->trace_element),
            panel[0]->scalar_value, panel[3]->scalar_value);
    } // end case where panels are scalar
    else 
    {
        // rpanel has matrix pointers
        matns_ptr_t *panel = (matns_ptr_t *)rpanel;

        // Exception checking
        if( (panel[0]->row_level != panel[1]->row_level) || 
            (panel[0]->row_level != panel[2]->row_level) ||
            (panel[0]->row_level != panel[3]->row_level) )
        {
            fprintf(stderr,"ERROR: In %s ",__func__);
            fprintf(stderr," trying to insert matrix with non-equal row");
            fprintf(stderr," sub-levels %d %d %d %d\n",
		panel[0]->row_level, panel[1]->row_level,
		panel[2]->row_level, panel[3]->row_level );
            exit(1);
        }

        if( (panel[0]->col_level != panel[1]->col_level) || 
            (panel[0]->col_level != panel[2]->col_level) ||
            (panel[0]->col_level != panel[3]->col_level) )
        {
            fprintf(stderr,"ERROR: In %s ",__func__);
            fprintf(stderr,"trying to insert matrix with non-equal col");
            fprintf(stderr," sub-levels %d %d %d %d\n",
                panel[0]->col_level, panel[1]->col_level,
                panel[2]->col_level, panel[3]->col_level );
            exit(1);
        }

        // define new matrix in terms of its four panels
        for (int i = 0; i < 4; i++) {
          if (panel[i] == NULL) {
	    fprintf(stderr,"Error looking up submatrices\n");
	    exit(1);
          }
          // each panel is a submatrix of the new matrix
          m_ptr->subMatList[i] = panel[i]->packedID;
          matrix_appears_as_sub_count_increment_ptr(panel[i]);
        }
        // testing for zero or identity matrices
        if (panel[1]->iszero && panel[2]->iszero) {
          if (panel[0]->iszero && panel[3]->iszero) {
	      // [0, 0; 0, 0] 
              m_ptr->iszero = 1;
          }
          else if (panel[0]->isid && panel[3]->isid) {
              // [I, 0; 0, I]; only possible for square matrices
              m_ptr->isid = 1;
          }
        }

        // Inside of the larc_smatrix_t there is a union which holds
        // either the trace_element (for nonscalars) or scalar_value
        // special SQUARE MATRIX stuff
        if (row_level == col_level) {
            // trace of matrix is sum of traces of the two diagonal blocks
            // we use trace_element because the submatrices are not scalars
            sca_init(&(m_ptr->trace_element));
            sca_add(&(m_ptr->trace_element),
                    panel[0]->trace_element, panel[3]->trace_element);
        } // end square matrix tests
    } // end case where panels are nonscalar
  }  // end MATRIX type

  else if (mat_type == COL_VECTOR) {

    if (row_level==1)
    {
       // panel has scalar pointers 
       mats_ptr_t *panel = (mats_ptr_t *)rpanel;
       m_ptr->subMatList[0] = panel[0]->packedID;
       m_ptr->subMatList[2] = panel[2]->packedID;
       m_ptr->subMatList[1] = m_ptr->subMatList[3] = MATRIX_ID_INVALID;
       scalar_appears_as_sub_count_increment_ptr(panel[0]);
       scalar_appears_as_sub_count_increment_ptr(panel[2]);

       if (panel[0]->iszero && panel[2]->iszero) m_ptr->iszero = 1;
    } // end case submatrices are scalar
    else
    {
       // panel has matrix pointers
       matns_ptr_t *panel = (matns_ptr_t *)rpanel;

       // Exception checking
       if( panel[0]->col_level != panel[2]->col_level )
       {
	   fprintf(stderr,"ERROR: In %s trying to insert column vector with non-equal col sub-levels %d %d\n", __func__,
		panel[0]->col_level, panel[2]->col_level );
	   exit(1);
       }
       if ((panel[0] == NULL) || (panel[2] == NULL)) {
         fprintf(stderr,"Error looking up submatrices\n");
         exit(1);
            // return MATRIX_PTR_INVALID;
       }

       m_ptr->subMatList[0] = panel[0]->packedID;
       m_ptr->subMatList[2] = panel[2]->packedID;
       m_ptr->subMatList[1] = m_ptr->subMatList[3] = MATRIX_ID_INVALID;
       matrix_appears_as_sub_count_increment_ptr(panel[0]);
       matrix_appears_as_sub_count_increment_ptr(panel[2]);
    
    if (panel[0]->iszero && panel[2]->iszero) m_ptr->iszero = 1;
   } // end case submatrices are not scalar
  }  // end of COL_VECTOR type
 
  else if (mat_type == ROW_VECTOR) {

    if (col_level==1)
    {
       // panel has scalar pointers 
       mats_ptr_t *panel = (mats_ptr_t*)rpanel;

       m_ptr->subMatList[0] = panel[0]->packedID;
       m_ptr->subMatList[1] = panel[1]->packedID;
       m_ptr->subMatList[2] = m_ptr->subMatList[3] = MATRIX_ID_INVALID;
       scalar_appears_as_sub_count_increment_ptr(panel[0]);
       scalar_appears_as_sub_count_increment_ptr(panel[1]);

       if (panel[0]->iszero && panel[1]->iszero) m_ptr->iszero = 1;
    } // end submatrices are scalar
    else
    {
       matns_ptr_t *panel = (matns_ptr_t* )rpanel;

       // exception checking
       if( panel[0]->row_level != panel[1]->row_level )
       {
	   fprintf(stderr,"ERROR: In %s",__func__);
           fprintf(stderr,"trying to insert row vector with non-equal ");
           fprintf(stderr,"row sub-levels %d %d\n",
		   panel[0]->row_level, panel[1]->row_level );
	   exit(1);
       }
       if ((panel[0] == NULL) || (panel[1] == NULL)) {
	   fprintf(stderr,"ERROR: In %s",__func__);
           fprintf(stderr,"Error looking up submatrices\n");
           exit(1);
           // return MATRIX_PTR_INVALID;
       }
   
       m_ptr->subMatList[0] = panel[0]->packedID;
       m_ptr->subMatList[1] = panel[1]->packedID;
       m_ptr->subMatList[2] = m_ptr->subMatList[3] = MATRIX_ID_INVALID;
       matrix_appears_as_sub_count_increment_ptr(panel[0]);
       matrix_appears_as_sub_count_increment_ptr(panel[1]);

       if (panel[0]->iszero && panel[1]->iszero) m_ptr->iszero = 1;
   } // end case submatrices are nonscalar
  }  // end of ROW_VECTOR type

  if (matrix_store.num_nonscalars+1<matrix_store.num_nonscalars)
  {
    fprintf(stderr,"ERROR: overflow in matrix_store.num_nonscalars\n");
  }

  matrix_store.num_nonscalars++;
  if( matrix_store.num_nonscalars > matrix_store.max_num_nonscalars ) {
      matrix_store.max_num_nonscalars = matrix_store.num_nonscalars;
  }
  
  // Record the matrix pointer in table indexed by matrixID
  add_recordPTR_to_indexTable((record_ptr_t)m_ptr,m_ptr->packedID);
  
  return packedID;
}


/************************************************************************
 *               preload_matrix_and_scalar_stores                                   *
 *   This should not be done until the op_store has been initialized.   *
 *   This function loads all the basic 1 and 2 bit matrices.            *
 *   It then loads all the identity matrices and rectangular zero       *
 *   matrices up to the maximum level size store.largest_level. It      *
 *   also loads all the integer Hadamard matrices up to largest_level.  *
 ***********************************************************************/   
int preload_matrix_and_scalar_stores(void) {

  int verbose = 0;
 
  int top_level = matrix_store.largest_level;
  if (verbose) {
    printf("Inside routine %s with:\n",__func__);
    printf("   retrieved largest level %d\n", top_level);
  }
  
  // 0, 1, -1 are constants needed for preload of zero, identity, integer Hadamard
  // NOTE: scalar ops must be instantiated before preloading matrix store!
  scalarType scalar0;
  scalarType scalar1;
  scalarType scalarM1;
  sca_init(&scalar0);
  sca_init(&scalar1);
  sca_init(&scalarM1);

  sca_set_2ldoubles(&scalar0, 0.0L, 0.0L);
  sca_set_2ldoubles(&scalar1, 1.0L, 0.0L);
  sca_set_2ldoubles(&scalarM1, -1.0L, 0.0L);
  if (verbose) {
    printf("Inside routine %s with:\n",__func__);
    char *scalar0_string = sca_get_readable_approx_str(scalar0);
    char *scalar1_string = sca_get_readable_approx_str(scalar1);
    char *scalarM1_string = sca_get_readable_approx_str(scalarM1);
    printf("loaded %s, %s, %s.\n", scalar0_string, scalar1_string, scalarM1_string);
    free(scalarM1_string);
    free(scalar1_string);
    free(scalar0_string);
  }

  // panels of sub matrices
  int64_t sm[4];

  // LOAD ALL THE ZERO MATRICES
  matrix_store.zero[0][0] = get_scalarPTR_for_scalarVal(scalar0)->packedID;
  if (!lock_matrix(matrix_store.zero[0][0])) {
      fprintf(stderr,"FAIL: scalar zero preload tried to lock a matrix which does not exist\n");
      exit(1);
    }
  if (verbose) {
    printf("Inside routine %s:\n",__func__);
    printf("   preloaded store.zero[0][0]\n");
  }
  
  for (int k = 1; k <= top_level; k++) {

    // Case: ROW_VECTOR  (0,2**k)
    sm[0] = matrix_store.zero[0][k-1];
    sm[1] = matrix_store.zero[0][k-1];
    sm[2] = MATRIX_ID_INVALID;
    sm[3] = MATRIX_ID_INVALID;
    matrix_store.zero[0][k] = get_pID_from_array_of_four_sub_pIDs(sm, 0, k);
    if (!lock_matrix(matrix_store.zero[0][k])) {
      fprintf(stderr,"FAIL: row Z preload tried to lock a matrix which does not exist\n");
      exit(1);
    }

    // Case: COL_VECTOR (2**k,0)
    sm[0] = matrix_store.zero[k-1][0];
    sm[1] = MATRIX_ID_INVALID;
    sm[2] = matrix_store.zero[k-1][0];
    sm[3] = MATRIX_ID_INVALID;
    matrix_store.zero[k][0] = get_pID_from_array_of_four_sub_pIDs(sm, k, 0);
    if (!lock_matrix(matrix_store.zero[k][0])) {
      fprintf(stderr,"FAIL: col Z preload tried to lock a matrix which does not exist\n");
      exit(1);
    }

  }
  if (verbose) {
    printf("Inside routine %s:\n",__func__);
    printf("   preloaded row and column zero vectors into matrix_store.zero\n");
  }


  // Case: MATRIX  (2**i,2**j) with (0 < i,j <= k) 
  for (int i = 1; i <= top_level; i++) {
    for (int j = 1; j <= top_level; j++) {
      sm[0] = matrix_store.zero[i-1][j-1];
      sm[1] = matrix_store.zero[i-1][j-1];
      sm[2] = matrix_store.zero[i-1][j-1];
      sm[3] = matrix_store.zero[i-1][j-1];
      matrix_store.zero[i][j] = get_pID_from_array_of_four_sub_pIDs(sm,i,j);
      if (!lock_matrix(matrix_store.zero[i][j])) {
	fprintf(stderr,"FAIL: matrix Z preload tried to lock a matrix which does not exist\n");
	exit(1);
      }
    }
  }

  // Prestore all the IDENTITY matrices (square)
  matrix_store.identity[0] = get_scalarPTR_for_scalarVal(scalar1)->packedID;
  if (!lock_matrix(matrix_store.identity[0])) {
      fprintf(stderr,"FAIL: scalar 1 preload tried to lock a matrix which does not exist\n");
      exit(1);
    }


  for (int i = 1; i <= top_level; i++) {
    sm[1] = sm[2] = matrix_store.zero[i-1][i-1];
    sm[0] = sm[3] = matrix_store.identity[i-1];
    matrix_store.identity[i] = get_pID_from_array_of_four_sub_pIDs(sm, i, i);
    if (!lock_matrix(matrix_store.identity[i])) {
      fprintf(stderr,"FAIL: I preload tried to lock a matrix which does not exist\n");
      exit(1);
    }
  }

  if (verbose) {
    printf("Inside routine %s:\n",__func__);
    printf("   preloaded all identity matrices into matrix_store.identity\n");
  }


  // Prestore all the integer HADAMARD matrices and scalar negative one.

  // Load scalar negative one into the otherwise un-used store.iHadamard[0]
  matrix_store.iHadamard[0] = get_scalarPTR_for_scalarVal(scalarM1)->packedID;
  // if -1 == 1 (eg. in boolean), avoid warning about locking a locked matrix. 
  if (sca_eq(scalarM1, scalar1) == 0){
      if (!lock_matrix(matrix_store.iHadamard[0])) {
          fprintf(stderr,"FAIL: scalar negative one preload tried to lock a nonexistent matrix\n");
          exit(1);
      }
  }

  sm[0] = sm[1] = sm[2] = matrix_store.identity[0];
  sm[3] = matrix_store.iHadamard[0];
  matrix_store.iHadamard[1] = get_pID_from_array_of_four_sub_pIDs(sm, 1, 1);
  if (!lock_matrix(matrix_store.iHadamard[1])) {
      fprintf(stderr,"FAIL: iHadamard[1] preload tried to lock a nonexistent matrix\n");
      exit(1);
    }

  if (verbose) {
    printf("Inside routine %s:\n",__func__);
    printf("   preloaded iHadamard[1]\n");
  }

  for (int i = 2; i <= top_level; i++) {
    sm[0] = sm[1] = sm[2] = matrix_store.iHadamard[i-1];
    sm[3] = scalar_mult(matrix_store.iHadamard[0], matrix_store.iHadamard[i-1]);
    matrix_store.iHadamard[i] = get_pID_from_array_of_four_sub_pIDs(sm,i,i);
    if (!lock_matrix(matrix_store.iHadamard[i])) {
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
    matns_ptr_t mat_ptr = va_arg(mat_ptrs, matns_ptr_t);
    if (mat_ptr == MATRIX_PTR_INVALID){
        at_least_one_invalid = 1;
        fprintf(stderr,"ERROR: In %s the %dth (0-up) matrix checked was not a valid matrix.\n", function, i);
    }
  }
  va_end(mat_ptrs);

  if (at_least_one_invalid)
    exit(EXIT_FAILURE);
}

int64_t get_zero_pID(mat_level_t row_level, mat_level_t col_level)
{
  if (row_level==0 && col_level == 0) return packedID_scalar0;
  mat_level_t highestInputLevel = MAX(row_level, col_level);
  if (highestInputLevel > matrix_store.largest_level) 
  {
    fprintf(stderr,"Error: requested zero matrix (levels %dx%d) larger than max_level input (%d) to create store.\n", row_level, col_level, matrix_store.largest_level);
    exit(1);
  }
  return matrix_store.zero[row_level][col_level];
}  

int64_t get_identity_pID(mat_level_t level)
{
  if (level==0) return packedID_scalar1;
  if (level > matrix_store.largest_level) 
    {
      fprintf(stderr,"Error: requested identity matrix has size %d\n", level);
      fprintf(stderr,"       which is larger than maximum matrix size %d.\n", matrix_store.largest_level);
      exit(1);
    }
  return matrix_store.identity[level];
}

int64_t get_iHadamard_pID(mat_level_t level) 
{
  if (level == 0)
  {
    fprintf(stderr,"Error in get_iHadamard_matrix_ptr: there is no sensible definition ");
    fprintf(stderr,"for a level 0 Hadamard matrix\n");
    exit(1);
  }

  if (level > matrix_store.largest_level) 
  {
    fprintf(stderr,"Error: requested integer Hadamard matrix larger than maximum size.\n");
    exit(1);
  }

  return matrix_store.iHadamard[level];
}

#if ALERT_ON_SNAP
extern uint64_t info_snap_occurred;
uint64_t info_snap_occurred = 0;

extern int print_on_snap;
int print_on_snap = 0;

void handle_snap(scalarType original_scalar, scalarType new_scalar, const char *tile_layout) {
    ++info_snap_occurred;
    if (print_on_snap) {
        char *from_str, *to_str;

        from_str = sca_get_readable_approx_str(original_scalar);
        to_str = sca_get_readable_approx_str(new_scalar);
        printf("===>>> %s SNAP  FROM: %s  TO: %s\n", tile_layout, from_str, to_str);
        free(to_str);
        free(from_str);

        from_str = sca_get_exact_str(original_scalar);
        to_str = sca_get_exact_str(new_scalar);
        printf("=EX=>> %s SNAP  FROM: %s  TO: %s\n", tile_layout, from_str, to_str);
        free(to_str);
        free(from_str);

        fflush(stdout);
    }
}
#endif // #if ALERT_ON_SNAP

#ifndef MAR // SPRmode
/* Returns the matrix PTR which matches scalar, starting from the hash */
/*!
 * \ingroup larc
 * \brief Finds the matrix pointer for a given scalar value (does not add it to the matrix store if it is not found)
 * \param hash A value indicating which hash chain would contain the scalar value
 * \param scalar A scalarType value
 * \return The scalar pointer to the stored scalar, or SCALAR_PTR_INVALID if not found
 */
mats_ptr_t find_scalar_SPRmode(uint64_t hash, scalarType scalar)
{
  int debug = 0;

  if (debug) {
    printf("In routine %s\n", __func__);
    printf(">>>calling hash_get_chain\n");
  }

  hash_node_t *n = hash_get_chain(matrix_store.scalar_hash_table, hash);
  uint32_t depth = 0;
#ifdef HASHSTATS
  (matrix_store.scalar_hash_table->num_accesses[hash])++;
#endif

  if (scratchVars.quick_use_in_use)
      fprintf(stderr,"%s reusing scratchVars.quick_use!\n",__func__);
  scratchVars.quick_use_in_use = 1;
  if (scratchVars.misc_in_use)
      fprintf(stderr,"%s reusing scratchVars.misc!\n",__func__);
  scratchVars.misc_in_use = 1;

  scalarType *input_center = &(scratchVars.quick_use);
  scalarType *stored_center = &(scratchVars.misc);

  return_SPR_region_center(input_center,scalar);

  while (n) {

    if (debug) printf(">>>in loop, comparing old and new SPR region centers\n");
    // find the scalar value stored in the scalarRecord for this node
    int64_t stored_packedID = n->packedID;
    mats_ptr_t stored_scaPTR = (mats_ptr_t)get_recordPTR_from_pID(
       stored_packedID, "", __func__, 0);

    // EXCEPTION CHECKING
    if (stored_scaPTR == SCALAR_PTR_INVALID) {
	fprintf(stderr,"Fail in %s, invalid matrixRecord in matrixStore!\n",
	        __func__);
	exit(1);
    }
 
    // confirm that we are in the correct region
    return_SPR_region_center(stored_center,stored_scaPTR->scalar_value);
    if (sca_eq(*stored_center,*input_center)) {
      matrix_store.scalar_hash_table->hits++;    
      if (!n->hits_maxxed) {
          n->record_hits++;
          if (n->record_hits == ((uint32_t) (-1))) {
              n->hits_maxxed = 1;
          }
      }

      if (depth > matrix_store.deepest) { matrix_store.deepest = depth; }

    // This code was a start at working on optimizing the hash chain order
    // It is probably obsolete at this point, but we are leaving it in place
    // in case it is a good start point for a future optimization
//    if (n != matrix_store.hash[hash]) {
//      struct node *head = matrix_store.hash[hash];
//      if (n->prev) n->prev->next = n->next;
//      if (n->next) n->next->prev = n->prev;
//      n->prev = NULL;
//      n->next = head;
//      if (head)
//        head->prev = n;
//      matrix_store.hash[hash] = n;
//    }

      if (debug) printf("Leaving %s with valid matrix pointer\n", __func__);

      mats_ptr_t final_answer = stored_scaPTR;
#if ALERT_ON_SNAP
      if (!sca_eq(scalar, final_answer->scalar_value)) {
          handle_snap(scalar, final_answer->scalar_value, "SPR");
      }
#endif // #if ALERT_ON_SNAP
      scratchVars.quick_use_in_use = 0;
      scratchVars.misc_in_use = 0;
      return final_answer;
    }  // if (equal region centers)
    n = n->next;
    depth++;
  } // end while
  scratchVars.quick_use_in_use = 0;
  scratchVars.misc_in_use = 0;

  matrix_store.scalar_hash_table->misses++;

  if (debug) printf("Leaving %s with invalid matrix pointer\n", __func__);

  return SCALAR_PTR_INVALID;
}

/* Inserts a new scalar matrix into the store and returns the new matrix PTR*/
// NOTE: a new copy of scalar is made for storage matrices (in scalar_value)
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
mats_ptr_t insert_scalar_SPRmode(uint64_t hash, scalarType scalar)
{
  
  int verbose = 0;
  if (verbose) {  
    printf("Inside %s: levels (0,0), hash %" PRIu64 "\n", __func__, hash);
  }

  // Create a new primal LARC matrix for this scalar
  
  // allocate space for the new LARC matrix
  mats_ptr_t s_ptr = calloc(1, sizeof(larc_smatrix_t));
  if (s_ptr == NULL) {
    ALLOCFAIL();
    return SCALAR_PTR_INVALID;
  }

  if (matrixID_tracker.matrixID_next+1==0)
  {
    fprintf(stderr,"ERROR: 64-bit matrixID_next about to overflow!\n");
    exit(0);
  }

  // set the basic parameters for the primal matrix for this scalar
  int64_t matrixID = matrixID_tracker.matrixID_next++;
  // since this is a scalar matrix, the matrixID will have the scalar bit set
  s_ptr->packedID = PID_FROM_SCALAR_MID(matrixID);
  if (VERBOSE>=DEBUG) {
    fprintf(stderr,"adding scalar with matrixID %ld", matrixID);
    fprintf(stderr," and packedID %ld", s_ptr->packedID);
    fprintf(stderr," to chain with hash %lu\n",hash);
  }

  // Initialize the tile_node list to be NULLs
  for (int i=0;i<4;++i) s_ptr->tile_node[i] = (hash_node_t*)NULL;

  if (verbose) { 
    printf("  allocated space, set row level, col level, and matrixID\n");
  }
  
  hash_node_t *scalar_node;
#ifndef HASH_CHAIN_GROWS_AT_HEAD
  // INSERT new items at tail of hash collision chain
  scalar_node = hash_insert_at_tail(matrix_store.scalar_hash_table,
           RECORD_PTR_INVALID, s_ptr->packedID, hash);
#endif
#ifdef HASH_CHAIN_GROWS_AT_HEAD
  // INSERT new items at head of hash collision chain
  scalar_node = hash_insert_at_head(matrix_store.scalar_hash_table,
           RECORD_PTR_INVALID, s_ptr->packedID, hash);
#endif

  // copy scalar_node into tile_node[0] - this may facilitate cleaning
  s_ptr->tile_node[0] = scalar_node;

  if (verbose) {
    printf("  appended hash chain with 0, 0levels\n");
  }
  
  // update matrix store parameters
  matrix_store.hist[0][0]++;
  s_ptr->info = 0;
  s_ptr->hold = 0;
  s_ptr->lock = 0;
  s_ptr->appears_as_sub_count = 0;
  
  if (verbose) { 
    printf("  set row level, col level, and matrixID\n");
  }
  
  if (verbose) {
    char *scalar_string = sca_get_readable_approx_str(scalar);
    printf("ADDING TO MATRIX SCALAR STORE: %ld -> %s\n", s_ptr->packedID, scalar_string);
    free(scalar_string);
  }

  // Inside of the larc_smatrix_t there is a union which holds
  // either the trace_element (for nonscalars) or scalar_value (for scalars)
  // scalar_value is a separate copy of the input scalar,
  // not a reference to it. 
  sca_init(&(s_ptr->scalar_value));
  sca_set(&(s_ptr->scalar_value), scalar);
    
  // test to see if scalars are zero or one
  if (scratchVars.quick_use_in_use)
    fprintf(stderr,"%s reusing scratchVars.quick_use!\n",__func__);
  scratchVars.quick_use_in_use = 1;

  scalarType *sca_id = &scratchVars.quick_use;
  sca_set_str(sca_id, "0");
  s_ptr->iszero = (0 != sca_eq(s_ptr->scalar_value, *sca_id));
  sca_set_str(sca_id, "1");
  s_ptr->isid = (0 != sca_eq(s_ptr->scalar_value, *sca_id));
  scratchVars.quick_use_in_use = 0;

  // check for overflow
  if (matrix_store.num_scalars+1 < matrix_store.num_scalars)
  {
    fprintf(stderr,"ERROR: overflow in matrix_store.num_scalars\n");
  }
  
  matrix_store.num_scalars++;
  if( matrix_store.num_scalars > matrix_store.max_num_scalars ) {
    matrix_store.max_num_scalars = matrix_store.num_scalars;
  }

  // Record the matrix pointer in table indexed by matrixID
  add_recordPTR_to_indexTable((record_ptr_t)s_ptr, s_ptr->packedID);
 
  return s_ptr;
}
#endif // SPRmode

int64_t get_pID_from_four_sub_pIDs(int64_t A_pID, int64_t B_pID,
    int64_t C_pID, int64_t D_pID, mat_level_t row_level, mat_level_t col_level) 
{
  int verbose = 0;
  // we do not check validity, since it is OK for some of these packedIDs
  // to be invalid

  // check consistency of row and column levels
  if (row_level==0 && col_level==0)
  {
    fprintf(stderr,"Inside %s:\n", __func__);
    fprintf(stderr," (0,0) level matrix has no submatrices\n");
    exit(1);
  }

  if (verbose) { 
    printf("Inside %s:\n", __func__);
    printf("  (%d,%d) level matrix \n",row_level,col_level);
  }

  // consistency checks
  if (A_pID == MATRIX_ID_INVALID)
  {
    fprintf(stderr,"in %s, submatrix in position 0 is invalid\n",__func__);
    exit(1);
  }
  if ((B_pID == MATRIX_ID_INVALID)&&(col_level>0))
  {
    fprintf(stderr,"in %s, submatrix in position 1 is invalid\n",__func__);
    fprintf(stderr,"and column level > 0\n");
    exit(1);
  }
  if ((col_level==0)&&(B_pID != MATRIX_ID_INVALID))
  {
    fprintf(stderr,"in %s, submatrix in position 1 is valid\n",__func__);
    fprintf(stderr,"and column level = 0\n");
    exit(1);
  }
  if ((C_pID == MATRIX_ID_INVALID)&&(row_level>0))
  {
    fprintf(stderr,"in %s, submatrix in position 2 is invalid\n",__func__);
    fprintf(stderr,"and row level > 0\n");
    exit(1);
  }
  if ((row_level==0)&&(C_pID != MATRIX_ID_INVALID))
  {
    fprintf(stderr,"in %s, submatrix in position 2 is valid\n",__func__);
    fprintf(stderr,"and row level = 0\n");
    exit(1);
  }
  if ((D_pID == MATRIX_ID_INVALID)&&(row_level>0)&&(col_level>0))
  {
    fprintf(stderr,"in %s, submatrix in position 3 is invalid\n",__func__);
    fprintf(stderr,"and both row and column level > 0\n");
    exit(1);
  }
  if (((row_level==0)||(col_level==0))&&(D_pID != MATRIX_ID_INVALID))
  {
    fprintf(stderr,"in %s, submatrix in position 3 is valid\n",__func__);
    fprintf(stderr,"and either row or column level is 0\n");
    exit(1);
  }

  if (verbose)
  {
    fprintf(stderr,"The four quadrant matrixIDs are\n");
    fprintf(stderr,"A_mID = %ld, B_mID = %ld, C_mID = %ld, D_mID = %ld\n",
                MID_FROM_PID(A_pID), MID_FROM_PID(B_pID),
                MID_FROM_PID(C_pID), MID_FROM_PID(D_pID));
  }

  matrix_type_t mat_type = matrix_type_from_levels(row_level,col_level);
  record_ptr_t  panel[4];
  panel[0] = get_recordPTR_from_pID(A_pID, "first", __func__,0);
  panel[1] = get_recordPTR_from_pID(B_pID, "second", __func__,0);
  panel[2] = get_recordPTR_from_pID(C_pID, "third", __func__,0);
  panel[3] = get_recordPTR_from_pID(D_pID, "fourth", __func__,0);

  // case 1: confirm all submatrices SCALAR
  if ((row_level<=1)&&(col_level<=1))
  {
    int errorFlag = 0;
    errorFlag |= !(IS_SCALAR(A_pID));
    errorFlag |= ((B_pID != MATRIX_ID_INVALID) && (!IS_SCALAR(B_pID)));
    errorFlag |= ((C_pID != MATRIX_ID_INVALID) && (!IS_SCALAR(C_pID)));
    errorFlag |= ((D_pID != MATRIX_ID_INVALID) && (!IS_SCALAR(D_pID)));
    if (errorFlag)
    {
      fprintf(stderr,"in %s, all submatrices are supposed\n",__func__);
      fprintf(stderr,"to be scalar, but at least one is not\n");
      fprintf(stderr,"A_mID = %ld, B_mID = %ld, C_mID = %ld, D_mID = %ld\n",
                MID_FROM_PID(A_pID), MID_FROM_PID(B_pID),
                MID_FROM_PID(C_pID), MID_FROM_PID(D_pID));
      exit(1);
    }
  } // end case scalar submatrices
  else // case 2: confirm all submatrices not SCALAR
  {
    int errorFlag = 0;
    errorFlag |= IS_SCALAR(A_pID);
    errorFlag |= ((B_pID != MATRIX_ID_INVALID) && IS_SCALAR(B_pID));
    errorFlag |= ((C_pID != MATRIX_ID_INVALID) && IS_SCALAR(C_pID));
    errorFlag |= ((D_pID != MATRIX_ID_INVALID) && IS_SCALAR(D_pID));
    if (errorFlag)
    {
      fprintf(stderr,"in %s, all submatrices are supposed\n",__func__);
      fprintf(stderr,"to be non-scalar, but at least one is not\n");
      fprintf(stderr,"A_mID = %ld, B_mID = %ld, C_mID = %ld, D_mID = %ld\n",
                MID_FROM_PID(A_pID), MID_FROM_PID(B_pID),
                MID_FROM_PID(C_pID), MID_FROM_PID(D_pID));
      exit(1);
    }
    // confirm dimensions of nonscalar submatrices
    for (int i=0; i<4; ++i)
    {
      matns_ptr_t sub_ptr = (matns_ptr_t)(panel[i]);
      if (sub_ptr == MATRIX_PTR_INVALID) continue;
      if ((row_level>1) && (sub_ptr->row_level + 1 != row_level)) {
          fprintf(stderr,"in %s, ",__func__);
          fprintf(stderr,"WARNING on panel %d: row levels inconsistent\n",i);
          fprintf(stderr,"Expected: %d, Found: %d\n", row_level-1,
    		  sub_ptr->row_level);
      }
      if ((col_level>1) && (sub_ptr->col_level + 1 != col_level)) {
          fprintf(stderr,"in %s, ",__func__);
          fprintf(stderr,"WARNING on panel %d: col levels inconsistent\n",i);
          fprintf(stderr,"Expected: %d, Found: %d\n", col_level-1,
    		  sub_ptr->col_level);
      }
    } // have confirmed submatrix types agree
  } // end case nonscalar submatrices

  if (verbose) { 
    printf("  About to calculate the matrix_type, hash, and see if matrix is stored\n");
  }
  int64_t subMatList[4] = {A_pID, B_pID, C_pID, D_pID};
  uint64_t hash = hash_from_matrix_subMatList(subMatList, mat_type,
                 matrix_store.nonscalar_hash_table->exponent);
  int64_t m_pID =  matrix_find_subMatList(hash, subMatList);

  // if unable to find, attempt to insert
  if (verbose) { 
    printf("  About to check if the matrix was not found\n");
  }

  if (m_pID == MATRIX_ID_INVALID)
  {
    if (verbose) { 
        printf("Inside %s, matrix index not found by matrix_find\n",__func__);
        printf("   CALLING  matrix_insert\n");
    }
    m_pID = matrix_insert_panel(hash, panel, row_level, col_level);
    if (verbose) {   
      printf("    RETURNED from matrix_insert with packedID %ld\n",m_pID);
    }
  }

  if (m_pID == MATRIX_ID_INVALID) {
    fprintf(stderr,"ERROR: inserting matrix from panels into store failed.\n");
    exit(EXIT_FAILURE);
  }

  return m_pID;
}

// python interface for scalar type
int64_t get_valID_from_valString(char *val)
{
  scalarType sca;
  sca_init(&sca);
  sca_set_str(&sca, val);
  mats_ptr_t s_ptr = get_scalarPTR_for_scalarVal(sca);
  sca_clear(&sca);
  return s_ptr->packedID;
}

/* The first argument of this routine is a scalarType (complex, double, int) */
mats_ptr_t get_scalarPTR_for_scalarVal(scalarType scalar)
{

  int verbose = 0;

  if (verbose) { 
    printf("Inside %s:\n", __func__);
    printf("  (0,0) level matrix");
    printf("  with value %s\n",sca_get_readable_approx_str(scalar));
  }

  if (verbose) { 
    printf("  About to calculate the matrix_type, hash, and see if matrix is stored\n");
  }

#ifndef MAR // SPRmode
mats_ptr_t s_ptr = get_PTR_scalar_record(scalar);
  if (verbose) {   
    printf("    RETURNED from get_PTR with scalar ptr %p\n", s_ptr);
  }
#endif // SPRmode
#ifdef MAR
  mats_ptr_t s_ptr = retrieve_PTR_scalar_record(scalar);
  if (verbose) {   
    printf("    RETURNED from retrieve_PTR with scalar ptr %p\n", s_ptr);
  }
#endif // MARmode

  
  // if couldn't insert, then bail out
  if (s_ptr == SCALAR_PTR_INVALID) {
    abort();
  }
  if (verbose) {
    printf("Leaving %s\n", __func__);
  }
  
  return s_ptr;
}

size_t scalar_store_count(void)
{
	return matrix_store.num_scalars;
}


size_t nonscalar_store_count(void)
{
	return matrix_store.num_nonscalars;
}

int64_t adjoint(int64_t m_pID)
{
  record_ptr_t r_ptr;
  check_validity_one_input(m_pID,__func__,&r_ptr);

  if (IS_SCALAR(m_pID))
  {
     // get the matrix pointers from the matrixIDs, and see if still in store
     uint64_t adj_hash = hash_from_op(m_pID,m_pID,ADJOINT);
     // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE
     int64_t adj_pID = op_get(ADJOINT,m_pID,m_pID,adj_hash);
     if (adj_pID == MATRIX_ID_INVALID) {
        mats_ptr_t s_ptr = (mats_ptr_t)r_ptr;
        if (sca_is_pure_real(s_ptr->scalar_value)) {
           // adjoint of real number is that number, guaranteed to be in store
           adj_pID = m_pID;
        } else {
           scalarType scalar;
           sca_init(&scalar); 
           // adjoint of complex number might not be in store
           sca_conj(&scalar, s_ptr->scalar_value);
           adj_pID = (get_scalarPTR_for_scalarVal(scalar))->packedID;
           sca_clear(&scalar);
        }
        op_set(ADJOINT, m_pID, m_pID, adj_pID, adj_hash);
     } // end if (not found in ops store)
     return adj_pID;
  } // end SCALAR

  // NON-SCALAR CASES

  // get the matrix pointers from the matrixIDs, and see if still in store
  matns_ptr_t m_ptr = (matns_ptr_t)r_ptr;
  uint64_t adj_hash = hash_from_op(m_pID,m_pID,ADJOINT);
  // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE
  int64_t adj_pID = op_get(ADJOINT,m_pID,m_pID,adj_hash);
  if (adj_pID == MATRIX_ID_INVALID) {
     matrix_type_t mat_type = matrix_type(m_ptr);
     matrix_type_t adjoint_type;
     if (mat_type == ROW_VECTOR) {adjoint_type = COL_VECTOR;}
     else if (mat_type == COL_VECTOR) {adjoint_type = ROW_VECTOR;}
     else  {adjoint_type = mat_type;}

     int64_t subMatList[4];
     subMatList[0] = adjoint(m_ptr->subMatList[0]);
     if (m_ptr->subMatList[2]==MATRIX_ID_INVALID)
        subMatList[1] = MATRIX_ID_INVALID;
     else subMatList[1] = adjoint(m_ptr->subMatList[2]);
     if (m_ptr->subMatList[1]==MATRIX_ID_INVALID)
        subMatList[2] = MATRIX_ID_INVALID;
     else subMatList[2] = adjoint(m_ptr->subMatList[1]);
     if (m_ptr->subMatList[3]==MATRIX_ID_INVALID)
        subMatList[3] = MATRIX_ID_INVALID;
     else subMatList[3] = adjoint(m_ptr->subMatList[3]);

     uint64_t hash = hash_from_matrix_subMatList(subMatList, adjoint_type,
            matrix_store.nonscalar_hash_table->exponent);
     adj_pID =  matrix_find_subMatList(hash, subMatList);
     if (adj_pID == MATRIX_ID_INVALID)
     {
         record_ptr_t panel[4];
         for (int i=0;i<4;++i)
            panel[i] = get_recordPTR_from_pID(subMatList[i],"",__func__,0);
         adj_pID = matrix_insert_panel(hash, panel,
                  m_ptr->col_level, m_ptr->row_level);
     }
     op_set(ADJOINT, m_pID, m_pID, adj_pID, adj_hash);
  } // end if (not found in ops store)   
  return adj_pID;
}

/* The total number of matrices created and stored in  */
/* the matrix store (the current number of stored     */
/* may be less if there was matrix store cleaning).     */
uint64_t num_matrices_created(void)
{
  uint64_t ret = matrixID_tracker.matrixID_next;
  return(ret);
}

record_ptr_t get_recordPTR_from_pID(int64_t packedID, const char* arg_no,
        const char* calling_routine, int print_flag)
{

  // confirm that passed packedID is not invalid
  if (packedID == MATRIX_ID_INVALID) return RECORD_PTR_INVALID;

  int64_t matrixID = MID_FROM_PID(packedID);
  uint64_t block_index, index;
  get_matrixID_table_indices(packedID, &block_index, &index);

  // grab highest matrixID
  uint64_t highest_matrixID = num_matrices_created()-1;

  // check bounds
  if (matrixID > highest_matrixID) {
    if (print_flag)
      fprintf(stderr,"ERROR: In %s, the %s argument, matrixID %" PRIu64 ", is out of range.\n",
          calling_routine, arg_no, matrixID);
    return(RECORD_PTR_INVALID);
  }

  // get matrix address
  record_ptr_t val_ptr =
    matrixID_tracker.ptrs_indexed_by_matrixID[block_index][index];

  // warn the user if the packedID doesn't refer to a valid matrix
  if ((val_ptr == RECORD_PTR_INVALID) && (print_flag==1)) {
    fprintf(stderr,"WARNING: In %s, problem with the %s argument:",
        calling_routine, arg_no);
    fprintf(stderr," matrix with matrixID %" PRIu64 " has been", matrixID);
    fprintf(stderr," removed from the Matrix Store.\n");
  }

  return val_ptr;
}

/*!
 * \ingroup larc
 * \brief Writes statistic information about a matrix in the matrix_store to a file
 * \param mat_pID The packedID of the matrix 
 * \param f A file pointer
 * \return 1 
 */
static int  
matrix_info_to_file(uint64_t mat_pID, FILE *f) {
  // EXCEPTIONS 
  if (f==NULL) { 
    fprintf(stderr,"WARN: %s called with bad file pointer\n", __func__); 
    return(0); 
  } 

  record_ptr_t r_ptr = get_recordPTR_from_pID(mat_pID,"",__func__,0);

  // HEADER: "matrixID   Levels    Value    Appears_As_Sub Count   Lock   Hold
  //   matrixID 
  fprintf(f,"      %" PRId64, MID_FROM_PID(mat_pID));
  //   Levels
  if (IS_SCALAR(mat_pID))
  {
    mats_ptr_t s_ptr = (mats_ptr_t)r_ptr;
    fprintf(f,"    \t( 0, 0)");
    //   Value
    char *trace_string = sca_get_readable_approx_str(s_ptr->scalar_value);
    fprintf(f," \t%s                   ", trace_string);
    free(trace_string);
    //   Appears_As_Sub Count = Reference Count 
    fprintf(f," \t%8u     %4d    %4d     ", matrix_appears_as_sub_count(s_ptr),
        s_ptr->lock, s_ptr->hold);
  } // end scalar branch
  else
  {
    matns_ptr_t m_ptr = (matns_ptr_t)r_ptr;
    mat_level_t row_level = m_ptr->row_level;
    mat_level_t col_level = m_ptr->col_level;
    fprintf(f,"    \t(%2d,%2d)", row_level, col_level);
    //   Value
    fprintf(f,"   [ ");
    for (int i = 0; i < 4; i++)
    {
      int64_t subMat_pID = m_ptr->subMatList[i];
      if (subMat_pID == MATRIX_ID_INVALID) fprintf(f,"  X ");
      else fprintf(f,"%" PRId64 " ", MID_FROM_PID(subMat_pID));
    } 
    fprintf(f,"]        ");
    //   Appears_As_Sub Count = Reference Count 
    fprintf(f," \t%8u     %4d    %4d     ", matrix_appears_as_sub_count(m_ptr),
        m_ptr->lock, m_ptr->hold);
  } // end nonscalar branch

  //   End of this line 
  fprintf(f,"\n");
  return(1);
}

/***************************************************************************
*                      fprint_store_info_for_matrixID_range                           *
*  This function prints information about matrices from matrixID      *
*  "start" to "end"  (0 <= start <= end < SIZE_MPT^2
*  The output file path and a user comment to be printed in the file are   *
*  also arguments.                                                         *
***************************************************************************/
int
fprint_store_info_for_matrixID_range(uint64_t start, uint64_t end,
                          char *outfilepath, char *comment)
{
  if (start > end) {
    fprintf(stderr,"ERROR in %s: impossible range %" PRIu64 " > %" PRIu64 "\n",
             __func__, start, end);
    return(0);
  }

  printf("Printing matrices with matrixIDs from %" PRIu64 " to %" PRIu64 " to file %s\n",
          start,end,outfilepath);
  
  int ret = 1;
  FILE *f = fopen(outfilepath, "w"); 
  fprintf(f,"Comment: %s\n",comment); 
  fprintf(f,"Total number of matrices that have been created is:           %" PRIu64 "\n", 
  	  num_matrices_created()); 

  int64_t max_num_matrices = ((int64_t)1)<<(2*LOG_SIZE_MPT);

  /* // fprintf(f,"Largest number of matrices ever stored simultaneously:        %zd\n",  */
  /* //	  store.max_num_nonscalars); */
  /* // fprintf(f,"Largest number of scalars ever stored simultaneously:         %zd\n",  */
  /* //	  store.max_num_scalars); */
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
  for (i = start; i <= end; ++i)
  { 
       // matrix_info_to_file expects a packedID, so we must recover one
       // determine whether pointer is to scalar or matrix structure
       int64_t scalar_pID = PID_FROM_SCALAR_MID(i);
       int64_t nonscalar_pID = PID_FROM_NONSCALAR_MID(i);
       // it doesn't matter which putative packedID is passed here
       // as the function only uses the matrixID part of it
       record_ptr_t r_ptr = get_recordPTR_from_pID(nonscalar_pID,"",__func__,0);
       if (r_ptr==RECORD_PTR_INVALID) continue;

       // Both scalar structure and matrix structure have packedID fields,
       // so the compiler will not complain and 'some' int64_t value will be
       // readable even when accessing the incorrect structure. (To minimize
       // problems, we made packedID the first thing in each structure.) We
       // check for the case where we get no correct answer or both answers
       // appear to be correct.
       int64_t mat_pID;
       int flag = 0;
       if (((matns_ptr_t)r_ptr)->packedID == nonscalar_pID)
          mat_pID = nonscalar_pID, flag++;
       if (((mats_ptr_t)r_ptr)->packedID == scalar_pID)
          mat_pID = scalar_pID, flag++;
       if (flag != 1)
       {
          fprintf(stderr,"ERROR in %s: cannot recover packedID\n",__func__);
          fprintf(stderr,"from given matrixID %" PRId64 "\n",i);
          exit(0);
       }
       // fprintf(f,"   --- call to matrix_info_to_file(%d,f) ---\n",i);
       ret = matrix_info_to_file(mat_pID,f);
     }
  fprintf(f,"\n======================================================================\n");

  fclose(f); 

  return ret;
}

/***************************************************************************
*                   fprint_nonscalar_hash_chain_info                         *
*  This function prints information about matrices to a file,			   *
*  given the hash value for that chain.                                    *
*  To get the appropriate hash value for the argument of this function     *
*     hashID = hash_pID(matrixID);				   *
*  The output file path and a user comment to be printed in the file are   *
*  also arguments.                                                         *
***************************************************************************/
int
fprint_nonscalar_hash_chain_info(uint64_t hash, char *outfilepath, char *comment)
{
  FILE *f; 

  if (strcmp(outfilepath,"stdout")) {
     printf("Printing matrix hash chain with hash value %" PRIu64 " to file %s\n",
          hash,outfilepath);
     f = fopen(outfilepath, "w"); 

  }
  else {
    printf("Printing matrix hash chain with hash value %" PRIu64 " to screen\n",hash);
    f = stdout;
  }
  
  int ret = 1;
  fprintf(f,"Comment: %s\n",comment); 
  fprintf(f,"This is the matrix hash chain info for hash value %" PRIu64 "\n", hash); 
  fprintf(f,"\n\n===================(Matrix Hash Chain)=============================\n");
  fprintf(f,"MatrixID        Levels           Value   ");
  fprintf(f,"        Appears_As_SubCount\n"); 

  hash_table_t *table_ptr = matrix_store.nonscalar_hash_table;
  hash_node_t *node_ptr = table_ptr->heads[hash];

  while (node_ptr)  {
    uint64_t packedID = node_ptr->packedID;
    ret = matrix_info_to_file(packedID,f);
    node_ptr = node_ptr->next;
  }
  fprintf(f,"\n======================================================================\n");

  if  (strcmp(outfilepath,"stdout")) fclose(f);   // if not stdout then close

  return ret;
} 

int
fprint_scalar_hash_chain_info(uint64_t hash, char *outfilepath, char *comment)
{
  FILE *f; 

  if (strcmp(outfilepath,"stdout")) {
     printf("Printing scalar hash chain with hash value %" PRIu64 " to file %s\n",
          hash,outfilepath);
     f = fopen(outfilepath, "w"); 

  }
  else {
    printf("Printing scalar hash chain with hash value %" PRIu64 " to screen\n",hash);
    f = stdout;
  }
  
  int ret = 1;
  fprintf(f,"Comment: %s\n",comment); 
  fprintf(f,"This is the scalar hash chain info for hash value %" PRIu64 "\n", hash); 
  fprintf(f,"\n\n===================(Scalar Hash Chain)=============================\n");
  fprintf(f,"MatrixID        Levels           Value   ");
  fprintf(f,"        Appears_As_SubCount\n"); 

  hash_table_t *table_ptr = matrix_store.scalar_hash_table;
  hash_node_t *node_ptr = table_ptr->heads[hash];

  while (node_ptr)  {
    uint64_t packedID = node_ptr->packedID;
    ret = matrix_info_to_file(packedID,f);
    node_ptr = node_ptr->next;
  }
  fprintf(f,"\n======================================================================\n");

  if  (strcmp(outfilepath,"stdout")) fclose(f);   // if not stdout then close

  return ret;
}


/***************************************************************************
*                   print_nonscalar_hash_chain_info                       *
*  This function prints information about matrices given the hash value    *
*  for that chain.                                                         *
*  To get the appropriate hash value for the argument of this function     *
*     hashID = hash_pID(matrixID);				   *
*  A user comment to be printed is also an argument   					   *
***************************************************************************/
int
print_nonscalar_hash_chain_info(uint64_t hash, char *comment)
{
  char *filename = "stdout";
  return fprint_nonscalar_hash_chain_info(hash, filename, comment);
} 

int
print_scalar_hash_chain_info(uint64_t hash, char *comment)
{
  char *filename = "stdout";
  return fprint_scalar_hash_chain_info(hash, filename, comment);
}

#ifdef MAR
static MAR_tile_index_t *get_tile_index_from_scalarPTR(mats_ptr_t s_ptr)
{
    // create local static variable to hold tile indices
    // on first call to routine
    static MAR_tile_index_t * temp_remove_tile_index_PTR;
    static int temp_is_initialized = 0;
    if (!temp_is_initialized) {
          temp_remove_tile_index_PTR = (MAR_tile_index_t *)malloc(sizeof(MAR_tile_index_t));

    if (temp_remove_tile_index_PTR == NULL) { ALLOCFAIL(); }
#ifdef IS_COMPLEX
          mpz_init(temp_remove_tile_index_PTR->real_index);
          mpz_init(temp_remove_tile_index_PTR->imag_index);
#else
          mpz_init(temp_remove_tile_index_PTR->index);
#endif // #ifdef IS_COMPLEX
          temp_is_initialized = 1;
    }

    int regionbitparam = get_regionbitparam();
    get_tile_index(temp_remove_tile_index_PTR, regionbitparam,
          s_ptr->scalar_value);
    return temp_remove_tile_index_PTR;
}
#endif // MAR

/******************************************************************
*      remove_matrix_from_nonscalar_store_by_recordPTR_andor_nodePTR
*******************************************************************/
/*!
 \ingroup larc
 * \brief A routine that tries to remove a matrix from the nonscalar store
 * \param r_ptr A pointer to the matrix to be removed
 * \param n The hash table node containing the matrix (may be set to NULL)
 * \param supplied_hash The index for the hash chain containing the matrix
 * \param m_pID The packedID for the matrix to be removed
 *
 * If the input n is not NULL, then n will be used; if it is NULL, the
 * function will search the hash table for the matrix.
 * The flag removeNbhrRecord should always be zero when this routine is called.
 * When in MAR mode, the routine will recursively call itself with the flag
 * equal to one to clean claimed neighbor regions when a primal record is
 * cleaned.
 *
 * \return NULL if n null, otherwise the pointer to the next hash node in the chain
 */
static
hash_node_t *remove_matrix_from_nonscalar_store_by_recordPTR_andor_nodePTR (
  record_ptr_t r_ptr, hash_node_t *n, uint64_t supplied_hash, int64_t m_pID)
{
  int verbose = 0;

  if (verbose) fprintf(stderr,"in %s, with input packedID = %ld (matrixID %ld)\n",
      __func__, m_pID, MID_FROM_PID(m_pID));

  if (IS_SCALAR(m_pID))
  {
    fprintf(stderr,"in %s: input packedID %ld is for a scalar!\n",
         __func__, m_pID);
    exit(0);
  }
  matns_ptr_t m_ptr;

#define CHECK_CONSISTENCY
#ifdef CHECK_CONSISTENCY
  // do a check that m_pID and r_ptr are consistent
  m_ptr = (matns_ptr_t)r_ptr;
  if (m_ptr->packedID != m_pID)
    fprintf(stderr,"in %s: m_pID is %ld, but value stored in record is %ld\n",
                    __func__,m_pID, m_ptr->packedID);
  m_ptr = (matns_ptr_t)get_recordPTR_from_pID(m_pID,"",__func__,0);
  if (m_ptr != (matns_ptr_t)r_ptr)
    fprintf(stderr,"in %s: input ptr is %p, but value obtained using pID is %p\n",
                    __func__,m_ptr, r_ptr);
#endif
#undef CHECK_CONSISTENCY

  // first, confirm that matrix should be removed
  // Can not remove the matrix in any of the following situation:
  //      matrix appears_as_subs for other matrices
  //      matrix has a lock or hold on it
  //      matrix is a neighbor Record in MAR associated with a primal record
  //           unless the removeNhbrRecord flag is set

  m_ptr = (matns_ptr_t)r_ptr;
  if (verbose) fprintf(stderr,"m_ptr->appears_as_sub_count = %d\n",
        m_ptr->appears_as_sub_count);
  if (verbose) fprintf(stderr,"m_ptr->lock = %d\n", m_ptr->lock);
  if (verbose) fprintf(stderr,"m_ptr->hold = %d\n", m_ptr->hold);
  if (m_ptr->appears_as_sub_count != 0)
    { if (n) return n->next; else return NULL; }
  if (m_ptr->lock) { if (n) return n->next; else return NULL; }
  if (m_ptr->hold) { if (n) return n->next; else return NULL; }
  if (verbose) fprintf(stderr,"passed tests - remove\n");

  // If the matrix is not a SCALAR, then recursively remove it
  //  *for each child of matrix, decrement appears_as_sub_count 
  //  *if child is now at zero counts, attempt to remove it
  //printf("DEBUG: record type not scalar in %s.\n", __func__); //debug
  for (int i=0;i<4;++i)
  {
    int64_t child_pID = m_ptr->subMatList[i];
    if (child_pID==MATRIX_ID_INVALID) continue;
    if (verbose) fprintf(stderr,"i=%d: child_pID = %ld\n",i,child_pID);
    record_ptr_t child_ptr =
        get_recordPTR_from_pID(child_pID, "", __func__, 0);
    if (child_ptr == RECORD_PTR_INVALID)
    {
      fprintf(stderr,"in %s, for packedID %ld:\n",__func__,m_pID);
      fprintf(stderr,"\tfound child with packedID %ld", child_pID);
      fprintf(stderr," that has invalid pointer\n");
      continue;
    }
    if (verbose) fprintf(stderr,"record not invalid\n");
    // decrement sub_count for child, then remove if allowed
    if (IS_SCALAR(child_pID)) {
         if (verbose) fprintf(stderr,"decrementing scalar subcount\n");
         scalar_appears_as_sub_count_decrement_ptr((mats_ptr_t)child_ptr);
         // scalars are removed by a different subroutine
    }
    else {
         if (verbose) fprintf(stderr,"decrementing matrix subcount\n");
         matrix_appears_as_sub_count_decrement_ptr((matns_ptr_t)child_ptr);
         uint64_t child_hash = hash_pID(child_pID);
         if (verbose) fprintf(stderr,"%s called from \n%s with pID %ld\n",
             __func__,__func__,child_pID);
         remove_matrix_from_nonscalar_store_by_recordPTR_andor_nodePTR(
             child_ptr, NULL, child_hash, child_pID);
    }
  } // loop over submatrices

  // First remove the hash node (and sew up the hash chain), but don't free
  // the memory for the matrix yet.
  // This is easy when we already know the hash node; otherwise we have to
  // use the record pointer and search the hash chain for the proper record.
  hash_node_t *next_n = NULL;
  if (n) {
    if (verbose) fprintf(stderr,"calling hash_node_remove_node\n");
    next_n = 
     hash_node_remove_node(matrix_store.nonscalar_hash_table, n, supplied_hash);
  }
  else {
    if (verbose) fprintf(stderr,"calling hash_node_remove\n");
    hash_node_remove(matrix_store.nonscalar_hash_table, RECORD_PTR_INVALID,
                       m_pID, supplied_hash);
  }

  // Remove the matrix pointer from the table indexed by matrixIDs
  if (verbose)
     fprintf(stderr,"calling invalidate_recordPTR_in_indexTable(%ld)\n",m_pID);
  invalidate_recordPTR_in_indexTable(m_pID);

  // Clean up required for Nonscalar Matrices

  // Get row level and col level for matrix.
  mat_level_t row_level;
  mat_level_t col_level;
  row_level = m_ptr->row_level;
  col_level = m_ptr->col_level;

  // Decrement the histogram and counts of matrices
  matrix_store.hist[row_level][col_level]--;
  matrix_store.num_nonscalars--;
  
  // special SQUARE MATRIX stuff: clear the matrix trace
  if (row_level == col_level) {
    // Free/clear the scalar
    sca_clear(&(m_ptr->trace_element));
  }

  // free the structure, note that no mallocs other than scalar occur in
  // the matrix
  free(m_ptr);

  return next_n;
}

/******************************************************************
*      remove_matrix_from_scalar_store_by_recordPTR_andor_nodePTR
*******************************************************************/
hash_node_t *remove_matrix_from_scalar_store_by_recordPTR_andor_nodePTR (
  record_ptr_t r_ptr, hash_node_t *n, uint64_t supplied_hash, int64_t m_pID)
{
  int verbose = 0;

  if (!IS_SCALAR(m_pID))
  {
    fprintf(stderr,"in %s: input packedID %ld is for a nonscalar!\n",
         __func__,m_pID);
    exit(0);
  }

  if (verbose)
    fprintf(stderr,"in %s, with input packedID = %ld (matrixID %ld)\n",
        __func__, m_pID, MID_FROM_PID(m_pID));

  mats_ptr_t s_ptr;

#define CHECK_CONSISTENCY
#ifdef CHECK_CONSISTENCY
  // do a check that m_pID and r_ptr are consistent
  s_ptr = (mats_ptr_t)r_ptr;
  if (s_ptr->packedID != m_pID)
    fprintf(stderr,"in %s: m_pID is %ld, but value stored in record is %ld\n",
                    __func__,m_pID, s_ptr->packedID);
  s_ptr = (mats_ptr_t)get_recordPTR_from_pID(m_pID,"",__func__,0);
  if (s_ptr != (mats_ptr_t)r_ptr)
    fprintf(stderr,"in %s: input ptr is %p, but value obtained using pID is %p\n",
                    __func__,s_ptr, r_ptr);
#endif
#undef CHECK_CONSISTENCY

  // first, confirm that matrix should be removed
  // Can not remove the matrix in any of the following situation:
  //      matrix appears_as_subs for other matrices
  //      matrix has a lock or hold on it
  //      matrix is a neighbor Record in MAR associated with a primal record
  //           unless the removeNhbrRecord flag is set
  if (verbose) fprintf(stderr,"scalar pID\n");
  s_ptr = (mats_ptr_t)r_ptr;
  if (verbose) fprintf(stderr,"s_ptr->appears_as_sub_count = %d\n",
        s_ptr->appears_as_sub_count);
  if (verbose) fprintf(stderr,"s_ptr->lock = %d\n", s_ptr->lock);
  if (verbose) fprintf(stderr,"s_ptr->hold = %d\n", s_ptr->hold);
  if (s_ptr->appears_as_sub_count != 0)
    { if (n) return n->next; else return NULL; }
  if (s_ptr->lock) { if (n) return n->next; else return NULL; }
  if (s_ptr->hold) { if (n) return n->next; else return NULL; }
  if (verbose) fprintf(stderr,"passed tests - remove\n");

  // Remove the matrix from the hash table for the scalar store

  // First remove the hash node (and sew up the hash chain), but don't free
  // the memory for the matrix yet.
  // This is easy when we already know the hash node; otherwise we have to
  // use the record pointer and search the hash chain for the proper record.
  hash_node_t *next_n = NULL;
  if (n) {
      if (verbose) fprintf(stderr,"calling hash_node_remove_node\n");
      next_n = 
        hash_node_remove_node(matrix_store.scalar_hash_table, n, supplied_hash);
  }
  else {
      if (verbose) fprintf(stderr,"calling hash_node_remove\n");
      hash_node_remove(matrix_store.scalar_hash_table, RECORD_PTR_INVALID,
                       m_pID, supplied_hash);
  }
  // Remove the matrix pointer from the table indexed by matrixIDs
  if (verbose)
     fprintf(stderr,"calling invalidate_recordPTR_in_indexTable(%ld)\n",m_pID);
  invalidate_recordPTR_in_indexTable(m_pID);

  // Clean up required for ScalarRecord

  // Decrement the histogram and counts of matrices and scalars
  matrix_store.hist[0][0]--;
  matrix_store.num_scalars--;
  sca_clear(&(s_ptr->scalar_value));

  // free the structure, note that no mallocs other than scalar occur in the matrix
  free(s_ptr);
  return next_n;
}

// python interface version
int remove_matrix_from_store(int64_t m_pID)
{
  int verbose = 0;
  // get the matrix pointers from the matrixIDs, and see if still in store
  record_ptr_t r_ptr = get_recordPTR_from_pID(m_pID, "", __func__,0);

  // Users could try to remove the same thing twice, but
  // we chose to not kill the program if they do this, just warn them.
  if (r_ptr == RECORD_PTR_INVALID)
  {
    fprintf(stderr,"in %s, packedID %ld has invalid record pointer\n",
        __func__, m_pID);
    fprintf(stderr,"no action taken\n");
    return 0;
  }
  uint64_t m_hash = hash_pID(m_pID);
  if (IS_SCALAR(m_pID))
  {
    if (verbose)
    {
      fprintf(stderr,
          "remove_matrix_from_scalar_store_by_recordPTR_andor_nodePTR\n");
      fprintf(stderr,"called from %s with pID %ld\n",__func__,m_pID);
    }
    remove_matrix_from_scalar_store_by_recordPTR_andor_nodePTR(
          r_ptr, NULL, m_hash, m_pID);
    return 1;
  }
  if (verbose)
  {
    fprintf(stderr,
          "remove_matrix_from_nonscalar_store_by_recordPTR_andor_nodePTR\n");
    fprintf(stderr,"called from %s with pID %ld\n",__func__,m_pID);
  }
  remove_matrix_from_nonscalar_store_by_recordPTR_andor_nodePTR(
          r_ptr, NULL, m_hash, m_pID);
  return 1;
}


/*****************************************************************
 *                hash_pID                      *
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
int64_t hash_pID(int64_t m_pID)
{
  record_ptr_t r_ptr = get_recordPTR_from_pID(m_pID, "", __func__, 0);
  if (r_ptr == RECORD_PTR_INVALID) { return(-1); }

  uint64_t hash;
  if (IS_SCALAR(m_pID))
  {
    mats_ptr_t s_ptr = (mats_ptr_t)r_ptr;
#ifndef MAR // SPRmode
    hash = region_hash_from_scalarType(s_ptr->scalar_value,
              matrix_store.scalar_hash_table->exponent);
#else // MARmode
#ifdef STORE_TILE_INDEX
      hash = hash_tile_index(s_ptr->tile);
#else
      hash = hash_tile_index(get_tile_index_from_scalarPTR(s_ptr));
#endif // #ifdef STORE_TILE_INDEX
#endif // MARmode
  }
  else
  {
    matns_ptr_t m_ptr = (matns_ptr_t)r_ptr;
    int64_t subMatList[4];
    for (int i=0;i<4;++i) { subMatList[i] = m_ptr->subMatList[i]; }
    // MATRIX is mat_type_t and will be removed from arguments eventually
//    fprintf(stderr,"in %s, call hash_from_matrix_subMatList with\n",__func__);
//    fprintf(stderr,"submatlist %ld, %ld, %ld, %ld\n",subMatList[0],
//        subMatList[1], subMatList[2], subMatList[3]);
    hash = hash_from_matrix_subMatList(subMatList, MATRIX,
           matrix_store.nonscalar_hash_table->exponent);
  }

  return hash;
}

/************************************************************************
 * Clean nonscalar matrices
 *                                                                      *
 ***********************************************************************/
/*!
 * \ingroup larc
 * \brief removes all eligible matrices from the nonscalar (matrix) store.
 * \result Returns 1.
 */
static int clean_nonscalar_matrices()
{
    int verbose = 0;
    uint64_t max = (uint64_t)1 << (matrix_store.nonscalar_hash_table->exponent);
    uint64_t hash_value;

    // loop though hash chains for each hashID without calling
    // clean_op_hash_chain so we don't have to repeat checks for invalid hashID
    hash_table_t *table_ptr = matrix_store.nonscalar_hash_table;
    for (hash_value = 0; hash_value < max; ++hash_value) {
        hash_node_t *node_ptr = table_ptr->heads[hash_value];

//        if (node_ptr) fprintf(stderr,"hash_value = %lu, head node_ptr = %p\n",
//                hash_value, node_ptr);

        while (node_ptr) {
            int64_t m_pID = node_ptr->packedID;
            if (m_pID == MATRIX_ID_INVALID)
            {
                fprintf(stderr,"ERROR in %s:\n",__func__);
                fprintf(stderr,"found matrixStore node with invalid packedID ");
                fprintf(stderr,"while traversing hash chain %lu\n",hash_value);
//                fprintf(stderr,"Removing the node, but ");
//                fprintf(stderr,"it's likely that cleaning is bugged.\n");
//                node_ptr = hash_node_remove_node(matrix_store.nonscalar_hash_table,
//                        node_ptr, hash_value);
//                continue;
         	fprintf(stderr,"The scalar store is messed up, exiting.\n");
		exit(1);

            }
            record_ptr_t r_ptr = get_recordPTR_from_pID(m_pID,"",__func__,0);
            if (r_ptr == RECORD_PTR_INVALID)
            {
                fprintf(stderr,"ERROR in %s:\n",__func__);
                fprintf(stderr,"packedID %lu with null record ",m_pID);
                fprintf(stderr,"found while traversing matrixStore ");
                fprintf(stderr,"hash chain %lu\n", hash_value);
//                fprintf(stderr,"Removing the node, but ");
//                fprintf(stderr,"it's likely that cleaning is bugged.\n");
//                node_ptr = hash_node_remove_node(matrix_store.nonscalar_hash_table,
//                        node_ptr, hash_value);
//                continue;
         	fprintf(stderr,"The scalar store is messed up, exiting.\n");
		exit(1);
            }
            else
            {
                if (verbose)
                {
                    fprintf(stderr,"remove_matrix_from_nonscalar_store_by_");
                    fprintf(stderr,"recordPTR_andor_nodePTR\n");
                    fprintf(stderr,"called from %s with pID %ld\n",
                            __func__, m_pID);
                    fprintf(stderr,"since this routine is recursive, the");
                    fprintf(stderr,"next node pointer is not determined\n");
                    fprintf(stderr,"until all submatrices have been removed\n");
                }
                node_ptr =
                remove_matrix_from_nonscalar_store_by_recordPTR_andor_nodePTR(
                    r_ptr, node_ptr, hash_value, m_pID);
            }
        }
    }
    return 1;
}

/************************************************************************
 * Clean matrix storage	                                                *
 *                                                                      *
 *  Removes all eligible matrices from the matrix store.                *
 *      The called function, remove_matrix_from_store, will         *
 *      recursively remove eligible children of any deleted matrix      *      	
 ***********************************************************************/
int clean_matrix_storage()
{
    clean_nonscalar_matrices();
#ifdef MAR
    clean_scalar_matrices_MAR(&matrix_store);
#else
    clean_scalar_matrices_SPR(&matrix_store);
#endif
    return 1;
}
 
int64_t get_counts_square_matrices_by_level(uint64_t i)
{
  int64_t diag_hist;

  if ((i < matrix_store.largest_level+1) && (i >= 0))
    {
      // allocate space
      // diag_hist = (int64_t *) calloc(matrix_store.largest_level+1,sizeof(int64_t));
      // // warning this may cause a memory leak
      diag_hist = (int64_t)matrix_store.hist[i][i];
      return(diag_hist);
    }
  else return(-1);

}


int64_t nonscalar_hash_chain_length(uint64_t i)
{
  
  int counter = 0;
  
  if ((i < matrix_store.nonscalar_hash_table->nentries) && (i >= 0))
    {
      hash_table_t *table_ptr = matrix_store.nonscalar_hash_table;
      hash_node_t *node_ptr = table_ptr->heads[i];
      // matns_ptr_t record_ptr;
      
      while (node_ptr)  {
        ++counter;
	// record_ptr = (matns_ptr_t) node_ptr->record_ptr;
	// uint64_t matrixID =  record_ptr->packedID;
	// ret = matrix_info_to_file(matrixID,f);
	node_ptr = node_ptr->next;
      }

      return counter;
    }  
  else return(-1);
}

int64_t scalar_hash_chain_length(uint64_t i)
{
  
  int counter = 0;
  
  if ((i < matrix_store.scalar_hash_table->nentries) && (i >= 0))
    {
      hash_table_t *table_ptr = matrix_store.scalar_hash_table;
      hash_node_t *node_ptr = table_ptr->heads[i];
      // matns_ptr_t record_ptr;
      
      while (node_ptr)  {
        ++counter;
	// record_ptr = (matns_ptr_t) node_ptr->record_ptr;
	// uint64_t matrixID =  record_ptr->packedID;
	// ret = matrix_info_to_file(matrixID,f);
	node_ptr = node_ptr->next;
      }

      return counter;
    }  
  else return(-1);
}



//TODO: needs updating once stores are separate
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
  hashstats_to_files(matrix_store.nonscalar_hash_table, accesses_file, nodes_file, report_file);
}
#endif

/****************************************************************** 
   The routine get_readableString_scalar_from_pID_and_coords returns a
   string version of the scalar located in (row,col) of the matrix with
   matrix packedID pID
   ALERT: It has a memory leak, because it is allocating space for
   the string and not deleting it.
********************************************************************/
char *get_readableString_scalar_from_pID_and_coords(int64_t pID,
    int64_t row, int64_t col)
{
    // call matrix pointer version of function with scalarType
    // scalarType val;
    // sca_init(&val);
    int64_t s_pID = get_scalarID_from_pID_and_coords(pID, row, col);
    mats_ptr_t sPTR = (mats_ptr_t)get_recordPTR_from_pID(s_pID,"",__func__,0);
    char *return_String = sca_get_readable_approx_str(sPTR->scalar_value);
    // sca_clear(&val);
    return return_String;
}


/****************************************************************** 
   The routine get_exactString_scalar_from_pID_and_coords returns a string version
   of the scalar located in (row,col) of the matrix with packedID pID
   ALERT: It has a memory leak, because it is allocating space for
   the string and not deleting it.
   TODO: perhaps have matrixID in the name.
********************************************************************/
char *get_exactString_scalar_from_pID_and_coords(int64_t pID, int64_t row, int64_t col)
{
    mats_ptr_t sPTR = get_scalarPTR_from_pID_and_coords(row, col, pID);
    return sca_get_exact_str(sPTR->scalar_value);
}

// Create a new matrix that is the same as m but with the entry at
// (row,col) set to v. This requires saving the panel values at each level
// of recursion and modifying the appropriate panel entry for the submatrix
// which has changed due to changing the scalar at (row,col)

/*!
 * \ingroup larc
 * \brief Set a value within a matrix.
 * \param pID The packedID for the input matrix
 * \param row The row within the matrix.
 * \param col The column within the matrix.
 * \param val_pID The packedID of the value to set.
 * \return The matrixID of a new matrix that is the same as the input matrix, except with the value at the specified row and column set to the new value.
 */
static int64_t replace_scalar_in_matrix_by_pID_and_coords(int64_t pID, int64_t row, int64_t col, int64_t val_pID)
{
  if (pID == MATRIX_ID_INVALID) return MATRIX_ID_INVALID;

  // handle the scalar case
  if (IS_SCALAR(pID))
  {
    if (row | col)
    {
      fprintf(stderr,"ERROR in %s: scalar matrixID %ld\n",__func__,pID);
      fprintf(stderr,"passed with row = %" PRId64 " and col = %" PRId64".\n",
                         row, col);
      exit(1);
    }
    return val_pID;
  }

  // get the matrix pointer from the matrixID
  matns_ptr_t m = (matns_ptr_t)get_recordPTR_from_pID(pID, "", __func__, 0);

  // calculate the number of rows and columns in the matrix
  uint64_t num_rows = (uint64_t)1 << (m->row_level);
  uint64_t num_cols = (uint64_t)1 << (m->col_level);

  // check to see if the indices for row, col are in span for the matrix
  if (row >= num_rows){
    fprintf(stderr,"ERROR in %s: the row %" PRId64 "\n", __func__ ,row);
    fprintf(stderr,"is too large for matrixID %ld\n", pID);
    exit(1);
  }

  if (col >= num_cols){
    fprintf(stderr,"ERROR in %s: the col %" PRId64 "\n", __func__, col);
    fprintf(stderr,"is too large for matrixID %ld\n", pID);
    exit(1);
  }

  int64_t sub_pID[4];

  // Handle the one-row case
  if (1 == num_rows){
      sub_pID[0] = m->subMatList[0];
      sub_pID[1] = m->subMatList[1];

      if (2*col < num_cols) {
            sub_pID[0] = replace_scalar_in_matrix_by_pID_and_coords(
                    sub_pID[0], row, col, val_pID);
      }
      else {
            sub_pID[1] = replace_scalar_in_matrix_by_pID_and_coords(
                    sub_pID[1], row, col - num_cols/2, val_pID);
      }

      return join(sub_pID[0], sub_pID[1]);
  }

  // Handle the one-column case
  if (1 == num_cols){
      sub_pID[0] = m->subMatList[0];
      sub_pID[1] = m->subMatList[2];

      if (2*row < num_rows)
            sub_pID[0] = replace_scalar_in_matrix_by_pID_and_coords(
                    sub_pID[0], row, col, val_pID);
      else
            sub_pID[1] = replace_scalar_in_matrix_by_pID_and_coords(
                    sub_pID[1], row - num_rows/2, col, val_pID);

      return stack(sub_pID[0], sub_pID[1]);
  }

  // At this point matrices have at least two columns and at least two rows
  for (int quad = 0; quad < 4; quad++) sub_pID[quad] = m->subMatList[quad];

  if (2*row < num_rows){
    // element to be changed is in top half of matrix
    if (2*col < num_cols)
        sub_pID[0] = replace_scalar_in_matrix_by_pID_and_coords(sub_pID[0],
            row         , col         , val_pID);
    else sub_pID[1] = replace_scalar_in_matrix_by_pID_and_coords(sub_pID[1],
            row         , col - num_cols/2, val_pID);
  } 
  else {
    // element to be changed is in bottom half of matrix
    if (2*col < num_cols)
        sub_pID[2] = replace_scalar_in_matrix_by_pID_and_coords(sub_pID[2],
            row - num_rows/2, col         , val_pID);
    else sub_pID[3] = replace_scalar_in_matrix_by_pID_and_coords(sub_pID[3],
            row - num_rows/2, col - num_cols/2, val_pID);
  }

  return get_pID_from_array_of_four_sub_pIDs(sub_pID,
                              m->row_level, m->col_level);
}

// MatrixID wrapper for replace_scalar_in_matrix_by_pID_and_coords. 
int64_t replace_scalar_in_matrix_by_string_and_coords(int64_t pID,
        int64_t row, int64_t col, char *val)
{
    if (pID == MATRIX_ID_INVALID) return MATRIX_ID_INVALID;

    // find matrixID of val
    int64_t val_pID = get_valID_from_valString(val);

    // handle scalar case
    if (IS_SCALAR(pID))
    {
       if ( row | col )
       {
          fprintf(stderr,"ERROR in %s: scalar with\n",__func__);
          fprintf(stderr,"matrixID %ld passed with row=%ld and col=%ld\n",
                  pID, row, col);
          exit(1);
       }
       return val_pID;
    } 

    // return input matrix if it already has val set at (row, col)
    int64_t orig_pID = get_scalarID_from_pID_and_coords(pID, row, col);
    if (orig_pID == val_pID) return pID;

    // call version of function with matrixID of value
    return replace_scalar_in_matrix_by_pID_and_coords(pID,
                row, col, val_pID);
}

mats_ptr_t get_scalarPTR_from_pID_and_coords( uint64_t row_i, uint64_t col_j,
	int64_t m_pID)
{
  if (m_pID == MATRIX_ID_INVALID) return NULL;

  if (IS_SCALAR(m_pID)) // matrix is scalar
  {
    if (row_i || col_j)
    {
      fprintf(stderr,"in %s: non-zero row or column index\n",__func__);
      fprintf(stderr,"into matrixID %ld: row=%ld, col=%ld\n",m_pID,row_i,col_j);
      return SCALAR_PTR_INVALID;
    }
    return (mats_ptr_t)get_recordPTR_from_pID(m_pID,"",__func__,0);
  }

  int64_t sca_pID =  get_scalarID_from_pID_and_coords(m_pID, row_i, col_j);
  return (mats_ptr_t)get_recordPTR_from_pID(sca_pID,"",__func__,0);
}

int64_t get_scalarID_from_pID_and_coords(int64_t m_pID, int64_t row_i,
            int64_t col_j) 
{
  if (m_pID == MATRIX_ID_INVALID) return MATRIX_ID_INVALID;

  if (IS_SCALAR(m_pID)) // matrix is scalar
  {
    if (row_i || col_j)
    {
      fprintf(stderr,"in %s: non-zero row or column index\n",__func__);
      fprintf(stderr,"into matrixID %ld: row=%ld, col=%ld\n",m_pID,row_i,col_j);
      return MATRIX_ID_INVALID;
    }
    return m_pID;
  }

  matns_ptr_t m_ptr = (matns_ptr_t)get_recordPTR_from_pID(m_pID,"",__func__,0);
  // calculate the total number of rows and columns in the matrix
  int row_level = m_ptr->row_level;
  int col_level = m_ptr->col_level;
  uint64_t num_rows = (uint64_t)1 << row_level;
  uint64_t num_cols = (uint64_t)1 << col_level;

  // check to see if the matrix big enough for input row and column indices
  if ( (row_i >= num_rows) || (col_j >= num_cols) )
  {
    fprintf(stderr,"ERROR in %s: matrix dimensions too small\n",__func__);
    fprintf(stderr,"for given input row/column\n");
    fprintf(stderr,"\tmatrix size is %lu X %lu\n",num_rows,num_cols);
    fprintf(stderr,"\trow_i = %lu, col_j = %lu\n",row_i, col_j);
    return MATRIX_ID_INVALID;
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
  uint64_t half_i_dim;  // half size i dimension
  uint64_t i_adjust;    // 0 or half_i_dim depending on size of i
  uint64_t i_top_bit = 0;
  uint64_t new_i = 0;  // mod off half dim
  uint64_t half_j_dim;  // half size j dimension
  uint64_t j_adjust; // (0 or half_j_dim) depending on size of j
  uint64_t j_top_bit = 0;
  uint64_t new_j = 0; // mod off half dim

  // Case: mat_type  MATRIX / COL_VECTOR, calculate new_i, and i_top_bit
  if (row_level > 0) {
    half_i_dim = (uint64_t)1 << (row_level-1);  // half size i dimension
    i_adjust = half_i_dim & row_i;   // 0 or half_i_dim depending on size of i
    i_top_bit = i_adjust >> (row_level-1);
    new_i = row_i - i_adjust; // mod off half dim
  } 

  // Case: mat_type  MATRIX / ROW_VECTOR, calculate new_j, and j_top_bit
  if  (col_level > 0) {
    half_j_dim = (uint64_t)1 << (col_level-1);  // half size j dimension
    j_adjust = half_j_dim & col_j;   // (0 or half_j_dim) depending on size of j
    j_top_bit = j_adjust >> (col_level-1);
    new_j = col_j - j_adjust; // mod off half dim
  }
  
  // Calculate the correct submatrix to use in recursive call
  uint64_t submat_index = i_top_bit * 2 + j_top_bit;
  return get_scalarID_from_pID_and_coords( m_ptr->subMatList[submat_index],
             new_i,new_j);
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
  int verbose = 0;
  uint64_t small_modulus = (uint64_t)1 << small_level;
  int64_t small_row = big_row % small_modulus;  /* small_level low bits */
  int64_t small_col = big_col % small_modulus; /* small_level low bits */
  int quad_index = 2*(small_row>>(small_level-1))+(small_col>>(small_level-1));

  if (verbose) {
    fprintf(stderr,"In %s:\n",__func__);
    fprintf(stderr,"big_row is %" PRId64 ", big_col is %" PRId64 ", small_level is %d\n",
	    big_row, big_col,small_level);
    fprintf(stderr,"small_row is set to %" PRId64 ", small_col is set to %" PRId64 "\n",
	    small_row, small_col);
    fprintf(stderr,"--> quad_index is %d\n", quad_index);
  }
  return quad_index;

}

int64_t get_matrix_with_single_nonzero_at_coords(mat_level_t level,
          int64_t row_i, int64_t col_j, int64_t scalarID)
{

  int verbose = 0;
  if (verbose) {
    fprintf(stderr,"In %s\n",__func__);
    fprintf(stderr,"level = %d, row_i = %" PRId64 ", col_j = %" PRId64 "\n",
	    level, row_i, col_j);
  }

  // If scalarID is to the zero scalar then do not need to do the recursion
  // and can return zero matrix
  if (matrix_is_zero(scalarID)) return get_zero_pID(level,level);
  //fprintf(stderr,"passed nonzero test\n");

  // when at scalar level, just return the passed matrixID
  if (level==0) return scalarID;

  // if level > zero, create a matrix from three zero matrixIDs, and one
  // matrixID for the matrix with the scalar in the appropriate position
  int64_t subMatList[4];

  // Get the matrixID for the zero matrix of one level down
  int64_t zero_matrixID = get_zero_pID(level-1,level-1);

  int i;

  // Loop through and intialize all the submatrixIDs to zero matrices
  for (i=0; i < 4 ; i++) {
    subMatList[i] = zero_matrixID;
  }
  //fprintf(stderr,"\tpanels initialized to zero\n");

  // Find out which submatrix should contain the scalar and set by a recursive
  // call
  int submat_number = get_quad_in_submatrix(row_i, col_j, level);

  subMatList[submat_number] = get_matrix_with_single_nonzero_at_coords
                 (level-1, row_i, col_j, scalarID);
  
  // build a matrix from the four panels
  int64_t result_ID = get_pID_from_array_of_four_sub_pIDs(subMatList,
                      level, level);

  return result_ID;
}

void free_matrix_store(void)
{
    // free nodes from matrix store hash table, along with the matrix
    // record pointed to by each node. This routine should only be called as
    // part of a LARC shutdown procedure.

    fprintf(stderr,"freeing memory for matrixStore structure\n");

    // some multiprecision types allocate memory for zerorealthresh

#ifndef MAR
#ifdef IS_RATIONAL
    mpq_clear(matrix_store.zerorealthresh);
#elif defined(IS_MP)
    mpfr_clear(matrix_store.zerorealthresh);
#endif
#endif
    int i;
    // free all arrays allocated in create_matrix_store

    // hist
    for (i=0;i<=matrix_store.largest_level;++i) free(matrix_store.hist[i]);
    free(matrix_store.hist);
    // identity, iHadamard
    free(matrix_store.identity);
    free(matrix_store.iHadamard);
    // zero
    for (i=0;i<=matrix_store.largest_level;++i) free(matrix_store.zero[i]);
    free(matrix_store.zero);
 
    // multiprecision scalarType global variables need to be cleared
    clear_globals();

    // MATRIX HASH TABLE
    fprintf(stderr,"\tfreeing matrix hash table\n");
    hash_table_t *table_ptr = matrix_store.nonscalar_hash_table;
    uint64_t table_size = (uint64_t)1 << table_ptr->exponent;
    uint64_t hashID;

    // hash_node_t contains a record_ptr_t, which for the matrix store
    // hides a pointer to larc_nsmatrix_t. The LARC matrix structure contains
    // a scalarType which may have had memory allocated, and if so it must
    // be cleared before the memory for the larc_nsmatrix_t is itself freed.
    for (hashID = 0; hashID < table_size; hashID ++){
        hash_node_t *node_ptr = table_ptr->heads[hashID];
        while (node_ptr)
        {
            hash_node_t *next_node_ptr = node_ptr->next;
            // free scalarType memory in matrix record
            // (this is a no-op unless scalarType is multiprecision)
            int64_t packedID = node_ptr->packedID;
            matns_ptr_t mat_ptr = (matns_ptr_t)get_recordPTR_from_pID(
                    packedID,"",__func__,0);
#ifdef MAR
            // MAR neighbor nodes point to the primary record, and 
            // we do not free the memory of the primary from them
            if (node_ptr->tile_offset_flag == 0) {
#endif
	    // we do not call sca_init() for the trace_element/scalar_value
	    // scalarType union unless the matrix is square (which includes
	    // scalars as both levels are zero). To avoid freeing an
	    // unallocated pointer, we must test for squareness.
	    if (mat_ptr->row_level == mat_ptr->col_level)
               sca_clear(&(mat_ptr->trace_element));
            // free the record pointer, then the node itself
            free(mat_ptr);
#ifdef MAR
            }
#endif
            free(node_ptr);
            // go to next node in chain
            node_ptr = next_node_ptr;
        }
    } // loop over all hash chains
    dealloc_hash(matrix_store.nonscalar_hash_table);

    // SCALAR HASH TABLE
    fprintf(stderr,"\tfreeing scalar hash table\n");
    table_ptr = matrix_store.scalar_hash_table;
    table_size = (uint64_t)1 << table_ptr->exponent;

    // hash_node_t contains a record_ptr_t, which for the matrix store
    // hides a pointer to larc_nsmatrix_t. The LARC matrix structure contains
    // a scalarType which may have had memory allocated, and if so it must
    // be cleared before the memory for the larc_nsmatrix_t is itself freed.
    for (hashID = 0; hashID < table_size; hashID ++){
        hash_node_t *node_ptr = table_ptr->heads[hashID];
        while (node_ptr)
        {
            hash_node_t *next_node_ptr = node_ptr->next;
            // free scalarType memory in matrix record
            // (this is a no-op unless scalarType is multiprecision)
            int64_t packedID = node_ptr->packedID;
            mats_ptr_t sca_ptr = (mats_ptr_t)get_recordPTR_from_pID(
                    packedID,"",__func__,0);
#ifdef MAR
            // MAR neighbor nodes point to the primary record, and 
            // we do not free the memory of the primary from them
            if (node_ptr->tile_offset_flag == 0) {
#endif
	    // we do not call sca_init() for the trace_element/scalar_value
	    // scalarType union unless the matrix is square (which includes
	    // scalars as both levels are zero). To avoid freeing an
	    // unallocated pointer, we must test for squareness.
            sca_clear(&(sca_ptr->scalar_value));
            // free the record pointer, then the node itself
            free(sca_ptr);
#ifdef MAR
            }
#endif
            free(node_ptr);
            // go to next node in chain
            node_ptr = next_node_ptr;
        }
    } // loop over all hash chains
    dealloc_hash(matrix_store.scalar_hash_table);

    fprintf(stderr,"\tfreeing matrixID tracker\n");
    // memory for ptrs_indexed_by_matrixID may be allocated during use of LARC
    // find number of allocated tables from the next assignible matrixID value
    int num_tables_alloc = 1 + ((matrixID_tracker.matrixID_next-1)>>LOG_SIZE_MPT);
    for (i=0;i<num_tables_alloc;++i)
    {
       free(matrixID_tracker.ptrs_indexed_by_matrixID[i]);
    }
    free(matrixID_tracker.ptrs_indexed_by_matrixID);
}


#ifdef MAR
#ifndef SWIG
mats_ptr_t find_record_from_tile_index(uint64_t hash_chain, MAR_tile_index_t *target_index_PTR)
{

  int verbose = 0;
  
  // We find the first node n in the hash chain
  hash_node_t *n = hash_get_chain(matrix_store.scalar_hash_table, hash_chain);

  //  Track the longest traversal of a hash chain  for matrix reports (at least for now)
  uint32_t depth = 0;

  // If HASHSTATS are defined we also track the number of accesses to each hash chain
#ifdef HASHSTATS
   (matrix_store.scalar_hash_table->num_accesses[hash_chain])++;
#endif
   
  // Compute hash filter of tile index we are looking for (i.e., the target tile index).
  unsigned int target_hashfilter =hashfilter_of_tileindex(target_index_PTR, 16);

  if (Gpti_in_use)
     fprintf(stderr,"in %s, reuse of GLOBAL_primary_tile_index_PTR!\n",
             __func__);
  Gpti_in_use = 1;

  // We traverse the hash chain nodes n looking for either a primary 
  // record or a nhbr record associated with this tile and some 
  // previously stored scalar (and its MAR region of one or more tiles).
  // If we find a previously stored tile we return the matrix PTR for the
  // primary tile. Otherwise we return NULL.
  while (n) {

    if (verbose) {
      printf("traversing hash chain looking for matrix record associated with tile index\n");
    }

    // We want the tile_index stored in the matrixRecord for this node
    // eventually a hash filter in the node will help cut down on 
    // memory retrieval however, at the moment we will have
    // to traverse into the matrix record until both:
    //   * this record is a scalar record AND
    //   * the indices are the same for this matrix record's MAR_tile_index_PTR
    //             as for the input target_index_PTR
    //  If we find such a record then (with a method dependent on whether this
    //  tile a primary or nhbr tile) we return the matrix PTR to the primary
    //  tile associated with this MAR region.

    // Get the pointer to the stored matrix Record
    int64_t primary_packedID = n->packedID;
//    mats_ptr_t primary_scaPTR = (mats_ptr_t)(n->record_ptr);
    mats_ptr_t primary_scaPTR = (mats_ptr_t)get_recordPTR_from_pID(
	primary_packedID, "", __func__, 0);

        // Do preliminary check to see if the hash filter is the same.
    if (n->hashfilter == target_hashfilter) {

        //  Check to see if target_index is the primary_index plus appropriate offset.
#ifdef STORE_TILE_INDEX
        if (tile_indices_equal_after_offset(primary_scaPTR->tile, target_index_PTR, n->tile_offset_flag))
#else
        int regionbitparam = get_regionbitparam();
        //printf("Region Bit Parameter = %d\n", regionbitparam);

        get_tile_index(GLOBAL_primary_tile_index_PTR, regionbitparam, primary_scaPTR->scalar_value);
        //printf("Stored version (correct): \n");
        //print_tile_index(primary_scaPTR->tile);
        //printf("Recomputed version (NOT A TRIGRAPH): \n");
        //print_tile_index(GLOBAL_primary_tile_index_PTR);
        //printf("\n");
        if (tile_indices_equal_after_offset(GLOBAL_primary_tile_index_PTR, target_index_PTR, n->tile_offset_flag))
#endif // #ifdef STORE_TILE_INDEX
	{

            // If the tile index for the matrix record matches the input tile index
            // then this is the record we want.

            // We increment the hit count for the hash table as a whole
            // and for this tile.
            // NOTE: the total number of hits for a MARregion will be the sum
            //             of the hits for its primary and nhbd regions.
            matrix_store.scalar_hash_table->hits++;    
            if (!n->hits_maxxed) {
                    n->record_hits++;
                    if (n->record_hits == ((uint32_t) (-1))) {
                        n->hits_maxxed = 1;
                    }
            }

            // We have been incrementing depth as we traversed the hash chain,
            // so we can now check to see whether the total depth of this
            // traversal breaks our record for the longest successful traversal.
            if (depth > matrix_store.deepest) {
                matrix_store.deepest = depth;
            }

            Gpti_in_use = 0;
            return primary_scaPTR;
  
	}  // indices are equal  (we have the right tile)
	    
    }  // hash filters are equal

    // if we have not found the right tile yet in the hash chain keep searching
    n = n->next;
    depth++;
    
  } // end while (hash nodes in chain)
  Gpti_in_use = 0;

  // At this point we have traversed the entire hash chain without finding the
  // sought indices. We increment the statistics on missed hash chain
  // traversals for the matrix store.
  matrix_store.scalar_hash_table->misses++;
  
  return NULL;
}  // find_record_from_tile_index



/*!
 * \ingroup larc
 * \brief Calculate and store the hash filter value.
 *
 * \param new_scalar_node Which node to update
 * \param tile_index_PTR Tile index of that node
 * \param neighbor_offset Neighbor offset relative to tile
 */
void calculate_and_store_hashfilter_and_flags(hash_node_t *new_scalar_node,
                                              const MAR_tile_index_t *tile_index_PTR,
                                              unsigned int neighbor_offset)
{
    new_scalar_node->hashfilter = hashfilter_of_tileindex(tile_index_PTR, 16);
    new_scalar_node->tile_offset_flag = neighbor_offset;
}


mats_ptr_t  insert_primary_scalar_record(const MAR_tile_index_t *target_tile_index_PTR,
				 const uint64_t hash_value,
                                 const scalarType target_scalar)
{
  int verbose = 0;
  if (verbose) printf("in insert_primary_scalar_record\n");

  // Allocate space for the new scalar record for this tile
  mats_ptr_t s_ptr = (mats_ptr_t)calloc(1, sizeof(larc_smatrix_t));
  if (s_ptr == NULL) {
    ALLOCFAIL();
    fprintf(stderr,"Failed allocation for matrix record for new scalar in %s\n",__func__);
    return SCALAR_PTR_INVALID;
  }

  if (verbose) { 
    printf("  allocated space for matrix record\n");
  }

  // Set the basic parameters for the primal tile in the matrix record for target_scalar

  if (matrixID_tracker.matrixID_next+1==0)
  {
    fprintf(stderr,"ERROR: 64-bit matrixID_next about to overflowed!\n");
    exit(0);
  }

  // the following line only works if routine is within file matrix_store.c
  int64_t nextID = matrixID_tracker.matrixID_next++;
  // since this is a scalar matrix, the matrixID will have the scalar bit set
  s_ptr->packedID = PID_FROM_SCALAR_MID(nextID);
  // s_ptr->MARnhbr = 0;   // right now we are creating the primal region
  s_ptr->info = 0;
  s_ptr->hold = 0;
  s_ptr->lock = 0;
  s_ptr->appears_as_sub_count = 0;

  if (verbose) { 
    printf("  set row level, col level, matrixID, initial hold, lock, and sub_count values\n");
    char *scalar_string = sca_get_readable_approx_str(target_scalar);
    printf("ADDING TO MATRIX SCALAR STORE: %" PRId64 " -> %s\n",
	       s_ptr->packedID, scalar_string);
    free(scalar_string);
  }

  // Inside of the larc_smatrix_t there is a union which holds
  // either the trace_element (for nonscalars) or scalar_value (for scalars)
  // scalar_value is a separate copy of the input scalar,
  // not a reference to it. 
  sca_init(&(s_ptr->scalar_value));
  sca_set(&(s_ptr->scalar_value), target_scalar);
    
  // test to see if scalars are zero or one, and set the flags
  if (scratchVars.quick_use_in_use)
    fprintf(stderr,"%s reusing scratchVars.quick_use!\n",__func__);
  scratchVars.quick_use_in_use = 1;

  scalarType *sca_id = &scratchVars.quick_use;
  sca_set_2ldoubles(sca_id, 0.0L, 0.0L);
  s_ptr->iszero = (0 != sca_eq(s_ptr->scalar_value, *sca_id));
  sca_set_2ldoubles(sca_id, 1.0L, 0.0L);
  s_ptr->isid = (0 != sca_eq(s_ptr->scalar_value, *sca_id));
  scratchVars.quick_use_in_use = 0;

#ifdef STORE_TILE_INDEX
  // We now create a MAR_tile_index structure and fill it with values from
  // the target_scalar
  s_ptr->tile = (MAR_tile_index_t *)malloc(sizeof(MAR_tile_index_t));
  if (s_ptr->tile == NULL) { ALLOCFAIL(); }
#ifdef IS_COMPLEX
  mpz_init(s_ptr->tile->real_index); 
  mpz_init(s_ptr->tile->imag_index); 
  mpz_set(s_ptr->tile->real_index, target_tile_index_PTR->real_index);
  mpz_set(s_ptr->tile->imag_index, target_tile_index_PTR->imag_index);
#else  
  mpz_init(s_ptr->tile->index);
  mpz_set(s_ptr->tile->index, target_tile_index_PTR->index);
#endif // #ifdef IS_COMPLEX
#endif // #ifdef STORE_TILE_INDEX

  // Initialize the tile_node list to be NULLs
  for (int i=0;i<4;++i) s_ptr->tile_node[i] = (hash_node_t*)NULL;

  hash_node_t  *primary_node_PTR;

  // At this point we have filled all values into the matrix record and now
  // create a node in the hash chain to load it into the matrix store
 #ifndef HASH_CHAIN_GROWS_AT_HEAD
  // INSERT new items at tail of hash collision chain
  primary_node_PTR = hash_insert_at_tail(matrix_store.scalar_hash_table,
         RECORD_PTR_INVALID, s_ptr->packedID, hash_value);
  calculate_and_store_hashfilter_and_flags(primary_node_PTR, target_tile_index_PTR, 0);
#endif
#ifdef HASH_CHAIN_GROWS_AT_HEAD
  // INSERT new items at head of hash collision chain
  primary_node_PTR = hash_insert_at_head(matrix_store.scalar_hash_table,
         RECORD_PTR_INVALID, s_ptr->packedID, hash_value);
  calculate_and_store_hashfilter_and_flags(primary_node_PTR, target_tile_index_PTR, 0);
#endif

  if (verbose) {
    printf("  successfully appended hash chain with new matrix record containing scalar\n");
  }

  // update matrix store parameters
  matrix_store.hist[0][0]++;

  if (matrix_store.num_scalars+1<matrix_store.num_scalars)
  {
    fprintf(stderr,"ERROR: overflow in matrix_store.num_scalars!\n");
  }

  matrix_store.num_scalars++;
  if( matrix_store.num_scalars > matrix_store.max_num_scalars ) {
    matrix_store.max_num_scalars = matrix_store.num_scalars;
  }

  // Record the matrix pointer in table indexed by matrixID
  // THIS FUNCTION IS LOCAL TO matrix_store.c
  add_recordPTR_to_indexTable((record_ptr_t)s_ptr, s_ptr->packedID);

  // we record the pointer to the hash node inside of the primary scalar record
  //  this facilitates cleaning, statistics keeping, etc.
  s_ptr->tile_node[0] = primary_node_PTR;

  if (verbose) printf("leaving insert_primary_scalar_record\n");
  return s_ptr;

}  // end of insert_primary_scalar_record



hash_node_t  *insert_nhbr_record(uint64_t nhbr_hash,
	                          MAR_tile_index_t *nhbr_tile_index_PTR,
	                          mats_ptr_t primary_ptr,
                                  unsigned int tile_offset)
{
  int verbose = 0;

//  // Allocate space for the new matrix record for this tile
//  matns_ptr_t nhbrMat_ptr = (matns_ptr_t)calloc(1, sizeof(larc_smatrix_t));
//  if (nhbrMat_ptr == NULL) {
//    ALLOCFAIL();
//    fprintf(stderr,"Failed allocation for matrix record for MAR nhbr in %s\n",__func__);
//    exit(1);
//  }

//  if (verbose) { 
//    printf("  allocated space for matrix record of nhbr tile\n");
//  }

//  // Set the basic parameters for the primal tile in the matrix record for target_scalar
//  // the following line only works if routine is within file matrix_store.c
//  nhbrMat_ptr->packedID = primary_ptr->packedID;
//  nhbrMat_ptr->row_level = 0;
//  nhbrMat_ptr->col_level = 0;
//  nhbrMat_ptr->MARnhbr = 1;   //  creating a nhbr tile
//  nhbrMat_ptr->info = 0;
//  nhbrMat_ptr->hold = 0;
//  nhbrMat_ptr->lock = 0;
//  nhbrMat_ptr->appears_as_sub_count = 0;

//  if (verbose) { 
//    printf("  set row level, col level, matrixID, initial hold, lock, and sub_count values\n");
//    char *scalar_string = sca_get_readable_approx_str(primary_ptr->scalar_value);
//    printf("ADDING NHBR RECORD TO MATRIX SCALAR STORE: %ld -> %s\n",
//	       nhbrMat_ptr->packedID, scalar_string);
//    free(scalar_string);
//  }

//  // Inside of the larc_smatrix_t there is a union which holds
//  // either the trace_element (for nonscalars) or scalar_value (for scalars)
//  // scalar_value is a separate copy of the input scalar,
//  // not a reference to it. 
//  sca_init(&(nhbrMat_ptr->scalar_value));
//  sca_set(&(nhbrMat_ptr->scalar_value), primary_ptr->scalar_value);
    
//  // copy flags from primary record
//  nhbrMat_ptr->iszero = primary_ptr->iszero;
//  nhbrMat_ptr->isid  = primary_ptr->isid;

//  // We now create a MAR_tile_index structure and fill it with values from
//  // the nhbr tile.
//  nhbrMat_ptr->tile = (MAR_tile_index_t *)malloc(sizeof(MAR_tile_index_t));
//#ifdef IS_COMPLEX
//  mpz_init(nhbrMat_ptr->tile->real_index); 
//  mpz_init(nhbrMat_ptr->tile->imag_index); 
//  mpz_set(nhbrMat_ptr->tile->real_index, nhbr_tile_index_PTR->real_index);
//  mpz_set(nhbrMat_ptr->tile->imag_index, nhbr_tile_index_PTR->imag_index);
//#else  
//  mpz_init(nhbrMat_ptr->tile->index);
//  mpz_set(nhbrMat_ptr->tile->index, nhbr_tile_index_PTR->index);
//#endif

//  // Put the PTR to the primary matrix in primaryMatrix[0]
//  nhbrMat_ptr->primaryMatrix[0] = primary_ptr;

//  // At this point we have filled all values into the matrix record and now
//  // we will insert the new matrix record into the matrix store by attaching
//  // to a new node in the hash chain, then we return the PTR to the node
//  // so that we can at the nodePTR into the primary record for cleaning purposes.

  // Create a hash node corresponding to the neighbor tile.
  hash_node_t  *nhbr_node_PTR;
  
  // create a node in the hash chain to load it into the matrix store
#ifndef HASH_CHAIN_GROWS_AT_HEAD
  // INSERT new items at tail of hash collision chain
  nhbr_node_PTR = hash_insert_at_tail(matrix_store.scalar_hash_table,
                     RECORD_PTR_INVALID, primary_ptr->packedID, nhbr_hash);
  calculate_and_store_hashfilter_and_flags(nhbr_node_PTR, nhbr_tile_index_PTR, tile_offset);
#endif
#ifdef HASH_CHAIN_GROWS_AT_HEAD
  // INSERT new items at head of hash collision chain
  nhbr_node_PTR = hash_insert_at_head(matrix_store.scalar_hash_table,
                     RECORD_PTR_INVALID, primary_ptr->packedID, nhbr_hash);
  calculate_and_store_hashfilter_and_flags(nhbr_node_PTR, nhbr_tile_index_PTR, tile_offset);
#endif

  if (verbose) {
//    printf("  successfully appended hash chain with new matrix record containing scalar\n");
    printf("  successfully appended hash chain with new hash node for neighbor tile\n");
  }

  // update matrix store parameter for number of nbhr tiles
  matrix_store.num_neighbors++;

  return nhbr_node_PTR;

}  // end of insert_nhbr_record


/* void   insert_nhbr_node_ptr_in_primary_record( */
/* 					      hash_node_t * nhbr_node_ptr, */
/* 					      matns_ptr_t primary_ptr) */
/* { */

/*   int i = 0; */
/*   while (primary_ptr->tile_node[i] != NULL) { */
/*     ++i; */
/*   } */
/*   primary_ptr->tile_node[i] = nhbr_node_ptr; */

/*   return; */

/* } // end of insert_nhbr_node_ptr_in_primary_record */



#endif // SWIG
#endif // MARmode

void check_validity_one_input(int64_t A_pID, const char *callingRoutine,
     record_ptr_t *A_Rptr)
{
  if (A_pID == MATRIX_ID_INVALID)
  {
     fprintf(stderr,"ERROR in %s: invalid matrixID passed\n",callingRoutine);
     fprintf(stderr,"to %s. matrixID is %ld \n",__func__,MID_FROM_PID(A_pID));
     exit(1);
  }
  record_ptr_t Arec_ptr = get_recordPTR_from_pID(A_pID, "",callingRoutine,0);
  if (Arec_ptr == RECORD_PTR_INVALID)
  {
     fprintf(stderr,"ERROR in %s: matrixID %ld is invalid\n",callingRoutine,
        MID_FROM_PID(A_pID));
     fprintf(stderr,"(found by %s)\n",__func__);
     exit(1);
  }
  *A_Rptr = Arec_ptr;
  return;
}

void check_validity_two_input(int64_t A_pID, int64_t B_pID,
     const char *callingRoutine, record_ptr_t *A_Rptr, record_ptr_t *B_Rptr)
{
  if ( A_pID == MATRIX_ID_INVALID || B_pID == MATRIX_ID_INVALID )
  {
     fprintf(stderr,"ERROR in %s: invalid matrixID passed\n",callingRoutine);
     fprintf(stderr,"to %s\n",__func__);
     fprintf(stderr,"matrixID of A is %ld \n",MID_FROM_PID(A_pID));
     fprintf(stderr,"matrixID of B is %ld \n",MID_FROM_PID(B_pID));
     exit(1);
  }

  record_ptr_t Arec_ptr=get_recordPTR_from_pID(A_pID, "first",callingRoutine,0);
  record_ptr_t Brec_ptr=get_recordPTR_from_pID(B_pID, "second",callingRoutine,0);

  int exitflag = 0;
  if (Arec_ptr == RECORD_PTR_INVALID)
  {
     fprintf(stderr,"ERROR in %s: matrixID %ld is invalid\n",callingRoutine,
        MID_FROM_PID(A_pID));
     fprintf(stderr,"(found by %s)\n",__func__);
     exitflag = 1;
  }
  if (Brec_ptr == RECORD_PTR_INVALID)
  {
     fprintf(stderr,"ERROR in %s: matrixID %ld is invalid\n",callingRoutine,
        MID_FROM_PID(B_pID));
     fprintf(stderr,"(found by %s)\n",__func__);
     exitflag = 1;
  }
  if (exitflag) exit(1);

  *A_Rptr = Arec_ptr;
  *B_Rptr = Brec_ptr;
  return;
}
