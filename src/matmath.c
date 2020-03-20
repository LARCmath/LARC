//                          matmath.c 
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

#include <complex.h>

#include "matmath.h"
#include "larc.h"
#include "scalars.h"
#include "global.h"


/* Python Interface for matrix_add(mat_ptr_t A_ptr, mat_ptr_t B_ptr):
   accepts matrixIDs for matrices to be added and converts
   inputs to pointers before calling. */

int64_t 
matrix_add_matrixID(int64_t A_mID, int64_t B_mID) {

  // get the matrix pointers from the matrixID, and see if still in store 
  mat_ptr_t A_ptr = get_matPTR_from_matID(A_mID, "first", __func__,0);
  mat_ptr_t B_ptr = get_matPTR_from_matID(B_mID, "second", __func__,0);

  // call matrix pointer version of function
  mat_ptr_t C_ptr = matrix_add(A_ptr, B_ptr);
  int64_t C_mID = get_matID_from_matPTR(C_ptr);
  return C_mID;
}


mat_ptr_t matrix_add_base(mat_ptr_t A_ptr, mat_ptr_t B_ptr)
{
  mat_ptr_t sum_ptr;

  // FOR SCALAR CASE
  if (matrix_type(A_ptr)==SCALAR)
  {
    scalarType *sum = &scratchVars.submit_to_store;
    sca_add(sum, matrix_trace(A_ptr), matrix_trace(B_ptr));
    
    // FIND OR CREATE APPROPRIATE MATRIX PTR
    sum_ptr = get_valMatPTR_from_val(*sum);
  }
  else
  {
    // ALL NON SCALAR CASES HAVE PANELS
    mat_ptr_t panel[4];
    
    for (int i = 0; i < 4; i++){
      mat_ptr_t Ai_ptr = matrix_sub(A_ptr, i);
      mat_ptr_t Bi_ptr = matrix_sub(B_ptr, i);
      if (Ai_ptr == MATRIX_PTR_INVALID || Bi_ptr == MATRIX_PTR_INVALID)
        panel[i] = MATRIX_PTR_INVALID;
      else
        panel[i] = matrix_add(Ai_ptr, Bi_ptr);
    }
    
    // FIND OR CREATE APPROPRIATE MATRIX PTR
    sum_ptr = get_matPTR_from_array_of_four_subMatPTRs(panel, matrix_row_level(A_ptr), matrix_col_level(A_ptr));
  }
  return sum_ptr;
}

void matrix_add_checks(mat_ptr_t A_ptr, mat_ptr_t B_ptr)
{
  if (matrix_row_level(A_ptr) != matrix_row_level(B_ptr)) {
    fprintf(stderr,"ERROR: attempting to add matrices of different row_level (%d and %d)\n",
           matrix_row_level(A_ptr), matrix_row_level(B_ptr));
    fprintf(stderr,"matrixIDs: for A = %ld, for B = %ld\n", A_ptr->matrixID,B_ptr->matrixID);
    exit(EXIT_FAILURE);
  }
  if (matrix_col_level(A_ptr) != matrix_col_level(B_ptr)) {
    fprintf(stderr,"ERROR: attempting to add matrices of different col_level (%d and %d)\n",
           matrix_col_level(A_ptr), matrix_col_level(B_ptr));
    fprintf(stderr,"matrixIDs: for A = %ld, for B = %ld\n", A_ptr->matrixID, B_ptr->matrixID);
    exit(EXIT_FAILURE);
  }
}

mat_ptr_t matrix_add_shortcuts(mat_ptr_t A_ptr, mat_ptr_t B_ptr)
{
  if (matrix_is_zero(A_ptr)) {
    return B_ptr;
  }
  if (matrix_is_zero(B_ptr)) {
    return A_ptr;
  }
  // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE
  return op_get(SUM, A_ptr, B_ptr);
}


mat_ptr_t matrix_add(mat_ptr_t A_ptr, mat_ptr_t B_ptr)
{
  exit_if_matrix_ptrs_invalid(__func__, 2, A_ptr, B_ptr);

  // VERIFY ADDITION IS PERMITTED
  matrix_add_checks(A_ptr, B_ptr);

  // This part sorts A and B locally so that A < B and you take advantage of commutativity
  if (A_ptr > B_ptr) {
    mat_ptr_t temp = A_ptr;
    A_ptr = B_ptr;
    B_ptr = temp;
  }

  // MATH RELATIONSHIP / IDENTITY SHORT CUTS / OP STORE CHECK
  mat_ptr_t sum_ptr = matrix_add_shortcuts(A_ptr, B_ptr);
  
  if (sum_ptr == MATRIX_PTR_INVALID)
    sum_ptr = matrix_add_base(A_ptr, B_ptr);

  // STORE RESULT IN OPERATIONS STORE 
  op_set(SUM, A_ptr, B_ptr, sum_ptr);

  return sum_ptr;
}

// Special version of matrix add that does not store op in operations store
// after calculating matrix sum. 
mat_ptr_t matrix_add_no_store(mat_ptr_t A_ptr, mat_ptr_t B_ptr)
{
  exit_if_matrix_ptrs_invalid(__func__, 2, A_ptr, B_ptr);

  // VERIFY ADDITION IS PERMITTED
  matrix_add_checks(A_ptr, B_ptr);

  // This part sorts A and B locally so that A < B and you take advantage of commutativity
  if (A_ptr > B_ptr) {
    mat_ptr_t temp = A_ptr;
    A_ptr = B_ptr;
    B_ptr = temp;
  }

  // MATH RELATIONSHIP / IDENTITY SHORT CUTS / OP STORE CHECK
  mat_ptr_t sum_ptr = matrix_add_shortcuts(A_ptr, B_ptr);

  if (sum_ptr == MATRIX_PTR_INVALID)
    sum_ptr = matrix_add_base(A_ptr, B_ptr);

  return sum_ptr;
}


/* Python Interface for matrix_diff(mat_ptr_t A_ptr, mat_ptr_t B_ptr):
   accepts matrixIDs for matrices to be subtracted and converts
   inputs to pointers before calling. */

int64_t 
matrix_diff_matrixID(int64_t A_mID, int64_t B_mID) {

  // get the matrix pointers from the matrixID, and see if still in store 
  mat_ptr_t A_ptr = get_matPTR_from_matID(A_mID, "first", __func__,0);
  mat_ptr_t B_ptr = get_matPTR_from_matID(B_mID, "second", __func__,0);

  // call matrix pointer version of function
  mat_ptr_t C_ptr = matrix_diff(A_ptr, B_ptr);
  int64_t C_mID = get_matID_from_matPTR(C_ptr);
  return C_mID;
}

mat_ptr_t
matrix_diff(mat_ptr_t A_ptr, mat_ptr_t B_ptr)
{
  // EXCEPTION CHECKING
  exit_if_matrix_ptrs_invalid(__func__, 2, A_ptr, B_ptr);
  if (matrix_row_level(A_ptr) != matrix_row_level(B_ptr)) {
    fprintf(stderr,"ERROR: attempting to subtract matrices of different row_level (%d and %d)\n",
           matrix_row_level(A_ptr), matrix_row_level(B_ptr));
    fprintf(stderr,"matrixIDs: for A = %ld, for B = %ld\n", A_ptr->matrixID, B_ptr->matrixID);
    exit(EXIT_FAILURE);
  }
  if (matrix_col_level(A_ptr) != matrix_col_level(B_ptr)) {
    fprintf(stderr,"ERROR: attempting to subtract matrices of different col_level (%d and %d)\n",
           matrix_col_level(A_ptr), matrix_col_level(B_ptr));
    fprintf(stderr,"matrixIDs: for A = %ld, for B = %ld\n", A_ptr->matrixID, B_ptr->matrixID);
    exit(EXIT_FAILURE);
  }

  // MATH RELATIONSHIP / IDENTITY SHORT CUTS
  // If we knowing that B-A has been calculated or that A is zero:                
  // any gain to substituting scalar multiplies with scalar
  // subtractions? 
  if (matrix_is_zero(B_ptr)) {
    return A_ptr;
  }

  // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  mat_ptr_t diff_ptr = op_get(DIFF, A_ptr, B_ptr);
  if (diff_ptr == MATRIX_PTR_INVALID) {     // WAS NOT FOUND
    
    // FOR SCALAR CASE
    //   CALCULATE RESULT OF OPERATION 
    if (matrix_type(A_ptr)==SCALAR)
    {
      scalarType *diff = &scratchVars.submit_to_store;
      sca_mult(diff, matrix_trace(B_ptr), scalarM1);
      sca_add(diff, matrix_trace(A_ptr), *diff);
      
      // FIND OR CREATE APPROPRIATE MATRIX PTR
      diff_ptr = get_valMatPTR_from_val(*diff);
    }
    else
    {
      
      // ALL NON SCALAR CASES HAVE PANELS
      mat_ptr_t panel[4];
      
      for (int i = 0; i < 4; i++) 
      {
        mat_ptr_t Ai_ptr = matrix_sub(A_ptr, i);
        mat_ptr_t Bi_ptr = matrix_sub(B_ptr, i);
        if (Ai_ptr == MATRIX_PTR_INVALID || Bi_ptr == MATRIX_PTR_INVALID)
          panel[i] = MATRIX_PTR_INVALID;
        else
          panel[i] = matrix_diff(Ai_ptr, Bi_ptr);
      }
      
      // FIND OR CREATE APPROPRIATE MATRIX PTR
      diff_ptr = get_matPTR_from_array_of_four_subMatPTRs(panel, matrix_row_level(A_ptr), matrix_col_level(A_ptr));
    }

    // STORE RESULT IN OPERATIONS STORE 
    op_set(DIFF, A_ptr, B_ptr, diff_ptr);
    
  }   // end matrix was not found in op store

  // RETURN VALUE
  return diff_ptr;
}

/* Python Interface for matrix_mult(mat_ptr_t A_ptr, mat_ptr_t B_ptr):
   accepts matrixIDs for matrices to be multiplied and converts
   inputs to pointers before calling. */
int64_t matrix_mult_matrixID(int64_t A_mID, int64_t B_mID)
{
  // get the matrix pointers from the matrixIDs, and see if still in store 
  mat_ptr_t A_ptr = get_matPTR_from_matID(A_mID, "first", __func__,0);
  mat_ptr_t B_ptr = get_matPTR_from_matID(B_mID, "second", __func__,0);

  // call matrix pointer version of function
  mat_ptr_t C_ptr = matrix_mult(A_ptr, B_ptr);
  int64_t C_mID = get_matID_from_matPTR(C_ptr);

  return C_mID;
}

/* Python Interface for matrix_mult(mat_ptr_t A_ptr, mat_ptr_t B_ptr, int cleanThresh):
   accepts matrixIDs for matrices to be multiplied and converts
   inputs to pointers before calling. */
int64_t matrix_mult_clean_matrixID(int64_t A_mID, int64_t B_mID, mat_level_t cleanThresh)
{

  // get the matrix pointers from the matrixIDs, and see if still in store 
  mat_ptr_t A_ptr = get_matPTR_from_matID(A_mID, "first", __func__,0);
  mat_ptr_t B_ptr = get_matPTR_from_matID(B_mID, "second", __func__,0);

  // call matrix pointer version of function
  mat_ptr_t C_ptr = matrix_mult_clean(A_ptr, B_ptr, cleanThresh);
  int64_t C_mID = get_matID_from_matPTR(C_ptr);

  return C_mID;
}

mat_ptr_t matrix_mult(mat_ptr_t A_ptr, mat_ptr_t B_ptr)
{
  //cheap way to be higher than max of levels
  mat_level_t noCleanThresh = matrix_row_level(A_ptr) + matrix_col_level(A_ptr)
    + matrix_row_level(B_ptr) + matrix_col_level(B_ptr) + 1;
  return matrix_mult_clean(A_ptr, B_ptr, noCleanThresh);
}

mat_ptr_t matrix_mult_clean(mat_ptr_t A_ptr, mat_ptr_t B_ptr, mat_level_t cleanThresh)
{
  int verbose = VERBOSE;
  int reporting = 0; // updates user on periodic cleaning and progress
  mat_level_t min_reporting_level = 29; 

  // EXCEPTION CHECKING
  exit_if_matrix_ptrs_invalid(__func__, 2, A_ptr, B_ptr);
  
  if (matrix_col_level(A_ptr) != matrix_row_level(B_ptr)) {
    fprintf(stderr,"ERROR: attempting to multiply matrix with col_level %d by matrix with row_level %d\n", 
          matrix_col_level(A_ptr), matrix_row_level(B_ptr));
    fprintf(stderr,"matrixIDs: for A = %ld, for B = %ld\n", 
          A_ptr->matrixID, B_ptr->matrixID);

    if (verbose>BASIC) {    
      matrix_store_report("stdout"); 
      op_store_report("stdout"); 
      rusage_report(0,"stdout");
    }
    // if you want a hash report add
    // hash_report(store.hash_table, stdout, 0);
    // or hash_report(hash_table_t *table, FILE *fp, int verbose)
    exit(EXIT_FAILURE);
  }

  matrix_type_t mat_type_A = matrix_type(A_ptr);
  matrix_type_t mat_type_B = matrix_type(B_ptr); 

  if (verbose==DEBUG) { // DEBUGGING
    printf("\n Matrix types for A (%ld) and B (%ld) are %d %d\n", A_ptr->matrixID, B_ptr->matrixID, mat_type_A, mat_type_B);
    printf("   Levels for A are %d %d\n",matrix_row_level(A_ptr), matrix_col_level(A_ptr));
    printf("   Levels for B are %d %d\n",matrix_row_level(B_ptr), matrix_col_level(B_ptr));
  }

  mat_level_t max_level = matrix_pair_max_level(A_ptr, B_ptr);
  reporting = reporting && (max_level >= min_reporting_level);

  // MATH IDENTITY SHORT CUTS
  if (matrix_is_id(A_ptr)) {
    if (reporting)
      printf("  L%d: shortcut (first matrix is identity).\n", max_level);
    return B_ptr;
  }
  if (matrix_is_id(B_ptr)) {
    if (reporting)
      printf("  L%d: shortcut (second matrix is identity).\n", max_level);
    return A_ptr;
  }
  if (matrix_is_zero(A_ptr) || matrix_is_zero(B_ptr)) {
    if (reporting)
      printf("  L%d: shortcut (zero matrix).\n", max_level);
    return get_zero_matrix_ptr(matrix_row_level(A_ptr), matrix_col_level(B_ptr));
  }

  // If A and/or B is SCALAR, then problem reduces to scalar_mult
  if (mat_type_A == SCALAR) {
    if (verbose==DEBUG) {
      // DEBUGGING
      printf("In %s calling scalar mult A scalar\n", __func__);
    }
    return scalar_mult(A_ptr, B_ptr);
  }
  if (mat_type_B == SCALAR) {
    if (verbose==DEBUG) {
      // DEBUGGING
      printf("In %s calling scalar mult B scalar\n", __func__);
    }
    return scalar_mult(B_ptr, A_ptr);
  }

  // If A is COL_VECTOR and B is ROW_VECTOR, then problem reduces
  // to kronecker_product.
  if ((mat_type_A == COL_VECTOR) && (mat_type_B == ROW_VECTOR)) {
    if (verbose==DEBUG) {
      // DEBUGGING
      printf("In %s calling kronecker product A col B row\n", __func__);
    }
    return kronecker_product(A_ptr, B_ptr);
  }

  // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  mat_ptr_t product_ptr = op_get(PRODUCT, A_ptr, B_ptr);
  
  if (product_ptr == MATRIX_PTR_INVALID){     // WAS NOT FOUND
 
    if (scratchVars.top_level == -1)
      scratchVars.top_level = max_level;

    if (max_level >= cleanThresh){
      //NOTE: holds are actually counters so when we internally holding
      //and releasing here won't change whether a matrix is held outside
      //of this routine. 
      set_hold_matrix(A_ptr);
      set_hold_matrix(B_ptr);
    }

    mat_ptr_t panel[4];

    // A = [A0, A1]     B = [B0, B1]    C = [C0, C1] = [A0B0 + A1B2, A0B1 + A1B3]
    //     [A2, A3]         [B2, B3]        [C2, C3] = [A2B0 + A3B2, A2B1 + A3B3]
    // 2*i + j = ind | Ai0 | Ai1 | B0j | B1j | Ai0B0j + Ai1B1j = Cind
    //   0   0    0     A0    A1    B0    B2     A0B0 + A1B2   = C0
    //   0   1    1     A0    A1    B1    B3     A0B1 + A1B3   = C1
    //   1   0    2     A2    A3    B0    B2     A2B0 + A3B2   = C2
    //   1   1    3     A2    A3    B1    B3     A2B1 + A3B3   = C3
    for (int i = 0; i < 2; i++){
      for (int j = 0; j < 2; j++){
        // calculate quadrant (i,j) of the new matrix:
        int ind = 2*i+j;
        mat_ptr_t Ai0 = matrix_sub(A_ptr, 2*i);     
        mat_ptr_t Ai1 = matrix_sub(A_ptr, 2*i+1);
        mat_ptr_t B0j = matrix_sub(B_ptr, j);
        mat_ptr_t B1j = matrix_sub(B_ptr, j+2);
        // A_ptr and B_ptr are held at this point, so submatrices 
        // Ai0, Ai1, B0j, B1j cannot be removed. 

        // this replaces the ROW/COL VECTOR vs MATRIX comparisons
        if ((Ai0 == MATRIX_PTR_INVALID) || (Ai1 == MATRIX_PTR_INVALID)
              || (B0j == MATRIX_PTR_INVALID) || (B1j == MATRIX_PTR_INVALID)){
          panel[ind] = MATRIX_PTR_INVALID;
          continue;
        }

        int64_t mat_store_size1, mat_store_size2;
        mat_store_size1 = (int64_t) (matrix_store_matrixCount() + matrix_store_scalarCount());

        mat_ptr_t Product0 = matrix_mult_clean(Ai0, B0j, cleanThresh);

        mat_store_size2 = (int64_t) (matrix_store_matrixCount() + matrix_store_scalarCount());
        if (reporting)
          printf("  L%d: mult0 adds %ld matrix records to store (now total %ld).\n", max_level, mat_store_size2-mat_store_size1, mat_store_size2);

        mat_ptr_t Product1 = matrix_mult_clean(Ai1, B1j, cleanThresh);

        mat_store_size1 = (int64_t) (matrix_store_matrixCount() + matrix_store_scalarCount());
        if (reporting)
          printf("  L%d: mult1 adds %ld matrix records to store (now total %ld).\n", max_level, mat_store_size1-mat_store_size2, mat_store_size1);

        // NOTE: when at zero threshold, there's little point in storing operations because the
        // matrices quickly become invalid. On the other hand, the ops store quickly fills up
        // (at lower levels where there is usually a limit on the number of matrices
        // we don't have a limit because matrices are readded under a different name all the time)
        // Until individual op records can be removed, I don't think we should store. 
        if (cleanThresh == 0)
            panel[ind] = matrix_add_no_store(Product0, Product1);
        else
            panel[ind] = matrix_add(Product0, Product1);

        mat_store_size2 = (int64_t) (matrix_store_matrixCount() + matrix_store_scalarCount());
        if (reporting)
          printf("  L%d: add adds %ld matrix records to store (now total %ld).\n", max_level, mat_store_size2-mat_store_size1, mat_store_size2);

        if (max_level == cleanThresh){
          set_hold_matrix(panel[ind]); 

          clean_matrix_store();  
          empty_op_store();

          mat_store_size1 = (int64_t) (matrix_store_matrixCount() + matrix_store_scalarCount());
          if (reporting)
            printf("  L%d: clean(g) removes %ld matrix records from store (now total %ld).\n", max_level, mat_store_size2-mat_store_size1, mat_store_size1);
        }
        else if (max_level > cleanThresh){
          set_hold_matrix(panel[ind]); 

          // since we just (recursively) cleaned at a lower level, all matrices involved in 
          // the calculations of P0 and P1 were cleaned EXCEPT P0 and P1. 
          // so we'll do a targeted removal of P0 and P1 instead of a general
          // clean - to save time/effort.
          //
          // NOTE: this handles the rare case where Product0==Product1
          // We hold Product0 until after we have tried to remove Product1;
          // if they are the same, attempting to remove Product1 will correctly
          // not occur.
          set_hold_matrix(Product0);
          remove_matrix_from_mat_store_by_matrix(Product1);
          release_hold_matrix(Product0); 
          remove_matrix_from_mat_store_by_matrix(Product0);

          mat_store_size1 = (int64_t) (matrix_store_matrixCount() + matrix_store_scalarCount());
          if (reporting)
            printf("  L%d: clean(t) removes %ld matrix records from store (now total %ld).\n", max_level, mat_store_size2-mat_store_size1, mat_store_size1);

          // the op_store is rubust in that it won't trip over missing matrices
          // so I think we should skip repairing/cleaning the op_store.
        }
      }
    }

    // DOT PRODUCT: A is ROW_VECTOR and B is COL_VECTOR and result is SCALAR
    if ((mat_type_A == ROW_VECTOR) && (mat_type_B == COL_VECTOR)) {
      product_ptr = panel[0];
    }
    else {
      // FIND OR CREATE APPROPRIATE MATRIX PTR
      product_ptr = get_matPTR_from_array_of_four_subMatPTRs(panel, matrix_row_level(A_ptr), matrix_col_level(B_ptr));
    }

    // We missed a case????????
    if (product_ptr == MATRIX_PTR_INVALID) {
      fprintf(stderr,"ERROR: something went very wrong in %s\n", __func__);
      exit(1);
    }

    // STORE RESULT IN OPERATIONS STORE 
    if ((cleanThresh > 0) || (max_level == scratchVars.top_level))
        op_set(PRODUCT, A_ptr, B_ptr, product_ptr);  

    // RELEASE ALL HOLDS MADE LOCALLY
    if (max_level >= cleanThresh) {
      release_hold_matrix(A_ptr);
      release_hold_matrix(B_ptr);
      for (int i = 0; i < 4; i++){
        if (panel[i] != MATRIX_PTR_INVALID){
            release_hold_matrix(panel[i]);
        }
      }
    }

    scratchVars_exitroutine(max_level);

  }      // end matrix-was-not-found-in-op-store
  
  return product_ptr;
}

/* Python Interface for scalar_mult(mat_ptr_t A_ptr, mat_ptr_t B_ptr):
   accepts matrixIDs for matrices to be multiplied and converts
   inputs to pointers before calling. */

int64_t scalar_mult_matrixID(int64_t A_mID, int64_t B_mID) {

  // get the matrix pointers from the matrixIDs, and see if still in store 
  mat_ptr_t A_ptr = get_matPTR_from_matID(A_mID, "first", __func__,0);
  mat_ptr_t B_ptr = get_matPTR_from_matID(B_mID, "second", __func__,0);

  // call matrix pointer version of function
  mat_ptr_t C_ptr = scalar_mult(A_ptr, B_ptr);
  int64_t C_mID = get_matID_from_matPTR(C_ptr);

  return C_mID;
}


mat_ptr_t scalar_mult(mat_ptr_t A_ptr, mat_ptr_t B_ptr)
{
  // EXCEPTION CHECKING
  exit_if_matrix_ptrs_invalid(__func__, 2, A_ptr, B_ptr);
  if (matrix_type(A_ptr) != SCALAR) {
    fprintf(stderr,"ERROR: attempted scalar multiply with nonscalar of type %d\n", matrix_type(A_ptr));
    fprintf(stderr,"matrixIDs: for A = %ld, for B = %ld\n", A_ptr->matrixID, B_ptr->matrixID);
    exit(EXIT_FAILURE);
  }
  
  matrix_type_t mat_type_B = matrix_type(B_ptr);

  // MATH IDENTITY SHORT CUTS
  if (matrix_is_zero(A_ptr) || matrix_is_zero(B_ptr)) {
    return get_zero_matrix_ptr(matrix_row_level(B_ptr), matrix_col_level(B_ptr));
  }
  
  // NO MEMOIZATION OF SCALAR-SCALAR MULTIPLY (faster to recompute value)
  if (mat_type_B == SCALAR) {
      
    scalarType *product = &scratchVars.submit_to_store;
    sca_mult(product, matrix_trace(A_ptr), matrix_trace(B_ptr));

    // FIND OR CREATE APPROPRIATE MATRIX PTR
    mat_ptr_t product_ptr = get_valMatPTR_from_val(*product);
    return (product_ptr);
  }

  // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  mat_ptr_t product_ptr = op_get(KRONECKER, A_ptr, B_ptr);
  if (product_ptr == MATRIX_PTR_INVALID) {   // WAS NOT FOUND

    // CALCULATE RESULT OF OPERATION 
    // ALL NON SCALAR CASES HAVE PANELS
    mat_ptr_t panel[4];

    for (int i = 0; i < 4; i++) {
        mat_ptr_t Bi_ptr = matrix_sub(B_ptr, i);
        if (Bi_ptr == MATRIX_PTR_INVALID)
          panel[i] = MATRIX_PTR_INVALID;
        else
          panel[i] = scalar_mult(A_ptr, Bi_ptr);
    }
	
    //   FIND OR CREATE APPROPRIATE MATRIX ID
    product_ptr = get_matPTR_from_array_of_four_subMatPTRs(panel, matrix_row_level(B_ptr), matrix_col_level(B_ptr));
    
    //   STORE RESULT IN OPERATIONS STORE 
    op_set(KRONECKER, A_ptr, B_ptr, product_ptr);
  }
  
  // RETURN VALUE
  return product_ptr;
}


/* Python Interface for scalar_divide(mat_ptr_t A_ptr, mat_ptr_t B_ptr):
   accepts matrixIDs for matrices to be processed and converts
   inputs to pointers before calling. */

int64_t scalar_divide_matrixID(int64_t A_mID, int64_t B_mID) {

  // get the matrix pointers from the matrixIDs, and see if still in store 
  mat_ptr_t A_ptr = get_matPTR_from_matID(A_mID, "first", __func__,0);
  mat_ptr_t B_ptr = get_matPTR_from_matID(B_mID, "second", __func__,0);

  // call matrix pointer version of function
  mat_ptr_t C_ptr = scalar_divide(A_ptr, B_ptr);
  int64_t C_mID = get_matID_from_matPTR(C_ptr);

  return C_mID;
}


mat_ptr_t scalar_divide(mat_ptr_t A_ptr, mat_ptr_t B_ptr)
{
  // EXCEPTION CHECKING
  exit_if_matrix_ptrs_invalid(__func__, 2, A_ptr, B_ptr);
  if (matrix_type(B_ptr) != SCALAR) {
    fprintf(stderr,"ERROR: attempted scalar divide with nonscalar of type %d\n", matrix_type(B_ptr));
    fprintf(stderr,"matrixIDs: for A = %ld, for B = %ld\n", A_ptr->matrixID, B_ptr->matrixID);
    exit(EXIT_FAILURE);
  }
  
  matrix_type_t mat_type_A = matrix_type(A_ptr);

  // CAN'T DIVIDE BY ZERO
  if (matrix_is_zero(B_ptr)) {
    fprintf(stderr,"ERROR: attempted scalar divide by zero\n");
    fprintf(stderr,"matrixIDs: for A = %ld, for B = %ld\n", A_ptr->matrixID, B_ptr->matrixID);
    exit(EXIT_FAILURE);
  }

  // MATH IDENTITY SHORT CUTS
  if (matrix_is_zero(A_ptr)) {
    return get_zero_matrix_ptr(matrix_row_level(A_ptr), matrix_col_level(A_ptr));
  }
  
  // NO MEMOIZATION OF SCALAR-SCALAR MULTIPLY (faster to recompute value)
  if (mat_type_A == SCALAR) {
      
    scalarType *quotient = &scratchVars.submit_to_store;
    sca_divide(quotient, matrix_trace(A_ptr), matrix_trace(B_ptr));

    // FIND OR CREATE APPROPRIATE MATRIX PTR
    mat_ptr_t quotient_ptr = get_valMatPTR_from_val(*quotient);
    return (quotient_ptr);
  }

  // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  mat_ptr_t quotient_ptr = op_get(QUOTIENT_SCALAR, A_ptr, B_ptr);
  if (quotient_ptr == MATRIX_PTR_INVALID) {   // WAS NOT FOUND

    // CALCULATE RESULT OF OPERATION 
    // ALL NON SCALAR CASES HAVE PANELS
    mat_ptr_t panel[4];

    for (int i = 0; i < 4; i++) {
        mat_ptr_t Ai_ptr = matrix_sub(A_ptr, i);
        if (Ai_ptr == MATRIX_PTR_INVALID)
          panel[i] = MATRIX_PTR_INVALID;
        else
          panel[i] = scalar_divide(Ai_ptr, B_ptr);
    }
	
    //   FIND OR CREATE APPROPRIATE MATRIX ID
    quotient_ptr = get_matPTR_from_array_of_four_subMatPTRs(panel, matrix_row_level(A_ptr), matrix_col_level(A_ptr));
    
    //   STORE RESULT IN OPERATIONS STORE 
    op_set(QUOTIENT_SCALAR, A_ptr, B_ptr, quotient_ptr);
  }

  // RETURN VALUE
  return quotient_ptr;
}


/* rewritten by Jenny Z and Steve C on 2016-02-11  */ 
/* Fast routine for calculating Kronecker product of 
   matrix A with an identity matrix I of level levelI */
/*!
 * \ingroup larc
 * \brief An efficient routine for Kronecker product when the matrix on the left is an identity matrix with power-of-two dimension
 * \param A_ptr The pointer to the matrix on the right
 * \param levelI The log base 2 of the dimension of the identity matrix
 * \return The Kronecker product of I_{level} otimes A
 */
static mat_ptr_t
tensor_with_identity_on_left(mat_ptr_t A_ptr, mat_level_t levelI)
{

  exit_if_matrix_ptrs_invalid(__func__, 1, A_ptr);

  // MATH IDENTITY SHORT CUTS
  //   If Identity matrix I is SCALAR, then return A
  if (levelI == 0) { return A_ptr; }

  // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  mat_ptr_t identity_ptr = get_identity_matrix_ptr(levelI);
  mat_ptr_t product_ptr = op_get(KRONECKER, identity_ptr, A_ptr);

  if (product_ptr == MATRIX_PTR_INVALID) {   // WAS NOT FOUND
    // CALCULATE RESULT OF OPERATION 
    
    // ALL CASES:
    //    calculate Kronecker product recursively
    mat_level_t row_levelA = matrix_row_level(A_ptr);
    mat_level_t col_levelA = matrix_col_level(A_ptr);
    
    mat_ptr_t panel[4] = {A_ptr, 
			  get_zero_matrix_ptr(row_levelA, col_levelA), 
			  get_zero_matrix_ptr(row_levelA, col_levelA),
			  A_ptr};
    
    product_ptr = tensor_with_identity_on_left(
	     get_matPTR_from_array_of_four_subMatPTRs(panel, row_levelA+1, col_levelA+1), levelI-1);
  }   // end create new matrix
  
  // STORE RESULT IN OPERATIONS STORE 
  op_set(KRONECKER, identity_ptr, A_ptr, product_ptr);
  
  // RETURN VALUE
  return product_ptr;
}

/* Python Interface for kronecker_product(mat_ptr_t A_ptr, mat_ptr_t B_ptr):
   accepts matrixIDs for matrices to be added and converts
   inputs to pointers before calling. */

int64_t kronecker_product_matrixID(int64_t A_mID, int64_t B_mID) {

  // get the matrix pointers from the matrixIDs, and see if still in store 
  mat_ptr_t A_ptr = get_matPTR_from_matID(A_mID, "first", __func__,0);
  mat_ptr_t B_ptr = get_matPTR_from_matID(B_mID, "second", __func__,0);

  // call matrix pointer version of function
  mat_ptr_t C_ptr = kronecker_product(A_ptr, B_ptr);
  int64_t C_mID = get_matID_from_matPTR(C_ptr);

  return C_mID;
}

mat_ptr_t
kronecker_product(mat_ptr_t A_ptr, mat_ptr_t B_ptr)
{

  // EXCEPTIONS (more below)
  exit_if_matrix_ptrs_invalid(__func__, 2, A_ptr, B_ptr);

  // MATH IDENTITY SHORT CUTS
  if (matrix_is_zero(A_ptr) || matrix_is_zero(B_ptr)) {
    return get_zero_matrix_ptr(matrix_row_level(A_ptr) + matrix_row_level(B_ptr),
                             matrix_col_level(A_ptr) + matrix_col_level(B_ptr));
  }
  // Use special function in the case when matrix on left is identity matrix
  //  (This function puts the result in the KRONECKER op store)
  if (matrix_is_id(A_ptr)) {
    return (tensor_with_identity_on_left(B_ptr, matrix_row_level(A_ptr)));    
  }
  matrix_type_t mat_type_A = matrix_type(A_ptr);
  matrix_type_t mat_type_B = matrix_type(B_ptr);
  // When either A or B is scalar, then the kronecker product is the scalar
  // product (the result will be stored in the op_store KRONECKER unless both
  // A and B are scalars)
  if (mat_type_A == SCALAR) {
    return scalar_mult(A_ptr, B_ptr);
  }
  if (mat_type_B == SCALAR) {
    return scalar_mult(B_ptr, A_ptr);
  }


  //   CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  mat_ptr_t product_ptr = op_get(KRONECKER, A_ptr, B_ptr);
  if (product_ptr == MATRIX_PTR_INVALID) {   // WAS NOT FOUND
    
    // CALCULATE RESULT OF OPERATION 
    mat_ptr_t panel[4];
    mat_level_t new_row_level = matrix_row_level(A_ptr) + matrix_row_level(B_ptr);
    mat_level_t new_col_level = matrix_col_level(A_ptr) + matrix_col_level(B_ptr);
    
    // CASE: A and B are both ROW_VECTOR 
    // former version had join outside of kronecker product
    if ((mat_type_A == ROW_VECTOR) && (mat_type_B == ROW_VECTOR)) {
      panel[0] = kronecker_product(matrix_sub(A_ptr,0),B_ptr);
      panel[1] = kronecker_product(matrix_sub(A_ptr,1),B_ptr);
      panel[2] = MATRIX_PTR_INVALID;
      panel[3] = MATRIX_PTR_INVALID;
    }
    
    // CASE: A and B are both COL_VECTOR
    // former version did kronecker product with submatrices of B rather than
    // the entire B matrix - fixed
    if ((mat_type_A == COL_VECTOR) && (mat_type_B == COL_VECTOR)) {
      panel[0] = kronecker_product(matrix_sub(A_ptr,0), B_ptr);
      panel[1] = MATRIX_PTR_INVALID;
      panel[2] = kronecker_product(matrix_sub(A_ptr,2), B_ptr);
      panel[3] = MATRIX_PTR_INVALID;
    }
    
    // CASE: A is COL_VECTOR and B is ROW_VECTOR (correct)
    // [a_0,a_2]^T \otimes [b_0,b_1] -> 
    // [[a_0b_0,a_0b_1],[a_2b_0,a_2b_1]
    if ((mat_type_A == COL_VECTOR) && (mat_type_B == ROW_VECTOR)) {
      panel[0] = kronecker_product(matrix_sub(A_ptr, 0), matrix_sub(B_ptr,0));
      panel[1] = kronecker_product(matrix_sub(A_ptr, 0), matrix_sub(B_ptr,1));
      panel[2] = kronecker_product(matrix_sub(A_ptr, 2), matrix_sub(B_ptr,0));
      panel[3] = kronecker_product(matrix_sub(A_ptr, 2), matrix_sub(B_ptr,1));
    }
    
    // CASE: A is ROW_VECTOR and B is COL_VECTOR (correct)
    // [a_0,a_1] \otimes [b_0,b_2]^T -> 
    // [[a_0b_0,a_1b_0],[a_0b_2,a_1b_2]
    if ((mat_type_A == ROW_VECTOR) && (mat_type_B == COL_VECTOR)) {
      panel[0] = kronecker_product(matrix_sub(A_ptr, 0), matrix_sub(B_ptr,0));
      panel[1] = kronecker_product(matrix_sub(A_ptr, 1), matrix_sub(B_ptr,0));
      panel[2] = kronecker_product(matrix_sub(A_ptr, 0), matrix_sub(B_ptr,2));
      panel[3] = kronecker_product(matrix_sub(A_ptr, 1), matrix_sub(B_ptr,2));
    }
    
    // CASE: A is ROW_VECTOR and B is MATRIX (fix by Matt Calef)
    // former version had join outside of kronecker product
    // [a_0,a_1] \otimes [[b_0,b_1],[b_2,b_3]] ->
    // [[a_0 join(b_0,b_1),a_1 join(b_0,b_1)],
    //    [a_0 join(b_2,b_3),a_1 join(b_2,b_3)]
    if ((mat_type_A == ROW_VECTOR) && (mat_type_B == MATRIX)) {
      panel[0] = kronecker_product(matrix_sub(A_ptr,0),join(matrix_sub(B_ptr,0),matrix_sub(B_ptr,1)));
      panel[1] = kronecker_product(matrix_sub(A_ptr,1),join(matrix_sub(B_ptr,0),matrix_sub(B_ptr,1)));
      panel[2] = kronecker_product(matrix_sub(A_ptr,0),join(matrix_sub(B_ptr,2),matrix_sub(B_ptr,3)));
      panel[3] = kronecker_product(matrix_sub(A_ptr,1),join(matrix_sub(B_ptr,2),matrix_sub(B_ptr,3)));
    }
    
    // CASE: A is COL_VECTOR and B is MATRIX
    // [a_0,a_2]^T \otimes B -> [a_0 B, a_2 B]^T
    // -> [a_0 stack(b_0,b_2),a_0 stack(b_1,b_3)],
    //          [a_2 stack(b_0,b_2),a_2 stack(b_1,b_3)]]
    // former version had stack outside of kronecker product
    if ((mat_type_A == COL_VECTOR) && (mat_type_B == MATRIX)) {
        panel[0] = kronecker_product(matrix_sub(A_ptr,0),stack(matrix_sub(B_ptr,0),matrix_sub(B_ptr,2)));
        panel[1] = kronecker_product(matrix_sub(A_ptr,0),stack(matrix_sub(B_ptr,1),matrix_sub(B_ptr,3)));
        panel[2] = kronecker_product(matrix_sub(A_ptr,2),stack(matrix_sub(B_ptr,0),matrix_sub(B_ptr,2)));
        panel[3] = kronecker_product(matrix_sub(A_ptr,2),stack(matrix_sub(B_ptr,1),matrix_sub(B_ptr,3)));
    }
    
    // CASE: A is MATRIX - not necessary to break down B 
    // (simplification by Matt Calef)
    if ((mat_type_A == MATRIX)) {
      for (int i = 0; i < 4 ; i++) {
	panel[i] = kronecker_product(matrix_sub(A_ptr, i), B_ptr);
      }
    }

    // find or create the matrix id of the result
    product_ptr = get_matPTR_from_array_of_four_subMatPTRs(panel, new_row_level, new_col_level);
    
    // ERROR CHECKING (since there are so many cases)
    if (product_ptr == MATRIX_PTR_INVALID) {
      fprintf(stderr,"ERROR: in %s something went bad\n", __func__);

      exit(1);
    }
    
    //   STORE RESULT IN OPERATIONS STORE 
    op_set(KRONECKER, A_ptr, B_ptr, product_ptr);
  }
  
  //   RETURN VALUE
  return product_ptr;
  
}


/********************************************************************
 *                         join()                                   *
 *  Returns the matrix in which two identical sized input matrices  *
 *  have been appended side to side                                 *
 *  Uses the operation store JOIN                                   *
 *******************************************************************/

/* Python Interface for join(mat_ptr_t A_ptr, mat_ptr_t B_ptr):
   accepts matrixIDs for matrices to be added and converts
   inputs to pointers before calling. */

int64_t 
join_matrixID(int64_t A_mID, int64_t B_mID) {

  // get the matrix pointers from the matrixIDs, and see if still in store 
  mat_ptr_t A_ptr = get_matPTR_from_matID(A_mID, "first", __func__,0);
  mat_ptr_t B_ptr = get_matPTR_from_matID(B_mID, "second", __func__,0);

  // call matrix pointer version of function
  mat_ptr_t C_ptr = join(A_ptr, B_ptr);
  int64_t C_mID = get_matID_from_matPTR(C_ptr);

  return C_mID;
}

mat_ptr_t join(mat_ptr_t A_ptr, mat_ptr_t B_ptr) {

  // EXCEPTION CHECKING
  exit_if_matrix_ptrs_invalid(__func__, 2, A_ptr, B_ptr);

  if ( (matrix_row_level(A_ptr) != matrix_row_level(B_ptr)) ||
       (matrix_col_level(A_ptr) != matrix_col_level(B_ptr)) ) {
    fprintf(stderr,"In function %s, input matrices were not the same size.\n",__func__);
    exit(1);
  }

  matrix_type_t mat_type_A = matrix_type(A_ptr);


  //   CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  mat_ptr_t out_ptr = op_get(JOIN, A_ptr, B_ptr);
  if (out_ptr == MATRIX_PTR_INVALID) {   // WAS NOT FOUND
    
    //   CALCULATE RESULT OF OPERATION
    mat_ptr_t panel[4];
    mat_level_t row_level = matrix_row_level(A_ptr);
    mat_level_t new_col_level = matrix_col_level(A_ptr)+1;
    
    // SCALAR CASE
    if (mat_type_A == SCALAR) {
      panel[0] = A_ptr;
      panel[1] = B_ptr;
      panel[2] = MATRIX_PTR_INVALID;
      panel[3] = MATRIX_PTR_INVALID;
    }
    
    // ROW_VECTOR CASE
    // note: with A a row vector, the result of
    //          join(matrix_sub(A_ptr,0),matrix_sub(A_ptr,1));
    // is A_ptr. This case should be caught by memoization, but
    // we could make it explicit
    if (mat_type_A == ROW_VECTOR) {
      panel[0] = join(matrix_sub(A_ptr,0),matrix_sub(A_ptr,1));
      panel[1] = join(matrix_sub(B_ptr,0),matrix_sub(B_ptr,1));
      panel[2] = MATRIX_PTR_INVALID;
      panel[3] = MATRIX_PTR_INVALID;
    }
    
    // COL_VECTOR CASE
    if (mat_type_A == COL_VECTOR) {
      panel[0] = matrix_sub(A_ptr,0);
      panel[1] = matrix_sub(B_ptr,0);
      panel[2] = matrix_sub(A_ptr,2);
      panel[3] = matrix_sub(B_ptr,2);
    }
    
    // MATRIX CASE
    if (mat_type_A == MATRIX) {
      panel[0] = join(matrix_sub(A_ptr,0),matrix_sub(A_ptr,1));
      panel[1] = join(matrix_sub(B_ptr,0),matrix_sub(B_ptr,1));
      panel[2] = join(matrix_sub(A_ptr,2),matrix_sub(A_ptr,3));
      panel[3] = join(matrix_sub(B_ptr,2),matrix_sub(B_ptr,3));
    }

    //   FIND OR CREATE MATRIX PTR
    out_ptr = get_matPTR_from_array_of_four_subMatPTRs(panel,row_level,new_col_level);
    
    //   STORE RESULT IN OPERATIONS STORE 
    op_set(JOIN, A_ptr, B_ptr, out_ptr);

  } // end calculate a new matrix

  //   RETURN VALUE
  return out_ptr;
  
}


/********************************************************************
 *                         stack()                                  *
 *  Returns the matrix in which two identical sized input matrices  *
 *  have been stacked with the first on top of the second.          *
 *  Uses the operation store STACK                                  *
 *******************************************************************/

/* Python Interface for stack(mat_ptr_t A_ptr, mat_ptr_t B_ptr):
   accepts matrixIDs for matrices to be added and converts
   inputs to pointers before calling. */

int64_t 
stack_matrixID(int64_t A_mID, int64_t B_mID) {

  // get the matrix pointers from the matrixIDs, and see if still in store 
  mat_ptr_t A_ptr = get_matPTR_from_matID(A_mID, "first", __func__,0);
  mat_ptr_t B_ptr = get_matPTR_from_matID(B_mID, "second", __func__,0);

  // call matrix pointer version of function
  mat_ptr_t C_ptr = stack(A_ptr, B_ptr);
  int64_t C_mID = get_matID_from_matPTR(C_ptr);

  return C_mID;
}

mat_ptr_t stack(mat_ptr_t A_ptr, mat_ptr_t B_ptr) {

  // EXCEPTION CHECKING
  exit_if_matrix_ptrs_invalid(__func__, 2, A_ptr, B_ptr);
  if ( (matrix_row_level(A_ptr) != matrix_row_level(B_ptr)) ||
       (matrix_col_level(A_ptr) != matrix_col_level(B_ptr)) ) {
    fprintf(stderr,"In function %s, input matrices were not the same size.\n",__func__);
    exit(1);
  }

  matrix_type_t mat_type_A = matrix_type(A_ptr);

  //   CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  mat_ptr_t out_ptr = op_get(STACK, A_ptr, B_ptr);
  if (out_ptr == MATRIX_PTR_INVALID) {  // WAS NOT FOUND
    
    //   CALCULATE RESULT OF OPERATION
    mat_ptr_t panel[4];
    mat_level_t new_row_level = matrix_row_level(A_ptr)+1;
    mat_level_t col_level = matrix_col_level(A_ptr);
    
    // SCALAR CASE
    if (mat_type_A == SCALAR) {
      panel[0] = A_ptr;
      panel[1] = MATRIX_PTR_INVALID;
      panel[2] = B_ptr;
      panel[3] = MATRIX_PTR_INVALID;
    }
    
    // COL_VECTOR CASE
    // note: with A a column vector, the result of
    //          stack(matrix_sub(A_ptr,0),matrix_sub(A_ptr,2));
    // is A_ptr. This case should be caught by memoization, but
    // we could make it explicit
    if (mat_type_A == COL_VECTOR) {
      panel[0] = stack(matrix_sub(A_ptr,0),matrix_sub(A_ptr,2));
      panel[1] = MATRIX_PTR_INVALID;
      panel[2] = stack(matrix_sub(B_ptr,0),matrix_sub(B_ptr,2));
      panel[3] = MATRIX_PTR_INVALID;
    }
    
    // ROW_VECTOR CASE
    if (mat_type_A == ROW_VECTOR) {
      panel[0] = matrix_sub(A_ptr,0);
      panel[1] = matrix_sub(A_ptr,1);
      panel[2] = matrix_sub(B_ptr,0);
      panel[3] = matrix_sub(B_ptr,1);
    }
    
    // MATRIX CASE
    if (mat_type_A == MATRIX) {
      panel[0] = stack(matrix_sub(A_ptr,0),matrix_sub(A_ptr,2));
      panel[1] = stack(matrix_sub(A_ptr,1),matrix_sub(A_ptr,3));
      panel[2] = stack(matrix_sub(B_ptr,0),matrix_sub(B_ptr,2));
      panel[3] = stack(matrix_sub(B_ptr,1),matrix_sub(B_ptr,3));
    }

    // FIND OR CREATE MATRIX ID
    out_ptr = get_matPTR_from_array_of_four_subMatPTRs(panel,new_row_level,col_level);
    
    //   STORE RESULT IN OPERATIONS STORE 
    op_set(STACK, A_ptr, B_ptr, out_ptr);

  } // end calculate new matrix

  //   RETURN VALUE
  return out_ptr;
  
}


/***********************************************************************
 *                    matrix_entrySquared                              *
 *  Returns the matrix in which each scalar has been squared           *
 *     and then multiplied by the given scale factor                   *
 *  Uses the operation store ENTRYSQUARE                               *
 *  This routine should recursively generate the matrix which has      *
 *  for its (i,j) element the norm squared of the (i,j) element of     *
 *  the input matrix multiplied by a scale_factor.                     *
 ***********************************************************************/

/* Python Interface for matrix_entrySquared(mat_ptr_t A_ptr, double scale_factor):
   accepts matrixIDs for matrices to be added and converts
   inputs to pointers before calling. */

int64_t matrix_entrySquared_matrixID(int64_t m_mID, char *scale_factor) {

  // get the matrix pointer from the matrixIDs, and see if still in store 
  mat_ptr_t m_ptr = get_matPTR_from_matID(m_mID, "", __func__,0);

  // convert char* scale_factor into scalarType scale. 
  scalarType scale;
  sca_init(&scale);
  sca_set_str(&scale, scale_factor);

  // put scale_factor into the matrix store and get its pointer *here*
  // so that we don't have to do it an exponential number of times
  mat_ptr_t scale_ptr = get_valMatPTR_from_val(scale);
  sca_clear(&scale);

  // call matrix pointer version of function
  mat_ptr_t C_ptr = matrix_entrySquared(m_ptr, scale_ptr);
  int64_t C_mID = get_matID_from_matPTR(C_ptr);

  return C_mID;
}

mat_ptr_t matrix_entrySquared(mat_ptr_t m_ptr, mat_ptr_t scale_ptr)
{
  // EXCEPTION CHECKING
  exit_if_matrix_ptrs_invalid(__func__, 2, m_ptr, scale_ptr);

  //   this routine is only designed to work on square matrices
  //   (so no submatrix will be MATRIX_PTR_INVALID)
  if (matrix_row_level(m_ptr) != matrix_col_level(m_ptr))  
  {
    fprintf(stderr,"Function %s requires a square matrix\n", __func__);
    exit(1);
  }
  
  // MATH IDENTITY SHORT CUTS
  // we may know squaring elements will leave matrix unchanged
  if (matrix_is_zero(m_ptr)) 
      return m_ptr;
  if (matrix_is_id(m_ptr)) 
      return scalar_mult(scale_ptr, m_ptr);

  // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  mat_ptr_t mSq_ptr = op_get(ENTRYSQUARE, m_ptr, scale_ptr);

  if (mSq_ptr == MATRIX_PTR_INVALID) {  

    // CASE: SCALAR
    if (matrix_type(m_ptr) == SCALAR) 
      { 
	scalarType *element = &scratchVars.submit_to_store;
	scalarType *conj_elt = &scratchVars.calc_conj;
	sca_conj(conj_elt, matrix_trace(m_ptr));
	sca_mult(element, *conj_elt, matrix_trace(m_ptr));
	sca_mult(element, *element, matrix_trace(scale_ptr));
	mSq_ptr = get_valMatPTR_from_val(*element);
      }

    // CASE: NON-SCALAR
    else
      { 
	mat_ptr_t new_sub[4];
	new_sub[0] = matrix_entrySquared(matrix_sub(m_ptr,0), scale_ptr);
	new_sub[1] = matrix_entrySquared(matrix_sub(m_ptr,1), scale_ptr);
	new_sub[2] = matrix_entrySquared(matrix_sub(m_ptr,2), scale_ptr);
	new_sub[3] = matrix_entrySquared(matrix_sub(m_ptr,3), scale_ptr);
	mSq_ptr = get_matPTR_from_array_of_four_subMatPTRs(new_sub,
                matrix_row_level(m_ptr), matrix_col_level(m_ptr));
      }
    
    // STORE RESULT IN OPERATIONS STORE 
    op_set(ENTRYSQUARE, m_ptr, scale_ptr, mSq_ptr);
  }
  
  return mSq_ptr;
}

int64_t iHadamard_times_matrix_matrixID(int64_t A_mID) {

  mat_ptr_t A_ptr = get_matPTR_from_matID(A_mID, "", __func__, 0);
  mat_ptr_t result_ptr = iHadamard_times_matrix(A_ptr);

  return get_matID_from_matPTR(result_ptr);
}

// We have decided for now that default Hadamard will be the integer
// (unnormalized) version
// routine takes the matrix index of matrix A and returns matrix index of iHadamard * A
mat_ptr_t iHadamard_times_matrix(mat_ptr_t A_ptr) {

  // EXCEPTION CHECKING
  exit_if_matrix_ptrs_invalid(__func__, 1, A_ptr);

  mat_level_t row_level = matrix_row_level(A_ptr);
  mat_level_t col_level = matrix_col_level(A_ptr);

  if (row_level == 0) {
    fprintf(stderr,"ERROR: in %s, there is no sensible definition for level 0 Hadamard\n",__func__);
    exit(1);
  }

  // find the index of the pre-stored iHadamard matrix of this level
  mat_ptr_t iHad_ptr = get_iHadamard_matrix_ptr(row_level);

  mat_ptr_t iHad_times_A_ptr = op_get(PRODUCT, iHad_ptr, A_ptr);

  if (iHad_times_A_ptr == MATRIX_PTR_INVALID){ 
    // if no stored value and the matrix level is 1, calculate the product HH1 * A
    // using the preloaded value for the integer iHadamard
    if (row_level == 1) 
        return matrix_mult(get_iHadamard_matrix_ptr(1), A_ptr);
  
    mat_ptr_t matptr_scalarM1 = get_valMatPTR_from_val(scalarM1);
  
    // if no stored value and the matrix larger than 2x2, calculate the iHadamard
    // product recursively as follows (this should use one less multiply per
    // submatrix than matrix_mult(iHad_ptr,A_ptr)
    mat_ptr_t submatrix[4];
    submatrix[0] = iHadamard_times_matrix(matrix_add(matrix_sub(A_ptr, 0), matrix_sub(A_ptr,2)));
    submatrix[1] = iHadamard_times_matrix(matrix_add(matrix_sub(A_ptr, 1), matrix_sub(A_ptr,3)));
    submatrix[2] = iHadamard_times_matrix(matrix_add(matrix_sub(A_ptr, 0), scalar_mult(matptr_scalarM1, matrix_sub(A_ptr, 2))));
    submatrix[3] = iHadamard_times_matrix(matrix_add(matrix_sub(A_ptr, 1), scalar_mult(matptr_scalarM1, matrix_sub(A_ptr, 3))));
    iHad_times_A_ptr = get_matPTR_from_array_of_four_subMatPTRs(submatrix, row_level, col_level);
  
    // STORE IN OP STORE
    op_set(PRODUCT, iHad_ptr, A_ptr, iHad_times_A_ptr);
  }

  return iHad_times_A_ptr;
}

int64_t matrix_times_iHadamard_matrixID(int64_t A_mID) {

  mat_ptr_t A_ptr = get_matPTR_from_matID(A_mID, "", __func__, 0);
  mat_ptr_t result_ptr = matrix_times_iHadamard(A_ptr);

  return get_matID_from_matPTR(result_ptr);
}

// We have decided for now that default iHadmard will be the integer (unnormalized) version
// routine takes the matrix index of matrix A and returns matrix index of A * iHadmard
mat_ptr_t matrix_times_iHadamard(mat_ptr_t A_ptr) {

  // EXCEPTION CHECKING
  exit_if_matrix_ptrs_invalid(__func__, 1, A_ptr);

  mat_level_t row_level = matrix_row_level(A_ptr);
  mat_level_t col_level = matrix_col_level(A_ptr);

  if (col_level == 0) {
    fprintf(stderr,"ERROR: in %s, there is no sensible definition for level 0 Hadamard\n",__func__);
    exit(1);
  }

  // find the index of the iHadamard matrix of this level
  mat_ptr_t iHad_ptr = get_iHadamard_matrix_ptr(col_level);

  mat_ptr_t A_times_iHad_ptr = op_get(PRODUCT, A_ptr, iHad_ptr); 
  
  if (A_times_iHad_ptr == MATRIX_PTR_INVALID){ 
    // if no stored value and the matrix level is 1, calculate the product A * HH1
    // using the preloaded value for the integer iHadamard
    if (col_level == 1) 
        return matrix_mult(A_ptr, get_iHadamard_matrix_ptr(1));
  
    mat_ptr_t matptr_scalarM1 = get_valMatPTR_from_val(scalarM1);
  
    // if no stored value and the matrix larger than 2x2,
    // calculate the iHadamard product recursively as follows
    // this should be cheaper than matrix_mult(A_ptr,iHad_ptr);
    mat_ptr_t  submatrix[4];
    submatrix[0] = matrix_times_iHadamard(matrix_add(matrix_sub(A_ptr, 0), matrix_sub(A_ptr, 1)));
    submatrix[1] = matrix_times_iHadamard(matrix_add(matrix_sub(A_ptr, 0), scalar_mult(matptr_scalarM1, matrix_sub(A_ptr, 1))));
    submatrix[2] = matrix_times_iHadamard(matrix_add(matrix_sub(A_ptr, 2), matrix_sub(A_ptr, 3)));
    submatrix[3] = matrix_times_iHadamard(matrix_add(matrix_sub(A_ptr, 2), scalar_mult(matptr_scalarM1, matrix_sub(A_ptr, 3))));
    A_times_iHad_ptr = get_matPTR_from_array_of_four_subMatPTRs(submatrix, row_level, col_level);
  
    // STORE IN OP STORE
    op_set(PRODUCT, A_ptr, iHad_ptr, A_times_iHad_ptr);
  }

  return A_times_iHad_ptr;
} 
 

void matrix_tracenorm(scalarType *ret, mat_ptr_t mat_ptr, scalarType scale_factor)
{
  exit_if_matrix_ptrs_invalid(__func__, 1, mat_ptr);

  mat_ptr_t matadj_ptr = matrix_adjoint(mat_ptr);

  // Recall: A * adj(A) = det(A) * I
  mat_ptr_t det_times_id_mat = matrix_mult(matadj_ptr, mat_ptr);

  sca_mult(ret, scale_factor, matrix_trace(det_times_id_mat));
}

char *matrix_tracenorm_matrixID(int64_t m_mID, char *scale_factor)
{
  mat_ptr_t mat_ptr = get_matPTR_from_matID(m_mID, "", __func__, 0);

  scalarType scalar;
  sca_init(&scalar);
  sca_set_str(&scalar, scale_factor);

  scalarType mat_trace; 
  sca_init(&mat_trace);
  matrix_tracenorm(&mat_trace, mat_ptr, scalar);

  char *mat_trace_str = sca_get_str(mat_trace);
  sca_clear(&scalar);
  sca_clear(&mat_trace);

  return mat_trace_str;
}


char *scalar_string_from_matrixID(int64_t m_mID)
{
   // This function checks to see if the matrixID passed to it refers to a
   // matrix of SCALAR type. If it is not a SCALAR, it returns an empty string.
   // If it is a SCALAR, it calls the function matrix_trace_from_matrixID()
   // since the scalar value is currently stored in the matrix store record as
   // "scalarType trace_element". If that ever changes, this function will need
   // to be revised.
   // The string returned is stored in malloc'd memory. If this function
   // is to be called frequently enough to lock up a substantial amount of
   // memory, the returned strings should be freed.

   // do basic error checking
   mat_ptr_t mat_ptr = get_matPTR_from_matID(m_mID, "",__func__,0);     
   exit_if_matrix_ptrs_invalid(__func__, 1, mat_ptr);

   // confirm that the matrix is a SCALAR
   matrix_type_t matrixType = matrix_type(mat_ptr);
   if (matrixType != SCALAR)
   {
      fprintf(stderr,"ERROR: matrixID %zd passed to %s",m_mID,__func__);
      fprintf(stderr,"is not a level (0,0) matrix (SCALAR).\n");
      fprintf(stderr,"Returning empty string.\n");
      char *ret_string = malloc(1*sizeof(char));
      ret_string = "";
      return ret_string;
   }

   // get the scalar value in string format
   return matrix_trace_from_matrixID(m_mID);
}

char *matrix_trace_from_matrixID(int64_t m_mID)
{
   // This routine returns a string giving the trace of the matrix m_mID.
   // If the matrix is a SCALAR, the trace of the matrix is the value of
   // the scalar. The trace is not defined for a non-square MATRIX (or any
   // ROW_VECTOR or COL_VECTOR); in these cases an empty string is returned.
   // In any case, the string returned is malloc'd; if this function is to
   // be called frequently enough to lock up substantial amounts of memory,
   // the returned strings should be freed.

   // do basic error checking
   mat_ptr_t mat_ptr = get_matPTR_from_matID(m_mID, "",__func__,0);     
   exit_if_matrix_ptrs_invalid(__func__, 1, mat_ptr);

   // confirm that the matrix has a well-defined trace
   matrix_type_t matrixType = matrix_type(mat_ptr);
   if (matrixType != SCALAR)
   {
      int call_error = 0;
      if ( (matrixType==ROW_VECTOR) || (matrixType==COL_VECTOR) || 
         (matrix_row_level(mat_ptr) != matrix_col_level(mat_ptr)) )
         { call_error = 1; }
      if (call_error)
      {
         fprintf(stderr,"ERROR: matrixID %zd passed to %s",m_mID,__func__);
         fprintf(stderr,"is for a non-square matrix.\n");
         fprintf(stderr,"Returning empty string.\n");
         char *ret_string = malloc(1*sizeof(char));
         ret_string = "";
         return ret_string;
      }
   }

   // return the trace of the matrix in string format
   char *mat_trace_str = sca_get_str(matrix_trace(mat_ptr));
   return mat_trace_str;
}


/* Python Interface for matrix_add(mat_ptr_t A_ptr, mat_ptr_t B_ptr):
   accepts matrixIDs for matrices to be added and converts
   inputs to pointers before calling. */
char *matrix_maxnorm_matrixID(int64_t mat_mID){
  // get the matrix pointers from the matrixID, and see if still in store
  mat_ptr_t mat_ptr = get_matPTR_from_matID(mat_mID, "", __func__, 0);

  // call matrix pointer version of function with scalarType
  scalarType max_norm;
  sca_init(&max_norm);
  matrix_maxnorm(&max_norm, mat_ptr);

  char *max_norm_str = sca_get_str(max_norm);
  sca_clear(&max_norm);

  return max_norm_str;
}

/********************************************************************
 * This routine uses a user defined norm on scalars.                *
 * It returns the maximum of this norm on all scalars in the matrix.*
 * If there is no norm function defined, then the routine returns   *
 * the largest scalar in the matrix.                                *
 ********************************************************************/
void matrix_maxnorm_custom_norm(scalarType *max, mat_ptr_t mat_ptr, void (*norm)(scalarType*, const scalarType))
{
    exit_if_matrix_ptrs_invalid(__func__, 1, mat_ptr);

    if (matrix_is_zero(mat_ptr)){
        // Nothing needs to be done.
    }

    else if (matrix_type(mat_ptr)==SCALAR)
    {
        scalarType *scalar_temp = &scratchVars.misc;

	// scalar_temp is set to the norm of the scalar value
        if (norm != NULL) {
            norm(scalar_temp, matrix_trace(mat_ptr));
        } else {
            sca_set(scalar_temp, matrix_trace(mat_ptr));
        }

	// if this is the max value seen so far, then keep it
        // sca_cmp(a,b) > 0 iff a > b. 
        if (sca_cmp(*max, *scalar_temp) < 0)
        {
            sca_set(max, *scalar_temp);
        }
    }

    else
    {
        for (int i = 0; i < 4; i++)
        {
            mat_ptr_t panel = matrix_sub(mat_ptr, i);
            if (panel != MATRIX_PTR_INVALID)
            {
                matrix_maxnorm_custom_norm(max, panel, norm);
            }
        }   // end loop through four sub panels
	// Error: if somethingToCompareWith is still 0 then we have never calculated any max
    } // else case (have a panel)
}

void matrix_maxnorm(scalarType *max, mat_ptr_t mat_ptr)
{
    exit_if_matrix_ptrs_invalid(__func__, 1, mat_ptr);

    // initialize max to zero
    sca_set_2ldoubles(max, 0.0L, 0.0L);

    matrix_maxnorm_custom_norm(max, mat_ptr, sca_norm);
}


/* Python Interface for matrix_l2norm(mat_ptr_t A_ptr):
   accepts matrixID for the matrix to compute the L2 norm of,
   and converts inputs to pointers before calling
   the internal routine. */
char *matrix_l2norm_matrixID(int64_t mat_mID){
  // get the matrix pointers from the matrixID, and see if still in store
  mat_ptr_t mat_ptr = get_matPTR_from_matID(mat_mID, "", __func__, 0);

  // call matrix pointer version of function with scalarType
  scalarType l2_norm;
  sca_init(&l2_norm);
  matrix_l2norm(&l2_norm, mat_ptr);

  char *l2_norm_str = sca_get_str(l2_norm);
  sca_clear(&l2_norm);

  return l2_norm_str;
}

void matrix_l2norm(scalarType *l2norm, mat_ptr_t mat_ptr)
{
    exit_if_matrix_ptrs_invalid(__func__, 1, mat_ptr);

    // initialize running sum to zero
    sca_set_2ldoubles(l2norm, 0.0L, 0.0L);

    // accumulate the squares of the L2 norm of all the elements
    matrix_sumnorm_custom_norm(l2norm, mat_ptr, NULL);

    // take the square root of the sum of squares of the element norms
    sca_sqrt(l2norm, *l2norm);
}

/********************************************************************
 * This routine uses a user defined norm on scalars.                *
 * It adds the sum of this norm on all scalars in the matrix        *
 * to the value in running_sum.                                     *
 * If there is no norm function defined, then the routine adds      *
 * the sum of the L2 norm squared for each scalar in the matrix.    *
 * Note: the routine uses a scratch variable to perform             *
 *       computations at the scalar level.                          *
 ********************************************************************/
void matrix_sumnorm_custom_norm(scalarType *running_sum, mat_ptr_t mat_ptr, void (*norm)(scalarType*, const scalarType))
{
    exit_if_matrix_ptrs_invalid(__func__, 1, mat_ptr);

    if (matrix_is_zero(mat_ptr)){
        // Nothing needs to be done.
    }

    else if (matrix_type(mat_ptr)==SCALAR)
    {
        scalarType *scalar_temp = &scratchVars.misc;

	// scalar_temp is set to the norm of the scalar value
        if (norm != NULL) {
            sca_set(scalar_temp, matrix_trace(mat_ptr));
            norm(scalar_temp, *scalar_temp);
        } else {
            // if there is no norm function given, then the norm is
            // the L2 norm squared of the scalar value
            sca_conj(scalar_temp, matrix_trace(mat_ptr));
            sca_mult(scalar_temp, *scalar_temp, matrix_trace(mat_ptr));
        }

	// now add this into the running sum
        sca_add(running_sum, *running_sum, *scalar_temp);
    }

    else
    {
        for (int i = 0; i < 4; i++)
        {
            mat_ptr_t panel = matrix_sub(mat_ptr, i);
            if (panel != MATRIX_PTR_INVALID)
            {
                matrix_sumnorm_custom_norm(running_sum, panel, norm);
            }
        }   // end loop through four sub panels
	// Error: if somethingToCompareWith is still 0 then we have never calculated any max
    } // else case (have a panel)
}

char *matrix_count_entries_matrixID(int64_t mat_mID, char *scalar_str)
{
    scalarType scalar;
    mpz_t ret;

    // interpret scalar_str as scalarType
    sca_init(&scalar);
    sca_set_str(&scalar, scalar_str);

    // get the matrix pointers from the matrixID, and see if still in store
    mat_ptr_t mat_ptr = get_matPTR_from_matID(mat_mID, "", __func__, 0);

    // call matrix pointer version of function
    mpz_init(ret);
    matrix_count_entries(ret, mat_ptr, scalar);

    // prepare string output from mpz_t ret
    size_t out_str_size = mpz_sizeinbase(ret, 10) + 2;
    char *out_str = calloc(out_str_size, sizeof(char));
    if (NULL == out_str){
        fprintf(stderr,"ERROR: allocation failure in %s.\n", __func__);
        exit(1);
    }
    gmp_snprintf(out_str, out_str_size, "%Zd", ret);

    // clean scalarType and mpz_t variables
    sca_clear(&scalar);
    mpz_clear(ret);

    return out_str;
}

void matrix_count_entries_base(mpz_t count, mat_ptr_t mat_ptr, scalarType scalar)
{
    // IDENTITY SHORT CUT
    if (matrix_is_zero(mat_ptr)){
        if (sca_eq(scalar, scalar0)){
            mpz_t *matSize = &scratchVars.counter;
            mpz_ui_pow_ui(*matSize, 2, matrix_row_level(mat_ptr) + matrix_col_level(mat_ptr));
            mpz_add(count, count, *matSize);
        }
        return;
    }

    // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE
    //TODO: how to store in op store? 
    //In the past we would do something like:
    //op_set(ZEROCOUNT, mat_ptr, mat_ptr, cnt_ptr);
    //but cnt_ptr holds scalarType which is not always mpz_t. 

    if (matrix_type(mat_ptr) == SCALAR){
        // if scalar matches matrix entry, increment count by one. 
        if (sca_eq(matrix_trace(mat_ptr), scalar))
            mpz_add_ui(count, count, 1);
        return;
    }
    else
    {
        // increase count in each panel
        for (int i = 0; i < 4; i++)
        {
            mat_ptr_t panel = matrix_sub(mat_ptr, i);
            if (panel != MATRIX_PTR_INVALID)
                matrix_count_entries_base(count, panel, scalar);
        }
    }

    // STORE RESULT IN OPERATIONS STORE
    //op_set(ZEROCOUNT, mat_ptr, mat_ptr, cnt_ptr);

}

void matrix_count_entries(mpz_t count, mat_ptr_t mat_ptr, scalarType scalar)
{
    // EXCEPTION CHECKING
    exit_if_matrix_ptrs_invalid(__func__, 1, mat_ptr);

    // initialize count to 0.
    mpz_set_ui(count, 0);

    matrix_count_entries_base(count, mat_ptr, scalar);
}

char *matrix_list_scalars_larcMatrix(char *path)
{
    scalarType *scalars;
    int64_t numScalars;
    matrix_get_scalars_larcMatrix(&scalars, &numScalars, path);

    int64_t string_size = numScalars + 1;
    for (int64_t i = 0; i < numScalars; i++) {
        char *scalar_string = sca_get_str(scalars[i]);
        string_size += strlen(scalar_string);
        free(scalar_string);
    }

    char *scalarList = calloc(string_size, sizeof(char));
    for (int64_t i = 0; i < numScalars; i++){
        char *scalar_string = sca_get_str(scalars[i]);
        scalarList = strcat(scalarList, scalar_string);
        free(scalar_string);
        scalarList = strcat(scalarList, ",");
    }

    scalarList[strlen(scalarList)-1] = 0;

    return scalarList;
}


// Allocates and fills a list of scalars from a LARCMatrix file, in order of appearance. 
void matrix_get_scalars_larcMatrix(scalarType **scalars_ptr, int64_t *numScalars, char *path)
{
  scalarType *scalars = NULL;
  *numScalars = 0;
  printf("ready to start: num set to %ld\n", *numScalars);
  
  int verbose = VERBOSE;
  FILE *f = fopen(path, "r");
  if (f == NULL)
  {
    fprintf(stderr,"%s: no file found at\n\t%s:\n\texiting\n", __func__, path);
    exit(1);
  }
  json_t *j = j_parse_file(f);

  int64_t max_matrixID = j_lookup_num64(j, "matrixID_max");

  json_t *t = j_key_lookup(j, "table");

  // check that we can allocate an array this large
  // (i.e., total number of bytes does not exceed SIZE_MAX)
  if (SIZE_MAX/(max_matrixID + 1)/sizeof(mat_ptr_t) == 0){
      fprintf(stderr,"Error is %s: matrixID_max too large - try renumbering.\n", __func__);
      exit(1);
  }

  mat_ptr_t *map = calloc(max_matrixID+1, sizeof(mat_ptr_t));
  if (!map){
    fprintf(stderr,"Error in %s: failed to allocate array.\n", __func__);
    exit(1);
  }
  
  int64_t len = j_key_count(t);

  if (verbose>BASIC)
    {printf("In %s, before the table line reader\n", __func__);}
  
  for (int64_t i=0; i<len; i++)
  {
    json_t *p = j_key_index(t, i);
    // index == current matrixID
    int64_t index = (int64_t) atoll(p->name); // maybe atol suffices? 
    if (verbose>BASIC)
      {printf("For i=%"PRId64", and index=%"PRId64"\n",i,index);}
    if (j_is_array(p))
    {
      int len1 = j_array_get_length(p);
      if (verbose>BASIC)
        {printf("The length of this line is %d\n",len1);}
      switch(len1)
      {
        case 3:			/* expect a scalar link */
          if (verbose>BASIC) {printf("In scalar case \n");}
    	  int row_level = j_get_num64(j_array_index(p, 0));
    	  int col_level = j_get_num64(j_array_index(p, 1));
          if (!j_is_string(j_array_index(p,2))){
              fprintf(stderr,"ERROR in %s: expected scalar values to be given as strings. Perhaps try 'legacy' version.\n", __func__);
              exit(1);
          }
          scalars = (scalarType *) realloc(scalars, (*numScalars+1)*sizeof(scalarType));
          if (scalars == NULL){
            fprintf(stderr,"ERROR: Could not grow list of scalars to size %ld in %s.\n", *numScalars+1, __func__);
            exit(0);
          }
          sca_init(&(scalars[*numScalars]));
          const char *val_str = j_get_string(j_array_index(p,2));
          sca_set_str(&(scalars[*numScalars]), val_str);
          *numScalars += 1;
          if (row_level || col_level) {
                  fprintf(stderr,"error in %s:\n\texpected ", __func__ );
                  fprintf(stderr,"scalar, but levels greater than zero!\n");
          }
          break;

        case 6:			/* expect a 4-tuple */
          break;
          
        default:
          fprintf(stderr,"Error in %s(%s): unexpected number of entries per line\n",
            __func__, path);
          fprintf(stderr,"Expected num entries to be 3 or 6, but had %d entries\n",len1);
          exit(1);
	  break;
      }
    }
  }
  *scalars_ptr = scalars;
  j_set_null(j); free(j);
}


gmp_randstate_t state;


/***********************************************************************************
 * random_bool_matrix_from_count, random_bool_matrixID_from_count
 *
 * Note added by JZ 9-2019:
 *
 *   Andy Wills wrote this clever recursive algorithm for producing a 
 *   "random" Boolean matrix with exactly k 1's.  
 *   It is based on a combinatorics counting problem that is sometimes called
 *   "stars and bars", which associates a partition of k into four parts with
 *   a random sequence of k *'s and three divider bars |. 
 *   For example, the partition of 6 into four parts (3,1,0,2) 
 *   is associated with sequence of stars and bars   ***|*||**.
 *   Andy implements the stars and bars decider, by choosing the positions 
 *   of the three bars in this k+3 long string of stars and bars with a 
 *   random number generator.  
 *
 *   This method is recursive and it can return all of the possible matrices 
 *   with exactly k 1's.  
 *
 *   However, one can see that if you choose partitions into four parts with 
 *   equal probability, then it does not produce all matrices with exactly k 1's 
 *   with uniform probability. 
 *     For example, consider all four by four matrices with exactly two 1's 
 *     and fourteen 0's.  The first step of the recursion is to allocate 
 *     how many of the two 1's are in each of the four two by two quadrant 
 *     submatrices.  If we select the allocation (1,1,0,0) as frequently as 
 *     we select the allocation (2,0,0,0), then this does not produce a uniform
 *     distribution of all four by four matrices with exactly two 1's, because
 *     - there are 4 x 4 x 1 x 1 = 16  four by four matrices with 
 *       a single 1 in the first quadrant, a single 1 in the second quadrant, 
 *       and no ones in the third and fourth quadrants, whereas,
 *     - there are 6 x 1 x 1 x 1 = 6  four by four matrices with 
 *       two 1's in the first quadrant, and no ones in the other quadrants.
 *  
 *   TODO: 
 *     - we should decide if we want to warn the user that we are picking one
 *       kind of uniform distribution, but that it is not the same as the more 
 *       standard idea of an uniform distribution of matrices with exactly k 1's. 
 *
 **********************************************************************************/
mat_ptr_t random_bool_matrix_from_count(mat_level_t row_level, mat_level_t col_level,
        mpz_t scalarNum, mat_ptr_t scalarPtr)
{
    int reporting = 0;
    mat_level_t min_reporting_level = 30; 

    if (0 > mpz_cmp_ui(scalarNum, 0)){
        char *trace_string = sca_get_str(matrix_trace(scalarPtr));
        gmp_fprintf(stderr,"ERROR: the number of entries with value %s must be positive in %s, not %Zd.\n", 
                trace_string, __func__, scalarNum);
        free(trace_string);
        exit(0);
    }

    // reporting, if reporting. 
    if ((row_level >= min_reporting_level) || (col_level >= min_reporting_level)){
        if (reporting) {
            char *trace_string = sca_get_str(matrix_trace(scalarPtr));
            gmp_fprintf(stderr,"  %s: requesting a dim level %dx%d with %Zd entries with value %s.\n", 
                    __func__, row_level, col_level, scalarNum, trace_string);
            free(trace_string);
        }
    }

    // MATH SHORTCUTS
    // if no <scalar>s, return zero matrix
    if (0 == mpz_cmp_ui(scalarNum, 0))
        return get_zero_matrix_ptr(row_level, col_level);

    mpz_t total_entries_in_matrix_size;
    mpz_init(total_entries_in_matrix_size); 
    mpz_ui_pow_ui(total_entries_in_matrix_size, 2, row_level);
    mpz_mul_2exp(total_entries_in_matrix_size, total_entries_in_matrix_size, col_level);
    // check that matrix can contain the requested number of 1s
    int requested_cmpTo_total = mpz_cmp(scalarNum, total_entries_in_matrix_size);
    if (0 == requested_cmpTo_total){
        ;
        // do stuff here if there is a shortcut for creating matrix with all 
        // the same entry
    }
    if (0 < requested_cmpTo_total){
        // note that we won't fail this during recursive calls, because 
        // we check first. 
        fprintf(stderr,"ERROR: function %s failed to create matrix: too many entries requested for matrix of this size!\n", __func__);
        exit(0);
    }

    if ((row_level == 0) && (col_level == 0)){
        mpz_clear(total_entries_in_matrix_size);
        return scalarPtr;
    }

    // decide randomly how the scalarNum <scalar>s should be distributed between the 
    // available quadrants: 
    // initialize random state and set seed based on system time
    // Note: using time to set seed is adequate for our purposes, but not 
    // cryptographically secure. 
    static int first = 1;
    if (first){
        gmp_randinit_default(state);
        gmp_randseed_ui(state, time(NULL));
        first = 0;
    }

    int rands_are_good = 0;
    int i;
    // the number of quadrant submatrices for row and col vectors is 2
    // for other matrices it is 4.   The stars and bars mentioned above
    // needs the number of bars = quads-1.
    int quads = ((col_level * row_level) == 0) ? 2 : 4;
    mpz_t *rands = malloc(quads * sizeof(mpz_t));
    for (i = 0; i < quads; i++)
        mpz_init(rands[i]);
    // Note: mpz_divexact* best choice when exact division is known. 
    mpz_divexact_ui(total_entries_in_matrix_size, total_entries_in_matrix_size, quads);
    // for stars (<scalar>s) and bars (quadrants) method, need #stars + #bars-1. we'll 
    // store it in scalarNum since it already exists and revert it afterwards 
    // (just in case).
    mpz_add_ui(scalarNum, scalarNum, quads-1);
    while (!rands_are_good){
        for (i = 0; i < quads-1; i++){
            // there are scalarNum ones to distribute between 4 quadrants so with the 
            // stars and bars method, we need 3 random numbers between 0 and 
            // scalarNum + 3 (exclusive). 
            // GMP routine mpz_urandomm(rand,state,n) generates a random integer 
            // in the range 0 to n-1, inclusive. Earlier, we added 3 to scalarNum
            // for this purpose so
            //      mpz_urandomm(rands[i], state, scalarNum); 
            // is correct while
            //      mpz_urandomm(rands[i], state, scalarNum+3); 
            // is not (and cannot be done with mpz_t types anyway).
	    // If it is a row or col vector we added 1 to scalarNum instead of 3.
            mpz_urandomm(rands[i], state, scalarNum); 
            // check that rands[i] is not already in the list
            for (int j = 0; j < i; j++){
                if (0 == mpz_cmp(rands[j], rands[i])){
                    i --;
                    break;
                }
            }
        }
        // fill in the last bin number
        mpz_set(rands[quads-1], scalarNum); //recall scalarNum is temporarily 3 higher
        // sort the random values so we can interpret as bins
        if (quads == 4){
            if (0 < mpz_cmp(rands[0], rands[1]))
                mpz_swap(rands[0], rands[1]);
            if (0 < mpz_cmp(rands[1], rands[2]))
                mpz_swap(rands[1], rands[2]);
            if (0 < mpz_cmp(rands[0], rands[1]))
                mpz_swap(rands[0], rands[1]);
        }
        else if (quads == 2){
            // no sort needed in this case. 
            ;
        }

        //printf("Sorted rands are \n");
        //for (i = 0; i < quads; i++)
        //    gmp_printf("%Zd \n", rands[i]);

        // differences in sorted rands is the numbers of stars per bin
        for (i = quads - 1; i > 0; i--){
            mpz_sub(rands[i], rands[i], rands[i-1]);
            mpz_sub_ui(rands[i], rands[i], 1);
        }
        // verify that we aren't putting too many <scalar>s for a quadrant to handle
        // this should be rare if 1s are sparse!
        rands_are_good = 1; 
        for (i = 0; i < quads; i++){
            if (0 < mpz_cmp(rands[i], total_entries_in_matrix_size)){
                rands_are_good = 0; 
                break;
            }
        }

    }
    // restore scalarNum to its original value
    mpz_sub_ui(scalarNum, scalarNum, quads-1);

    mat_ptr_t panel[4];
    // COL VECTOR
    if (col_level == 0){
        panel[0] = random_bool_matrix_from_count(row_level-1, 0, rands[0], scalarPtr);
        panel[1] = MATRIX_PTR_INVALID;
        panel[2] = random_bool_matrix_from_count(row_level-1, 0, rands[1], scalarPtr);
        panel[3] = MATRIX_PTR_INVALID;
    }
    // ROW VECTOR
    else if (row_level == 0){
        panel[0] = random_bool_matrix_from_count(0, col_level-1, rands[0], scalarPtr);
        panel[1] = random_bool_matrix_from_count(0, col_level-1, rands[1], scalarPtr);
        panel[2] = MATRIX_PTR_INVALID;
        panel[3] = MATRIX_PTR_INVALID;
    }
    // MATRIX
    else {
        for (i = 0; i < 4; i++){
            panel[i] = random_bool_matrix_from_count(row_level-1, col_level-1, rands[i], scalarPtr);
        }
    }

    // cleanup
    mpz_clear(total_entries_in_matrix_size);
    for (i = 0; i < quads; i++)
        mpz_clear(rands[i]);
    free(rands);

    // reporting, if reporting. 
    if ((row_level >= min_reporting_level) || (col_level >= min_reporting_level)){
        if (reporting) {
            char *trace_string = sca_get_str(matrix_trace(scalarPtr));
            gmp_printf("  %s: returning a dim level %dx%d with %Zd entries with value %s.\n",
                    __func__, row_level, col_level, scalarNum, trace_string);
            free(trace_string);
        }
    }

    // FIND OR CREATE APPROPRIATE MATRIX PTR
    return get_matPTR_from_array_of_four_subMatPTRs(panel, row_level, col_level); 
}

int64_t random_bool_matrixID_from_count(mat_level_t row_level, mat_level_t col_level, char* numOnes)
{
    mpz_t num;
    mpz_init(num);

    // scan number of desired ones into an mpz_t struct
    gmp_sscanf(numOnes, "%Zd", num);

    // call matrix pointer version of function 
    mat_ptr_t onePtr = get_valMatPTR_from_val(scalar1);

    mat_ptr_t mat_ptr = random_bool_matrix_from_count(row_level, col_level, num, onePtr);
    mpz_clear(num);

    // get the matrix pointers from the matrixID, and see if still in store
    return get_matID_from_matPTR(mat_ptr);
}

int64_t matrix_basischange_A_by_B_matrixID(int64_t B_mID, int64_t A_mID)
{
    // get the matrix pointers from the matrixID, and see if still in store    
    mat_ptr_t B_ptr = get_matPTR_from_matID(B_mID, "first", __func__, 0);
    mat_ptr_t A_ptr = get_matPTR_from_matID(A_mID, "second", __func__, 0);

    // call matrix pointer version of function
    mat_ptr_t C_ptr = matrix_basischange_A_by_B(B_ptr, A_ptr);
    return get_matID_from_matPTR(C_ptr);
}

mat_ptr_t matrix_basischange_A_by_B(mat_ptr_t B_ptr, mat_ptr_t A_ptr)
{
    // this routine should be used when the goal is to return the matrix which
    // is A conjugated by B. For this operation, there are two matrix
    // multiplies, and we do not care about the intermediate result, so we 
    // clean it from the matrix store when it is no longer needed (with the
    // proper protections so that we don't clean any submatrices needed for
    // other stored matrices).

    // EXCEPTION CHECKING
    exit_if_matrix_ptrs_invalid(__func__, 2, A_ptr, B_ptr);

    // matrix A must be square for conjugation to make sense
    if (matrix_col_level(A_ptr)!= matrix_row_level(A_ptr)) {
       fprintf(stderr,"ERROR: in %s, matrix A_ptr",__func__);
       fprintf(stderr," (with matrixID %" PRIu64 ") is not square!\n",
                A_ptr->matrixID);
       exit(EXIT_FAILURE);
    }

    // mismatch in row level of B and col level of A will be detected 
    // in matrix_mult()

    mat_ptr_t result_ptr = op_get(BASIS_CHANGE,B_ptr,A_ptr);
    if (result_ptr != MATRIX_PTR_INVALID) { // was found in op store
       return result_ptr;
    }

    // result was not found in operations store - calculation proceeds

    // ensure that cleaning doesn't take out the input matrices
    // NOTE: holds are actually counters, so internally holding and releasing
    // here won't change whether a matrix is held outside of this routine
    set_hold_matrix(A_ptr);
    set_hold_matrix(B_ptr);
    // calcluate B*A*B^{\dagger}
    // (for complex types, \dagger is complex conjugate transpose;
    //  for real types, it is simply transpose)
    mat_ptr_t temp1_ptr = matrix_mult(B_ptr,A_ptr);
    mat_ptr_t conjAdj_ptr = matrix_adjoint(B_ptr);
    result_ptr = matrix_mult(temp1_ptr,conjAdj_ptr);
    set_hold_matrix(result_ptr);
    // safely remove all intermediate valued matrices
    set_hold_matrix(conjAdj_ptr);
    remove_matrix_from_mat_store_by_matrix(temp1_ptr);
    release_hold_matrix(conjAdj_ptr);
    remove_matrix_from_mat_store_by_matrix(conjAdj_ptr);
    // release held inputs and output
    release_hold_matrix(A_ptr);
    release_hold_matrix(B_ptr);
    release_hold_matrix(result_ptr);

    // store result of calculation in operations store
    op_set(BASIS_CHANGE, B_ptr, A_ptr, result_ptr);

    return result_ptr;
}

int64_t matrix_saxpy_matrixID(int64_t A_mID, int64_t B_mID,
        int64_t scalar_a_mID, int64_t scalar_b_mID)
{
    // get the matrix pointers from the matrixID, and see if still in store    
    mat_ptr_t A_ptr = get_matPTR_from_matID(A_mID, "first", __func__, 0);
    mat_ptr_t B_ptr = get_matPTR_from_matID(B_mID, "second", __func__, 0);
    mat_ptr_t scalar_a_ptr = get_matPTR_from_matID(scalar_a_mID,
        "third", __func__, 0);
    mat_ptr_t scalar_b_ptr = get_matPTR_from_matID(scalar_b_mID,
        "fourth", __func__, 0);

    // call matrix pointer version of function
    mat_ptr_t C_ptr = matrix_saxpy(A_ptr, B_ptr, scalar_a_ptr, scalar_b_ptr);
    return get_matID_from_matPTR(C_ptr);
}

mat_ptr_t matrix_saxpy(mat_ptr_t A_ptr, mat_ptr_t B_ptr, mat_ptr_t scalar_a_ptr,
        mat_ptr_t scalar_b_ptr)
{
    // this routine performs the saxpy operation that is common in many
    // linear algebra applications. The result C = a*A + b*B, where a and b
    // are scalar values and A and B are similarly-dimensioned matrices
    // (usually row or column vectors).

    // EXCEPTION CHECKING
    exit_if_matrix_ptrs_invalid(__func__, 4, A_ptr, B_ptr, scalar_a_ptr,
        scalar_b_ptr);

    // confirm dimensions of inputs
    if (matrix_type(scalar_a_ptr) != SCALAR)
    {
       fprintf(stderr,"ERROR in %s: 3rd argument should be a scalar,",__func__);
       fprintf(stderr,"but\n\thas row level %d and column level %d!\n",
               matrix_row_level(scalar_a_ptr),matrix_col_level(scalar_a_ptr));
        exit(EXIT_FAILURE);
    }
    if (matrix_type(scalar_b_ptr) != SCALAR)
    {
       fprintf(stderr,"ERROR in %s: 4th argument should be a scalar,",__func__);
       fprintf(stderr,"but\n\thas row level %d and column level %d!\n",
               matrix_row_level(scalar_b_ptr),matrix_col_level(scalar_b_ptr));
        exit(EXIT_FAILURE);
    }

    if (matrix_col_level(A_ptr) != matrix_col_level(B_ptr)
        || matrix_row_level(A_ptr) != matrix_row_level(B_ptr))
    {
       fprintf(stderr,"ERROR: in %s, ",__func__);
       fprintf(stderr,"input matrices do not have the same dimensions!\n");
       fprintf(stderr,"\tmatrix A (with matrixID %" PRIu64 ")",A_ptr->matrixID);
       fprintf(stderr," has row level %d and column level %d",
                matrix_row_level(A_ptr),matrix_col_level(A_ptr));
       fprintf(stderr,"\tmatrix B (with matrixID %" PRIu64 ")",B_ptr->matrixID);
       fprintf(stderr," has row level %d and column level %d",
                matrix_row_level(B_ptr),matrix_col_level(B_ptr));
       exit(EXIT_FAILURE);
    }

    // scalar_mult will check KRONECKER op_store for a*A
    mat_ptr_t result1_ptr = scalar_mult(scalar_a_ptr,A_ptr);
    // scalar_mult will check KRONECKER op_store for b*B
    mat_ptr_t result2_ptr = scalar_mult(scalar_b_ptr,B_ptr);

    // matrix_add will check SUM op store for (a*A) + (b*B)
    mat_ptr_t result_ptr = matrix_add(result1_ptr, result2_ptr);

    return result_ptr;
}

int64_t vector_dot_product_matrixID(int64_t A_mID, int64_t B_mID, int verbose)
{
    // get the matrix pointers from the matrixID, and see if still in store    
    mat_ptr_t A_ptr = get_matPTR_from_matID(A_mID, "first", __func__, 0);
    mat_ptr_t B_ptr = get_matPTR_from_matID(B_mID, "second", __func__, 0);

    // call matrix pointer version of function
    mat_ptr_t C_ptr = vector_dot_product(A_ptr, B_ptr, verbose);
    return get_matID_from_matPTR(C_ptr);
}

mat_ptr_t vector_dot_product(mat_ptr_t A_ptr, mat_ptr_t B_ptr, int verbose)
{
    // This routine returns a matrix pointer for the scalar product
    // of two vectors. MATRIX_PTR_INVALID is returned if either matrix
    // passed is a MATRIX (both row level and col level > 0) or the
    // dimensions of the vectors passed do not agree.
    // If verbose is higher than BASIC, a lot of information about what
    // the routine is doing is printed to stderr.

    // EXCEPTION CHECKING
    exit_if_matrix_ptrs_invalid(__func__, 2, A_ptr, B_ptr);
    
    matrix_type_t A_type = matrix_type(A_ptr);
    matrix_type_t B_type = matrix_type(B_ptr);

    // shortcut - if both "vectors" are 1x1 
    if (A_type==SCALAR && B_type==SCALAR) {
       if (verbose>BASIC)  {
          int64_t A_mID = get_matID_from_matPTR(A_ptr);
          int64_t B_mID = get_matID_from_matPTR(B_ptr);
          fprintf(stderr,"%s doing 'dot product' of",__func__);
          fprintf(stderr,"two scalar values with matrixIDs ");
          fprintf(stderr,"%zd and %zd", A_mID, B_mID);
       }
       return scalar_mult(A_ptr,B_ptr);
    }

    // various tests for cases that don't fit dot product paradigm
    // including the case where only one of A or B is a scalar
    if ( (A_type==MATRIX) || (B_type==MATRIX) ||
         (A_type==SCALAR) || (B_type==SCALAR) )
    {
       fprintf(stderr,"ERROR: %s called with\n",__func__);
       int64_t A_mID = get_matID_from_matPTR(A_ptr);
       int64_t B_mID = get_matID_from_matPTR(B_ptr);
       fprintf(stderr,"matrix A with matrixID %zd has level (%d,%d) and\n",
          A_mID,matrix_row_level(A_ptr),matrix_col_level(A_ptr));
       fprintf(stderr,"matrix B with matrixID %zd has level (%d,%d)\n",
          B_mID,matrix_row_level(B_ptr),matrix_col_level(B_ptr));
       return MATRIX_PTR_INVALID;
    }

    // at this point, A and B are either ROW_VECTOR or COL_VECTOR types
    // check to make sure dimensions of the vectors agree
    int64_t rdimA = matrix_row_level(A_ptr);
    int64_t cdimA = matrix_col_level(A_ptr);
    int64_t rdimB = matrix_row_level(B_ptr);
    int64_t cdimB = matrix_col_level(B_ptr);
    int64_t dimA = (rdimA) ? rdimA : cdimA;
    int64_t dimB = (rdimB) ? rdimB : cdimB;
    if (dimA != dimB) {
       fprintf(stderr,"ERROR: vectors passed to %s\n",__func__);
       fprintf(stderr,"do not have the same nonzero level values\n");
       int64_t A_mID = get_matID_from_matPTR(A_ptr);
       int64_t B_mID = get_matID_from_matPTR(B_ptr);
       fprintf(stderr,"matrix A with matrixID %zd has level (%d,%d) and\n",
          A_mID,matrix_row_level(A_ptr),matrix_col_level(A_ptr));
       fprintf(stderr,"matrix B with matrixID %zd has level (%d,%d)\n",
          B_mID,matrix_row_level(B_ptr),matrix_col_level(B_ptr));
       return MATRIX_PTR_INVALID;
    }

    // if we get this far, assume that the user knows what they are doing, and
    // do a row_vector*col_vector matrix_mult call. Depending on the inputs,
    // we may need to reorder the vectors or apply the adjoint function to one
    // of them.
    if (A_type==ROW_VECTOR)
    {
       if (B_type==COL_VECTOR) {
          if (verbose>BASIC) {
             fprintf(stderr,"in %s, A is a row vector and\n",__func__);
             fprintf(stderr,"B is a column vector - returning A*B\n");
          }
          return matrix_mult(A_ptr,B_ptr);
       }
       // B_type==ROW_VECTOR
       if (verbose>BASIC) {
          fprintf(stderr,"in %s, A and B are both row vectors\n",__func__);
          fprintf(stderr,"returning A*B^{dagger}\n");
       }
       mat_ptr_t conj_ptr = matrix_adjoint(B_ptr);
       return matrix_mult(A_ptr,conj_ptr);
    }

    // A_type == COL_VECTOR
    if (B_type==ROW_VECTOR) {
       if (verbose>BASIC) {
             fprintf(stderr,"in %s, A is a column vector and\n",__func__);
             fprintf(stderr,"B is a row vector - returning B*A\n");
       }
          return matrix_mult(B_ptr,A_ptr);
    }
    // B_type==COL_VECTOR
    if (verbose>BASIC) {
       fprintf(stderr,"in %s, A and B are both column vectors\n",__func__);
       fprintf(stderr,"returning A^{dagger}*B\n");
    }
    mat_ptr_t conj_ptr = matrix_adjoint(A_ptr);
    return matrix_mult(conj_ptr,B_ptr);
}
