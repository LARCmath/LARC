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

#include <complex.h>

#include "matmath.h"
#include "larc.h"
#include "scalars.h"
#include "global.h"

/*!
 * \file matmath.c
 * \brief This file contains most of the linear algebraic routines of
 * LARC, along with other non-linear routines and routines that manipulate
 * matrices in other ways.
 */

int64_t matrix_add(int64_t A_pID, int64_t B_pID)
{

  // The following checks to ensure that both A_pID and B_ID 
  // have valid record pointers (eg, not RECORD_PTR_INVALID)
  record_ptr_t A_Rptr, B_Rptr;
  check_validity_two_input(A_pID, B_pID,__func__, &A_Rptr, &B_Rptr);

  if (IS_SCALAR(A_pID) && IS_SCALAR(B_pID))
  {
    mats_ptr_t A_ptr = (mats_ptr_t)A_Rptr;
    mats_ptr_t B_ptr = (mats_ptr_t)B_Rptr;
    if (scratchVars.submit_to_store_in_use)
        fprintf(stderr,"%s reusing scratchVars.submit_to_store!\n",__func__);
    scratchVars.submit_to_store_in_use = 1;

    scalarType *sum = &scratchVars.submit_to_store;
    sca_add(sum, A_ptr->scalar_value, B_ptr->scalar_value);
    
    // FIND OR CREATE APPROPRIATE MATRIX
    int64_t ret_pID =  get_scalarPTR_for_scalarVal(*sum)->packedID;
    scratchVars.submit_to_store_in_use = 0;
    return ret_pID;
  }

  if ( IS_SCALAR(A_pID) || IS_SCALAR(B_pID) ) // one SCALAR, one not
  {
    fprintf(stderr,"ERROR in %s: inputs different dimensions!\n",__func__);
    fprintf(stderr,"matrixIDs: for A = %" PRId64 ", for B = %" PRId64 "\n",
         MID_FROM_PID(A_pID), MID_FROM_PID(B_pID));
    fprintf(stderr,"one is a scalar, the other is not.\n");
    exit(EXIT_FAILURE);
  }

  // A and B are both not SCALARs
  matns_ptr_t A_ptr = (matns_ptr_t)A_Rptr;
  matns_ptr_t B_ptr = (matns_ptr_t)B_Rptr;

  // VERIFY ADDITION IS PERMITTED
  if (A_ptr->row_level != B_ptr->row_level)
  {
    fprintf(stderr,"ERROR: attempting to add matrices of different row_level (%d and %d)\n",
           A_ptr->row_level, B_ptr->row_level);
    fprintf(stderr,"additional information:\n");
    fprintf(stderr,"matrixIDs: for A = %" PRId64 ", for B = %" PRId64 "\n",
        MID_FROM_PID(A_pID), MID_FROM_PID(B_pID));
    fprintf(stderr,"pointers: for A = %p, for B = %p\n", A_ptr,B_ptr);
    fprintf(stderr,"column level: for A = %d, for B = %d\n",
           A_ptr->col_level, B_ptr->col_level);
    fprintf(stderr,"scalar or trace value: for A = %s, for B = %s\n",
           sca_get_readable_approx_str(A_ptr->trace_element),
           sca_get_readable_approx_str(B_ptr->trace_element));
    exit(EXIT_FAILURE);
  }

  if (A_ptr->col_level != B_ptr->col_level)
  {
    fprintf(stderr,"ERROR: attempting to add matrices of different col_level (%d and %d)\n",
           A_ptr->col_level, B_ptr->col_level);
    fprintf(stderr,"additional information:\n");
    fprintf(stderr,"matrixIDs: for A = %" PRId64 ", for B = %" PRId64 "\n",
           MID_FROM_PID(A_pID),MID_FROM_PID(B_pID));
    fprintf(stderr,"pointers: for A = %p, for B = %p\n", A_ptr,B_ptr);
    fprintf(stderr,"row level: for A = %d, for B = %d\n",
           A_ptr->row_level, B_ptr->row_level);
    fprintf(stderr,"scalar or trace value: for A = %s, for B = %s\n",
           sca_get_readable_approx_str(A_ptr->trace_element),
           sca_get_readable_approx_str(B_ptr->trace_element));
    exit(EXIT_FAILURE);
  }

  // This part sorts A and B locally so that A < B and you take
  // advantage of commutativity
  if (A_pID > B_pID) {
    int64_t temp = A_pID;
    A_pID = B_pID;
    B_pID = temp;
    matns_ptr_t tempPTR = A_ptr;
    A_ptr = B_ptr;
    B_ptr = tempPTR;
  }

  // MATH RELATIONSHIP / IDENTITY SHORT CUTS
  if (A_ptr->iszero) return B_pID;
  if (B_ptr->iszero) return A_pID;

  // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE
  uint64_t hash = hash_from_op(A_pID, B_pID, SUM);
  int64_t sum_pID = op_get(SUM, A_pID, B_pID, hash);
  if (sum_pID != MATRIX_ID_INVALID) return sum_pID;
  
  int64_t sub_pID[4];
    
  for (int i = 0; i < 4; i++)
  {
    int64_t Ai_pID = A_ptr->subMatList[i];
    int64_t Bi_pID = B_ptr->subMatList[i];
    if (Ai_pID == MATRIX_ID_INVALID || Bi_pID == MATRIX_ID_INVALID)
      sub_pID[i] = MATRIX_ID_INVALID;
    else
      sub_pID[i] = matrix_add(Ai_pID, Bi_pID);
  }
    
  // FIND OR CREATE APPROPRIATE MATRIX PTR
  sum_pID = get_pID_from_array_of_four_sub_pIDs(sub_pID, 
            A_ptr->row_level, A_ptr->col_level);

  // STORE RESULT IN OPERATIONS STORE 
//  printf("calling op_set from matrix_add\n");
  op_set(SUM, A_pID, B_pID, sum_pID, hash);

  return sum_pID;
}

int64_t matrix_diff(int64_t A_pID, int64_t B_pID)
{

  // The following checks to ensure that both A_pID and B_ID 
  // have valid record pointers (eg, not RECORD_PTR_INVALID)
  record_ptr_t A_Rptr, B_Rptr;
  check_validity_two_input(A_pID, B_pID,__func__, &A_Rptr, &B_Rptr);

  if ( IS_SCALAR(A_pID) && IS_SCALAR(B_pID) ) // both SCALAR
  {
     mats_ptr_t A_ptr = (mats_ptr_t)A_Rptr;
     mats_ptr_t B_ptr = (mats_ptr_t)B_Rptr;
     if (scratchVars.submit_to_store_in_use)
         fprintf(stderr,"%s reusing scratchVars.submit_to_store!\n",__func__);
     scratchVars.submit_to_store_in_use = 1;

     scalarType *diff = &scratchVars.submit_to_store;
     sca_mult(diff, B_ptr->scalar_value, scalarM1);
     sca_add(diff, A_ptr->scalar_value, *diff);
 
     // FIND OR CREATE APPROPRIATE MATRIX PTR
     int64_t ret_pID =  get_scalarPTR_for_scalarVal(*diff)->packedID;
     scratchVars.submit_to_store_in_use = 0;
     return ret_pID;
  }

  if ( IS_SCALAR(A_pID) || IS_SCALAR(B_pID) ) // one SCALAR, one not
  {
    fprintf(stderr,"ERROR in %s: inputs different dimensions!\n",__func__);
    fprintf(stderr,"matrixIDs: for A = %" PRId64 ", for B = %" PRId64 "\n",
         MID_FROM_PID(A_pID), MID_FROM_PID(B_pID));
    fprintf(stderr,"one is a scalar, the other is not.\n");
    exit(EXIT_FAILURE);
  }

  // both matrices non-SCALAR
  matns_ptr_t A_ptr = (matns_ptr_t)A_Rptr;
  matns_ptr_t B_ptr = (matns_ptr_t)B_Rptr;

  // EXCEPTION CHECKING
  if (A_ptr->row_level != B_ptr->row_level) {
    fprintf(stderr,"ERROR: attempting to subtract matrices of different row_level (%d and %d)\n",
           A_ptr->row_level, B_ptr->row_level);
    fprintf(stderr,"matrixIDs: for A = %" PRId64 ", for B = %" PRId64 "\n",
        MID_FROM_PID(A_ptr->packedID), MID_FROM_PID(B_ptr->packedID));
    exit(EXIT_FAILURE);
  }
  if (A_ptr->col_level != B_ptr->col_level) {
    fprintf(stderr,"ERROR: attempting to subtract matrices of different col_level (%d and %d)\n",
           A_ptr->col_level, B_ptr->col_level);
    fprintf(stderr,"matrixIDs: for A = %" PRId64 ", for B = %" PRId64 "\n",
        MID_FROM_PID(A_ptr->packedID), MID_FROM_PID(B_ptr->packedID));
    exit(EXIT_FAILURE);
  }

  // MATH IDENTITIES: subtract zero matrix
  if (B_ptr->iszero) return A_pID;

  // MATH IDENTITIES: subtract from zero matrix
  if (A_ptr->iszero) return scalar_mult(packedID_scalarM1, B_pID);

  // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  uint64_t hash = hash_from_op(A_pID, B_pID, DIFF);
  int64_t diff_pID = op_get(DIFF, A_pID, B_pID, hash);
  if (diff_pID != MATRIX_ID_INVALID) return diff_pID;

  int64_t sub_pID[4];
    
  for (int i = 0; i < 4; i++) 
  {
      int64_t Ai_pID = A_ptr->subMatList[i];
      int64_t Bi_pID = B_ptr->subMatList[i];
      if (Ai_pID == MATRIX_ID_INVALID || Bi_pID == MATRIX_ID_INVALID)
        sub_pID[i] = MATRIX_ID_INVALID;
      else
        sub_pID[i] = matrix_diff(Ai_pID, Bi_pID);
  }
    
  // FIND OR CREATE APPROPRIATE MATRIX
  diff_pID = get_pID_from_array_of_four_sub_pIDs(sub_pID, 
             A_ptr->row_level, A_ptr->col_level);

  // STORE RESULT IN OPERATIONS STORE 
//    printf("calling op_set from matrix_diff_PTR\n");
  op_set(DIFF, A_pID, B_pID, diff_pID, hash);

  return diff_pID;
}

int64_t matrix_mult(int64_t A_pID, int64_t B_pID)
{

  // The following checks to ensure that both A_pID and B_ID 
  // have valid record pointers (eg, not RECORD_PTR_INVALID)
  record_ptr_t A_Rptr, B_Rptr;
  check_validity_two_input(A_pID, B_pID,__func__, &A_Rptr, &B_Rptr);

  // If A and/or B is SCALAR, then problem reduces to scalar_mult
  if (IS_SCALAR(A_pID)) return scalar_mult(A_pID, B_pID);
  if (IS_SCALAR(B_pID)) return scalar_mult(B_pID, A_pID);

  // both A and B are non-SCALAR
  matns_ptr_t A_ptr = (matns_ptr_t)A_Rptr;
  matns_ptr_t B_ptr = (matns_ptr_t)B_Rptr;

  //cheap way to be higher than max of levels (and turn off cleaning)
  mat_level_t noCleanThresh = A_ptr->row_level + A_ptr->col_level
    + B_ptr->row_level + B_ptr->col_level + 1;

  return matrix_mult_clean(A_pID, B_pID, noCleanThresh);
}

int64_t matrix_mult_clean(int64_t A_pID, int64_t B_pID,
           mat_level_t cleanThresh)
{
  // This recursive routine cleans the matrix store if the maximum level of the
  // input matrices (that is, the largest of row_level(A), col_level(B), and
  // col_level(A)==row_level(B)) is at or larger than the cleanThresh value. It
  // performs a full cleaning and also deletes the OperationsStore when the
  // maximum level is equal to cleanThresh; when it is greater than
  // cleanThresh, to save time only the components that are summed to get each
  // quadrant submatrix are cleaned.
  //
  // This routine skips validity checking which is automatically done when
  // using the Python-callable versions of the code (that take matrixIDs as
  // inputs rather than matrix pointers). As such, the routine has been written 
  // to ensure that any data needed for the multiply is not unintentionally 
  // cleaned before use. (In this case, the important matrices are the inputs
  // to the routine, those which form the quadrant submatrices of the result,
  // and the two submatrix products which are added to form each quadrant
  // submatrix of the result.) The calls to set_hold_matrix and
  // release_hold_matrix protect each matrix from cleaning until all
  // computation which requires it has been completed. (As always, no matrix
  // is cleaned if it has a hold or lock set, or if it is a quadrant submatrix
  // of a larger matrix currently in the matrixStore.)

  // The following checks to ensure that both A_pID and B_ID 
  // have valid record pointers (eg, not RECORD_PTR_INVALID)
  record_ptr_t A_Rptr, B_Rptr;
  check_validity_two_input(A_pID, B_pID,__func__, &A_Rptr, &B_Rptr);

  int verbose = VERBOSE;

  // If A and/or B is SCALAR, then problem reduces to scalar_mult
  if (IS_SCALAR(A_pID)) {
    if (verbose==DEBUG) {
      // DEBUGGING
      printf("In %s calling scalar mult A scalar\n", __func__);
    }
    return scalar_mult(A_pID, B_pID);
  }

  if (IS_SCALAR(B_pID)) {
    if (verbose==DEBUG) {
      // DEBUGGING
      printf("In %s calling scalar mult B scalar\n", __func__);
    }
    return scalar_mult(B_pID, A_pID);
  }

  // both A and B are non-scalar
  // get the matrix pointers from the matrixIDs
  matns_ptr_t A_ptr = (matns_ptr_t)A_Rptr;
  matns_ptr_t B_ptr = (matns_ptr_t)B_Rptr;

  // Confirm that matrices have proper dimension for multiply
  if (A_ptr->col_level != B_ptr->row_level) {
    fprintf(stderr,"ERROR: attempting to multiply matrix with col_level %d by matrix with row_level %d\n", 
          A_ptr->col_level, B_ptr->row_level);
    fprintf(stderr,"matrixIDs: for A = %" PRId64 ", for B = %" PRId64 "\n", 
          MID_FROM_PID(A_pID), MID_FROM_PID(B_pID));

    if (verbose>BASIC) {    
      matrix_store_report("stdout"); 
      op_store_report("stdout"); 
      memory_and_time_report(0,"stdout");
    }
    // if you want a hash report add
    // hash_report(store.hash_table, stdout, 0);
    // or hash_report(hash_table_t *table, FILE *fp, int verbose)
    exit(EXIT_FAILURE);
  }

  matrix_type_t mat_type_A = matrix_type(A_ptr);
  matrix_type_t mat_type_B = matrix_type(B_ptr); 

  if (verbose==DEBUG) { // DEBUGGING
    printf("\n Matrix types for A (%" PRId64 ") and B (%" PRId64 ") are %d %d\n", A_pID, B_pID, mat_type_A, mat_type_B);
    printf("   Levels for A are %d %d\n",A_ptr->row_level, A_ptr->col_level);
    printf("   Levels for B are %d %d\n",B_ptr->row_level, B_ptr->col_level);
  }

  // This replaces routine matrix_pair_max_level, called nowhere else
  mat_level_t max_level = A_ptr->row_level;
  max_level = MAX(max_level,A_ptr->col_level);
  max_level = MAX(max_level,B_ptr->col_level);

  int reporting = 0; // updates user on periodic cleaning and progress
  mat_level_t min_reporting_level = 29; 

  reporting = reporting && (max_level >= min_reporting_level);

  // MATH IDENTITY SHORT CUTS
  if (A_ptr->isid) {
    if (reporting)
      printf("  L%d: shortcut (first matrix is identity).\n", max_level);
    return B_pID;
  }

  if (B_ptr->isid) {
    if (reporting)
      printf("  L%d: shortcut (second matrix is identity).\n", max_level);
    return A_pID;
  }

  if ( A_ptr->iszero || B_ptr->iszero )
  {
    if (reporting)
      printf("  L%d: shortcut (zero matrix).\n", max_level);
    return get_zero_pID(A_ptr->row_level, B_ptr->col_level);
  }

  // If A is COL_VECTOR and B is ROW_VECTOR, then problem reduces
  // to kronecker_product.
  if ((mat_type_A == COL_VECTOR) && (mat_type_B == ROW_VECTOR)) {
    if (verbose==DEBUG) {
      // DEBUGGING
      printf("In %s calling kronecker product A col B row\n", __func__);
    }
    return kronecker_product(A_pID, B_pID);
  }

  // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  uint64_t hash = hash_from_op(A_pID, B_pID, PRODUCT);
  int64_t product_pID = op_get(PRODUCT, A_pID, B_pID, hash);

  if (product_pID != MATRIX_ID_INVALID) return product_pID;

  if (scratchVars.top_level == -1)
    scratchVars.top_level = max_level;

  if (max_level >= cleanThresh){
    //NOTE: holds are actually counters so when we internally holding
    //and releasing here won't change whether a matrix is held outside
    //of this routine. 
    set_hold_matrix(A_pID);
    set_hold_matrix(B_pID);
  }

  int64_t sub_pID[4];

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
      int64_t Ai0 = A_ptr->subMatList[2*i];     
      int64_t Ai1 = A_ptr->subMatList[2*i+1];
      int64_t B0j = B_ptr->subMatList[j];
      int64_t B1j = B_ptr->subMatList[j+2];
      // A and B are held at this point, so submatrices 
      // Ai0, Ai1, B0j, B1j cannot be removed. 

      // this replaces the ROW/COL VECTOR vs MATRIX comparisons
      if ((Ai0 == MATRIX_ID_INVALID) || (Ai1 == MATRIX_ID_INVALID)
            || (B0j == MATRIX_ID_INVALID) || (B1j == MATRIX_ID_INVALID)){
        sub_pID[ind] = MATRIX_ID_INVALID;
        continue;
      }

      int64_t mat_store_size1, mat_store_size2;
      mat_store_size1 = (int64_t) (nonscalar_store_count() + scalar_store_count());

      int64_t Product0 = matrix_mult_clean(Ai0, B0j, cleanThresh);

      mat_store_size2 = (int64_t) (nonscalar_store_count() + scalar_store_count());
      if (reporting)
        printf("  L%d: mult0 adds %" PRId64 " matrix records to store (now total %" PRId64 ").\n", max_level, mat_store_size2-mat_store_size1, mat_store_size2);

      // set hold on Product0, so that the call to matrix_mult_clean
      // that calculates Product1 can't clean it
      set_hold_matrix(Product0);

      int64_t Product1 = matrix_mult_clean(Ai1, B1j, cleanThresh);

      mat_store_size1 = (int64_t) (nonscalar_store_count() + scalar_store_count());
      if (reporting)
        printf("  L%d: mult1 adds %" PRId64 " matrix records to store (now total %" PRId64 ").\n", max_level, mat_store_size1-mat_store_size2, mat_store_size1);

      // release hold on Product0
      release_hold_matrix(Product0);

      // Older versions of this routine had code which handled the case
      // (cleanThresh == 0) seperately, using a routine which did the matrix
      // addition without storing the results in the operations store. This
      // was not carried over into newer versions.
      sub_pID[ind] = matrix_add(Product0, Product1);

      mat_store_size2 = (int64_t) (nonscalar_store_count() + scalar_store_count());
      if (reporting)
        printf("  L%d: add adds %" PRId64 " matrix records to store (now total %" PRId64 ").\n", max_level, mat_store_size2-mat_store_size1, mat_store_size2);

      if (max_level == cleanThresh){
        set_hold_matrix(sub_pID[ind]); 

        clean_matrix_storage();  
        empty_op_store();

        mat_store_size1 = (int64_t) (nonscalar_store_count() + scalar_store_count());
        if (reporting)
          printf("  L%d: clean(g) removes %" PRId64 " matrix records from store (now total %" PRId64 ").\n", max_level, mat_store_size2-mat_store_size1, mat_store_size1);
      }
      else if (max_level > cleanThresh){
        set_hold_matrix(sub_pID[ind]); 

        // since we just (recursively) cleaned at a lower level, all matrices
        // involved in  the calculations of P0 and P1 were cleaned EXCEPT P0
        // and P1. So we'll do a targeted removal of P0 and P1 instead of a
        // general clean - to save time/effort.
        //
        // NOTE: this handles the rare case where Product0==Product1
        if (Product0 != Product1)
               remove_matrix_from_store(Product0);
        remove_matrix_from_store(Product1);

        mat_store_size1 = (int64_t) (nonscalar_store_count() + scalar_store_count());
        if (reporting)
          printf("  L%d: clean(t) removes %" PRId64 " matrix records from store (now total %" PRId64 ").\n", max_level, mat_store_size2-mat_store_size1, mat_store_size1);

        // the op_store is robust in that it won't trip over missing matrices
        // so I think we should skip repairing/cleaning the op_store.
      }
    } // loop over j
  } // loop over i

  // DOT PRODUCT: A is ROW_VECTOR and B is COL_VECTOR and result is SCALAR
  if ((mat_type_A == ROW_VECTOR) && (mat_type_B == COL_VECTOR)) {
    product_pID = sub_pID[0];
  }
  else {
    // FIND OR CREATE APPROPRIATE MATRIX PTR
    product_pID = get_pID_from_array_of_four_sub_pIDs(sub_pID,
                  A_ptr->row_level, B_ptr->col_level);
  }

  // We missed a case????????
  if (product_pID == MATRIX_ID_INVALID) {
    fprintf(stderr,"ERROR: something went very wrong in %s\n", __func__);
    exit(1);
  }

  // STORE RESULT IN OPERATIONS STORE 
  if ((cleanThresh > 0) || (max_level == scratchVars.top_level))
  {
//        printf("calling op_set from matrix_mult_clean\n");
      op_set(PRODUCT, A_pID, B_pID, product_pID, hash);
  }

  // RELEASE ALL HOLDS MADE LOCALLY
  if (max_level >= cleanThresh)
  {
    release_hold_matrix(A_pID);
    release_hold_matrix(B_pID);
    for (int i = 0; i < 4; i++)
      if (sub_pID[i] != MATRIX_ID_INVALID)
        release_hold_matrix(sub_pID[i]);
  }

  scratchVars_exitroutine(max_level);

  return product_pID;
}

int64_t scalar_mult(int64_t A_pID, int64_t B_pID)
{

  // The following checks to ensure that both A_pID and B_ID 
  // have valid record pointers (eg, not RECORD_PTR_INVALID)
  record_ptr_t A_Rptr, B_Rptr;
  check_validity_two_input(A_pID, B_pID,__func__, &A_Rptr, &B_Rptr);

  // EXCEPTION CHECKING
  if (IS_SCALAR(A_pID) == 0)
  {
    fprintf(stderr,"ERROR: attempted scalar multiply with nonscalar\n");
    fprintf(stderr,"matrixIDs: for A = %" PRId64 ", for B = %" PRId64 "\n",
              MID_FROM_PID(A_pID), MID_FROM_PID(B_pID));
    exit(EXIT_FAILURE);
  }

  // NO MEMOIZATION OF SCALAR-SCALAR MULTIPLY (faster to recompute value)
  if (IS_SCALAR(B_pID)) // SCALAR
  {
    if (A_pID == packedID_scalar0) return A_pID;
    if (B_pID == packedID_scalar0) return B_pID;
    mats_ptr_t A_ptr = (mats_ptr_t)A_Rptr;
    mats_ptr_t B_ptr = (mats_ptr_t)B_Rptr;
    if (scratchVars.submit_to_store_in_use)
        fprintf(stderr,"%s reusing scratchVars.submit_to_store!\n",__func__);
    scratchVars.submit_to_store_in_use = 1;

    scalarType *product = &scratchVars.submit_to_store;
    sca_mult(product, A_ptr->scalar_value, B_ptr->scalar_value);

    // FIND OR CREATE APPROPRIATE MATRIX PTR
    int64_t ret_pID = get_scalarPTR_for_scalarVal(*product)->packedID;
    scratchVars.submit_to_store_in_use = 0;
    return ret_pID;
  }

  // get the matrix pointer from the matrixID
  matns_ptr_t B_ptr = (matns_ptr_t)B_Rptr;

  // MATH IDENTITY SHORT CUTS
  if (A_pID == packedID_scalar0 || B_ptr->iszero)
    return get_zero_pID(B_ptr->row_level, B_ptr->col_level);
  
  // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  uint64_t hash = hash_from_op(A_pID, B_pID, KRONECKER);
  int64_t product_pID = op_get(KRONECKER, A_pID,B_pID, hash);
  if (product_pID != MATRIX_ID_INVALID) return product_pID;

  // CALCULATE RESULT OF OPERATION when B is a matrix
  int64_t sub_pID[4];

  for (int i = 0; i < 4; i++) {
      int64_t Bi_pID = B_ptr->subMatList[i];
      if (Bi_pID == MATRIX_ID_INVALID)
        sub_pID[i] = MATRIX_ID_INVALID;
      else
        sub_pID[i] = scalar_mult(A_pID, Bi_pID);
  }

  //   FIND OR CREATE APPROPRIATE MATRIX ID
  product_pID = get_pID_from_array_of_four_sub_pIDs(sub_pID,
                B_ptr->row_level, B_ptr->col_level);
  
  //   STORE RESULT IN OPERATIONS STORE 
//    printf("calling op_set from scalar_mult\n");
  op_set(KRONECKER, A_pID, B_pID, product_pID, hash);
  
  return product_pID;
}

int64_t scalar_divide(int64_t A_pID, int64_t B_pID)
{

  // The following checks to ensure that both A_pID and B_ID 
  // have valid record pointers (eg, not RECORD_PTR_INVALID)
  record_ptr_t A_Rptr, B_Rptr;
  check_validity_two_input(A_pID, B_pID,__func__, &A_Rptr, &B_Rptr);

  if (IS_SCALAR(B_pID) == 0)
  {
    fprintf(stderr,"ERROR: attempted scalar divide with nonscalar!\n");
    fprintf(stderr,"matrixIDs: for A = %" PRId64 ", for B = %" PRId64 "\n",
              MID_FROM_PID(A_pID), MID_FROM_PID(B_pID));
    exit(EXIT_FAILURE);
  }

  // EXCEPTION CHECKING
  
  // CAN'T DIVIDE BY ZERO
  if (B_pID==packedID_scalar0) {
    fprintf(stderr,"ERROR: attempted scalar divide by zero\n");
    fprintf(stderr,"matrixIDs: for A = %" PRId64 ", for B = %" PRId64 "\n",
              MID_FROM_PID(A_pID), MID_FROM_PID(B_pID));
    exit(EXIT_FAILURE);
  }

  // MATH IDENTITY SHORT CUTS
  if (matrix_is_zero(A_pID)) return A_pID;
  
  // NO MEMOIZATION OF SCALAR-SCALAR MULTIPLY (faster to recompute value)
  if (IS_SCALAR(A_pID)) // SCALAR
  {
    mats_ptr_t A_ptr = (mats_ptr_t)A_Rptr;
    mats_ptr_t B_ptr = (mats_ptr_t)B_Rptr;
    if (scratchVars.submit_to_store_in_use)
        fprintf(stderr,"%s reusing scratchVars.submit_to_store!\n",__func__);
    scratchVars.submit_to_store_in_use = 1;

    scalarType *quotient = &scratchVars.submit_to_store;
    sca_divide(quotient, A_ptr->scalar_value, B_ptr->scalar_value);

    // FIND OR CREATE APPROPRIATE MATRIX PTR
    int64_t ret_pID = get_scalarPTR_for_scalarVal(*quotient)->packedID;
    scratchVars.submit_to_store_in_use = 0;
    return ret_pID;
  }

  // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  uint64_t hash = hash_from_op(A_pID,B_pID,QUOTIENT_SCALAR);
  int64_t quotient_pID = op_get(QUOTIENT_SCALAR, A_pID, B_pID, hash);
  if (quotient_pID != MATRIX_ID_INVALID) return quotient_pID;

  // get the matrix pointer from the matrixID
  matns_ptr_t A_ptr = (matns_ptr_t)A_Rptr;

  // CALCULATE RESULT OF OPERATION FOR NON_SCALAR A MATRICES
  int64_t sub_pID[4];

  for (int i = 0; i < 4; i++) {
      int64_t Ai_pID = A_ptr->subMatList[i];
      if (Ai_pID == MATRIX_ID_INVALID)
        sub_pID[i] = MATRIX_ID_INVALID;
      else
        sub_pID[i] = scalar_divide(Ai_pID, B_pID);
  }

  //   FIND OR CREATE APPROPRIATE MATRIX ID
  quotient_pID = get_pID_from_array_of_four_sub_pIDs(sub_pID,
                 A_ptr->row_level, A_ptr->col_level);
  
  //   STORE RESULT IN OPERATIONS STORE 
//    printf("calling op_set from scalar_divide\n");
  op_set(QUOTIENT_SCALAR, A_pID, B_pID, quotient_pID, hash);

  return quotient_pID;
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
static int64_t
tensor_with_identity_on_left(int64_t A_pID, mat_level_t levelI)
{

  // MATH IDENTITY SHORT CUTS
  //   If Identity matrix I is SCALAR, then return A
  if (levelI == 0) { return A_pID; }

  // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  int64_t identity_pID = get_identity_pID(levelI);
  uint64_t hash = hash_from_op(identity_pID,A_pID,KRONECKER);
  int64_t product_pID = op_get(KRONECKER, identity_pID, A_pID, hash);

  if (product_pID != MATRIX_ID_INVALID) return product_pID;

  matns_ptr_t A_ptr = (matns_ptr_t)get_recordPTR_from_pID(A_pID,"",__func__,0);

  // CALCULATE RESULT OF OPERATION 
    
  // ALL CASES:
  //    calculate Kronecker product recursively
  mat_level_t row_levelA = A_ptr->row_level;
  mat_level_t col_levelA = A_ptr->col_level;
  
  int64_t sub_pID[4] = {A_pID, get_zero_pID(row_levelA, col_levelA),
                          get_zero_pID(row_levelA, col_levelA), A_pID};
  
  product_pID = tensor_with_identity_on_left(
    get_pID_from_array_of_four_sub_pIDs(sub_pID,row_levelA+1,col_levelA+1),
    levelI-1);

  // STORE RESULT IN OPERATIONS STORE 
//  printf("calling op_set from tensor_with_identity_from_left\n");
  op_set(KRONECKER, identity_pID, A_ptr->packedID, product_pID, hash);

  return product_pID;
}

int64_t kronecker_product(int64_t A_pID, int64_t B_pID)
{

  // The following checks to ensure that both A_pID and B_ID 
  // have valid record pointers (eg, not RECORD_PTR_INVALID)
  record_ptr_t A_Rptr, B_Rptr;
  check_validity_two_input(A_pID, B_pID,__func__, &A_Rptr, &B_Rptr);

  // When either A or B is scalar, then the kronecker product is the scalar
  // product (the result will be stored in the op_store KRONECKER unless both
  // A and B are scalars)
  if (IS_SCALAR(A_pID)) return scalar_mult(A_pID,B_pID);
  if (IS_SCALAR(B_pID)) return scalar_mult(B_pID,A_pID);

  // get the matrix pointers from the matrixIDs
  matns_ptr_t A_ptr = (matns_ptr_t)A_Rptr;
  matns_ptr_t B_ptr = (matns_ptr_t)B_Rptr;

  // MATH IDENTITY SHORT CUTS
  if (A_ptr->iszero || B_ptr->iszero)
  {
    return get_zero_pID(
                A_ptr->row_level + B_ptr->row_level,
                A_ptr->col_level + B_ptr->col_level);
  }
  // Use special function in the case when matrix on left is identity matrix
  //  (This function puts the result in the KRONECKER op store)
  if (A_ptr->isid)
    return tensor_with_identity_on_left(B_pID, A_ptr->row_level);    


  matrix_type_t mat_type_A = matrix_type(A_ptr);
  matrix_type_t mat_type_B = matrix_type(B_ptr);


  //   CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  uint64_t hash = hash_from_op(A_pID, B_pID, KRONECKER);
  int64_t product_pID = op_get(KRONECKER, A_pID, B_pID, hash);

  if (product_pID != MATRIX_ID_INVALID) return product_pID;

  // CALCULATE RESULT OF OPERATION 
//  matns_ptr_t panel[4];
  int64_t sub_pID[4];
  mat_level_t new_row_level = A_ptr->row_level + B_ptr->row_level;
  mat_level_t new_col_level = A_ptr->col_level + B_ptr->col_level;
    
  // CASE: A and B are both ROW_VECTOR 
  // former version had join() outside of kronecker product
  if ((mat_type_A == ROW_VECTOR) && (mat_type_B == ROW_VECTOR)) {
    sub_pID[0] = kronecker_product(A_ptr->subMatList[0],B_pID);
    sub_pID[1] = kronecker_product(A_ptr->subMatList[1],B_pID);
    sub_pID[2] = MATRIX_ID_INVALID;
    sub_pID[3] = MATRIX_ID_INVALID;
  }
    
  // CASE: A and B are both COL_VECTOR
  // former version did kronecker product with submatrices of B rather than
  // the entire B matrix - fixed
  if ((mat_type_A == COL_VECTOR) && (mat_type_B == COL_VECTOR)) {
    sub_pID[0] = kronecker_product(A_ptr->subMatList[0], B_pID);
    sub_pID[1] = MATRIX_ID_INVALID;
    sub_pID[2] = kronecker_product(A_ptr->subMatList[2], B_pID);
    sub_pID[3] = MATRIX_ID_INVALID;
  }
  
  // CASE: A is COL_VECTOR and B is ROW_VECTOR (correct)
  // [a_0,a_2]^T \otimes [b_0,b_1] -> 
  // [[a_0b_0,a_0b_1],[a_2b_0,a_2b_1]
  if ((mat_type_A == COL_VECTOR) && (mat_type_B == ROW_VECTOR)) {
    sub_pID[0] = kronecker_product(A_ptr->subMatList[0],
	B_ptr->subMatList[0]);
    sub_pID[1] = kronecker_product(A_ptr->subMatList[0],
	B_ptr->subMatList[1]);
    sub_pID[2] = kronecker_product(A_ptr->subMatList[2],
	B_ptr->subMatList[0]);
    sub_pID[3] = kronecker_product(A_ptr->subMatList[2],
	B_ptr->subMatList[1]);
  }
  
  // CASE: A is ROW_VECTOR and B is COL_VECTOR (correct)
  // [a_0,a_1] \otimes [b_0,b_2]^T -> 
  // [[a_0b_0,a_1b_0],[a_0b_2,a_1b_2]
  if ((mat_type_A == ROW_VECTOR) && (mat_type_B == COL_VECTOR)) {
    sub_pID[0] = kronecker_product(A_ptr->subMatList[0],
        B_ptr->subMatList[0]);
    sub_pID[1] = kronecker_product(A_ptr->subMatList[1],
        B_ptr->subMatList[0]);
    sub_pID[2] = kronecker_product(A_ptr->subMatList[0],
        B_ptr->subMatList[2]);
    sub_pID[3] = kronecker_product(A_ptr->subMatList[1],
        B_ptr->subMatList[2]);
  }
    
  // CASE: A is ROW_VECTOR and B is MATRIX (fix by Matt Calef)
  // former version had join() outside of kronecker product
  // [a_0,a_1] \otimes [[b_0,b_1],[b_2,b_3]] ->
  // [[a_0 join(b_0,b_1),a_1 join(b_0,b_1)],
  //    [a_0 join(b_2,b_3),a_1 join(b_2,b_3)]
  if ((mat_type_A == ROW_VECTOR) && (mat_type_B == MATRIX)) {
     int64_t temp_pID =join(B_ptr->subMatList[0],B_ptr->subMatList[1]);
     sub_pID[0] = kronecker_product(A_ptr->subMatList[0],temp_pID);
     sub_pID[1] = kronecker_product(A_ptr->subMatList[1],temp_pID);
     temp_pID = join(B_ptr->subMatList[2],B_ptr->subMatList[3]);
     sub_pID[2] = kronecker_product(A_ptr->subMatList[0],temp_pID);
     sub_pID[3] = kronecker_product(A_ptr->subMatList[1],temp_pID);
  }
  
  // CASE: A is COL_VECTOR and B is MATRIX
  // [a_0,a_2]^T \otimes B -> [a_0 B, a_2 B]^T
  // -> [a_0 stack(b_0,b_2),a_0 stack(b_1,b_3)],
  //          [a_2 stack(b_0,b_2),a_2 stack(b_1,b_3)]]
  // former version had stack() outside of kronecker product
  if ((mat_type_A == COL_VECTOR) && (mat_type_B == MATRIX)) {
     int64_t temp_pID=stack(B_ptr->subMatList[0],B_ptr->subMatList[2]);
     sub_pID[0] = kronecker_product(A_ptr->subMatList[0],temp_pID);
     sub_pID[2] = kronecker_product(A_ptr->subMatList[2],temp_pID);
     temp_pID = stack(B_ptr->subMatList[1],B_ptr->subMatList[3]);
     sub_pID[1] = kronecker_product(A_ptr->subMatList[0],temp_pID);
     sub_pID[3] = kronecker_product(A_ptr->subMatList[2],temp_pID);
  }
  
  // CASE: A is MATRIX - not necessary to break down B 
  // (simplification by Matt Calef)
  if (mat_type_A == MATRIX) {
    for (int i = 0; i < 4 ; i++) {
       sub_pID[i] = kronecker_product(A_ptr->subMatList[i], B_pID);
    }
  }

  // find or create the matrix id of the result
  product_pID = get_pID_from_array_of_four_sub_pIDs(sub_pID,
                new_row_level, new_col_level);

  // ERROR CHECKING (since there are so many cases)
  if (product_pID == MATRIX_ID_INVALID) {
      fprintf(stderr,"ERROR: in %s something went bad\n", __func__);
      exit(1);
  }
    
  //   STORE RESULT IN OPERATIONS STORE 
//    printf("calling op_set from kronecker_product\n");
  op_set(KRONECKER, A_pID, B_pID, product_pID, hash);

  return product_pID;
}

/********************************************************************
 *                         join()                          *
 *  Returns the matrix in which two identical sized input matrices  *
 *  have been appended side to side                                 *
 *  Uses the operation store JOIN                                   *
 *******************************************************************/

int64_t join(int64_t A_pID, int64_t B_pID)
{

  // The following checks to ensure that both A_pID and B_ID 
  // have valid record pointers (eg, not RECORD_PTR_INVALID)
  record_ptr_t A_Rptr, B_Rptr;
  check_validity_two_input(A_pID, B_pID,__func__, &A_Rptr, &B_Rptr);

  if (IS_SCALAR(A_pID) && IS_SCALAR(B_pID)) // both SCALARS
  {
    uint64_t hash = hash_from_op(A_pID, B_pID, JOIN);
    int64_t out_pID = op_get(JOIN, A_pID, B_pID, hash);
    if (out_pID != MATRIX_ID_INVALID) return out_pID;

    mat_level_t row_level = 0;
    mat_level_t new_col_level = 1;
    out_pID = get_pID_from_four_sub_pIDs(A_pID, B_pID, MATRIX_ID_INVALID,
       MATRIX_ID_INVALID, row_level, new_col_level);
    op_set(JOIN, A_pID, B_pID, out_pID, hash);
    return out_pID;
  }

  if (IS_SCALAR(A_pID) || IS_SCALAR(B_pID)) // one SCALAR, one not
  {
    fprintf(stderr,"In function %s, input matrices were not the same size.\n",__func__);
    exit(1);
  }

  // both not SCALARS
  uint64_t hash = hash_from_op(A_pID, B_pID, JOIN);
  int64_t out_pID = op_get(JOIN, A_pID, B_pID, hash);
  if (out_pID != MATRIX_ID_INVALID) return out_pID;

  // get the matrix pointers from the matrixIDs
  matns_ptr_t A_ptr = (matns_ptr_t)A_Rptr;
  matns_ptr_t B_ptr = (matns_ptr_t)B_Rptr;

  if ( (A_ptr->row_level != B_ptr->row_level) ||
       (A_ptr->col_level != B_ptr->col_level) ) {
    fprintf(stderr,"In function %s, input matrices were not the same size.\n",__func__);
    exit(1);
  }

  matrix_type_t mat_type_A = matrix_type(A_ptr);

  //   CALCULATE RESULT OF OPERATION
  int64_t sub_pID[4];
  mat_level_t row_level = A_ptr->row_level;
  mat_level_t new_col_level = A_ptr->col_level + 1;
  
  // SCALAR CASE already handled above

  // ROW_VECTOR CASE
  // note: with A a row vector, the result of
  //          join(A_ptr->subMatList[0],A_ptr->subMatList[1]);
  // is A_pID. This case should be caught by memoization, but
  // we could make it explicit
  if (mat_type_A == ROW_VECTOR) {
    sub_pID[0] = join(A_ptr->subMatList[0],A_ptr->subMatList[1]);
    sub_pID[1] = join(B_ptr->subMatList[0],B_ptr->subMatList[1]);
    sub_pID[2] = MATRIX_ID_INVALID;
    sub_pID[3] = MATRIX_ID_INVALID;
  }
  
  // COL_VECTOR CASE
  else if (mat_type_A == COL_VECTOR) {
    sub_pID[0] = A_ptr->subMatList[0];
    sub_pID[1] = B_ptr->subMatList[0];
    sub_pID[2] = A_ptr->subMatList[2];
    sub_pID[3] = B_ptr->subMatList[2];
  }
  
  // MATRIX CASE
  else {
    sub_pID[0] = join(A_ptr->subMatList[0],A_ptr->subMatList[1]);
    sub_pID[1] = join(B_ptr->subMatList[0],B_ptr->subMatList[1]);
    sub_pID[2] = join(A_ptr->subMatList[2],A_ptr->subMatList[3]);
    sub_pID[3] = join(B_ptr->subMatList[2],B_ptr->subMatList[3]);
  }

//   FIND OR CREATE MATRIX PTR
  out_pID = get_pID_from_array_of_four_sub_pIDs(sub_pID,
            row_level,new_col_level);
  
  //   STORE RESULT IN OPERATIONS STORE 
//    printf("calling op_set from join\n");
  op_set(JOIN, A_pID, B_pID, out_pID, hash);

  return out_pID;
}

/********************************************************************
 *                         stack()                         *
 *  Returns the matrix in which two identical sized input matrices  *
 *  have been stacked with the first on top of the second.          *
 *  Uses the operation store STACK                                  *
 *******************************************************************/

int64_t stack(int64_t A_pID, int64_t B_pID)
{

  // The following checks to ensure that both A_pID and B_ID 
  // have valid record pointers (eg, not RECORD_PTR_INVALID)
  record_ptr_t A_Rptr, B_Rptr;
  check_validity_two_input(A_pID, B_pID,__func__, &A_Rptr, &B_Rptr);

  if (IS_SCALAR(A_pID) && IS_SCALAR(B_pID)) // both SCALARS
  {
    uint64_t hash = hash_from_op(A_pID, B_pID, STACK);
    int64_t out_pID = op_get(STACK, A_pID, B_pID, hash);
    if (out_pID != MATRIX_ID_INVALID) return out_pID;

    mat_level_t new_row_level = 1;
    mat_level_t col_level = 0;
    out_pID = get_pID_from_four_sub_pIDs(A_pID, MATRIX_ID_INVALID, B_pID,
       MATRIX_ID_INVALID, new_row_level, col_level);
    op_set(STACK, A_pID, B_pID, out_pID, hash);
    return out_pID;
  }

  if (IS_SCALAR(A_pID) || IS_SCALAR(B_pID)) // one SCALAR, one not
  {
    fprintf(stderr,"In function %s, input matrices were not the same size.\n",__func__);
    exit(1);
  }

  // both not SCALARS
  uint64_t hash = hash_from_op(A_pID, B_pID, STACK);
  int64_t out_pID = op_get(STACK, A_pID, B_pID, hash);
  if (out_pID != MATRIX_ID_INVALID) return out_pID;

  // get the matrix pointers from the matrixIDs
  matns_ptr_t A_ptr = (matns_ptr_t)A_Rptr;
  matns_ptr_t B_ptr = (matns_ptr_t)B_Rptr;

  if ( (A_ptr->row_level != B_ptr->row_level) ||
       (A_ptr->col_level != B_ptr->col_level) ) {
    fprintf(stderr,"In function %s, input matrices were not the same size.\n",__func__);
    exit(1);
  }

  matrix_type_t mat_type_A = matrix_type(A_ptr);

  //   CALCULATE RESULT OF OPERATION
  int64_t sub_pID[4];
  mat_level_t new_row_level = A_ptr->row_level + 1;
  mat_level_t col_level = A_ptr->col_level;
  
  // SCALAR CASE already handled above
  
  // COL_VECTOR CASE
  // note: with A a column vector, the result of
  //          stack(A_ptr->subMatList[0],A_ptr->subMatList[2]);
  // is A_pID. This case should be caught by memoization, but
  // we could make it explicit
  if (mat_type_A == COL_VECTOR) {
    sub_pID[0] = stack(A_ptr->subMatList[0],A_ptr->subMatList[2]);
    sub_pID[1] = MATRIX_ID_INVALID;
    sub_pID[2] = stack(B_ptr->subMatList[0],B_ptr->subMatList[2]);
    sub_pID[3] = MATRIX_ID_INVALID;
  }
  
  // ROW_VECTOR CASE
  else if (mat_type_A == ROW_VECTOR) {
    sub_pID[0] = A_ptr->subMatList[0];
    sub_pID[1] = A_ptr->subMatList[1];
    sub_pID[2] = B_ptr->subMatList[0];
    sub_pID[3] = B_ptr->subMatList[1];
  }
  
  // MATRIX CASE
  else
  {
    sub_pID[0] = stack(A_ptr->subMatList[0],A_ptr->subMatList[2]);
    sub_pID[1] = stack(A_ptr->subMatList[1],A_ptr->subMatList[3]);
    sub_pID[2] = stack(B_ptr->subMatList[0],B_ptr->subMatList[2]);
    sub_pID[3] = stack(B_ptr->subMatList[1],B_ptr->subMatList[3]);
  }

  // FIND OR CREATE MATRIX ID
  out_pID = get_pID_from_array_of_four_sub_pIDs(sub_pID, 
            new_row_level, col_level);
  
  //   STORE RESULT IN OPERATIONS STORE 
//    printf("calling op_set from stack\n");
  op_set(STACK, A_pID, B_pID, out_pID, hash);

  return out_pID;
}

/***********************************************************************
 *                matrix_entrySquared                         *
 *  Returns the matrix in which each scalar has been squared           *
 *     and then multiplied by the given scale factor                   *
 *  Uses the operation store ENTRYSQUARE                               *
 *  This routine should recursively generate the matrix which has      *
 *  for its (i,j) element the norm squared of the (i,j) element of     *
 *  the input matrix multiplied by a scale_factor.                     *
 ***********************************************************************/

/*!
 * \ingroup larc
 *
 * \brief A recursive utility function for finding the norm squared of every element of a matrix multiplied by some scaling factor
 *
 * \param m_pID The packedID of the input matrix
 * \param scale_ptr A pointer to a scalar matrix holding the value by which we will scale the squared element
 *
 * \return The packedID of the output matrix
 */
static int64_t recursive_matrix_entrySquared(int64_t m_pID, mats_ptr_t scale_ptr)
{

  // MATH IDENTITY SHORT CUTS
  // we know squaring and scaling elements of zero matrix does not change it
  if (matrix_is_zero(m_pID)) return m_pID;

  int64_t scale_pID = scale_ptr->packedID;

  // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  uint64_t hash = hash_from_op(m_pID, scale_pID, ENTRYSQUARE);
  int64_t mSq_pID = op_get(ENTRYSQUARE, m_pID, scale_pID, hash);

  if (mSq_pID != MATRIX_ID_INVALID) return mSq_pID;

  if (IS_SCALAR(m_pID)) // SCALAR
  { 
     mats_ptr_t s_ptr = (mats_ptr_t)get_recordPTR_from_pID(m_pID, "", __func__,0);
     if (scratchVars.submit_to_store_in_use)
         fprintf(stderr,"%s reusing scratchVars.submit_to_store!\n",__func__);
     scratchVars.submit_to_store_in_use = 1;
     if (scratchVars.calc_conj_in_use)
         fprintf(stderr,"%s reusing scratchVars.calc_conj!\n",__func__);
     scratchVars.calc_conj_in_use = 1;

     scalarType *element = &scratchVars.submit_to_store;
     scalarType *conj_elt = &scratchVars.calc_conj;
     sca_conj(conj_elt, s_ptr->scalar_value);
     sca_mult(element, *conj_elt, s_ptr->scalar_value);
     scratchVars.calc_conj_in_use = 0;
     sca_mult(element, *element, scale_ptr->scalar_value);
     mSq_pID = (get_scalarPTR_for_scalarVal(*element))->packedID;
     scratchVars.submit_to_store_in_use = 0;
     op_set(ENTRYSQUARE, m_pID, scale_pID, mSq_pID, hash);
     return mSq_pID;
  }

  // CASE: NON-SCALAR
  matns_ptr_t m_ptr = (matns_ptr_t)get_recordPTR_from_pID(m_pID, "", __func__,0);

  if (m_ptr->isid) return scalar_mult(scale_pID, m_pID);

  int64_t sub_pID[4];
  sub_pID[0] = recursive_matrix_entrySquared(m_ptr->subMatList[0], scale_ptr);
  sub_pID[1] = recursive_matrix_entrySquared(m_ptr->subMatList[1], scale_ptr);
  sub_pID[2] = recursive_matrix_entrySquared(m_ptr->subMatList[2], scale_ptr);
  sub_pID[3] = recursive_matrix_entrySquared(m_ptr->subMatList[3], scale_ptr);
  mSq_pID = get_pID_from_array_of_four_sub_pIDs(sub_pID,
            m_ptr->row_level, m_ptr->col_level);
    
  // STORE RESULT IN OPERATIONS STORE 
//    printf("calling op_set from matrix_entrySquared_PTR\n");
  op_set(ENTRYSQUARE, m_pID, scale_pID, mSq_pID, hash);

  return mSq_pID;
}

int64_t matrix_entrySquared(int64_t m_pID, char *scale_factor)
{

  // The following checks to ensure that m_pID
  // has valid record pointer (eg, not RECORD_PTR_INVALID)
  record_ptr_t m_Rptr;
  check_validity_one_input(m_pID, __func__, &m_Rptr);

  if (IS_SCALAR(m_pID)==0) // non-SCALAR
  {
    matns_ptr_t m_ptr = (matns_ptr_t)m_Rptr;

    //   this routine is only designed to work on square matrices
    //   (so no submatrix will be MATRIX_ID_INVALID)
    if (m_ptr->row_level != m_ptr->col_level)  
    {
      fprintf(stderr,"Function %s requires a square matrix\n", __func__);
      exit(1);
    }
  }

  // convert char* scale_factor into scalarType scale. 
  scalarType scale;
  sca_init(&scale);
  sca_set_str(&scale, scale_factor);

  // put scale_factor into the matrix store and get its pointer *here*
  // so that we don't have to do it an exponential number of times
  mats_ptr_t scale_ptr = get_scalarPTR_for_scalarVal(scale);
  sca_clear(&scale);

  return recursive_matrix_entrySquared(m_pID, scale_ptr);
}

int64_t iHadamard_times_matrix(int64_t A_pID)
{

  // The following checks to ensure that A_pID
  // has valid record pointer (eg, not RECORD_PTR_INVALID)
  record_ptr_t A_Rptr;
  check_validity_one_input(A_pID, __func__, &A_Rptr);

  if (IS_SCALAR(A_pID)) // SCALAR
  {
    fprintf(stderr,"ERROR: in %s, there is no sensible definition for level 0 Hadamard\n",__func__);
    exit(1);
  }

  // find the index of the pre-stored iHadamard matrix of this level
  matns_ptr_t A_ptr = (matns_ptr_t)A_Rptr;

  mat_level_t row_level = A_ptr->row_level;
  mat_level_t col_level = A_ptr->col_level;

  if (row_level == 0) {
    fprintf(stderr,"ERROR: in %s, there is no sensible definition for level 0 Hadamard\n",__func__);
    exit(1);
  }

  int64_t iHad_pID = get_iHadamard_pID(row_level);

  uint64_t hash = hash_from_op(iHad_pID, A_pID, PRODUCT);
  int64_t iHad_times_A_pID = op_get(PRODUCT, iHad_pID, A_pID, hash);

  if (iHad_times_A_pID != MATRIX_ID_INVALID) return iHad_times_A_pID;

  // if no stored value and the matrix level is 1, calculate the product HH1 * A
  // using the preloaded value for the integer iHadamard
  if (row_level == 1) return matrix_mult(iHad_pID, A_pID);

  // if no stored value and the matrix larger than 2x2, calculate the iHadamard
  // product recursively as follows (this should use one less multiply per
  // submatrix than matrix_mult(iHad_pID,A_pID)
  int64_t sub_pID[4];
  sub_pID[0] = iHadamard_times_matrix(matrix_add(
      A_ptr->subMatList[0], A_ptr->subMatList[2]));
  sub_pID[1] = iHadamard_times_matrix(matrix_add(
      A_ptr->subMatList[1], A_ptr->subMatList[3]));
  sub_pID[2] = iHadamard_times_matrix(matrix_add(
      A_ptr->subMatList[0], scalar_mult(
         packedID_scalarM1, A_ptr->subMatList[2])));
  sub_pID[3] = iHadamard_times_matrix(matrix_add(
      A_ptr->subMatList[1], scalar_mult(
         packedID_scalarM1, A_ptr->subMatList[3])));
  iHad_times_A_pID = get_pID_from_array_of_four_sub_pIDs(sub_pID,
                     row_level, col_level);

  // STORE IN OP STORE
//    printf("calling op_set from iHadamard_timex_matrix_PTR\n");
  op_set(PRODUCT, iHad_pID, A_pID, iHad_times_A_pID, hash);

  return iHad_times_A_pID;
}

int64_t matrix_times_iHadamard(int64_t A_pID) 
{

  // The following checks to ensure that A_pID
  // has valid record pointer (eg, not RECORD_PTR_INVALID)
  record_ptr_t A_Rptr;
  check_validity_one_input(A_pID, __func__, &A_Rptr);

  if (IS_SCALAR(A_pID)) // SCALAR
  {
    fprintf(stderr,"ERROR: in %s, there is no sensible definition for level 0 Hadamard\n",__func__);
    exit(1);
  }

// We have decided for now that default iHadmard will be the integer
// (unnormalized) version
// routine takes the matrix index of matrix A and returns matrix index of
// A * iHadmard
  matns_ptr_t A_ptr = (matns_ptr_t)A_Rptr;

  mat_level_t row_level = A_ptr->row_level;
  mat_level_t col_level = A_ptr->col_level;

  if (col_level == 0) {
    fprintf(stderr,"ERROR: in %s, there is no sensible definition for level 0 Hadamard\n",__func__);
    exit(1);
  }

  // find the index of the iHadamard matrix of this level
  int64_t iHad_pID = get_iHadamard_pID(col_level);

  uint64_t hash = hash_from_op(A_pID, iHad_pID, PRODUCT);
  int64_t A_times_iHad_pID = op_get(PRODUCT, A_pID, iHad_pID, hash);

  if (A_times_iHad_pID != MATRIX_ID_INVALID) return A_times_iHad_pID;

  // if no stored value and the matrix level is 1, calculate the product A * HH1
  // using the preloaded value for the integer iHadamard
  if (col_level == 1) return matrix_mult(A_pID, iHad_pID);
  
  // if no stored value and the matrix larger than 2x2,
  // calculate the iHadamard product recursively as follows
  // this should be cheaper than matrix_mult(A_pID,iHad_pID);
  int64_t sub_pID[4];
  sub_pID[0] = matrix_times_iHadamard(matrix_add(
      A_ptr->subMatList[0], A_ptr->subMatList[1]));
  sub_pID[1] = matrix_times_iHadamard(matrix_add(
      A_ptr->subMatList[0], scalar_mult(packedID_scalarM1,
                                                 A_ptr->subMatList[1])));
  sub_pID[2] = matrix_times_iHadamard(matrix_add(
      A_ptr->subMatList[2], A_ptr->subMatList[3]));
  sub_pID[3] = matrix_times_iHadamard(matrix_add(
      A_ptr->subMatList[2], scalar_mult(packedID_scalarM1,
                                                 A_ptr->subMatList[3])));
  A_times_iHad_pID = get_pID_from_array_of_four_sub_pIDs(sub_pID,
                     row_level, col_level);

  // STORE IN OP STORE
//    printf("calling op_set from matrix_times_iHadamard_PTR\n");
  op_set(PRODUCT, A_pID, iHad_pID, A_times_iHad_pID, hash);

  return A_times_iHad_pID;
}

char *tracenorm(int64_t m_pID, char *scale_factor)
{

  // The following checks to ensure that m_pID
  // has valid record pointer (eg, not RECORD_PTR_INVALID)
  record_ptr_t m_Rptr;
  check_validity_one_input(m_pID, __func__, &m_Rptr);

  // do basic error checking
  if (m_pID == MATRIX_ID_INVALID)
  {
     fprintf(stderr,"ERROR in %s: invalid packedID passed\n",__func__);
     exit(1);
  }

  scalarType scalar;
  sca_init(&scalar);
  sca_set_str(&scalar, scale_factor);

  scalarType mat_trace; 
  sca_init(&mat_trace);

  int64_t matadj_pID = adjoint(m_pID);
  int64_t det_times_id_mat = matrix_mult(matadj_pID,m_pID);

  if (IS_SCALAR(m_pID)) // SCALAR
  {
     mats_ptr_t sca_ptr =
       (mats_ptr_t)get_recordPTR_from_pID(det_times_id_mat,"",__func__,0);
     sca_mult(&mat_trace, scalar, sca_ptr->scalar_value);
  }
  else
  {
     matns_ptr_t mat_ptr =
          (matns_ptr_t)get_recordPTR_from_pID(det_times_id_mat,"",__func__,0);
     sca_mult(&mat_trace, scalar, mat_ptr->trace_element);
  }

  char *mat_trace_str = sca_get_readable_approx_str(mat_trace);
  sca_clear(&scalar);
  sca_clear(&mat_trace);

  return mat_trace_str;
}


char *get_scalar_value_string(int64_t m_pID)
{
   // This function checks to see if the packedID passed to it refers to a
   // matrix of SCALAR type. If it is not a SCALAR, it returns an empty string.
   // If it is a SCALAR, it calls the function traceID()
   // since the scalar value is currently stored in the matrix store record
   // from a union which stores either the 'scalarType trace_element' or the
   // "scalarType scalar_value".  
   // The string returned is stored in malloc'd memory. If this function
   // is to be called frequently enough to lock up a substantial amount of
   // memory, the returned strings should be freed.

   // The following checks to ensure that m_pID
   // has valid record pointer (eg, not RECORD_PTR_INVALID)
   record_ptr_t m_Rptr;
   check_validity_one_input(m_pID, __func__, &m_Rptr);

   // if scalar, return scalar_value
   if (IS_SCALAR(m_pID)) return sca_get_readable_approx_str(
         ((mats_ptr_t)m_Rptr)->scalar_value);

   // otherwise print error message and return empty string
   fprintf(stderr,"ERROR: matrixID %" PRId64 " passed to %s",
           MID_FROM_PID(m_pID),__func__);
   fprintf(stderr,"is not a level (0,0) matrix (SCALAR).\n");
   fprintf(stderr,"Returning empty string.\n");
   char *ret_string = malloc(1*sizeof(char));
   if (ret_string == NULL) { ALLOCFAIL(); }
   ret_string = "";
   return ret_string;
}

char *traceID(int64_t m_pID)
{
   // This routine returns a string giving the trace of the matrix m_pID.
   // If the matrix is a SCALAR, the trace of the matrix is the value of
   // the scalar. The trace is not defined for a non-square MATRIX (or any
   // ROW_VECTOR or COL_VECTOR); in these cases an empty string is returned.
   // In any case, the string returned is malloc'd; if this function is to
   // be called frequently enough to lock up substantial amounts of memory,
   // the returned strings should be freed.

   // The following checks to ensure that m_pID
   // has valid record pointer (eg, not RECORD_PTR_INVALID)
   record_ptr_t m_Rptr;
   check_validity_one_input(m_pID, __func__, &m_Rptr);

   if (IS_SCALAR(m_pID)) // SCALAR
   {
      mats_ptr_t sca_ptr = (mats_ptr_t)m_Rptr;
      return sca_get_readable_approx_str(sca_ptr->scalar_value);
   }

   matns_ptr_t mat_ptr = (matns_ptr_t)m_Rptr;

   // confirm that the matrix has a well-defined trace
   matrix_type_t matrixType = matrix_type(mat_ptr);

   if ( (matrixType==ROW_VECTOR) || (matrixType==COL_VECTOR) || 
      (mat_ptr->row_level != mat_ptr->col_level) )
   {
      fprintf(stderr,"ERROR: matrixID %" PRId64 " passed to %s",
              MID_FROM_PID(m_pID),__func__);
      fprintf(stderr,"is for a non-square matrix.\n");
      fprintf(stderr,"Returning empty string.\n");
      char *ret_string = malloc(1*sizeof(char));
      if (ret_string == NULL) { ALLOCFAIL(); }
      ret_string = "";
      return ret_string;
   }

   // return the trace of the matrix in string format
   return sca_get_readable_approx_str(mat_ptr->trace_element);
}

/*!
 * \ingroup larc
 * \brief Utility for finding the number of entries with a given value in a matrix.
 *
 * \param count   The count is returned through this parameter. 
 * \param mat_pID A packedID for a matrix in LARC format. 
 * \param scalar  A scalar value to count. 
 */
static void matrix_count_entries_base(mpz_t count, int64_t mat_pID, scalarType scalar)
{
    // IDENTITY SHORT CUT
    if (matrix_is_zero(mat_pID)){
        // if scalar is zero, add number of elements in matrix to count;
        // if scalar is not zero, leave count unchanged.
        if (sca_eq(scalar, scalar0)){
            if (IS_SCALAR(mat_pID)) // scalar
                mpz_add_ui(count, count, 1);
            else
            {
                if (scratchVars.counter_in_use)
                   fprintf(stderr,"%s reusing scratchVars.counter!\n",__func__);
                scratchVars.counter_in_use = 1;
                mpz_t *matSize = &scratchVars.counter;
                matns_ptr_t mat_ptr = (matns_ptr_t)get_recordPTR_from_pID(mat_pID,
                           "",__func__,0);
                mpz_ui_pow_ui(*matSize, 2,
                    mat_ptr->row_level + mat_ptr->col_level);
                mpz_add(count, count, *matSize);
                scratchVars.counter_in_use = 0;
            }
        }
        return;
    }

    // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE
    //TODO: how to store in op store? 
    //In the past we would do something like:
    //op_set(ZEROCOUNT, mat_ptr, mat_ptr, cnt_ptr);
    //but cnt_ptr holds scalarType which is not always mpz_t. 

    if (IS_SCALAR(mat_pID)) // SCALAR
    {
        mats_ptr_t sca_ptr = (mats_ptr_t)get_recordPTR_from_pID(mat_pID,
                   "",__func__,0);
        // if scalar matches matrix entry, increment count by one. 
        if (sca_eq(sca_ptr->scalar_value, scalar))
            mpz_add_ui(count, count, 1);
        return;
    }
    else
    {
        // increase count in each submatrix
        matns_ptr_t mat_ptr = (matns_ptr_t)get_recordPTR_from_pID(mat_pID,
                   "",__func__,0);
        for (int i = 0; i < 4; i++)
        {
            int64_t subMat = mat_ptr->subMatList[i];
            if (subMat != MATRIX_ID_INVALID)
                matrix_count_entries_base(count, subMat, scalar);
        }
    }
}

char *matrix_count_entries(int64_t mat_pID, char *scalar_str)
{

   // The following checks to ensure that mat_pID
   // has valid record pointer (eg, not RECORD_PTR_INVALID)
   record_ptr_t mat_Rptr;
   check_validity_one_input(mat_pID, __func__, &mat_Rptr);

    scalarType scalar;
    mpz_t count;

    // interpret scalar_str as scalarType
    sca_init(&scalar);
    sca_set_str(&scalar, scalar_str);

    // call matrix pointer version of function
    mpz_init(count);
    mpz_set_ui(count, 0);

    matrix_count_entries_base(count, mat_pID, scalar);

    // prepare string output from mpz_t count
    size_t out_str_size = mpz_sizeinbase(count, 10) + 2;
    char *out_str = calloc(out_str_size, sizeof(char));
    if (NULL == out_str){
        fprintf(stderr,"ERROR: allocation failure in %s.\n", __func__);
        exit(1);
    }
    gmp_snprintf(out_str, out_str_size, "%Zd", count);

    // clean scalarType and mpz_t variables
    sca_clear(&scalar);
    mpz_clear(count);

    return out_str;
}

void get_array_of_scalars_in_larcMatrixFile(scalarType **scalars_ptr, int64_t *numScalars, char *path)
{
  scalarType *scalars = NULL;
  *numScalars = 0;
  printf("ready to start: num set to %" PRId64 "\n", *numScalars);
  
  int verbose = VERBOSE;
  FILE *f = fopen(path, "r");
  if (f == NULL)
  {
    fprintf(stderr,"%s: no file found at\n\t%s:\n\texiting\n", __func__, path);
    exit(1);
  }
  json_t *j = j_parse_file(f);

  json_t *t = j_key_lookup(j, "table");

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
        case 3:                        /* expect a scalar link */
          if (verbose>BASIC) {printf("In scalar case \n");}
              int row_level = j_get_num64(j_array_index(p, 0));
              int col_level = j_get_num64(j_array_index(p, 1));
          if (!j_is_string(j_array_index(p,2))){
              fprintf(stderr,"ERROR in %s: expected scalar values to be given as strings. Perhaps try 'legacy' version.\n", __func__);
              exit(1);
          }
          scalars = (scalarType *) realloc(scalars, (*numScalars+1)*sizeof(scalarType));
          if (scalars == NULL){
            fprintf(stderr,"ERROR: Could not grow list of scalars to size %" PRId64 " in %s.\n", *numScalars+1, __func__);
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

        case 6:                        /* expect a 4-tuple */
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

char *get_list_of_scalars_in_larcMatrixFile(char *path)
{
    scalarType *scalars;
    int64_t numScalars;
    get_array_of_scalars_in_larcMatrixFile(&scalars, &numScalars, path);

    int64_t string_size = numScalars + 1;
    for (int64_t i = 0; i < numScalars; i++) {
        char *scalar_string = sca_get_readable_approx_str(scalars[i]);
        string_size += strlen(scalar_string);
        free(scalar_string);
    }

    char *scalarList = calloc(string_size, sizeof(char));
    if (scalarList == NULL) { ALLOCFAIL(); }
    for (int64_t i = 0; i < numScalars; i++){
        char *scalar_string = sca_get_readable_approx_str(scalars[i]);
        scalarList = strcat(scalarList, scalar_string);
        free(scalar_string);
        scalarList = strcat(scalarList, ",");
    }

    scalarList[strlen(scalarList)-1] = 0;

    return scalarList;
}

// This routine works for small matrices, but according to the TODO in the
// code, doesn't always work.
static void mpz_locate_entries_larcMatrixFile(int64_t ***locations_ptr, mpz_t count, char *path, scalarType scalar, unsigned maxno)
{
    FILE *file = fopen(path, "r");
    int64_t i;

    struct count_info {
        mpz_t count;
        int64_t subMatIDs[4];
    };

    struct sca_info {
        int64_t matID;
        int64_t row;
        int64_t col;
    };

    if (file == NULL)
    {
        fprintf(stderr,"%s: no file found at\n\t%s:\n\texiting\n", __func__, path);
        exit(1);
    }
    json_t *jData = j_parse_file(file);
    //printf("  parsed json\n");

    int64_t matrixID = j_lookup_num64(jData, "matid");
    //printf("  matrixID is %ld\n", matrixID);
    json_t *jTable = j_key_lookup(jData, "table");

    // check that we can allocate an array this large
    // (i.e., total number of bytes does not exceed SIZE_MAX)
    if (SIZE_MAX/(matrixID + 1)/sizeof(mpz_t) == 0) {
        fprintf(stderr,"Error is %s: matrixID_max too large - try renumbering.\n", __func__);
        exit(1);
    }

    struct count_info *records = malloc((matrixID + 1) * sizeof(struct count_info));
    if (records == NULL) { ALLOCFAIL(); }
    // matrixID + 1 because we'll be indexing by [0, matrixID].
    for (i = 0; i < matrixID + 1; i ++){
        mpz_init(records[i].count);
        mpz_set_ui(records[i].count, 0);
    }

    // for each key/matrixID in table, count number of zeros
    int64_t jTableLen = j_key_count(jTable);
    int64_t scaID = -1; //to record whether we have found table entry for scalar yet.
    scalarType temp;
    sca_init(&temp);
    for (i = 0; i < jTableLen; i++){
        json_t *entry = j_key_index(jTable, i);
        int64_t matID = (int64_t) atoll(entry->name);

        if (j_is_array(entry)){
            // If we have found the entry in the matrix holding scalar, then
            // we don't need to consider SCALAR entries anymore.
            if (scaID == -1){
                // SCALAR entries have length 3
                // if scalar value == scalar, record 1 for entry matID
                if (3 == j_array_get_length(entry)){
                    // get string value of table entry
                    if (!j_is_string(j_array_index(entry,2))){
                        fprintf(stderr,"ERROR in %s(%s): entry 2 of line not string (%" PRId64 " of table)\n",
                                __func__, path, i);
                        exit(1);
                    }
                    // set temp scalarType with entry string
                    sca_set_str(&temp, j_get_string(j_array_index(entry, 2)));
                    //printf("  entry %ld: scalar %s\n", i, sca_get_readable_approx_str(temp));
                    // if temp == scalar, we found the entry so set count to one
                    // for this matID.
                    if (0 != sca_eq(temp, scalar)){
                        // record that we don't need to check scalars anymore.
                        // by recording matrixID of scalar
                        scaID = matID;
                        // record count of one
                        mpz_set_ui(records[matID].count, 1);
                    }
                }
            }
            // If we have not found the entry in the matrix holding scalar,
            // then the counts of any matrices from table entries must be zero
            // so not worth looking at them.
            else {
                // SUBMATRIX entries have length 6
                if (6 == j_array_get_length(entry)){
                    for (int j = 0; j < 4; j++){
                        // get the number of zeros for each submatrix and
                        // add to the total count for current matID.
                        int64_t subMatID = j_get_num64(j_array_index(entry, j+2));
                        if (subMatID != MATRIX_ID_INVALID)
                            mpz_add(records[matID].count, records[matID].count, records[subMatID].count);
                        // record the submatrices that contributed
                        records[matID].subMatIDs[j] = subMatID;
                    }
                }
            }
        } // end (if j_is_array(entry))
    } // end loop over len
    sca_clear(&temp);
    j_set_null(jData); free(jData);


//    //Debug
//    for (i = 0; i < matrixID + 1; i ++){
//        gmp_printf("MatrixID %ld has %Zd entries\n", i, numZeros[i]);
//    }


    // if there were no entries with scalar, return NULL for empty list
    if (mpz_cmp_ui(records[matrixID].count, 0) == 0){
        mpz_set_ui(count, 0);
        *locations_ptr = NULL;
        free(records);
        return;
    }

    // (~arbitrarily) we'll say that if there are more than 2*32 appearances
    // of a number, we'll refuse to locate them all. In particular, we'll
    // protect against overflowing <<total>>.
    if (mpz_cmp_ui(records[matrixID].count, maxno) > 0){
        fprintf(stderr,"error in %s: too many appearances of scalar to return.\n", __func__);
        exit(0);
    }

    // allocate an array of info for each scalar occurence
    int64_t total = mpz_get_si(records[matrixID].count);
    //TODO: in my example, this total value is coming out WAY high. for rd11.json
    // error in mpz_locate_entries_larcMatrixFile: failed to allocate loc array with total -172825635700342774.
    // create and initialize array of length numZeros[matrixID] and width 3
    // [current matrixID, partial row bin expansion, partial col bin expansion]
    struct sca_info *locs = calloc(total, sizeof(struct sca_info));
    if (locs == NULL){
        fprintf(stderr,"error in %s: failed to allocate location info array with size %" PRId64 ".\n", __func__, total);
        exit(0);
    }
    for (i = 0; i < total; i++){
        // set current matrixID to final matrix matrixID for each scalar
        locs[i].matID = matrixID;
    }

    // until we hit the scalar level, make passes through the list of scalar
    // locations
    //int p = 0; // temp var for debug statement below
    while (locs[0].matID != scaID){
        for (i = 0; i < total; ){ // incrementing happens inside loop
            // for a given scalar location, consider the matrixID to do next
            int64_t curID = locs[i].matID;
            int64_t j,k;
            // go through each of the 4 quadrants
            for (j = 0; j < 4; j++){
                // make edits to / replace matrixID-to-do-next with quadrant
                int64_t subMatID = records[curID].subMatIDs[j];
                // and do this to exactly the number of scalars that the
                // quadrant contributes. Note that since total fits in an
                // int64_t, records[subMatID].count must be able to as well.
                if (subMatID != MATRIX_ID_INVALID){
                    for (k = 0; mpz_cmp_ui(records[subMatID].count, k) > 0; k++){
                        // reset next matrixID
                        locs[i+k].matID = subMatID;
                        // include row/col location info (add and shift)
                        locs[i+k].row = (locs[i+k].row << 1) + (j>1);
                        locs[i+k].col = (locs[i+k].col << 1) + (j%2);
                        // debug:
                        //printf("touched ind %ld on pass %d: subMatID = %ld, curID = %ld: j=%ld -> add %d to r, add %ld to c\n", i+k, p, subMatID, curID, j, j>1, j%2);
                    }
                    i += k; //increment i by the progress we've made, namely k=records[subMatID]
                }
            }
        }
        //p += 1; //part of debug print statement
    }

//    // debug
//    for (i = 0; i < total; i++){
//        printf("scalar occ %ld: row %ld col %ld\n", i, locs[i].row, locs[i].col);
//    }

    //printf("Preparing output\n");
    mpz_set(count, records[matrixID].count);
    *locations_ptr = malloc(total * sizeof(int64_t*));
    if (locations_ptr == NULL) { ALLOCFAIL(); }
    for (i = 0; i < total; i++){
        (*locations_ptr)[i] = malloc(2 * sizeof(int64_t));
        (*locations_ptr)[i][0] = locs[i].row;
        (*locations_ptr)[i][1] = locs[i].col;
    }

    for(i = 0; i < matrixID + 1; i ++)
        mpz_clear(records[i].count);
    free(records);
    free(locs);
}

int64_t **locate_entries_larcMatrixFile(char *path, char *scalar_str, unsigned maxno)
{
    scalarType scalar;
    mpz_t ret;
    int64_t **locations;
    int64_t **locations2;

    // interpret scalar_str as scalarType
    sca_init(&scalar);
    sca_set_str(&scalar, scalar_str);

    // run count entries routine
    mpz_init(ret);
    mpz_locate_entries_larcMatrixFile(&locations, ret, path, scalar, maxno);

    // prepare for swig interface
    // as of now, the ends of the r/c lists are marked with -1
    // and the end of the array is marked with NULL
    int64_t i;
    locations2 = malloc((mpz_get_ui(ret) + 1) * sizeof(int64_t *));
    if (locations2 == NULL) { ALLOCFAIL(); }
    for (i = 0; mpz_cmp_ui(ret, i) > 0; i++){
        locations2[i] = malloc(3*sizeof(int64_t));
        if (locations2[i] == NULL) { ALLOCFAIL(); }
        locations2[i][0] = locations[i][0];
        locations2[i][1] = locations[i][1];
        locations2[i][2] = -1;
    }
    locations2[i] = NULL;

    // clean scalarType and mpz_t variables
    sca_clear(&scalar);
    mpz_clear(ret);
    free(locations);

    return locations2;
}

gmp_randstate_t state;


/***********************************************************************************
 * random_bool_matrix_from_count, random_Boolean_matrix_from_count
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

/*!
 * \ingroup larc
 *
 * \brief Generate a random matrix by randomly distributing a given 
 * number of scalars into quadrants of zeros. 
 *
 * \param row_level row level of matrix
 * \param col_level col level of matrix
 * \param scalarNum  number of entries with value one desired in matrix
 * \param scalarPtr  pointer to scalar value - increases efficiency.
 *
 * \return packedID of matrix generated
 */
static int64_t recursive_random_bool_matrix_from_count(
        mat_level_t row_level, mat_level_t col_level,
        mpz_t scalarNum, mats_ptr_t scalarPtr)
{
    int reporting = 0;
    mat_level_t min_reporting_level = 30; 

    if (0 > mpz_cmp_ui(scalarNum, 0)){
        char *trace_string = sca_get_readable_approx_str(scalarPtr->scalar_value);
        gmp_fprintf(stderr,"ERROR: the number of entries with value %s must be positive in %s, not %Zd.\n", 
                trace_string, __func__, scalarNum);
        free(trace_string);
        exit(0);
    }

    // reporting, if reporting. 
    if ((row_level >= min_reporting_level) || (col_level >= min_reporting_level)){
        if (reporting) {
            char *trace_string = sca_get_readable_approx_str(scalarPtr->scalar_value);
            gmp_fprintf(stderr,"  %s: requesting a dim level %dx%d with %Zd entries with value %s.\n", 
                    __func__, row_level, col_level, scalarNum, trace_string);
            free(trace_string);
        }
    }

    // MATH SHORTCUTS
    // if no <scalar>s, return zero matrix
    if (0 == mpz_cmp_ui(scalarNum, 0))
        return get_zero_pID(row_level, col_level);

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
        return scalarPtr->packedID;
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
    if (rands == NULL) { ALLOCFAIL(); }
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

    int64_t sub_pID[4];
    // COL VECTOR
    if (col_level == 0){
        sub_pID[0] = recursive_random_bool_matrix_from_count(row_level-1, 0, rands[0], scalarPtr);
        sub_pID[1] = MATRIX_ID_INVALID;
        sub_pID[2] = recursive_random_bool_matrix_from_count(row_level-1, 0, rands[1], scalarPtr);
        sub_pID[3] = MATRIX_ID_INVALID;
    }
    // ROW VECTOR
    else if (row_level == 0){
        sub_pID[0] = recursive_random_bool_matrix_from_count(0, col_level-1, rands[0], scalarPtr);
        sub_pID[1] = recursive_random_bool_matrix_from_count(0, col_level-1, rands[1], scalarPtr);
        sub_pID[2] = MATRIX_ID_INVALID;
        sub_pID[3] = MATRIX_ID_INVALID;
    }
    // MATRIX
    else {
        for (i = 0; i < 4; i++){
            sub_pID[i] = recursive_random_bool_matrix_from_count(
                row_level-1, col_level-1, rands[i], scalarPtr);
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
            char *trace_string = sca_get_readable_approx_str(scalarPtr->scalar_value);
            gmp_printf("  %s: returning a dim level %dx%d with %Zd entries with value %s.\n",
                    __func__, row_level, col_level, scalarNum, trace_string);
            free(trace_string);
        }
    }

    // FIND OR CREATE APPROPRIATE MATRIX PTR
    return get_pID_from_array_of_four_sub_pIDs(sub_pID,row_level,col_level); 
}

int64_t random_Boolean_matrix_from_count(mat_level_t row_level, mat_level_t col_level, char* numOnes)
{
    mpz_t num;
    mpz_init(num);

    // scan number of desired ones into an mpz_t struct
    gmp_sscanf(numOnes, "%Zd", num);

    // call matrix pointer version of function 
    mats_ptr_t onePtr = get_scalarPTR_for_scalarVal(scalar1);

    int64_t mat_pID = recursive_random_bool_matrix_from_count(
        row_level, col_level, num, onePtr);
    mpz_clear(num);
    return mat_pID;
}

int64_t matrix_basischange_A_by_B(int64_t B_pID, int64_t A_pID)
{
    // this routine should be used when the goal is to return the matrix which
    // is A conjugated by B. For this operation, there are two matrix
    // multiplies, and we do not care about the intermediate result, so we 
    // clean it from the matrix store when it is no longer needed (with the
    // proper protections so that we don't clean any submatrices needed for
    // other stored matrices).

    // The following checks to ensure that both A_pID and B_ID 
    // have valid record pointers (eg, not RECORD_PTR_INVALID)
    record_ptr_t A_Rptr, B_Rptr;
    check_validity_two_input(A_pID, B_pID,__func__, &A_Rptr, &B_Rptr);

    if (IS_SCALAR(A_pID)== 0) // non-SCALAR
    {
        matns_ptr_t A_ptr = (matns_ptr_t)A_Rptr;

        // matrix A must be square for conjugation to make sense
        if (A_ptr->col_level != A_ptr->row_level) {
           fprintf(stderr,"ERROR: in %s, matrix A_ptr",__func__);
           fprintf(stderr," (with matrixID %" PRIu64 ") is not square!\n",
                       MID_FROM_PID(A_pID));
           exit(EXIT_FAILURE);
        }
    }

    // mismatch in row level of B and col level of A will be detected 
    // in matrix_mult()

    uint64_t hash = hash_from_op(B_pID, A_pID, BASIS_CHANGE);
    int64_t result_pID = op_get(BASIS_CHANGE,B_pID,A_pID, hash);
    if (result_pID != MATRIX_ID_INVALID) // was found in op store
       return result_pID;

    // result was not found in operations store - calculation proceeds

    // ensure that cleaning doesn't take out the input matrices
    // NOTE: holds are actually counters, so internally holding and releasing
    // here won't change whether a matrix is held outside of this routine
    set_hold_matrix(A_pID);
    set_hold_matrix(B_pID);
    // calcluate B*A*B^{\dagger}
    // (for complex types, \dagger is complex conjugate transpose;
    //  for real types, it is simply transpose)
    int64_t temp1_pID = matrix_mult(B_pID,A_pID);

    int64_t conjAdj_pID = adjoint(B_pID);
    result_pID = matrix_mult(temp1_pID,conjAdj_pID);
    set_hold_matrix(result_pID);
    // safely remove all intermediate valued matrices
    set_hold_matrix(conjAdj_pID);
    remove_matrix_from_store(temp1_pID);
    release_hold_matrix(conjAdj_pID);
    remove_matrix_from_store(conjAdj_pID);
    // release held inputs and output
    release_hold_matrix(A_pID);
    release_hold_matrix(B_pID);
    release_hold_matrix(result_pID);

    // store result of calculation in operations store
//    printf("calling op_set from matrix_basischange_A_by_B\n");
    op_set(BASIS_CHANGE, B_pID, A_pID, result_pID, hash);

    return result_pID;
}

int64_t matrix_saxpy(int64_t A_pID, int64_t B_pID,
        int64_t scalar_a_pID, int64_t scalar_b_pID)
{
    // this routine performs the saxpy operation that is common in many
    // linear algebra applications. The result C = a*A + b*B, where a and b
    // are scalar values and A and B are similarly-dimensioned matrices
    // (usually row or column vectors).

    // The following checks to ensure that both A_pID and B_ID 
    // have valid record pointers (eg, not RECORD_PTR_INVALID)
    record_ptr_t A_Rptr, B_Rptr;
    check_validity_two_input(A_pID, B_pID,__func__, &A_Rptr, &B_Rptr);

    // The following checks to ensure that both scalar_a_pID and scalar_b_pID 
    // have valid record pointers (eg, not RECORD_PTR_INVALID)
    record_ptr_t a_Rptr, b_Rptr;
    check_validity_two_input(scalar_a_pID, scalar_b_pID,__func__,
          &a_Rptr, &b_Rptr);
    
    // confirm that scalar inputs are in fact scalars
    if (0 == IS_SCALAR(scalar_a_pID)) // not SCALAR
    {
       matns_ptr_t x = (matns_ptr_t)a_Rptr;
       fprintf(stderr,"ERROR in %s: 3rd argument should be a scalar,",__func__);
       fprintf(stderr,"but\n\thas row level %d and column level %d!\n",
               x->row_level,x->col_level);
        exit(EXIT_FAILURE);
    }
    if (0 == IS_SCALAR(scalar_b_pID)) // not SCALAR
    {
       matns_ptr_t x = (matns_ptr_t)b_Rptr;
       fprintf(stderr,"ERROR in %s: 4th argument should be a scalar,",__func__);
       fprintf(stderr,"but\n\thas row level %d and column level %d!\n",
               x->row_level,x->col_level);
        exit(EXIT_FAILURE);
    }

    // check to see if A and B are also scalars
    if  (IS_SCALAR(A_pID) && IS_SCALAR(B_pID))
    {
       int64_t result1_pID = scalar_mult(A_pID,scalar_a_pID);
       int64_t result2_pID = scalar_mult(B_pID,scalar_b_pID);
       return matrix_add(result1_pID,result2_pID);
    }
    else if  (IS_SCALAR(A_pID) || IS_SCALAR(B_pID))
    {
       // size disagreement - exit
       fprintf(stderr,"ERROR: in %s, ",__func__);
       fprintf(stderr,"input matrices do not have the same dimensions!\n");
       if (IS_SCALAR(A_pID))
       {
          fprintf(stderr,"\tmatrix A (with matrixID %" PRIu64 ")",
              MID_FROM_PID(A_pID));
          fprintf(stderr,"is a scalar, and matrrix B is not.\n");
       }
       else
       {
          fprintf(stderr,"\tmatrix B (with matrixID %" PRIu64 ")",
              MID_FROM_PID(B_pID));
          fprintf(stderr,"is a scalar, and matrrix A is not.\n");
       }
       exit(EXIT_FAILURE);
    }

    // neither A nor B is a SCALAR matrix
    matns_ptr_t A_ptr = (matns_ptr_t)A_Rptr;
    matns_ptr_t B_ptr = (matns_ptr_t)B_Rptr;

    // check for size agreement
    if (A_ptr->col_level != B_ptr->col_level
        || A_ptr->row_level != B_ptr->row_level)
    {
       fprintf(stderr,"ERROR: in %s, ",__func__);
       fprintf(stderr,"input matrices do not have the same dimensions!\n");
       fprintf(stderr,"\tmatrix A (with matrixID %" PRIu64 ")",MID_FROM_PID(A_pID));
       fprintf(stderr," has row level %d and column level %d",
                A_ptr->row_level,A_ptr->col_level);
       fprintf(stderr,"\tmatrix B (with matrixID %" PRIu64 ")",MID_FROM_PID(B_pID));
       fprintf(stderr," has row level %d and column level %d",
                B_ptr->row_level,B_ptr->col_level);
       exit(EXIT_FAILURE);
    }

    // scalar_mult will check KRONECKER op_store for a*A
    int64_t result1_pID = scalar_mult(scalar_a_pID,A_pID);
    // scalar_mult will check KRONECKER op_store for b*B
    int64_t result2_pID = scalar_mult(scalar_b_pID,B_pID);

    // matrix_add will check SUM op store for (a*A) + (b*B)
    int64_t result_pID = matrix_add(result1_pID, result2_pID);

    return result_pID;
}

int64_t vector_dot_product(int64_t A_pID, int64_t B_pID, int verbose)
{
    // This routine returns a packedID for the scalar product
    // of two vectors. MATRIX_ID_INVALID is returned if either matrix
    // passed is a MATRIX (both row level and col level > 0) or the
    // dimensions of the vectors passed do not agree.
    // If verbose is higher than BASIC, a lot of information about what
    // the routine is doing is printed to stderr.
    // get the matrix pointers from the packedID, and see if still in store    

    // The following checks to ensure that both A_pID and B_ID 
    // have valid record pointers (eg, not RECORD_PTR_INVALID)
    record_ptr_t A_Rptr, B_Rptr;
    check_validity_two_input(A_pID, B_pID,__func__, &A_Rptr, &B_Rptr);

    // shortcut - if both "vectors" are 1x1 
    if ( IS_SCALAR(A_pID) && IS_SCALAR(B_pID) ) {
       if (verbose>BASIC)  {
          fprintf(stderr,"%s doing 'dot product' of",__func__);
          fprintf(stderr,"two scalar values with matrixIDs ");
          fprintf(stderr,"%" PRId64 " and %" PRId64 "\n",
              MID_FROM_PID(A_pID), MID_FROM_PID(B_pID));
       }
       return scalar_mult(A_pID,B_pID);
    }

    // if only one of A or B is SCALAR, problem
    if ( IS_SCALAR(A_pID) || IS_SCALAR(B_pID) ) {
       fprintf(stderr,"ERROR: %s called with\n",__func__);
       fprintf(stderr,"matrix A with matrixID %" PRId64 " and\n",
          MID_FROM_PID(A_pID));
       fprintf(stderr,"matrix B with matrixID %" PRId64 ")\n",
          MID_FROM_PID(B_pID));
       matns_ptr_t x;
       if (IS_SCALAR(A_pID))
       {
           x = (matns_ptr_t)B_Rptr;
           fprintf(stderr,"matrix A is scalar and matrix B has levels ");
       }
       else
       {
           x = (matns_ptr_t)A_Rptr;
           fprintf(stderr,"matrix B is scalar and matrix A has levels ");
       }
       fprintf(stderr,"(%d,%d)\n",x->row_level,x->col_level);
       return MATRIX_ID_INVALID;
    }

    // both A and B now known not to be SCALAR
    matns_ptr_t A_ptr = (matns_ptr_t)A_Rptr;
    matns_ptr_t B_ptr = (matns_ptr_t)B_Rptr;
    matrix_type_t A_type = matrix_type(A_ptr);
    matrix_type_t B_type = matrix_type(B_ptr);

    // test for cases that don't fit dot product paradigm
    if ( (A_type==MATRIX) || (B_type==MATRIX) )
    {
       fprintf(stderr,"ERROR: %s called with\n",__func__);
       fprintf(stderr,"matrix A with packedID %" PRId64 " has level (%d,%d) and\n",
          A_pID,A_ptr->row_level,A_ptr->col_level);
       fprintf(stderr,"matrix B with packedID %" PRId64 " has level (%d,%d)\n",
          B_pID,B_ptr->row_level,B_ptr->col_level);
       return MATRIX_ID_INVALID;
    }

    // at this point, A and B are either ROW_VECTOR or COL_VECTOR types
    // check to make sure dimensions of the vectors agree
    int rdimA = A_ptr->row_level;
    int cdimA = A_ptr->col_level;
    int rdimB = B_ptr->row_level;
    int cdimB = B_ptr->col_level;
    int dimA = (rdimA) ? rdimA : cdimA;
    int dimB = (rdimB) ? rdimB : cdimB;
    if (dimA != dimB) {
       fprintf(stderr,"ERROR: vectors passed to %s\n",__func__);
       fprintf(stderr,"do not have the same nonzero level values\n");
       fprintf(stderr,"matrix A with packedID %" PRId64 " has level (%d,%d) and\n",
          A_pID,rdimA,cdimA);
       fprintf(stderr,"matrix B with packedID %" PRId64 " has level (%d,%d)\n",
          B_pID,rdimB,cdimB);
       return MATRIX_ID_INVALID;
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
          return matrix_mult(A_pID,B_pID);
       }
       // B_type==ROW_VECTOR
       if (verbose>BASIC) {
          fprintf(stderr,"in %s, A and B are both row vectors\n",__func__);
          fprintf(stderr,"returning A*B^{dagger}\n");
       }
       int64_t conj_pID = adjoint(B_pID);
       return matrix_mult(A_pID,conj_pID);
    }

    // A_type == COL_VECTOR
    if (B_type==ROW_VECTOR) {
       if (verbose>BASIC) {
             fprintf(stderr,"in %s, A is a column vector and\n",__func__);
             fprintf(stderr,"B is a row vector - returning B*A\n");
       }
          return matrix_mult(B_pID,A_pID);
    }
    // B_type==COL_VECTOR
    if (verbose>BASIC) {
       fprintf(stderr,"in %s, A and B are both column vectors\n",__func__);
       fprintf(stderr,"returning A^{dagger}*B\n");
    }
    int64_t conj_ID = adjoint(A_pID);
    return matrix_mult(conj_ID,B_pID);
}

/****************** 
 *                *
 * NORM FUNCTIONS *
 *                *
 ******************/


/*!
 * \ingroup larc
 * \brief A recursive routine to both calculate the matrix norm and find the matrix element with maximal norm
 *
 * \param matID The input packedID
 * \param norm On call, should be zero; on return, is the value of the norm of the input matrix
 * \param norm_pIDptr A pointer to the packedID of the norm value
 * \param maxScalarIDptr On call, should point to scalar zero; on return, a pointer to the packedID of the element of the input matrix with largest norm value
 * \param whichNorm An enum value; one of \{ L_infty, L_1, L_2 \}
 */
static void get_norm_and_maxElement(int64_t mat_pID, scalarType *norm,        
        int64_t *norm_pIDptr, int64_t *maxScalarIDptr, norm_type_t whichNorm)
{
        // The following checks to ensure that mat_pID
        // has valid record pointer (eg, not RECORD_PTR_INVALID)
        record_ptr_t m_Rptr;
        check_validity_one_input(mat_pID, __func__, &m_Rptr);

        // these values should be passed into this routine
        // sca_set_str(norm,"0");
        // *maxScalarIDPTR = packedID_scalar0;

        // we may need some local scalarType variables, but won't 
        // initialize them unless needed
        scalarType oldval, newval;

        // We will return *maxScalarPTR, *norm, and *norm_pIDptr to the calling
        // routine

        // SCALAR case: matrix norm is just norm
        //         and pointer to scalar must be pointer to max in matrix
        if (IS_SCALAR(mat_pID))
        {
                *maxScalarIDptr = mat_pID;
                mats_ptr_t scaPTR = (mats_ptr_t)m_Rptr;
//                printf("scalar - value is %s\n",
//                        sca_get_readable_approx_str(scaPTR->scalar_value));
                // the value we need in localNorm is either the norm
                // of the scalar or (for L2 norm) its norm squared
                sca_init(&newval);
                sca_norm(&newval, scaPTR->scalar_value);
                if (whichNorm == L_2) {
                        sca_mult(norm, newval, newval);
                }
                else { sca_set(norm, newval); }
                // put the localNorm value into the MatrixStore
                *norm_pIDptr = (get_scalarPTR_for_scalarVal(*norm))->packedID;
                // clean up local scalarType memory
                sca_clear(&newval);
                return;
        }

        // check the OperationsStore to see if we already calculated
        // the norm for this matrix
        uint64_t normHash, maxHash1, maxHash2;
        normHash = hash_from_op(mat_pID,(uint64_t)whichNorm,NORM_VAL);

        int64_t norm_pID = op_get(NORM_VAL,mat_pID,
                     (uint64_t)whichNorm,normHash);

        if ( norm_pID != MATRIX_ID_INVALID )
        {
                *norm_pIDptr = norm_pID;
                mats_ptr_t normPTR = (mats_ptr_t)get_recordPTR_from_pID(
                    norm_pID,"",__func__,0);
                // we have already calculated the value we are 
                // looking for; use the pointer to set the value
                sca_set(norm, normPTR->scalar_value);
//                printf("found norm, its value is %s\n",
//                        sca_get_readable_approx_str(norm));

                // the pointer to the element with max norm is also
                // in the OperationStore
                maxHash1 = hash_from_op(mat_pID,mat_pID,MAX_ELEMENT);
                *maxScalarIDptr = op_get(MAX_ELEMENT,mat_pID,mat_pID,maxHash1);
                return;
        }

        // need a scalarType variable to hold the norm value for the matrix
        // for L\infty norm, it is the norm of the maximum scalar in the matrix
        // for L1 norm, it is the sum of the norms of all scalars in the matrix
        // for L2 norm, it is the sum of the squares of all the norms
        //         (square root taken later)
        scalarType quadNorm;
        sca_init(&quadNorm);

        // possible efficiencies that would complicate the code with branches:
        // * could check to see if maxScalarPTR already in store
        //         and not calculate it again
        // * oldval and newval not needed unless calculating maxScalarPTR for
        //         L1 and L2 norm, so initialization could sometimes be skipped

        // need extra scalarType variables for L1 and L2 norms
        sca_init(&oldval);
        sca_init(&newval);

        // set quadrant norm to zero, local pointer to
        // point at stored zero (true on entry)
        // sca_set_str(norm,"0");
        // *maxScalarIDPTR = packedID_scalar0;
        
        // set the value of the norm of maxScalarPTR (currently zero)
        sca_set_str(&oldval,"0");

        matns_ptr_t matPTR =(matns_ptr_t)m_Rptr;
        for (int i=0;i<4;++i) 
        {
                int64_t quadNorm_pID;
                // set the passed variables to initial zero values
                sca_set_str(&quadNorm,"0");
                int64_t quadMaxScalarID = packedID_scalar0;
                // call recursion on each of the four quadrant submatrices
//                printf("i=%d: recursion\n",i);
                int64_t quad_pID = matPTR->subMatList[i];
                if (quad_pID != MATRIX_ID_INVALID)
                    get_norm_and_maxElement(quad_pID, &quadNorm,
                        &quadNorm_pID, &quadMaxScalarID, whichNorm);
//                printf("returned from recursion\n");
//                printf("i=%d; quadnorm is %s\n",i,sca_get_readable_approx_str(quadNorm));

                // for L_infty, keep maximum norm value
                //         and pointer to this value
                // for L_1 and L_2, accumulate sums of norm values
                //         and keep track of whether the max element changes
                if ( whichNorm == L_infty) 
                {
                        if (sca_cmp(quadNorm, *norm) > 0 )
                        {
//                                printf("updating norm value (L_infty)\n");
                                sca_set(norm,quadNorm);
                                *maxScalarIDptr = quadMaxScalarID;
                        }
                }
                else 
                {
                        // norm += quadNorm
                        sca_set(&newval, *norm);
                        sca_add(norm, newval, quadNorm);
                        // determine whether to change maxScalarPTR
                        mats_ptr_t newptr = (mats_ptr_t)get_recordPTR_from_pID(
                             quadMaxScalarID,"",__func__,0);
                        sca_norm(&newval, newptr->scalar_value);
                        if (sca_cmp(newval,oldval) > 0)
                        {
//                                printf("updating norm value (L_1 or L_2)\n");
                                *maxScalarIDptr = quadMaxScalarID;
                                sca_set(&oldval, newval);
                        }
                }
        } // loop over quadrants

        // check to see if the maximum element pointer is already in
        // the OperationStore, and if not add it
        // (this will happen if a user calculates two different norms for
        // the same matrix)
        maxHash2 = hash_from_op(mat_pID,mat_pID,MAX_ELEMENT);
        if (op_get(MAX_ELEMENT,mat_pID,mat_pID,maxHash2) == MATRIX_ID_INVALID)
        {
//           printf("calling op_set from get_norm_and_maxElement, maxHash2\n");
           op_set(MAX_ELEMENT,mat_pID,mat_pID,*maxScalarIDptr,maxHash2);
        }

        // put the norm value into the MatrixStore
        *norm_pIDptr = (get_scalarPTR_for_scalarVal(*norm))->packedID;
        // record the norm calculation in the OperationStore
//        printf("calling op_set from get_norm_and_maxElement, normHash\n");
//        normHash = hash_from_op(matID,(uint64_t)whichNorm,NORM_VAL);
        op_set(NORM_VAL,mat_pID,(uint64_t)whichNorm,*norm_pIDptr,normHash);

        // clean up local scalarType variables before exiting routine
        sca_clear(&oldval);
        sca_clear(&newval);
        sca_clear(&quadNorm);
        return;
}

int64_t normID(int64_t mat_pID, int whichNorm)
{
        // convert norm indicator from string to enum value
        norm_type_t whichNormEnum;
        switch (whichNorm)
        {
                case 0:
                        whichNormEnum = L_infty;
                        break;

                case 1:
                        whichNormEnum = L_1;
                        break;

                case 2:
                        whichNormEnum = L_2;
                        break;

                default:
                        printf("The whichNorm argument to %s ",__func__);
                        printf("must be one of 0, 1, or 2\n");
                        printf("Your input was %d\n", whichNorm);
                        printf("...proceeding as if you input 0\n");
                        whichNormEnum = L_infty;
        }

        int64_t norm_pID;
        mats_ptr_t scalarPTR = SCALAR_PTR_INVALID;
        // We create a scalarType variable to hold the maximum norm value
        scalarType norm;
        sca_init(&norm);
        uint64_t normHash, maxHash;


        // SHORTCUT: if matPTR a scalar, we just get its norm
        // (which is the same for all three matrix norms)
        // SHORTCUT: if looking for Linfty norm and max element is stored,
        // just take the norm of the max element
        if (IS_SCALAR(mat_pID))
           scalarPTR = (mats_ptr_t)get_recordPTR_from_pID(mat_pID, "",__func__,0);
        else if (whichNormEnum == L_infty) 
        {
                maxHash = hash_from_op(mat_pID,mat_pID,MAX_ELEMENT);
                int64_t scalarID = op_get(MAX_ELEMENT,mat_pID,mat_pID,maxHash);
                scalarPTR = (mats_ptr_t)get_recordPTR_from_pID(scalarID,
                                "",__func__,0);
        }
        if (scalarPTR != SCALAR_PTR_INVALID)
        {
                sca_norm(&norm,scalarPTR->scalar_value);
                norm_pID = (get_scalarPTR_for_scalarVal(norm))->packedID;
                // clean up local scalarType variable before exiting routine
                sca_clear(&norm);
                return norm_pID;
        }

        // check to see if norm has already been stored
        normHash = hash_from_op(mat_pID,(uint64_t)whichNormEnum,NORM_VAL);
        norm_pID = op_get(NORM_VAL,mat_pID,(uint64_t)whichNormEnum,normHash);
        if (norm_pID == MATRIX_ID_INVALID)
        {

                // the minimum norm value should be zero, so we start there
                sca_set_str(&norm,"0");
                // Initialize a value to scalar zero, updated by recursion
                int64_t maxScalarID = packedID_scalar0;
        
                // call the recursive function to set both the scalar
                // norm value and the pointer to a scalar in the matrix
                // which has this maximal value
                get_norm_and_maxElement(mat_pID, &norm, &norm_pID,
                        &maxScalarID, whichNormEnum);
        }
        else
        {
                mats_ptr_t normPTR = (mats_ptr_t)get_recordPTR_from_pID(
                      norm_pID,"",__func__,0);
                sca_set(&norm,normPTR->scalar_value);
        }

        // the value stored for L2 norm is the square of the norm
        if (whichNormEnum == L_2)
        {
                // set normPTR to point to the square root 
                if (scratchVars.submit_to_store_in_use)
                    fprintf(stderr,"%s reusing scratchVars.submit_to_store!\n",__func__);
                scratchVars.submit_to_store_in_use = 1;

                scalarType *rootNorm = &scratchVars.submit_to_store;
                sca_sqrt(rootNorm, norm);
                norm_pID = (get_scalarPTR_for_scalarVal(*rootNorm))->packedID;
                scratchVars.submit_to_store_in_use = 0;
        }

        // clean up local scalarType variable before exiting routine
        sca_clear(&norm);
        return norm_pID;
}

int64_t matrix_element_with_maxNorm(int64_t mat_pID)
{
        // The following checks to ensure that mat_pID
        // has valid record pointer (eg, not RECORD_PTR_INVALID)
        record_ptr_t m_Rptr;
        check_validity_one_input(mat_pID, __func__, &m_Rptr);

        if (mat_pID == MATRIX_ID_INVALID) exit(1);

        // SHORTCUT: if matID a scalar, it is the max element in its matrix
        if (IS_SCALAR(mat_pID)) return mat_pID;

        // If we have already calculated any norm for this matrix, our
        // answer is in the op store. If not, we call the norm function
        // to add the correct value to the op store, then get it again
        uint64_t maxHash = hash_from_op(mat_pID,mat_pID,MAX_ELEMENT);
        int64_t scalarID = op_get(MAX_ELEMENT,mat_pID,mat_pID,maxHash);
        if (scalarID == MATRIX_ID_INVALID)
        {
                // this points to *the norm of* the maximum element
                int64_t tempID = normID(mat_pID, L_infty);
                // line added to avoid compiler warning
                if (tempID == MATRIX_ID_INVALID) exit(1);
                // this (now) points to the maximum element
                scalarID = op_get(MAX_ELEMENT,mat_pID,mat_pID,maxHash);
        }
        return scalarID;
}

int64_t create_const_matrix(int64_t const_pID, mat_level_t row_level,
        mat_level_t col_level)
{
        // When the matrix is a square matrix, its smallest submatrix with
        // four equal quadrant submatrices is level-(1,1), and the constant
        // scalar pointer is used to build it. The new matrix is similarly
        // the building block for the next larger matrix, and so on until
        // the desired matrix is constructed.
        //
        // For non_square matrices, there are two possibilities. If either
        // row_level or col_level is zero, the desired vector/matrix has
        // MAX(0,row_level-col_level) and column level
        // MAX(0,col_level-row_level), and we must build this vector.
        // If both are nonzero, then four copies of the vector of this size
        // are the quadrant submatrices of the next larger matrix, and we
        // recurse on successively larger matrices until the desired matrix
        // is constructed.

        // The following checks to ensure that const_pID
        // has valid record pointer (eg, not RECORD_PTR_INVALID)
        record_ptr_t m_Rptr;
        check_validity_one_input(const_pID, __func__, &m_Rptr);

        // SHORTCUT: if the constant is zero, the matrix is already stored
        if (const_pID==packedID_scalar0)
            return get_zero_pID(row_level, col_level);

        // We recursively build our constant matrix from the bottom up.
       
        mat_level_t level = MAX(row_level,col_level);
        mat_level_t level_diff = row_level - col_level;
        mat_level_t abs_diff = labs(level_diff);

        int64_t temp_pID, sub_pID[4];
//        matns_ptr_t tempPTR, panel[4];
        int64_t i;

//        // initialize our pointer to the input pointer
        temp_pID = const_pID;
        mat_level_t row_offset = 0, col_offset = 0;

        // if level_diff = 0, then this next block of code does nothing.
        //
        // Otherwise, we generate the base row or column vector which we
        // will use to build the desired matrix. We start from the input
        // scalar and recursively build successively larger vectors.
        if (level_diff>0)
        {
            // more rows: generate column vector
            sub_pID[1] = sub_pID[3] = MATRIX_ID_INVALID;
            for (i=1;i<=abs_diff;++i)
            {
                sub_pID[0] = sub_pID[2] = temp_pID;
                temp_pID = get_pID_from_array_of_four_sub_pIDs(sub_pID,i,0);
            }
            col_offset = abs_diff;
        }
        else if (level_diff<0)
        {
            // more columns: generate row vector
            sub_pID[2] = sub_pID[3] = MATRIX_ID_INVALID;
            for (i=1;i<=abs_diff;++i)
            {
                sub_pID[0] = sub_pID[1] = temp_pID;
                temp_pID = get_pID_from_array_of_four_sub_pIDs(sub_pID,0,i);
            }
            row_offset = abs_diff;
        }

        // If the desired matrix is a vector, abs_diff==level and this loop
        // does nothing.
        //
        // Otherwise, we use our base matrix pointer (row or column vector,
        // or scalar) to generate a matrix for which it is the four quadrant
        // submatrices. This process is repeated for larger matrices until
        // the desired level-(row_level,col_level) matrix is constructed.
        for (i=abs_diff+1;i<=level;++i)
        {
            sub_pID[0] = sub_pID[1] =
                    sub_pID[2] = sub_pID[3] = temp_pID;
            temp_pID = get_pID_from_array_of_four_sub_pIDs(sub_pID,
                       i-row_offset, i-col_offset);
        }

        return temp_pID;
}

int64_t apply_function_to_matrix_values(int64_t m_pID,
        void (*func)(scalarType*, const scalarType), op_type_t op_memz)
{
    // This function produces a new matrix by applying a function to each
    // scalar in an input matrix. The expectation is that *func will produce a
    // result that cannot be obtained via linear algebra, for example, squaring
    // each element of a matrix. It is possible to use this routine to perform
    // linear algebraic operations such as multiplication by a scalar, but for
    // these cases memoization is improved if the predefined LARC linear
    // algebra function is used.

    // Memoization is enabled for this function unless the input op_memz is set
    // to INVALID_OP. Because BAD THINGS WILL HAPPEN if the op_type_t for a
    // different operation is chosen, we protect against this possibility for
    // all operations predefined for use by LARC. However, the user is
    // responsible for picking a different user-defined operation for each
    // different value of *func, e.g. FUNC_A for squaring each element and
    // FUNC_B for setting to zero all elements below some threshold.

    // To repeat: EACH different function passed MUST have a different
    // op_type_t value assigned, or this function will return incorrect
    // answers!


    // The following checks to ensure that m_pID
    // has valid record pointer (eg, not RECORD_PTR_INVALID)
    record_ptr_t m_Rptr;
    check_validity_one_input(m_pID, __func__, &m_Rptr);

    if ((op_memz < FUNC_A) || (op_memz > INVALID_OP))
    {
        fprintf(stderr,"In %s, invalid operation type given for\n", __func__);
        fprintf(stderr,"memoization purposes; exiting.\n");
        exit(0);
    }

    uint64_t hash = hash_from_op(m_pID, m_pID, op_memz);
    int64_t new_pID = op_get(op_memz,m_pID,m_pID,hash);
    if (new_pID != MATRIX_ID_INVALID) return new_pID;

    if (IS_SCALAR(m_pID))
    {
        mats_ptr_t s_ptr = (mats_ptr_t)m_Rptr;
        if (scratchVars.submit_to_store_in_use)
           fprintf(stderr,"%s reusing scratchVars.submit_to_store!\n",__func__);
        scratchVars.submit_to_store_in_use = 1;

        // apply func() to the scalar and put into the matrixStore
        scalarType *new = &scratchVars.submit_to_store;
        func(new,s_ptr->scalar_value);
        new_pID = (get_scalarPTR_for_scalarVal(*new))->packedID;
        scratchVars.submit_to_store_in_use = 0;
    }
    else
    {
        // recursively apply this function to each (valid) quadrant of the
        // input matrix
        matns_ptr_t m_ptr = (matns_ptr_t)m_Rptr;
        int64_t subMatList[4];

        for (int i = 0; i < 4; i++)
        {
            if (m_ptr->subMatList[i] != MATRIX_ID_INVALID)
                subMatList[i] = apply_function_to_matrix_values(
                     m_ptr->subMatList[i], func, op_memz);
            else subMatList[i] = MATRIX_ID_INVALID;
        }
        // add matrix with these submatrices to the matrix store
        new_pID = get_pID_from_array_of_four_sub_pIDs(subMatList,
                  m_ptr->row_level, m_ptr->col_level);
    }

//    printf("calling op_set from apply_function_to_matrix_values\n");
    op_set(op_memz, m_pID, m_pID, new_pID, hash);
    return new_pID;
}
