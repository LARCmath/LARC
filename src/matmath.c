//                          matmath.c 
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

#include <complex.h>

#include "matmath.h"


// inline int64_t get_matrixID_from_ptr(mat_add_t m_ptr) {
// return ((m_ptr==MATRIX_PTR_INVALID) ? MATRIX_ID_INVALID: m_ptr->matrixID);
//static int64_t
//record_mat_ptr_by_matrixID( mat_add_t m_ptr) {
//panel[0] = mat_ptr_from_matrixID(A_id);
// m_ptr->matrixID = store.matrixID_next++;
//uint64_t  
//num_matrices_created(void)
//{
//  uint64_t ret = store.matrixID_next;
//  return(ret);
//}

/* Python Interface for matrix_add(mat_add_t A_ptr, mat_add_t B_ptr):
   accepts matrixIDs for matrices to be added and converts
   inputs to pointers before calling. */

int64_t 
matrix_add_matrixID(int64_t A_mID, int64_t B_mID) {

  // get the matrix pointers from the matrixID, and see if still in store 
  mat_add_t A_ptr = mat_ptr_from_matrixID(A_mID, "first", __func__,0);
  mat_add_t B_ptr = mat_ptr_from_matrixID(B_mID, "second", __func__,0);

  // call matrix pointer version of function
  mat_add_t C_ptr = matrix_add(A_ptr, B_ptr);
  int64_t C_mID = get_matrixID_from_ptr(C_ptr);
  return C_mID;
}

mat_add_t
matrix_add(mat_add_t A_ptr, mat_add_t B_ptr)
{
  // EXCEPTION CHECKING
  if ((A_ptr == MATRIX_PTR_INVALID) ||  (B_ptr == MATRIX_PTR_INVALID)){
    printf("ERROR: In %s an argument was not a valid matrix\n", __func__);
    exit(1);
  }
  if (matrix_row_level(A_ptr) != matrix_row_level(B_ptr)) {
    printf("ERROR: attempting to add matrices of different row_level (%d and %d)\n",
           matrix_row_level(A_ptr), matrix_row_level(B_ptr));
    printf("matrixIDs: for A = %ld, for B = %ld\n", A_ptr->matrixID,B_ptr->matrixID);
    exit(EXIT_FAILURE);
  }
  if (matrix_col_level(A_ptr) != matrix_col_level(B_ptr)) {
    printf("ERROR: attempting to add matrices of different col_level (%d and %d)\n",
           matrix_col_level(A_ptr), matrix_col_level(B_ptr));
    printf("matrixIDs: for A = %ld, for B = %ld\n", A_ptr->matrixID, B_ptr->matrixID);
    exit(EXIT_FAILURE);
  }

  // MATH RELATIONSHIP / IDENTITY SHORT CUTS
  // This part sorts A and B so that A < B and you take advantage of commutativity
  if (A_ptr > B_ptr) {
    mat_add_t temp;
    temp = A_ptr;
    A_ptr = B_ptr;
    B_ptr = temp;
  }
  if (matrix_is_zero(A_ptr)) {
    return B_ptr;
  }
  if (matrix_is_zero(B_ptr)) {
    return A_ptr;
  }


  // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  mat_add_t sum_ptr = op_get(SUM, A_ptr, B_ptr);
  if (sum_ptr == MATRIX_PTR_INVALID) {     // WAS NOT FOUND
    
    // FOR SCALAR CASE
    //   CALCULATE RESULT OF OPERATION 
    if (matrix_type(A_ptr)==SCALAR)
    {

      ScalarType sum = matrix_trace(A_ptr)+ matrix_trace(B_ptr);
      
      // FIND OR CREATE APPROPRIATE MATRIX PTR
      sum_ptr = matrix_get_ptr_scalar(sum);
    }
    else
    {
      
      // ALL NON SCALAR CASES HAVE PANELS
      mat_add_t panel[4];
      
      // CASE: ROW_VECTOR
      if (matrix_type(A_ptr)==ROW_VECTOR)
	{
	  // panels 2 and 3 are MATRIX_PTR_INVALID
	  panel[0] = matrix_add(matrix_sub(A_ptr, 0), matrix_sub(B_ptr, 0));
	  panel[1] = matrix_add(matrix_sub(A_ptr, 1), matrix_sub(B_ptr, 1));
	  panel[2] = panel[3] = MATRIX_PTR_INVALID;
	}
      
      // CASE: COL_VECTOR
      if (matrix_type(A_ptr)==COL_VECTOR)
	{
	  // panels 1 and 3 are MATRIX_PTR_INVALID
	  panel[0] = matrix_add(matrix_sub(A_ptr, 0), matrix_sub(B_ptr, 0));
	  panel[2] = matrix_add(matrix_sub(A_ptr, 2), matrix_sub(B_ptr, 2));
	  panel[1] = panel[3] = MATRIX_PTR_INVALID;
	}
      
      // CASE: MATRIX
      if (matrix_type(A_ptr)==MATRIX)
	{
	  for (int i = 0; i < 4; i++) 
	    {
	      panel[i] = matrix_add(matrix_sub(A_ptr, i), matrix_sub(B_ptr, i));
	    }
	}
      
      // FIND OR CREATE APPROPRIATE MATRIX PTR
      sum_ptr = matrix_get_ptr_panel(panel,
                  matrix_row_level(A_ptr), matrix_col_level(A_ptr));
    }

    // STORE RESULT IN OPERATIONS STORE 
    op_set(SUM, A_ptr, B_ptr, sum_ptr);
    
  }   // end matrix was not found in op store

  // RETURN VALUE
  return sum_ptr;
}

/* Python Interface for matrix_diff(mat_add_t A_ptr, mat_add_t B_ptr):
   accepts matrixIDs for matrices to be subtracted and converts
   inputs to pointers before calling. */

int64_t 
matrix_diff_matrixID(int64_t A_mID, int64_t B_mID) {

  // get the matrix pointers from the matrixID, and see if still in store 
  mat_add_t A_ptr = mat_ptr_from_matrixID(A_mID, "first", __func__,0);
  mat_add_t B_ptr = mat_ptr_from_matrixID(B_mID, "second", __func__,0);

  // call matrix pointer version of function
  mat_add_t C_ptr = matrix_diff(A_ptr, B_ptr);
  int64_t C_mID = get_matrixID_from_ptr(C_ptr);
  return C_mID;
}

mat_add_t
matrix_diff(mat_add_t A_ptr, mat_add_t B_ptr)
{
  // EXCEPTION CHECKING
  if ((A_ptr == MATRIX_PTR_INVALID) ||  (B_ptr == MATRIX_PTR_INVALID)){
    printf("ERROR: In %s an argument was not a valid matrix\n", __func__);
    exit(1);
  }
  if (matrix_row_level(A_ptr) != matrix_row_level(B_ptr)) {
    printf("ERROR: attempting to subtract matrices of different row_level (%d and %d)\n",
           matrix_row_level(A_ptr), matrix_row_level(B_ptr));
    printf("matrixIDs: for A = %ld, for B = %ld\n", A_ptr->matrixID, B_ptr->matrixID);
    exit(EXIT_FAILURE);
  }
  if (matrix_col_level(A_ptr) != matrix_col_level(B_ptr)) {
    printf("ERROR: attempting to subtract matrices of different col_level (%d and %d)\n",
           matrix_col_level(A_ptr), matrix_col_level(B_ptr));
    printf("matrixIDs: for A = %ld, for B = %ld\n", A_ptr->matrixID, B_ptr->matrixID);
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
  mat_add_t diff_ptr = op_get(DIFF, A_ptr, B_ptr);
  if (diff_ptr == MATRIX_PTR_INVALID) {     // WAS NOT FOUND
    
    // FOR SCALAR CASE
    //   CALCULATE RESULT OF OPERATION 
    if (matrix_type(A_ptr)==SCALAR)
    {

      ScalarType diff = matrix_trace(A_ptr) - matrix_trace(B_ptr);
      
      // FIND OR CREATE APPROPRIATE MATRIX PTR
      diff_ptr = matrix_get_ptr_scalar(diff);
    }
    else
    {
      
      // ALL NON SCALAR CASES HAVE PANELS
      mat_add_t panel[4];
      
      // CASE: ROW_VECTOR
      if (matrix_type(A_ptr)==ROW_VECTOR)
	{
	  // panels 2 and 3 are MATRIX_PTR_INVALID
	  panel[0] = matrix_diff(matrix_sub(A_ptr, 0), matrix_sub(B_ptr, 0));
	  panel[1] = matrix_diff(matrix_sub(A_ptr, 1), matrix_sub(B_ptr, 1));
	  panel[2] = panel[3] = MATRIX_PTR_INVALID;
	}
      
      // CASE: COL_VECTOR
      if (matrix_type(A_ptr)==COL_VECTOR)
	{
	  // panels 1 and 3 are MATRIX_PTR_INVALID
	  panel[0] = matrix_diff(matrix_sub(A_ptr, 0), matrix_sub(B_ptr, 0));
	  panel[2] = matrix_diff(matrix_sub(A_ptr, 2), matrix_sub(B_ptr, 2));
	  panel[1] = panel[3] = MATRIX_PTR_INVALID;
	}
      
      // CASE: MATRIX
      if (matrix_type(A_ptr)==MATRIX)
	{
	  for (int i = 0; i < 4; i++) 
	    {
	      panel[i] = matrix_diff(matrix_sub(A_ptr, i), matrix_sub(B_ptr, i));
	    }
	}
      
      // FIND OR CREATE APPROPRIATE MATRIX PTR
      diff_ptr = matrix_get_ptr_panel(panel,
                  matrix_row_level(A_ptr), matrix_col_level(A_ptr));
    }

    // STORE RESULT IN OPERATIONS STORE 
    op_set(DIFF, A_ptr, B_ptr, diff_ptr);
    
  }   // end matrix was not found in op store

  // RETURN VALUE
  return diff_ptr;
}

/* Python Interface for matrix_mult(mat_add_t A_ptr, mat_add_t B_ptr):
   accepts matrixIDs for matrices to be multiplied and converts
   inputs to pointers before calling. */

int64_t 
matrix_mult_matrixID(int64_t A_mID, int64_t B_mID) {

  // get the matrix pointers from the matrixIDs, and see if still in store 
  mat_add_t A_ptr = mat_ptr_from_matrixID(A_mID, "first", __func__,0);
  mat_add_t B_ptr = mat_ptr_from_matrixID(B_mID, "second", __func__,0);
  // the following check is performed inside matrix_mult():
  // if ((A_ptr == MATRIX_PTR_INVALID) || (B_ptr == MATRIX_PTR_INVALID)) { exit(1); }

  // call matrix pointer version of function
  mat_add_t C_ptr = matrix_mult(A_ptr, B_ptr);
  int64_t C_mID = get_matrixID_from_ptr(C_ptr);
  return C_mID;
}

mat_add_t
matrix_mult(mat_add_t A_ptr, mat_add_t B_ptr)
{
  int verbose = 0;

  // EXCEPTION CHECKING
  if ((A_ptr == MATRIX_PTR_INVALID) ||  (B_ptr == MATRIX_PTR_INVALID)){
    printf("ERROR: In %s an argument was not a valid matrix\n", __func__);
    exit(1);
  }
  
  if (matrix_col_level(A_ptr) != matrix_row_level(B_ptr)) {
    printf("ERROR: attempting to multiply matrix with col_level %d by matrix with row_level %d\n", 
           matrix_col_level(A_ptr), matrix_row_level(B_ptr));
    printf("matrixIDs: for A = %ld, for B = %ld\n", 
           A_ptr->matrixID, B_ptr->matrixID);

    if (verbose) {    
    matrix_store_report("stdout"); 
    op_store_report("stdout"); 
    rusage_report(0,"stdout");
    }
    //hash_report(store.hash_table, stdout, 0);  or hash_report(hash_table_t *table, FILE *fp, int verbose)
    exit(EXIT_FAILURE);
  }
  
  matrix_type_t mat_type_A = matrix_type(A_ptr);
  matrix_type_t mat_type_B = matrix_type(B_ptr);

  if (verbose) {
    // DEBUGGING
    printf("Matrix types for A and B are %d %d\n",mat_type_A, mat_type_B);
    printf("   Levels for A are %d %d\n",matrix_row_level(A_ptr), matrix_col_level(A_ptr));
    printf("   Levels for B are %d %d\n",matrix_row_level(B_ptr), matrix_col_level(B_ptr));
  }

  // MATH IDENTITY SHORT CUTS
  if (matrix_is_id(A_ptr)) {
    return B_ptr;
  }
  if (matrix_is_id(B_ptr)) {
    return A_ptr;
  }
  if (matrix_is_zero(A_ptr) || matrix_is_zero(B_ptr)) {
    return get_zero_matrix_ptr(matrix_row_level(A_ptr),matrix_col_level(B_ptr));
  }
  // If A and/or B is SCALAR, then problem reduces to scalar_mult
  // and the operation is stored in the KRONECKER op store.
  if (mat_type_A == SCALAR) {
    
    if (verbose) {
      // DEBUGGING
      printf("In %s calling scalar mult A scalar\n", __func__);
    }
    
    return scalar_mult(A_ptr, B_ptr);
  }
  if (mat_type_B == SCALAR) {

    if (verbose) {
      // DEBUGGING
      printf("In %s calling scalar mult B scalar\n", __func__);
    }

    return scalar_mult(B_ptr, A_ptr);
  }
  
  
  // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  mat_add_t product_ptr = op_get(PRODUCT, A_ptr, B_ptr);
  if (product_ptr == MATRIX_PTR_INVALID)     // WAS NOT FOUND
    {

      // DOT PRODUCT: A is ROW_VECTOR and B is COL_VECTOR and result is SCALAR
      //   FIND OR CREATE APPROPRIATE MATRIX PTR
      if ((mat_type_A == ROW_VECTOR) && (mat_type_B == COL_VECTOR)) 
	{
	  product_ptr = matrix_add(matrix_mult(matrix_sub(A_ptr,0),matrix_sub(B_ptr,0)),
				   matrix_mult(matrix_sub(A_ptr,1),matrix_sub(B_ptr,2)));
	}
      else
	{
	  
	  //  For all remaining cases we use a panel
	  mat_add_t panel[4];
	  
	  // CASE: ROW_VECTOR times MATRIX
	  if ((mat_type_A == ROW_VECTOR) && (mat_type_B == MATRIX)) 
	    {
	      panel[0] = matrix_add(matrix_mult(matrix_sub(A_ptr, 0),matrix_sub(B_ptr, 0)),
				    matrix_mult(matrix_sub(A_ptr, 1),matrix_sub(B_ptr, 2)));
	      panel[1] = matrix_add(matrix_mult(matrix_sub(A_ptr, 0),matrix_sub(B_ptr, 1)),
				    matrix_mult(matrix_sub(A_ptr, 1),matrix_sub(B_ptr, 3)));
	      panel[2] = MATRIX_PTR_INVALID;
	      panel[3] = MATRIX_PTR_INVALID;
	    }
	  
	  // CASE: COL_VECTOR times ROW_VECTOR  (there is no COL_VECTOR times MATRIX)
	  if (mat_type_A == COL_VECTOR) // in this case B is a ROW_VECTOR
	    {
	      panel[0] = matrix_mult(matrix_sub(A_ptr, 0),matrix_sub(B_ptr, 0));
	      panel[1] = matrix_mult(matrix_sub(A_ptr, 0),matrix_sub(B_ptr, 1));
	      panel[2] = matrix_mult(matrix_sub(A_ptr, 2),matrix_sub(B_ptr, 0));
	      panel[3] = matrix_mult(matrix_sub(A_ptr, 2),matrix_sub(B_ptr, 1));
	    }
	  
	  // CASE: MATRIX times COL_VECTOR
	  if ((mat_type_A == MATRIX) && (mat_type_B == COL_VECTOR)) 
	    {
	      panel[0] = matrix_add(matrix_mult(matrix_sub(A_ptr, 0), matrix_sub(B_ptr, 0)),
				    matrix_mult(matrix_sub(A_ptr, 1), matrix_sub(B_ptr, 2)));
	      panel[1] = MATRIX_PTR_INVALID;
	      panel[2] = matrix_add(matrix_mult(matrix_sub(A_ptr, 2), matrix_sub(B_ptr, 0)),
				    matrix_mult(matrix_sub(A_ptr, 3), matrix_sub(B_ptr, 2)));
	      panel[3] = MATRIX_PTR_INVALID;
	    }
	  
	  // CASE: MATRIX times MATRIX
	  if ((mat_type_A == MATRIX) && (mat_type_B == MATRIX)) 
	    {
	      panel[0] = matrix_add(matrix_mult(matrix_sub(A_ptr, 0), matrix_sub(B_ptr, 0)),
				    matrix_mult(matrix_sub(A_ptr, 1), matrix_sub(B_ptr, 2)));
	      panel[1] = matrix_add(matrix_mult(matrix_sub(A_ptr, 0), matrix_sub(B_ptr, 1)),
				    matrix_mult(matrix_sub(A_ptr, 1), matrix_sub(B_ptr, 3)));
	      panel[2] = matrix_add(matrix_mult(matrix_sub(A_ptr, 2), matrix_sub(B_ptr, 0)),
				    matrix_mult(matrix_sub(A_ptr, 3), matrix_sub(B_ptr, 2)));
	      panel[3] = matrix_add(matrix_mult(matrix_sub(A_ptr, 2), matrix_sub(B_ptr, 1)),
				    matrix_mult(matrix_sub(A_ptr, 3), matrix_sub(B_ptr, 3)));
	    }
	  
	  //   FIND OR CREATE APPROPRIATE MATRIX PTR
	  product_ptr = matrix_get_ptr_panel(panel, 
                matrix_row_level(A_ptr), matrix_col_level(B_ptr));

          // We missed a case????????
	  if (product_ptr == MATRIX_PTR_INVALID) {
	    printf ("ERROR: something went very wrong in %s\n", __func__);
	    exit(1);
	  }

	}
      
      // STORE RESULT IN OPERATIONS STORE 
      op_set(PRODUCT, A_ptr, B_ptr, product_ptr);
    }      // end matrix was not found in op store
  
  // RETURN VALUE
  return product_ptr;
}

/* Python Interface for scalar_mult(mat_add_t A_ptr, mat_add_t B_ptr):
   accepts matrixIDs for matrices to be added and converts
   inputs to pointers before calling. */

int64_t scalar_mult_matrixID(int64_t A_mID, int64_t B_mID) {

  // get the matrix pointers from the matrixIDs, and see if still in store 
  mat_add_t A_ptr = mat_ptr_from_matrixID(A_mID, "first", __func__,0);
  mat_add_t B_ptr = mat_ptr_from_matrixID(B_mID, "second", __func__,0);
  // the following check is performed inside scalar_mult():
  // if ((A_ptr == MATRIX_PTR_INVALID) || (B_ptr == MATRIX_PTR_INVALID)) { exit(1); }

  // call matrix pointer version of function
  mat_add_t C_ptr = scalar_mult(A_ptr, B_ptr);
  int64_t C_mID = get_matrixID_from_ptr(C_ptr);
  return C_mID;
}

mat_add_t
scalar_mult(mat_add_t A_ptr, mat_add_t B_ptr)
{
  // EXCEPTION CHECKING
  if ((A_ptr == MATRIX_PTR_INVALID) ||  (B_ptr == MATRIX_PTR_INVALID)){
    printf("ERROR: In %s an argument was not a valid matrix\n", __func__);
    exit(1);
  }
  if (matrix_type(A_ptr) != SCALAR) {
    printf("ERROR: attempted scalar multiply with nonscalar of type %d\n", matrix_type(A_ptr));
    printf("matrixIDs: for A = %ld, for B = %ld\n", A_ptr->matrixID, B_ptr->matrixID);
    exit(EXIT_FAILURE);
  }
  
  matrix_type_t mat_type_B = matrix_type(B_ptr);

  // MATH IDENTITY SHORT CUTS
  if (matrix_is_zero(A_ptr) || matrix_is_zero(B_ptr)) {
    return get_zero_matrix_ptr(matrix_row_level(B_ptr),matrix_col_level(B_ptr));
  }
  
  // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  mat_add_t product_ptr = op_get(KRONECKER, A_ptr, B_ptr);
  if (product_ptr == MATRIX_PTR_INVALID) {   // WAS NOT FOUND
    
    
    // CALCULATE RESULT OF OPERATION 
    // CASE: B SCALAR
    if (mat_type_B == SCALAR) {
      
      ScalarType product = matrix_trace(A_ptr) * matrix_trace(B_ptr);

      product_ptr = matrix_get_ptr_scalar(product);
    }
    else
      {
	
	// non SCALAR all have panels
	mat_add_t panel[4];
	
	// CASE: B is ROW_VECTOR
	if (mat_type_B == ROW_VECTOR) {
	  panel[0] = scalar_mult(A_ptr, matrix_sub(B_ptr, 0));
	  panel[1] = scalar_mult(A_ptr, matrix_sub(B_ptr, 1));
	  panel[2] = MATRIX_PTR_INVALID;
	  panel[3] = MATRIX_PTR_INVALID;
	}
	
	// CASE: B is COL_VECTOR
	if (mat_type_B == COL_VECTOR) {
	  panel[0] = scalar_mult(A_ptr, matrix_sub(B_ptr, 0));
	  panel[1] = MATRIX_PTR_INVALID;
	  panel[2] = scalar_mult(A_ptr, matrix_sub(B_ptr, 2));
	  panel[3] = MATRIX_PTR_INVALID;
	}
	
	// CASE: B is MATRIX
	if (mat_type_B == MATRIX) {
	  for (int i = 0; i < 4; i++) {
	    panel[i] = scalar_mult(A_ptr, matrix_sub(B_ptr, i));
	  }
	}

      //   FIND OR CREATE APPROPRIATE MATRIX ID
      product_ptr = matrix_get_ptr_panel(panel, 
                matrix_row_level(B_ptr), matrix_col_level(B_ptr));

      }
    
    //   STORE RESULT IN OPERATIONS STORE 
    op_set(KRONECKER, A_ptr, B_ptr, product_ptr);
  }
  
  //   RETURN VALUE
  return product_ptr;
}



/* rewritten by Jenny Z and Steve C on 2016-02-11  */ 
/* Fast routine for calculating Kronecker product of 
   matrix A with an identity matrix I of level levelI */
static mat_add_t
tensor_with_identity_on_left(mat_add_t A_ptr, mat_level_t levelI)
{

  // MATH IDENTITY SHORT CUTS
  //   If Identity matrix I is SCALAR, then return A
  if (levelI == 0) { return A_ptr; }

  // CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  mat_add_t identity_ptr = get_identity_matrix_ptr(levelI);
  mat_add_t product_ptr = op_get(KRONECKER, identity_ptr, A_ptr);

  if (product_ptr == MATRIX_PTR_INVALID) {   // WAS NOT FOUND
    // CALCULATE RESULT OF OPERATION 
    
    // ALL CASES:
    //    calculate Kronecker product recursively
    mat_level_t row_levelA = matrix_row_level(A_ptr);
    mat_level_t col_levelA = matrix_col_level(A_ptr);
    
    mat_add_t panel[4] = {A_ptr, 
			  get_zero_matrix_ptr(row_levelA,col_levelA), 
			  get_zero_matrix_ptr(row_levelA,col_levelA),
			  A_ptr};
    
    product_ptr = tensor_with_identity_on_left(
	     matrix_get_ptr_panel(panel, row_levelA+1, col_levelA+1), levelI-1);
  }   // end create new matrix
  
  // STORE RESULT IN OPERATIONS STORE 
  op_set(KRONECKER, identity_ptr, A_ptr, product_ptr);
  
  // RETURN VALUE
  return product_ptr;
}

/* Python Interface for kronecker_product(mat_add_t A_ptr, mat_add_t B_ptr):
   accepts matrixIDs for matrices to be added and converts
   inputs to pointers before calling. */

int64_t 
kronecker_product_matrixID(int64_t A_mID, int64_t B_mID) {

  // get the matrix pointers from the matrixIDs, and see if still in store 
  mat_add_t A_ptr = mat_ptr_from_matrixID(A_mID, "first", __func__,0);
  mat_add_t B_ptr = mat_ptr_from_matrixID(B_mID, "second", __func__,0);
  // the following check is performed inside kronecker_product():
  // if ((A_ptr == MATRIX_PTR_INVALID) || (B_ptr == MATRIX_PTR_INVALID)) { exit(1); }

  // call matrix pointer version of function
  mat_add_t C_ptr = kronecker_product(A_ptr, B_ptr);
  int64_t C_mID = get_matrixID_from_ptr(C_ptr);
  return C_mID;
}

mat_add_t
kronecker_product(mat_add_t A_ptr, mat_add_t B_ptr)
{

  // EXCEPTIONS (more below)
  if ((A_ptr == MATRIX_PTR_INVALID) ||  (B_ptr == MATRIX_PTR_INVALID)){
    exit(1);
  }

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
  // When either A or B is scalar then the kronecker product is the scalar product
  //   (The result will be stored in the op_store KRONECKER)
  if (mat_type_A == SCALAR) {
    return scalar_mult(A_ptr, B_ptr);
  }
  if (mat_type_B == SCALAR) {
    return scalar_mult(B_ptr, A_ptr);
  }


  //   CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  mat_add_t product_ptr = op_get(KRONECKER, A_ptr, B_ptr);
  if (product_ptr == MATRIX_PTR_INVALID) {   // WAS NOT FOUND
    
    // CALCULATE RESULT OF OPERATION 
    mat_add_t panel[4];
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
    
    // CASE: A is ROW_VECTOR and B is MATRIX (fix by matt calef)
    // former version had join outside of kronecker product
    // [a_0,a_1] \otimes [[b_0,b_1],[b_2,b_3]] ->
    // [[a_0 join(b_0,b_1),a_1 join(b_0,b_1)],
    //    [a_0 join(b_2,b_3),a_1 join(b_2,b_3)]
    if ((mat_type_A == ROW_VECTOR) && (mat_type_B == MATRIX)) {
      panel[0] = kronecker_product(matrix_sub(A_ptr,0),join(matrix_sub(B_ptr,0),
                matrix_sub(B_ptr,1)));
      panel[1] = kronecker_product(matrix_sub(A_ptr,1),join(matrix_sub(B_ptr,0),
                matrix_sub(B_ptr,1)));
      panel[2] = kronecker_product(matrix_sub(A_ptr,0),join(matrix_sub(B_ptr,2),
                matrix_sub(B_ptr,3)));
      panel[3] = kronecker_product(matrix_sub(A_ptr,1),join(matrix_sub(B_ptr,2),                matrix_sub(B_ptr,3)));
    }
    
    // CASE: A is COL_VECTOR and B is MATRIX
    // (matt's fix is to use the adjoint - we need to think what to do instead)
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

#if 0    
    // CASE: A is MATRIX and B is ROW_VECTOR (simplification by matt calef)
    // (formerly broken up by joins, unnecessarily)
    if ((mat_type_A == MATRIX) && (mat_type_B == ROW_VECTOR)) {
      for (int i = 0; i < 4 ; i++) {
	panel[i] = kronecker_product(matrix_sub(A_ptr, i), B_ptr);
      }
    }
    
    // CASE: A is MATRIX and B is COL_VECTOR (simplification by matt calef)
    // (formerly broken up by stacks, unnecessarily)
    if ((mat_type_A == MATRIX) && (mat_type_B == COL_VECTOR)) {
      for (int i = 0; i < 4 ; i++) {
	panel[i] = kronecker_product(matrix_sub(A_ptr, i), B_ptr);
      }
    }
    
    // CASE: A and B are both MATRIX
    if ((mat_type_A == MATRIX) && (mat_type_B == MATRIX)) {
      for (int i = 0; i < 4 ; i++) {
	panel[i] = kronecker_product(matrix_sub(A_ptr, i), B_ptr);
      }
    }
#endif

    // find or create the matrix id of the result
    product_ptr = matrix_get_ptr_panel(panel, new_row_level, new_col_level);
    
    // ERROR CHECKING (since there are so many cases)
    if (product_ptr == MATRIX_PTR_INVALID) {
      printf("ERROR: in %s something went bad\n", __func__);

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

/* Python Interface for join(mat_add_t A_ptr, mat_add_t B_ptr):
   accepts matrixIDs for matrices to be added and converts
   inputs to pointers before calling. */

int64_t 
join_matrixID(int64_t A_mID, int64_t B_mID) {

  // get the matrix pointers from the matrixIDs, and see if still in store 
  mat_add_t A_ptr = mat_ptr_from_matrixID(A_mID, "first", __func__,0);
  mat_add_t B_ptr = mat_ptr_from_matrixID(B_mID, "second", __func__,0);
  // the following check is performed inside join():
  // if ((A_ptr == MATRIX_PTR_INVALID) || (B_ptr == MATRIX_PTR_INVALID)) { exit(1); }

  // call matrix pointer version of function
  mat_add_t C_ptr = join(A_ptr, B_ptr);
  int64_t C_mID = get_matrixID_from_ptr(C_ptr);
  return C_mID;
}

mat_add_t join(mat_add_t A_ptr, mat_add_t B_ptr) {

  // EXCEPTION CHECKING
  if ((A_ptr == MATRIX_PTR_INVALID) ||  (B_ptr == MATRIX_PTR_INVALID)){
    printf("ERROR: In %s an argument was not a valid matrix\n",__func__);
    exit(1);
  }

  if ( (matrix_row_level(A_ptr) != matrix_row_level(B_ptr)) ||
       (matrix_col_level(A_ptr) != matrix_col_level(B_ptr)) ) {
    printf("In function %s, input matrices were not the same size.\n",__func__);
    exit(1);
  }

  matrix_type_t mat_type_A = matrix_type(A_ptr);


  //   CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  mat_add_t out_ptr = op_get(JOIN, A_ptr, B_ptr);
  if (out_ptr == MATRIX_PTR_INVALID) {   // WAS NOT FOUND
    
    //   CALCULATE RESULT OF OPERATION
    mat_add_t panel[4];
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
    out_ptr = matrix_get_ptr_panel(panel,row_level,new_col_level);
    
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

/* Python Interface for stack(mat_add_t A_ptr, mat_add_t B_ptr):
   accepts matrixIDs for matrices to be added and converts
   inputs to pointers before calling. */

int64_t 
stack_matrixID(int64_t A_mID, int64_t B_mID) {

  // get the matrix pointers from the matrixIDs, and see if still in store 
  mat_add_t A_ptr = mat_ptr_from_matrixID(A_mID, "first", __func__,0);
  mat_add_t B_ptr = mat_ptr_from_matrixID(B_mID, "second", __func__,0);
  // the following check is performed inside stack():
  // if ((A_ptr == MATRIX_PTR_INVALID) || (B_ptr == MATRIX_PTR_INVALID)) { exit(1); }

  // call matrix pointer version of function
  mat_add_t C_ptr = stack(A_ptr, B_ptr);
  int64_t C_mID = get_matrixID_from_ptr(C_ptr);
  return C_mID;
}

mat_add_t stack(mat_add_t A_ptr, mat_add_t B_ptr) {

  // EXCEPTION CHECKING
  if ((A_ptr == MATRIX_PTR_INVALID) ||  (B_ptr == MATRIX_PTR_INVALID)){
    printf("ERROR: In %s an argument was not a valid matrix\n", __func__);
    exit(1);
  }  
  if ( (matrix_row_level(A_ptr) != matrix_row_level(B_ptr)) ||
       (matrix_col_level(A_ptr) != matrix_col_level(B_ptr)) ) {
    printf("In function %s, input matrices were not the same size.\n",__func__);
    exit(1);
  }

  matrix_type_t mat_type_A = matrix_type(A_ptr);

  //   CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  mat_add_t out_ptr = op_get(STACK, A_ptr, B_ptr);
  if (out_ptr == MATRIX_PTR_INVALID) {  // WAS NOT FOUND
    
    //   CALCULATE RESULT OF OPERATION
    mat_add_t panel[4];
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
    out_ptr = matrix_get_ptr_panel(panel,new_row_level,col_level);
    
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

/* Python Interface for matrix_entrySquared(mat_add_t A_ptr, double scale_factor):
   accepts matrixIDs for matrices to be added and converts
   inputs to pointers before calling. */

int64_t 
matrix_entrySquared_matrixID(int64_t m_mID, double scale_factor) {

  // get the matrix pointer from the matrixIDs, and see if still in store 
  mat_add_t m_ptr = mat_ptr_from_matrixID(m_mID, "first", __func__,0);
  // the following check is performed inside matrix_entrySquared():
  // if (m_ptr == MATRIX_PTR_INVALID) { exit(1); }

  // put scale_factor into the matrix store and get its pointer *here*
  // so that we don't have to do it an exponential number of times
  mat_add_t scale_ptr = matrix_get_ptr_scalar((ScalarType)scale_factor);

  // call matrix pointer version of function
  mat_add_t C_ptr = matrix_entrySquared(m_ptr, scale_factor, scale_ptr);
  int64_t C_mID = get_matrixID_from_ptr(C_ptr);
  return C_mID;
}

mat_add_t matrix_entrySquared(mat_add_t m_ptr, double scale_factor,
        mat_add_t scale_ptr)
{
  // EXCEPTION CHECKING
  if (m_ptr == MATRIX_PTR_INVALID) {
    printf("ERROR: In %s an argument was not a valid matrix\n",__func__);
    exit(1);
  }  
  //   this routine is only designed to work on square matrices
  //   (so no submatrix will be MATRIX_PTR_INVALID)
  if (matrix_row_level(m_ptr) != matrix_col_level(m_ptr))  
    {
      printf("Function %s requires a square matrix\n",__func__);
      exit(1);
    }
  
  // MATH IDENTITY SHORT CUTS
  // we may know squaring elements will leave matrix unchanged
  //if (m->iszero || (m->isid && scale_factor == 1.0)) { return m_ptr; } 
  if (matrix_is_zero(m_ptr)) return m_ptr;
  if (matrix_is_id(m_ptr)) return scalar_mult(scale_ptr, m_ptr);

  //   CHECK TO SEE IF SOLUTION IN OPERATIONS STORE,
  // mat_add_t scale_ptr = matrix_get_ptr_scalar(scale_factor);
  mat_add_t mSq_ptr = op_get(ENTRYSQUARE, m_ptr, scale_ptr);
  if (mSq_ptr == MATRIX_PTR_INVALID) {  // WAS NOT FOUND


    // CASE: SCALAR
    if (matrix_type(m_ptr) == SCALAR) 
      { 
	ScalarType element;
        //printf("m_ptr is %p, scale_factor is %g\n",m_ptr,scale_factor);
#ifdef USE_COMPLEX 
	element = conj(matrix_trace(m_ptr)) * matrix_trace(m_ptr); 
	//printf("squaring: started as %.25g+%.25gI, ended as %.25g+%.25gI\n", 
	//       creal(m_ptr->trace_element),cimag(m_ptr->trace_element),
	//       creal(element), cimag(element); 
#else
	element = matrix_trace(m_ptr) * matrix_trace(m_ptr); 
	//printf("squaring: started as %.25g, ended as %.25g\n", 
	//       m_ptr->trace_element, element);
#endif    
	element *= scale_factor;
#ifdef USE_COMPLEX 
	//printf("scaling: ended as %.25g+%.25gI\n", 
	//       creal(element),cimag(element));
#else
	//printf("scaling: ended as %.25g\n", element);
#endif    
	mSq_ptr = matrix_get_ptr_scalar(element);
      }

    // CASE: NON-SCALAR
    else
      { 
	mat_add_t new_sub[4];
	new_sub[0] = matrix_entrySquared(matrix_sub(m_ptr,0),scale_factor,scale_ptr);
	new_sub[1] = matrix_entrySquared(matrix_sub(m_ptr,1),scale_factor,scale_ptr);
	new_sub[2] = matrix_entrySquared(matrix_sub(m_ptr,2),scale_factor,scale_ptr);
	new_sub[3] = matrix_entrySquared(matrix_sub(m_ptr,3),scale_factor,scale_ptr);
	mSq_ptr = matrix_get_ptr_panel(new_sub,
                matrix_row_level(m_ptr),matrix_col_level(m_ptr));
      }
    
    //   STORE RESULT IN OPERATIONS STORE 
    op_set(ENTRYSQUARE, m_ptr, scale_ptr, mSq_ptr);
  }
  
  //   RETURN VALUE
  return mSq_ptr;
}

int64_t iHadamard_times_matrix_matrixID(int64_t A_mID) {
  mat_add_t A_ptr = mat_ptr_from_matrixID(A_mID, "first", __func__, 0);
  mat_add_t result_ptr = iHadamard_times_matrix(A_ptr);
  return get_matrixID_from_ptr(result_ptr);
}
// We have decided for now that default Hadamard will be the integer
// (unnormalized) version
mat_add_t iHadamard_times_matrix(mat_add_t A_ptr) {

  // EXCEPTION CHECKING
  if (A_ptr == MATRIX_PTR_INVALID) {
    printf("ERROR in %s:  has an argument that was not a valid matrix\n",
                __func__);
    exit(1);
  }

  // routine takes the matrix index of matrix A and returns matrix index of iHadamard * A
  
  mat_level_t row_level = matrix_row_level(A_ptr);
  mat_level_t col_level = matrix_col_level(A_ptr);

  // find the index of the iHadamard matrix of this level
  mat_add_t iHad_ptr = get_iHadamard_matrix_ptr(row_level);

  mat_add_t iHad_times_A_ptr = op_get(PRODUCT,iHad_ptr,A_ptr);

  if (row_level == 0) {
    printf("ERROR: in %s, there is no sensible definition for level 0 Hadamard\n",__func__);
    exit(1);
  }

  // return memoized (stored) value of iHad * A if it exists
  if (iHad_times_A_ptr != 0) { return iHad_times_A_ptr; }

  // if no stored value and the matrix level is 1, calculate the product HH1 * A
  // using the preloaded value for the integer iHadamard
  mat_add_t matptr_HH1 = get_iHadamard_matrix_ptr(1);
  mat_add_t matptr_scalarM1 = matrix_get_ptr_scalar((ScalarType)(-1));
  if (row_level == 1) return matrix_mult(matptr_HH1,A_ptr);

  // if no stored value and the matrix larger than 2x2, calculate the iHadamard
  // product recursively as follows (this should use one less multiply per
  // submatrix than matrix_mult(iHad_ptr,A_ptr)
  mat_add_t submatrix[4];
  submatrix[0] = iHadamard_times_matrix(matrix_add(matrix_sub(A_ptr,0),matrix_sub(A_ptr,2)));
  submatrix[1] = iHadamard_times_matrix(matrix_add(matrix_sub(A_ptr,1),matrix_sub(A_ptr,3)));
  submatrix[2] = iHadamard_times_matrix(matrix_add(matrix_sub(A_ptr,0),scalar_mult(matptr_scalarM1,matrix_sub(A_ptr,2))));
  submatrix[3] = iHadamard_times_matrix(matrix_add(matrix_sub(A_ptr,1),scalar_mult(matptr_scalarM1,matrix_sub(A_ptr,3))));
  iHad_times_A_ptr = matrix_get_ptr_panel(submatrix,row_level,col_level);

  // memoize the value of negative A
  op_set(PRODUCT,iHad_ptr,A_ptr,iHad_times_A_ptr);

  return iHad_times_A_ptr;
}

int64_t matrix_times_iHadamard_matrixID(int64_t A_mID) {
  mat_add_t A_ptr = mat_ptr_from_matrixID(A_mID, "first", __func__, 0);
  mat_add_t result_ptr = matrix_times_iHadamard(A_ptr);
  return get_matrixID_from_ptr(result_ptr);
}

// We have decided for now that default iHadmard will be the integer (unnormalized) version
mat_add_t matrix_times_iHadamard(mat_add_t A_ptr) {

  // EXCEPTION CHECKING
  if (A_ptr == MATRIX_PTR_INVALID) {
    printf("ERROR: In %s has an argument that was not a valid matrix\n",
                __func__);
    exit(1);
  }

  mat_level_t row_level = matrix_row_level(A_ptr);
  mat_level_t col_level = matrix_col_level(A_ptr);

  // find the index of the iHadamard matrix of this level
  mat_add_t iHad_ptr = get_iHadamard_matrix_ptr(col_level);

  // routine takes the matrix index of matrix A and returns matrix index of A * iHadmard
  mat_add_t A_times_iHad_ptr = op_get(PRODUCT,A_ptr,iHad_ptr); 

  if (col_level == 0) {
    printf("ERROR: in %s, there is no sensible definition for level 0 Hadamard\n",__func__);
    exit(1);
  }
  
  // return memoized (stored) value of iHad * A if it exists
  if (A_times_iHad_ptr != 0) { return A_times_iHad_ptr; }
  
  mat_add_t matptr_HH1 = get_iHadamard_matrix_ptr(1);
  mat_add_t matptr_scalarM1 = matrix_get_ptr_scalar((ScalarType)(-1));
  // if no stored value and the matrix level is 1, calculate the product A * HH1
  // using the preloaded value for the integer iHadamard
  if (col_level == 1) return matrix_mult(A_ptr,matptr_HH1);

  // if no stored value and the matrix larger than 2x2,
  // calculate the iHadamard product recursively as follows
  // this should be cheaper than matrix_mult(A_ptr,iHad_ptr);
  mat_add_t  submatrix[4];
  submatrix[0] = matrix_times_iHadamard(matrix_add(matrix_sub(A_ptr,0),matrix_sub(A_ptr,1)));
  submatrix[1] = matrix_times_iHadamard(matrix_add(matrix_sub(A_ptr,0),scalar_mult(matptr_scalarM1,matrix_sub(A_ptr,1))));
  submatrix[2] = matrix_times_iHadamard(matrix_add(matrix_sub(A_ptr,2),matrix_sub(A_ptr,3)));
  submatrix[3] = matrix_times_iHadamard(matrix_add(matrix_sub(A_ptr,2),scalar_mult(matptr_scalarM1,matrix_sub(A_ptr,3))));
  A_times_iHad_ptr = matrix_get_ptr_panel(submatrix,row_level,col_level);

  // memoize the value of negative A
  op_set(PRODUCT,A_ptr,iHad_ptr,A_times_iHad_ptr);

  return(A_times_iHad_ptr);
} 

ScalarType matrix_tracenorm(mat_add_t mat_ptr, ScalarType scale)
{
  mat_add_t matadj_ptr = matrix_adjoint(mat_ptr);
  mat_add_t test_ptr = matrix_mult(matadj_ptr,mat_ptr);
  return scale*matrix_trace(test_ptr);
}


/* Python Interface for matrix_add(mat_add_t A_ptr, mat_add_t B_ptr):
   accepts matrixIDs for matrices to be added and converts
   inputs to pointers before calling. */
double matrix_maxnorm_matrixID(int64_t mat_mID){
  // get the matrix pointers from the matrixID, and see if still in store
  mat_add_t mat_ptr = mat_ptr_from_matrixID(mat_mID, "first", __func__, 0);

  // call matrix pointer version of function
  return matrix_maxnorm(mat_ptr);
}

double matrix_maxnorm(mat_add_t mat_ptr)
{
    // EXCEPTION CHECKING
    if (mat_ptr == MATRIX_PTR_INVALID){
        printf("ERROR: In %s an argument was not a valid matrix\n", __func__);
        exit(1);
    }

    // IDENTITY SHORT CUTS
    if (matrix_is_zero(mat_ptr)){
        return 0.0;
    }
    if (matrix_is_id(mat_ptr)){
        return 1.0;
    }

    double maxnorm = 0.0;

    // FOR SCALAR CASE
    // CALCULATE RESULT OF OPERATION
    if (matrix_type(mat_ptr)==SCALAR)
    {
        // NOTE: fabs is type generic macro (>=C99) <tgmath.h>. should we use?
#ifdef USE_COMPLEX
        // can change to cabsl if we need a long double
        return cabs(matrix_trace(mat_ptr));
#elif defined(USE_REAL)
        // can change to fabsl if we need a long double
        return fabs(matrix_trace(mat_ptr));
#else
        // no long double equivalent here. 
        return fabs(matrix_trace(mat_ptr));
#endif
    }

    else
    {
        // ALL NON SCALAR CASES HAVE PANELS
        for (int i = 0; i < 4; i++)
        {
            mat_add_t panel = matrix_sub(mat_ptr, i);
            if (panel != MATRIX_PTR_INVALID)
            {
                double temp = matrix_maxnorm(panel);
                if (maxnorm < temp)
                {
                    maxnorm = temp;
                }
            }
        }
    }

    // TODO: STORE RESULT IN OPERATIONS STORE?
    // I don't think we want to think of the result as a matrix ptr
    // e.g. in the integer case, norm is a double? 
    // (not an integer, because in the complex case, norm is a real value)

    // RETURN VALUE
    return maxnorm;
}



/* Python Interface for matrix_add(mat_add_t A_ptr, mat_add_t B_ptr):
   accepts matrixIDs for matrices to be added and converts
   inputs to pointers before calling. */
double matrix_sparsity_matrixID(int64_t mat_mID){
  // get the matrix pointers from the matrixID, and see if still in store
  mat_add_t mat_ptr = mat_ptr_from_matrixID(mat_mID, "first", __func__, 0);

  // call matrix pointer version of function
  return matrix_sparsity(mat_ptr);
}

double matrix_sparsity(mat_add_t mat_ptr)
{
    // EXCEPTION CHECKING
    if (mat_ptr == MATRIX_PTR_INVALID){
        printf("ERROR: In %s an argument was not a valid matrix\n", __func__);
        exit(1);
    }

    // IDENTITY SHORT CUTS
    if (matrix_is_zero(mat_ptr)){
        return 1.0;
    }

    double sparsity_sum = 0.0;

    // FOR SCALAR CASE
    // because of the zero matrix shortcut, the only possible scalar case at this point is 
    // that it is nonzero -> return sparsity =  0/1. 
    if (matrix_type(mat_ptr) == SCALAR){
        return 0.0;
    }

    else
    {
        // ALL NON SCALAR CASES HAVE PANELS
        for (int i = 0; i < 4; i++)
        {
            mat_add_t panel = matrix_sub(mat_ptr, i);
            if (panel != MATRIX_PTR_INVALID)
            {
                double ms = matrix_sparsity(panel);
                //sparsity_sum += (double) matrix_sparsity(panel);
                sparsity_sum += ms;
            }
        }

        // CASE: ROW_VECTOR OR COL_VECTOR
        if ((matrix_type(mat_ptr) == ROW_VECTOR) || (matrix_type(mat_ptr) == COL_VECTOR))
        {
            sparsity_sum = sparsity_sum / 2.0;
        }

        // CASE: MATRIX
        else
        {
            sparsity_sum = sparsity_sum / 4.0;
        }
    }
    // RETURN VALUE
    return (double) sparsity_sum;
}

          


        
        




