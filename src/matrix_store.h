//                   matrix_store.h
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


#ifndef LARC_MATRIX_STORE_H
#define LARC_MATRIX_STORE_H

#ifdef __cplusplus
extern "C" {
#endif // #ifdef __cplusplus

#include <stddef.h>      // size_t
#include <stdarg.h>      // .../va_list/va_args/
#include "larc.h"

/*************************************************************************
 *         INLINE FUNCTIONS acting on matrix id (mat_ptr_t)              *
 ************************************************************************/
// functions return enum matrix_type_t 
/*!
 * \ingroup larc
 * \brief A utility function which determins the matrix_type enum value given the dimensions of the matrix
 * \param row_level The log base 2 of the number of rows
 * \param col_level The log base 2 of the number of columns
 * \return SCALAR, ROW_VECTOR, COL_VECTOR or MATRIX
 */
inline static matrix_type_t matrix_type_from_levels(mat_level_t row_level, mat_level_t col_level)
{
  if ((row_level == 0) && (col_level == 0)) { return SCALAR;}
  if (row_level == 0) { return ROW_VECTOR;}
  if (col_level == 0) { return COL_VECTOR;}
  return MATRIX;
}


/***********************************************************************
 *               MACRO FUNCTIONS acting on mat_ptr_t                   *
 ***********************************************************************/
// This value for matrix ptr is also for the fake panels in row and col vectors
#define MATRIX_PTR_INVALID ((mat_ptr_t)0)
// This value for matrix id 
#define MATRIX_ID_INVALID ((int64_t)-1)
// Boolean functions
/*!
 * \ingroup larc
 * \brief Determine if a matrix is invalid
 * \param m_ptr A handle to the matrix we are testing
 * \return 1 if m_ptr is MATRIX_PTR_INVALID, 0 otherwise
 */
#define matrix_is_invalid(m_ptr) (m_ptr == MATRIX_PTR_INVALID)
/*!
 * \ingroup larc
 * \brief Determine if a matrix is an identity matrix
 * \param m_ptr A handle to the matrix we are testing
 * \return 1 if m_ptr is an identity matrix, 0 otherwise
 */
#define matrix_is_id(m_ptr) (m_ptr && m_ptr->isid)
/*!
 * \ingroup larc
 * \brief Determine if a matrix is all-zeros
 * \param m_ptr A handle to the matrix we are testing
 * \return 1 if m_ptr is an all-zero matrix, 0 otherwise
 */
#define matrix_is_zero(m_ptr) (m_ptr && m_ptr->iszero)
/*!
 * \ingroup larc
 * \brief Determine if a matrix is locked (can never be cleaned from the matrix store)
 * \param m_ptr A handle to the matrix we are testing
 * \return 1 if m_ptr is locked, 0 otherwise
 */
#define matrix_has_lock(m_ptr) (m_ptr && m_ptr->lock)
/*!
 * \ingroup larc
 * \brief Determine if a matrix is held (can not currently be cleaned from the matrix store)
 * \param m_ptr A handle to the matrix we are testing
 * \return 0 if the matrix is not held, else the number of holds on the matrix
 */
#define matrix_has_hold(m_ptr) (m_ptr && m_ptr->hold)
// Functions retrieving values inside a larc_matrix structure

/*!
 * \ingroup larc
 * \brief Return the trace of a matrix.
 * \param m_ptr A handle to the matrix for which we want the trace.
 * \return The trace of the matrix.
 */
#define matrix_trace(m_ptr) (m_ptr->trace_element)
/*!
 * \ingroup larc
 * \brief Determine the row dimension of a matrix
 * \param m_ptr A handle to the matrix
 * \return The log base 2 of the row dimension
 */
#define matrix_row_level(m_ptr) (m_ptr->row_level)
/*!
 * \ingroup larc
 * \brief Determine the column dimension of a matrix
 * \param m_ptr A handle to the matrix
 * \return The log base 2 of the column dimension
 */
#define matrix_col_level(m_ptr) (m_ptr->col_level)
/*!
 * \ingroup larc
 * \brief Finds the pointer to a quadrant submatrix of a given matrix
 * \param m_ptr A handle to the matrix for which we want a submatrix
 * \param s An integer in [0,3] denoting which quadrant submatrix is desired, with 0 denoting the * upper left, 1 the upper right, 2 the lower left and 3 the lower right
 * \return A pointer to the desired submatrix
 */
#define matrix_sub(m_ptr,s) (m_ptr->submatrix[s])
/*!
 * \ingroup larc
 * \brief Determines whether the matrix is needed for storing a larger matrix
 * \param m_ptr A handle to the matrix
 * \return The count of the number of times a matrix is a submatrix of a larger stored matrix
 */
#define matrix_appears_as_sub_count(m_ptr) (m_ptr->appears_as_sub_count)

  

/******************************************************************
 *         ACCESSORS acting on matrix PTR (mat_ptr_t)              *
 *****************************************************************/

/*!
 * \ingroup larc
 *
 * \brief Get a value within a matrix at the location specified by row and col.
 *
 * \param retScalar Pointer for returning value.
 * \param mPTR A handle to the matrix from which to retrieve a value.
 * \param row The row within the matrix.
 * \param col The column within the matrix.
 */
void get_valPTR_from_matPTR_and_coords(scalarType *retScalar, mat_ptr_t mPTR, int64_t row, int64_t col);

#ifndef SWIG
/*!
 * \ingroup larc
 * \brief Set a value within a matrix.
 * \param m   A handle to the matrix for which we set a value.
 * \param row The row within the matrix.
 * \param col The column within the matrix.
 * \param val The value to set.
 * \return A new matrix that is the same as the input matrix, except
 * with the value at the specified row and column set to the value v.
 */
mat_ptr_t get_matPTR_from_oldMatPTR_newVal_and_coords(mat_ptr_t m, int64_t row, int64_t col, scalarType val);
#endif

/*!
 * \ingroup larc
 * \brief A utility function returning the matrix_type for a matrix pointer
 * \param m_ptr A pointer to a matrix in the matrix store
 * \returns SCALAR, ROW_VECTOR, COL_VECTOR, or MATRIX
 */
inline static matrix_type_t matrix_type(mat_ptr_t m_ptr) 
{
  mat_level_t row_level = matrix_row_level(m_ptr);
  mat_level_t col_level = matrix_col_level(m_ptr);
  return matrix_type_from_levels(row_level,col_level);
}

/*!
 * \ingroup larc
 * \brief A utility function for finding a matrixID
 * \param m_ptr A pointer to a matrix in the matrix store
 * \return The matrixID for the given matrix pointer
 */
inline static int64_t get_matID_from_matPTR(mat_ptr_t m_ptr) {
  return ((m_ptr==MATRIX_PTR_INVALID) ? MATRIX_ID_INVALID : m_ptr->matrixID);
}

/*************************************************************************
 ************************************************************************/

/*!
 * \ingroup larc
 * \brief Marks a stored matrix as never to be cleaned from the matrix store
 * \param m_ptr A pointer to the stored matrix
 * \returns 1 on success, 0 if an error occurred
 */
int lock_matrix(mat_ptr_t m_ptr);
/*!
 * \ingroup larc
 * \brief Increments the hold value for a matrix in the matrix store (if above zero, the matrix will not be cleaned from the store)
 * \param m_ptr A pointer to the stored matrix
 * \returns 1 on success, 0 if an error occurred
 */
int set_hold_matrix(mat_ptr_t m_ptr);
/*!
 * \ingroup larc
 * \brief Decrements the hold value for a matrix in the matrix store (if above zero, the matrix will not be cleaned from the store)
 * \param m_ptr A pointer to the stored matrix
 * \returns 1 on success, 0 if an error occurred
 */
int release_hold_matrix(mat_ptr_t m_ptr);

/*!
 * \ingroup larc
 * \brief Returns the largest row or column level of A_ptr or B_ptr
 * \param A_ptr A pointer to a stored matrix
 * \param B_ptr A pointer to a stored matrix
 * \return The log base 2 of the biggest dimension of either matrix
 */
mat_level_t matrix_pair_max_level(mat_ptr_t A_ptr, mat_ptr_t B_ptr);

/*!
 * \ingroup larc
 *
 * \brief Checks if each matrix ptr provided is a valid matrix pointer.
 * 
 * If any are not, prints ERROR message that includes given function name and 
 * exits. Accepts any number of matrix pointers. Examples: 
 * exit_if_matrix_ptrs_invalid(__func__, 1, A_ptr);
 * or 
 * exit_if_matrix_ptrs_invalid(__func__, 3, A_ptr, B_ptr, C_ptr);
 *
 * \param function      name of routine check occurs in for error message
 * \param mat_ptr_num   the number of mat_ptrs to check 
 * \param ...           matrix pointers to check if MATRIX_PTR_INVALID
 */
void exit_if_matrix_ptrs_invalid(const char *function, int mat_ptr_num, ...);

/* access to commonly used preloaded matrices */
/*!
 * \ingroup larc
 * \brief A utility function which returns a preloaded zero matrix
 * \param row_level The log base 2 of the row dimension
 * \param col_level The log base 2 of the column dimension
 * \return A pointer to the matrix of all zeros with the given dimensions
 */
mat_ptr_t get_zero_matrix_ptr(mat_level_t row_level, mat_level_t col_level);
/*!
 * \ingroup larc
 * \brief A utility function which returns a preloaded identity matrix
 * \param level The log base 2 of the row and column dimensions
 * \return A pointer to the identity matrix with the given dimensions
 */
mat_ptr_t get_identity_matrix_ptr(mat_level_t level);
/*!
 * \ingroup larc
 * \brief A utility function which returns a preloaded integer Hadamard matrix
 *
 * The normal Hadamard matrix is scaled so that all elements are +1 or -1
 *
 * \param level The log base 2 of the row and column dimensions
 * \return A pointer to the integer Hadamard matrix with the given dimensions
 */
mat_ptr_t get_iHadamard_matrix_ptr(mat_level_t level);

/*!
 * \ingroup larc
 * \brief A utility function which increments the number of times a matrix is a submatrix of a larger matrix (if greater than zero, the matrix will not be cleaned)
 * \param id A pointer to a matrix
 * \returns The number of times that the given matrix is a quadrant submatrix
 */
uint32_t matrix_appears_as_sub_count_increment(mat_ptr_t id);
/*!
 * \ingroup larc
 * \brief A utility function which decrements the number of times a matrix is a submatrix of a larger matrix (if greater than zero, the matrix will not be cleaned)
 * \param id A pointer to a matrix
 * \returns The number of times that the given matrix is a quadrant submatrix
 */
uint32_t matrix_appears_as_sub_count_decrement(mat_ptr_t id);

/*!
 * \ingroup larc
 *
 * \brief Forms a matrix from 4 panels of submatrices and returns the 
 * corresponding matrix ptr from the matrix store.
 *
 * If the matrix is not already in the matrix store,
 * the matrix is added to the store before the 
 * ptr is returned. If the matrix is a row vector (row_level is 0), panels 
 * 0 and 1 must be valid and panels 2 and 3 must be MATRIX_PTR_INVALID. If 
 * the matrix is a column vector (col_level is 0), panels 0 and 2 must be 
 * valid and panels 1 and 3 must be MATRIX_PTR_INVALID. Otherwise, all 
 * panels must be valid.
 *
 * \param panel     The 4 matrix panels of the new matrix.
 * \param row_level The row level of the new matrix.
 * \param col_level The column level of the new matrix.
 * \return A matrix ptr for the new matrix.
 */
mat_ptr_t get_matPTR_from_array_of_four_subMatPTRs(mat_ptr_t panel[4], mat_level_t row_level, mat_level_t col_level);

#ifndef SWIG
/*!
 * \ingroup larc
 * \brief Finds a scalar in the matrix store, or inserts it if not found
 * \param scalar A scalarType value
 * \return A matrix pointer pointing to the scalar in the store
 */
mat_ptr_t get_valMatPTR_from_val(scalarType scalar);
#endif

/*!
 * \ingroup larc
 * \brief Finds the adjoint (complex conjugate transpose) of a matrix
 * \param id A handle to the matrix for which we would like the adjoint.
 * \return A handle to the adjoint matrix.
 */
mat_ptr_t matrix_adjoint(mat_ptr_t id);

/*!
 * \ingroup larc
 * \brief Given a matrixID, returns the pointer to the matrix
 * \param matrixID The matrixID label for the input matrix
 * \param arg_no A user-provided clue (used in error reporting) - if a routine calls this function more than once, helps determine which call failed
 * \param calling_routine A user-provided clue (used in error reporting) - the routine which called this function
 * \param print_flag Set to 1 to get warnings, 0 to suppress
 * \return The pointer to that matrix, or MATRIX_PTR_INVALID
 */
mat_ptr_t get_matPTR_from_matID(int64_t matrixID, const char* arg_no,
        const char* calling_routine, int print_flag);

/*!
 * \ingroup larc
 * \brief A routine that tries to remove a matrix from the matrix store
 * \param m_ptr A pointer to the matrix to be removed
 * \param n The hash table node containing the matrix (may be set to NULL)
 *
 * If the input n is not NULL, then n and supplied_hash will be used.
 * If the input n is set to NULL, the function will calculate the hash of the
 * matrix and then search the hash table for the matrix.
 *
 * \param supplied_hash The hash determining the hash chain containing the matrix
 * \return 1 if the remove succeeded, 0 if not 
 */
int remove_matrix_from_mat_store_by_matrix_andor_node (mat_ptr_t  m_ptr, hash_node_t *n, uint64_t supplied_hash);
/*!
 * \ingroup larc
 * \brief A routine that tries to remove a matrix from the matrix store
 * \param m_ptr A pointer to the matrix to be removed
 * \return 1 if the remove succeeded, 0 if not 
 */
int remove_matrix_from_mat_store_by_matrix (mat_ptr_t  m_ptr);
/*!
 * \ingroup larc
 * \brief A routine that tries to remove a matrix from the matrix store
 * \param n The hash table node containing the matrix
 * \param supplied_hash The hash determining the hash chain containing the matrix
 * \return 1 if the remove succeeded, 0 if not 
 */
int remove_matrix_from_mat_store_by_node (hash_node_t *n, uint64_t supplied_hash);

/*!
 * \ingroup larc
 * \brief Removes a matrix pointer from the matrixID table (used when cleaning has removed a matrix from the matrix store)
 * \param m_ptr The pointer to be removed from the table
 * \return 1
 */
int table_mat_ptr_by_matrixID_Remove_entry(mat_ptr_t m_ptr);


/************************************************************************
 *               MACRO FUNCTIONS acting on matrix id (mat_ptr_t)       *
 ***********************************************************************/

/*!
 * \ingroup larc
 * \brief Utility function that finds the row level of a matrix
 * \param mat_mID The matrixID for a stored matrix
 * \return The row level of that matrix
 */
mat_level_t matrix_row_level_matrixID(int64_t mat_mID);
/*!
 * \ingroup larc
 * \brief Utility function that finds the column level of a matrix
 * \param mat_mID The matrixID for a stored matrix
 * \return The column level of that matrix
 */
mat_level_t matrix_col_level_matrixID(int64_t mat_mID);
/*!
 * \ingroup larc
 * \brief returns the matrixID for one of the quadrant submatrices of a matrix
 * \param mat_mID The matrixID for a stored matrix
 * \param s One of [0,3], identifies which quadrant matrixID will be returned
 * \return The matrixID for the chosen submatrix
 */
int64_t matrix_sub_matrixID(int64_t mat_mID, int s);

// formerly inlined
/*!
 * \ingroup larc
 * \brief A utility function that checks if a given matrix is all zeros
 * \param mat_mID The matrixID for a stored matrix
 * \return 1 if the matrix is all zeros, 0 otherwise
 */
int matrix_is_zero_matrixID(int64_t mat_mID);
/*!
 * \ingroup larc
 * \brief A utility function that checks if a given matrix is invalid
 * \param mat_mID The matrixID for a (possibly) stored matrix
 * \return 1 if the matrix is invalid, 0 if it is valid
 */
int matrix_is_invalid_matrixID(int64_t mat_mID);
    
  
/******************************************************************
 *         ACCESSORS acting on matrix id (int64_t)                *
 *****************************************************************/

/*!
 * \ingroup larc
 * \brief Finds the matrixID of a scalar element of a stored matrix
 * \param mID The matrixID for a stored matrix
 * \param row The row index of the desired scalar
 * \param col The column index of the desired scalar
 * \return The matrixID for the chosen scalar value
 */
int64_t get_valID_from_matID_and_coords(int64_t mID, int64_t row, int64_t col);

/*!
 * \ingroup larc
 * \brief Finds the matrixID assigned to a scalar value 
 * \param val The scalarType value in string format
 * \return The matrixID for the scalar value
 */
int64_t get_valID_from_valString(char *val);

/*!
 * \ingroup larc
 * \brief Finds the value of a scalar element of a stored matrix
 * \param mID The matrixID for a stored matrix
 * \param row The row index of the desired scalar
 * \param col The column index of the desired scalar
 * \return The scalarType value of the desired scalar, in string format
 */
char *get_valString_from_matID_and_coords(int64_t mID, int64_t row, int64_t col);

/*!
 * \ingroup larc
 * \brief Creates a new stored matrix which differs from a previously stored matrix by a single element
 * \param mID The matrixID for a stored matrix
 * \param row The row index of the desired scalar
 * \param col The column index of the desired scalar
 * \param val A scalarType value in string format
 * \return The matrixID of the newly created matrix
 */
int64_t get_matID_from_oldMatID_newValString_and_coords(int64_t mID,
        int64_t row, int64_t col, char *val);

/*!
 * \ingroup larc
 * \brief Retrieve the maximum level of row and columns allowed in the store 
 * \return The maximum level allowed
 */
mat_level_t maximum_level(void);

/*!
 * \ingroup larc
 * \brief Retrieve a locality-approximation parameter: the number of significant bits kept after rounding
 * \return The number of significant bits kept after rounding
 */
int get_sighash(void);

/*!
 * \ingroup larc
 * \brief Retreive the locality-approximation parameter: distance from zero under which the approximation is zero
 * \return The value of this distance
 */
double get_zerorealthresh(void);

/*!
 * \ingroup larc
 * \brief Prints a summary of the matrix store usage 
 * \param outfilepath The path to the file which will contain the summary
 */
void matrix_store_report(char *outfilepath);

/*!
 * \ingroup larc
 * \brief Creates the matrix store structure
 * \param exponent The log base 2 of the size of the store
 * \param max_level Sets the maximum size of any matrix to be stored
 * \param sighash A locality-approximation parameter
 * \param zerobitthresh A locality-approximation parameter
 * \return 1 on successfully creating the matrix store, 0 on failure
 */
int create_matrix_store(size_t exponent, mat_level_t max_level, int sighash, int zerobitthresh);

/* python interface versions using matrixID */

/*!
 * \ingroup larc
 * \brief Increments the hold value for a matrix in the matrix store (if above zero, the matrix will not be cleaned from the store)
 * \param m_mID The matrixID of the stored matrix
 * \returns 1 on success, 0 if an error occurred
 */
int set_hold_matrix_from_matrixID(int64_t m_mID);
/*!
 * \ingroup larc
 * \brief Decrements the hold value for a matrix in the matrix store (if above zero, the matrix will not be cleaned from the store)
 * \param m_mID The matrixID of the stored matrix
 * \returns 1 on success, 0 if an error occurred
 */
int release_hold_matrix_from_matrixID(int64_t m_mID);

/*!
 * \ingroup larc
 * \brief Preloads a set of matrices which are generally useful
 * \return 1
 */
int preload_matrix_store(void);

/* access to commonly used preloaded matrices */
/*!
 * \ingroup larc
 * \brief A utility function which returns a preloaded zero matrix
 * \param row_level The log base 2 of the row dimension
 * \param col_level The log base 2 of the column dimension
 * \return The matrixID for the matrix of all zeros with the given dimensions
 */
int64_t get_zero_matrixID(mat_level_t row_level, mat_level_t col_level);
/*!
 * \ingroup larc
 * \brief A utility function which returns a preloaded identity matrix
 * \param level The log base 2 of the row and column dimensions
 * \return The matrixID for the identity matrix with the given dimensions
 */
int64_t get_identity_matrixID(mat_level_t level);
/*!
 * \ingroup larc
 * \brief A utility function which returns a preloaded integer Hadamard matrix
 * \param level The log base 2 of the row and column dimensions
 * \return The matrixID for the integer Hadamard matrix with the given dimensions
 */
int64_t get_iHadamard_matrixID(mat_level_t level);


// python interface versions - needed to avoid passing void * type
/*!
 * \ingroup larc
 *
 * \brief Forms a matrix from 4 submatrix matrixIDs and returns the 
 * corresponding matrixID from the matrix store.
 *
 * If the matrix is not already in the matrix store,
 * the matrix is added to the store before the
 * matrixID is returned. If the matrix is a row vector (row_level is 0),
 * mIDs A and B must be valid and mIDs C and D must be MATRIX_ID_INVALID. 
 * If the matrix is a column vector (col_level is 0), mIDs A and C must be
 * valid and mIDs B and D must be MATRIX_ID_INVALID. Otherwise, all 
 * mIDs must be valid. 
 *
 * \param A_mID ID for submatrix 1/4 of the new matrix.
 * \param B_mID ID for submatrix 2/4 of the new matrix.
 * \param C_mID ID for submatrix 3/4 of the new matrix.
 * \param D_mID ID for submatrix 4/4 of the new matrix.
 * \param row_level The row level of the new matrix.
 * \param col_level The column level of the new matrix.
 *
 * \return The matrixID of the new matrix. 
 */
int64_t get_matID_from_four_subMatIDs(int64_t A_mID, int64_t B_mID,
        int64_t C_mID, int64_t D_mID, mat_level_t row_level, 
        mat_level_t col_level);

/*!
 * \ingroup larc
 * \brief Returns the number of scalars currently in memory
 * \return The desired value
 */
size_t matrix_store_scalarCount(void);

/*!
 * \ingroup larc
 * \brief Returns the number of matrices (excluding scalars) in memory 
 * \return The desired value
 */
size_t matrix_store_matrixCount(void);

/*!
 * \ingroup larc
 * \brief Finds the adjoint (complex conjugate transpose) of a matrix
 * \param m_mID A handle to the matrix for which we would like the adjoint.
 * \return The matrixID of the adjoint matrix.
 */
int64_t matrix_adjoint_matrixID(int64_t m_mID);

/* Useful for some linear algebra applications */
// int preload_iHadamard(void);

/*!
 * \ingroup larc
 * \brief The total number of matrices created in the matrix store
 * \return The desired number
 */
uint64_t  num_matrices_created(void);

/*!
 * \ingroup larc
 * \brief The total number of matrices stored in the matrix store
 *
 * The current number stored may be less than the number created if there was matrix store cleaning.
 * 
 * \return The desired number
 */
uint64_t num_matrices_in_store (void);

/*!
 * \ingroup larc
 * \brief Prints statistical information about matrices to a file
 * \param start The smallest matrixID for which info is output
 * \param end The largest matrixID for which info is output
 * \param outfilepath Identifies the file to which the info is printed
 * \param comment A string printed at the top of the file
 * \return 1 on success, 0 otherwise
 */
int matrix_store_info_to_file(uint64_t start, uint64_t end, 
			      char *outfilepath, char *comment);

/*!
 * \ingroup larc
 * \brief Prints statistical information about a matrix store hash chain to a file
 * \param hash The hash identifying which hash chain to report on
 * \param outfilepath Identifies the file to which the info is printed
 * \param comment A string printed at the top of the file
 * \return 1 on success, 0 otherwise
 */
int matrix_hash_chain_info_to_file(uint64_t hash, char *outfilepath, char *comment);

/*!
 * \ingroup larc
 * \brief Prints statistical information about a matrix store hash chain to stdout
 * \param hash The hash identifying which hash chain to report on
 * \param comment A string printed at the top of the file
 * \return 1 on success, 0 otherwise
 */
int matrix_hash_chain_info_to_screen(uint64_t hash, char *comment);

/*!
 * \ingroup larc
 * \brief A routine that tries to remove a matrix from the matrix store
 * \param m_mID The matrixID of the matrix to be removed
 * \return 1 if the remove succeeded, 0 if not 
 */
int remove_matrix_from_mat_store_by_matrixID (int64_t  m_mID);

/*!
 * \ingroup larc
 * \brief Returns the hash of the matrix pointer associated with the input matrixID
 * \param m_mID The matrixID
 * \return The desired hash value, or -1 if failed
 */
int64_t matrix_hashID_from_matrixID(int64_t m_mID);

/*! 
 * \ingroup larc
 * \brief Cleans unneeded matrices from the matrix store
 * \return 1
 */
int clean_matrix_store();

/*! 
 * \ingroup larc
 * \brief Cleans unneeded matrices from a particular hash chain of the matrix store
 * \param hash Identifies the hash chain to be cleaned
 * \return 1
 */
int clean_matrix_hash_chain(uint64_t hash);

/*!
 * \ingroup larc
 * \brief Gets the current number of square matrices of a given size
 * \param i The log base 2 of the row and column sizes
 * \return The number of matrices in the store of this size
 */
int64_t get_diag_hist(uint64_t i);

/*!
 * \ingroup larc
 * \brief Counts the number of entries in a specific hash chain in the matrix store
 * \param i Identifies the hash chain to be counted
 * \return The length of the chain
 */
int64_t matrix_hash_chain_length(uint64_t i);

/*!
 * \ingroup larc
 * \brief Returns the trace of a stored matrix in scalarType string format (if a scalar, the value of that scalar)
 * \param m_ID The matrixID of a stored matrix
 * \return The string representing the trace value
 */
char *matrix_trace_matrixID(int64_t m_ID);

#ifdef HASHSTATS
/*!
 * \ingroup larc
 * \brief Prints hash table statistics for the matrix store, and produces a standard hash_report to stdout or file
 * \param accesses_file A file which will hold an array of the number of access to each hash bucket
 * \param nodes_file A file which will hold an array of the number of nodes in each hash_chain
 * \param report_file A file which will receive a standard hash report
 */
void matrix_hashstats(char *accesses_file,  char *nodes_file, char *report_file);
#endif

/*!
 * \ingroup larc
 * \brief Finds the pointer to the (i,j) entry of a matrix
 * \param row_i The index i
 * \param col_j The index j
 * \param m_ptr The pointer to the full matrix M
 * \return The pointer to the scalar M(i,j)
 */
mat_ptr_t get_valMatPTR_from_matPTR_and_coords(long int row_i,
				    long int col_j, mat_ptr_t m_ptr);


/*!
 * \ingroup larc
 * \brief returns the quadrant location in a small_level sized submatrix of a location (big_row, big_col) in a big matrix.
 * \param big_row The index of the row in a big matrix (0 based)
 * \param big_col The index of the col in a big matrix (0 based)
 * \param small_level The level of the small submatrix of the big matrix.
 * \return The quadrant submatrix (0,1,2, or 3) inside the small_level submatrix in which the reduced coordinates of (big_row,big_col) lie.
 *
 * small_row The index of the row in the internal small_level size submatrix (0 based)
 * small_col The index of the col in the internal small_level size submatrix (0 based)
************************************************************************/
int get_quad_in_submatrix(int64_t big_row, int64_t big_col,
			  mat_level_t small_level);


/*!
 * \ingroup larc
 * \brief This routine returns a matrixPTR to a square matrix of the given level with a single nonzero entry one_scalar in position (row,col), where row and col are zero-indexed.
 * \param level The level of the square matrix to be created.
 * \param one_scalar A matrix PTR to a scalarType value (presumably nonzero).
 * \param row_i The zero-indexed row index for the nonzero element one_scalar.
 * \param col_j The zero-indexed col index for the nonzero element one_scalar.
 * \return Matrix pointer to a square matrix with the specified single nonzero entry.
 *
 *  Note: when the user makes this call, i and j need to be smaller than the dimension
*        of the matrix.  But, when this calls itself, i and j will not change even as
*        the level is reduced.  This is all handled gracefully in get_quad_in_submatrix.
************************************************************************/
mat_ptr_t get_matPTR_single_nonzero_using_valPTR_at_coords(mat_level_t level,
	  int64_t row_i, int64_t col_j, scalarType *one_scalar);

 

#ifdef __cplusplus
}
#endif // #ifdef __cplusplus

#endif    
