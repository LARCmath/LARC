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


#ifndef LARC_MATRIX_STORE_H
#define LARC_MATRIX_STORE_H

#ifdef __cplusplus
extern "C" {
#endif // #ifdef __cplusplus

#include <stddef.h>      // size_t
#include <stdarg.h>      // .../va_list/va_args/
#include "larc.h"

#define MATRIX_ID_INVALID -1

// The packedID is based on the matrixID but includes a bit which indicates
// whether a matrix is scalar or non-scalar. This is used throughout LARC to
// simplify code, but we need to be able to generate packedID from matrixID
// and to recover matrixID from packedID.
//
// There is more than one way to make this definition. The first method we
// tried is below, which makes the packedID twice the matrixID and then adds
// one when the matrix is scalar.
//
//#define SCALAR_BIT 1
//#define IS_SCALAR(packedID) (packedID & 1)
//#define MID_FROM_PID(packedID) ( packedID>>1 )
//#define PID_FROM_SCALAR_MID(matrixID) ( (matrixID<<1) ^ 1)
//#define PID_FROM_NONSCALAR_MID(matrixID) (matrixID<<1)
//
// An alternative is to make the lower bits of the packedID the same as the
// matrixID, and set a high bit to one when the matrix is scalar. Any bit
// position > 2*LOG_SIZE_MPT (currently this is 50) will work, since the
// matrixID can have no more than this many bits before LARC exits. We chose to
// use bit 60 below. In the future, we may wish to add additional flag bits to
// the packedID, which is why we are using this method rather than the first.
#define SCALAR_BIT 0x1000000000000000L
#define IS_SCALAR(packedID) (packedID & SCALAR_BIT)
#define MID_FROM_PID(packedID) \
  ( (packedID == MATRIX_ID_INVALID) ? packedID : \
                                      (packedID & ~SCALAR_BIT))
#define PID_FROM_SCALAR_MID(matrixID) ( matrixID | SCALAR_BIT)
#define PID_FROM_NONSCALAR_MID(matrixID) (matrixID)


/*!
 * \ingroup larc
 * \brief removes any flag bits in the packedID and returns a matrixID
 *
 * This routine is provided for Python access to the C macro which
 * clears any flag bits from the packedID and performs any other needed
 * operations to produce the matrixID that was assigned when a scalar
 * or matrix was added to the store.
 *
 * \param packedID A packedID
 * \return The matrixID corresponding to that packedID
 */
inline int64_t matrixID_from_packedID(int64_t packedID)
{
    return MID_FROM_PID(packedID);
}

#ifndef SWIG
/*!
 * \cond
 * Doxygen will ignore the macros below
 */

/***********************************************************************
 *               MACRO FUNCTIONS acting on record pointers             *
 ***********************************************************************/

#define RECORD_PTR_INVALID ((record_ptr_t)0)
#define MATRIX_PTR_INVALID ((matns_ptr_t)0)
#define SCALAR_PTR_INVALID ((mats_ptr_t)0)


// Boolean functions

/*!
 * \ingroup larc
 * \brief Determine if a record is invalid
 * \param r_ptr A handle to the record we are testing
 * \return 1 if r_ptr is RECORD_PTR_INVALID, 0 otherwise
 */
#define record_is_invalid(r_ptr) (r_ptr == RECORD_PTR_INVALID)

/*!
 * \ingroup larc
 * \brief Determine if a matrix is locked (can never be cleaned from the matrix store)
 *
 * This macro works for both scalar and nonscalar matrices
 *
 * \param m_ptr A handle to the matrix we are testing
 * \return 1 if m_ptr is locked, 0 otherwise
 */
#define matrix_has_lock(m_ptr) (m_ptr && (m_ptr)->lock)

/*!
 * \ingroup larc
 * \brief Determine if a matrix is held (can not currently be cleaned from the matrix store)
 *
 * This macro works for both scalar and nonscalar matrices
 *
 * \param m_ptr A handle to the matrix we are testing
 * \return 0 if the matrix is not held, else the number of holds on the matrix
 */
#define matrix_has_hold(m_ptr) ((( (m_ptr) ?  (m_ptr)->hold : 0 )))

// Functions retrieving values inside a larc_matrix structure

/*!
 * \ingroup larc
 * \brief Determines whether the matrix is needed for storing a larger matrix
 *
 * This works for both scalar and nonscalar matrices
 *
 * \param m_ptr A handle to the matrix
 * \return The count of the number of times a matrix is a submatrix of a larger stored matrix
 */
#define matrix_appears_as_sub_count(m_ptr) ((m_ptr)->appears_as_sub_count)

/*!
 * \endcond
 * Any other macros below will not be ignored by Doxygen
 */


/******************************************************************
 *         matrix_store_t                                         *
 *****************************************************************/


struct matrix_store_t {
  size_t num_nonscalars;    // the total number matrices in the nonscalar store currently
  size_t num_scalars;     // the total number of matrices in the scalar store currently
#ifdef MAR  
  size_t num_neighbors;     // the total number of neighbor tile records
#endif // #ifdef MAR  
  size_t max_num_nonscalars;  // the maximum number of nonscalar matrices that
			 // were ever in the store at the same time
  size_t max_num_scalars;   // the maximum number of scalar matrices that were
                         // ever in the store at the same time
  size_t **hist;         // histogram of matrices by row_level and col_level
  hash_table_t *nonscalar_hash_table;  
  hash_table_t *scalar_hash_table;  

  // PRELOADED VALUES IN MATRIX STORE THAT WE NEVER WANT TO ERASE AND USE IN MATH FUNCTIONS
  mat_level_t largest_level;    // largest row_level and col_level allowed in store, 
                               // also largest level of preloaded  matrices
  // LARC uses a locality sensitive hash which depends on separating
  // the space of scalars into regions, and then hashing all scalars in
  // a particular region to the same hash value.  This allows LARC to create
  // at most one matrix record per region.
  int regionbitparam;   // SPR parameter defining the width of most regions 
  int zeroregionbitparam;  // SPR parameter for regions near zero

#ifndef MAR
  // SPRmode
  int special_zeroregion;   // if 1 then have special regions near zero

  // the SPR zero threshold parameter is derived from zeroregionbitparam, and
  // determines which of the region centers defined through use of 
  // regionbitparam are collapsed to a zero region. Its type depends on the
  // chosen scalarType. 
#ifdef IS_RATIONAL
  mpq_t zerorealthresh;
#elif defined(IS_MP)
  mpfr_t zerorealthresh;
#else
  double zerorealthresh;
#endif // #ifdef IS_RATIONAL
#endif // #ifndef MAR

  int64_t **zero;     // doubly indexed matrix IDs for the preloaded zero matrices
  int64_t *identity;  // array of matrix IDs for the preloaded identity matrices (square)
  int64_t *iHadamard; // array of matrix IDs for the preloaded integer Hadamard matrices

  // POSSIBLY CAN BE REMOVED                      
  uint32_t deepest;    // longest successful traversal down a hash chain to retrieve a matrix or scalar

};



/******************************************************************
 *         ACCESSORS acting on matrix PTR (matns_ptr_t)              *
 *****************************************************************/

/*!
 * \ingroup larc
 * \brief A utility function returning the matrix_type for a matrix pointer
 *
 * Note that this function does not work for pointers to scalar matrices
 *
 * \param m_ptr A pointer to a nonscalar matrix (in the MatrixStore)
 * \returns SCALAR, ROW_VECTOR, COL_VECTOR, or MATRIX
 */
matrix_type_t matrix_type(matns_ptr_t m_ptr);

/*************************************************************************
 ************************************************************************/

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

/*!
 * \ingroup larc
 * \brief Finds a scalar in the matrix store, or inserts it if not found
 * \param scalar A scalarType value
 * \return A matrix pointer pointing to the scalar in the store
 */
mats_ptr_t get_scalarPTR_for_scalarVal(scalarType scalar);

/*!
 * \ingroup larc
 * \brief Given a packedID, returns a generic record pointer to the matrix or scalar associated to that packedID
 *
 * The record_ptr_t type is returned so that this function can be used for
 * matrices in both the MatrixStore and ScalarStore.
 *
 * \param packedID The packedID label for the input matrix
 * \param arg_no A user-provided clue (used in error reporting) - if a routine calls this function more than once, helps determine which call failed
 * \param calling_routine A user-provided clue (used in error reporting) - the routine which called this function
 * \param print_flag Set to 1 to get warnings, 0 to suppress
 * \return The pointer to that matrix or scalar, or RECORD_PTR_INVALID
 */
record_ptr_t get_recordPTR_from_pID(int64_t packedID, const char* arg_no,
        const char* calling_routine, int print_flag);

#endif


/************************************************************************
 *               MACRO FUNCTIONS acting on matrix id                   *
 ***********************************************************************/

/*!
 * \ingroup larc
 * \brief Utility function that finds the row level of a matrix
 * \param mat_pID The packedID for a stored matrix
 * \return The row level of that matrix
 */
mat_level_t matrix_row_level(int64_t mat_pID);

/*!
 * \ingroup larc
 * \brief Utility function that finds the column level of a matrix
 * \param mat_pID The packedID for a stored matrix
 * \return The column level of that matrix
 */
mat_level_t matrix_col_level(int64_t mat_pID);

/*!
 * \ingroup larc
 * \brief returns the packedID for one of the quadrant submatrices of a matrix
 * \param mat_pID The packedID for a stored matrix
 * \param s One of [0,3], identifies which quadrant packedID will be returned
 * \return The packedID for the chosen submatrix
 */
int64_t get_pID_of_indexed_submatrix(int64_t mat_pID, int s);

// formerly inlined
/*!
 * \ingroup larc
 * \brief A utility function that checks if a given matrix is all zeros
 * \param mat_pID The packedID for a stored matrix
 * \return 1 if the matrix is all zeros, 0 otherwise
 */
int matrix_is_zero(int64_t mat_pID);

/*!
 * \ingroup larc
 * \brief A utility function that checks if a given matrix is the identity
 * \param mat_pID The packedID for a stored matrix
 * \return 1 if the matrix is the identity, 0 otherwise
 */
int matrix_is_identity(int64_t mat_pID);

/*!
 * \ingroup larc
 * \brief A utility function that checks if a given matrix is invalid
 * \param mat_pID The packedID for a (possibly) stored matrix
 * \return 1 if the matrix is invalid, 0 if it is valid
 */
int matrix_is_invalid(int64_t mat_pID);
    
  
/******************************************************************
 *         ACCESSORS acting on matrix id (int64_t)                *
 *****************************************************************/

/*!
 * \ingroup larc
 * \brief Finds the packedID of a scalar element of a stored matrix
 * \param pID The packedID for a stored matrix
 * \param row The row index of the desired scalar
 * \param col The column index of the desired scalar
 * \return The packedID for the chosen scalar element
 */
int64_t get_scalarID_from_pID_and_coords(int64_t pID, int64_t row, int64_t col);

/*!
 * \ingroup larc
 * \brief Finds the packedID assigned to a scalar, or adds the scalar to the ScalarStore and returns the new packedID
 * \param val The scalarType value in string format
 * \return The packedID for the scalar value
 */
int64_t get_valID_from_valString(char *val);

/*!
 * \ingroup larc
 * \brief Finds the value of a scalar element of a stored matrix
 *
 * For REAL, COMPLEX, MPREAL, and MPCOMPLEX scalar types, the
 * returned string is in base 10 and thus may not be able to
 * exactly represent the internal binary value.
 *
 * \param pID The packedID for a stored matrix
 * \param row The row index of the desired scalar
 * \param col The column index of the desired scalar
 * \return The scalarType value of the desired scalar, in string format
 */
char *get_readableString_scalar_from_pID_and_coords(int64_t pID, int64_t row, int64_t col);

/*!
 * \ingroup larc
 * \brief Finds the exact value of a scalar element of a stored matrix
 *
 * For REAL, COMPLEX, MPREAL, and MPCOMPLEX scalar types, the returned
 * string is in hexadecimal format, and is guaranteed to be restored to
 * the exact same value if written out and read back in.
 *
 * \param pID The packedID for a stored matrix
 * \param row The row index of the desired scalar
 * \param col The column index of the desired scalar
 * \return The exact scalarType value of the desired scalar, in string format
 */
char *get_exactString_scalar_from_pID_and_coords(int64_t pID, int64_t row, int64_t col);

/*!
 * \ingroup larc
 * \brief Creates a new stored matrix which differs from a previously stored matrix by a single element
 * \param pID The packedID for a stored matrix
 * \param row The row index of the desired scalar
 * \param col The column index of the desired scalar
 * \param val A scalarType value in string format
 * \return The packedID of the modified matrix
 */
int64_t replace_scalar_in_matrix_by_string_and_coords(int64_t pID,
        int64_t row, int64_t col, char *val);

/*!
 * \ingroup larc
 * \brief Retrieve the exponent for the matrix store
 * \return The exponent used to create the matrix store hash table
 */
size_t get_nonscalar_store_exp(void);

/*!
 * \ingroup larc
 * \brief Retrieve the exponent for the scalar store
 * \return The exponent used to create the scalar store hash table
 */
size_t get_scalar_store_exp(void);

/*!
 * \ingroup larc
 * \brief Retrieve the maximum level of row and columns allowed in the store 
 * \return The maximum level allowed
 */
mat_level_t max_level_allowed_matrixStore(void);

/*!
 * \ingroup larc
 * \brief Retrieve a SPR parameter: the number of significant bits kept after rounding
 * \return The number of significant bits kept after rounding
 */
int get_regionbitparam(void);

/*!
 * \ingroup larc
 * \brief Retrieve a SPR parameter: the number of significant bits kept after rounding near zero
 * \return The number of significant bits kept after rounding
 */
int get_zeroregionbitparam(void);

#ifndef MAR // SPRmode  
/*!
 * \ingroup larc
 * \brief Retrieve the SPR parameter: distance which defines the region around zero
 * \return The value of this distance
 */
#ifdef IS_RATIONAL
const mpq_t* const get_zerorealthresh(void);
#elif defined(IS_MP)
const mpfr_t* const get_zerorealthresh(void);
#else
double get_zerorealthresh(void);
#endif  // IS_RATIONAL
#endif // SPRmode

#ifndef MAR //SPRmode  
/*!
 * \ingroup larc
 * \brief Retreive the SPR parameter determining whether zero regions are calculated differently
 * \return 1 if zero regions of SPR need to be calculated specially.
 */
int is_special_zeroregion(void);
#endif //SPRmode  

/*!
 * \ingroup larc
 * \brief Prints a summary of the matrix store usage 
 * \param outfilepath The path to the file which will contain the summary
 */
void matrix_store_report(char *outfilepath);

#ifndef SWIG
/*!
 * \ingroup larc
 * \brief Creates the matrix store structure
 * \param exponent The log base 2 of the size of the store
 * \param max_level Sets the maximum size of any matrix to be stored
 * \param regionbitparam A SPR parameter
 * \param zeroregionbitparam A SPR parameter
 * \return 1 on successfully creating the matrix store, 0 on failure
 */
int create_matrix_store(size_t exponent, mat_level_t max_level, int regionbitparam, int zeroregionbitparam);
#endif // #ifndef SWIG

/*!
 * \ingroup larc
 * \brief Sets the lock value for a matrix (if nonzero, the matrix will not be cleaned)
 * \param m_pID The packedID of the stored matrix
 * \returns 1 on success, 0 if an error occurred
 */
int lock_matrix(int64_t m_pID);

/*!
 * \ingroup larc
 * \brief Increments the hold value for a matrix (if above zero, the matrix will not be cleaned)
 * \param m_pID The packedID of the stored matrix
 * \returns 1 on success, 0 if an error occurred
 */
int set_hold_matrix(int64_t m_pID);

/*!
 * \ingroup larc
 * \brief Decrements the hold value for a matrix (if above zero, the matrix will not be cleaned)
 * \param m_pID The packedID of the stored matrix
 * \returns 1 on success, 0 if an error occurred
 */
int release_hold_matrix(int64_t m_pID);

#ifndef SWIG
/*!
 * \ingroup larc
 * \brief Preloads a set of matrices which are generally useful
 * \return 1
 */
int preload_matrix_and_scalar_stores(void);
#endif // #ifndef SWIG

/* access to commonly used preloaded matrices */
/*!
 * \ingroup larc
 * \brief A utility function which returns a preloaded zero matrix
 * \param row_level The log base 2 of the row dimension
 * \param col_level The log base 2 of the column dimension
 * \return The packedID for the matrix of all zeros with the given dimensions
 */
int64_t get_zero_pID(mat_level_t row_level, mat_level_t col_level);

/*!
 * \ingroup larc
 * \brief A utility function which returns a preloaded identity matrix
 * \param level The log base 2 of the row and column dimensions
 * \return The packedID for the identity matrix with the given dimensions
 */
int64_t get_identity_pID(mat_level_t level);

/*!
 * \ingroup larc
 * \brief A utility function which returns a preloaded integer Hadamard matrix
 * \param level The log base 2 of the row and column dimensions
 * \return The packedID for the integer Hadamard matrix with the given dimensions
 */
int64_t get_iHadamard_pID(mat_level_t level);

/*!
 * \ingroup larc
 *
 * \brief Forms a matrix from 4 submatrix packedIDs and returns the 
 * corresponding packedID from the matrix store.
 *
 * If the matrix is not already in the matrix store, the matrix is added
 * to the store before the packedID is returned. If the matrix is a row
 * vector (row_level is 0), pIDs A and B must be valid and pIDs C and D
 * must be MATRIX_ID_INVALID. If the matrix is a column vector (col_level
 * is 0), pIDs A and C must be valid and pIDs B and D must be
 * MATRIX_ID_INVALID. Otherwise, all pIDs must be valid. 
 *
 * \param A_pID ID for submatrix 1/4 of the new matrix.
 * \param B_pID ID for submatrix 2/4 of the new matrix.
 * \param C_pID ID for submatrix 3/4 of the new matrix.
 * \param D_pID ID for submatrix 4/4 of the new matrix.
 * \param row_level The row level of the new matrix.
 * \param col_level The column level of the new matrix.
 *
 * \return The packedID of the new matrix. 
 */
int64_t get_pID_from_four_sub_pIDs(int64_t A_pID, int64_t B_pID,
        int64_t C_pID, int64_t D_pID, mat_level_t row_level, 
        mat_level_t col_level);

/*!
 * \ingroup larc
 * \brief Forms a matrix from an array of 4 submatrix packedIDs and returns the 
 * corresponding packedID from the matrix store.
 *
 * \param pID An array of four packedIDs, to be the submatrices of the larger matrix
 * \param row_level The row level of the new matrix.
 * \param col_level The column level of the new matrix.
 *
 * \return The packedID of the new matrix. 
 */
inline int64_t get_pID_from_array_of_four_sub_pIDs(int64_t pID[4],
        mat_level_t row_level, mat_level_t col_level)
{

  return get_pID_from_four_sub_pIDs(pID[0], pID[1], pID[2], pID[3],
         row_level, col_level);
}

/*!
 * \ingroup larc
 * \brief Returns the number of scalars currently in memory
 * \return The desired value
 */
size_t scalar_store_count(void); /*!
 * \ingroup larc
 * \brief Returns the number of matrices (excluding scalars) in memory 
 * \return The desired value
 */
size_t nonscalar_store_count(void);

/*!
 * \ingroup larc
 * \brief Finds the adjoint (complex conjugate transpose) of a matrix
 * \param m_pID The packedID of the matrix for which we would like the adjoint
 * \return The packedID of the adjoint matrix.
 */
int64_t adjoint(int64_t m_pID);

/* Useful for some linear algebra applications */
// int preload_iHadamard(void);

/*!
 * \ingroup larc
 * \brief The total number of matrices created in the matrix store
 * \return The desired number
 */
uint64_t num_matrices_created(void);

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
 *
 * This function takes a range of matrixIDs (not packedIDs) and outputs
 * information about the matrices and scalars in this range.
 * 
 * \param start The smallest matrixID for which info is output
 * \param end The largest matrixID for which info is output
 * \param outfilepath Identifies the file to which the info is printed
 * \param comment A string printed at the top of the file
 * \return 1 on success, 0 otherwise
 */
int fprint_store_info_for_matrixID_range(uint64_t start, uint64_t end, 
			      char *outfilepath, char *comment);

/*!
 * \ingroup larc
 * \brief Prints statistical information about a matrix store hash chain to a file
 * \param hash The hash identifying which hash chain to report on
 * \param outfilepath Identifies the file to which the info is printed
 * \param comment A string printed at the top of the file
 * \return 1 on success, 0 otherwise
 */
int fprint_matrix_hash_chain_info(uint64_t hash, char *outfilepath, char *comment);

/*!
 * \ingroup larc
 * \brief Prints statistical information about a matrix store hash chain to stdout
 * \param hash The hash identifying which hash chain to report on
 * \param comment A string printed at the top of the file
 * \return 1 on success, 0 otherwise
 */
int print_matrix_hash_chain_info(uint64_t hash, char *comment);

/*!
 * \ingroup larc
 * \brief A routine that tries to remove a matrix from the matrix store
 * \param m_pID The packedID of the matrix to be removed
 * \return 1 if the remove succeeded, 0 if not 
 */
int remove_matrix_from_store(int64_t  m_pID);

/*!
 * \ingroup larc
 * \brief Returns the hash of the input packedID
 * \param m_pID The packedID to be hashed
 * \return The desired hash value, or -1 if failed
 */
int64_t hash_pID(int64_t m_pID);

/*! 
 * \ingroup larc
 * \brief Cleans unneeded matrices from both nonscalar store and scalar store
 * \return 1
 */
int clean_matrix_storage();

/*!
 * \ingroup larc
 * \brief Gets the current number of square matrices of a given size
 * \param i The log base 2 of the row and column sizes
 * \return The number of matrices in the store of this size
 */
int64_t get_counts_square_matrices_by_level(uint64_t i);

/*!
 * \ingroup larc
 * \brief Counts the number of entries in a specific hash chain in the matrix store
 * \param i Identifies the hash chain to be counted
 * \return The length of the chain
 */
int64_t matrix_hash_chain_length(uint64_t i);

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

#ifndef SWIG
/*!
 * \ingroup larc
 * \brief Finds the pointer to the (i,j) entry of a matrix
 *
 * This function is needed in io.c as well as matrix_store.c, so cannot be
 * made static.
 *
 * \param row_i The index i
 * \param col_j The index j
 * \param m_pID The packedID of the full matrix M
 * \return The pointer to the scalar at M(i,j)
 */
mats_ptr_t get_scalarPTR_from_pID_and_coords(uint64_t row_i,
				    uint64_t col_j, int64_t m_pID);
#endif

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
 * \brief Returns a packedID to a square matrix of the given level with a single nonzero entry in position (row,col), where row and col are zero-indexed.
 * \param level The level of the square matrix to be created.
 * \param scalarID A packedID for a scalarType value (presumably nonzero).
 * \param row_i The zero-indexed row index for the nonzero element one_scalar.
 * \param col_j The zero-indexed col index for the nonzero element one_scalar.
 * \return Matrix pointer to a square matrix with the specified single nonzero entry.
 *
 *  Note: when the user makes this call, i and j need to be smaller than the
 *  dimension of the matrix.  But, when this calls itself, i and j will not
 *  change even as the level is reduced.  This is all handled gracefully in
 *  get_quad_in_submatrix.
 */
int64_t get_matrix_with_single_nonzero_at_coords(mat_level_t level,
	  int64_t row_i, int64_t col_j, int64_t scalarID);
 
#ifndef SWIG
/*!
 * \brief removes all memory allocated by create_matrix_store() and some additional values allocated within initialize_larc()
 *
 * This routine is called by shutdown_larc(), and should never be called by
 * any other routine.
 */
void free_matrix_store(void);

/*!
 \ingroup larc
 * \brief A routine that tries to remove a matrix from the scalar store
 * \param r_ptr A pointer to the matrix to be removed
 * \param n The hash table node containing the matrix (may be set to NULL)
 * \param supplied_hash The index for the hash chain containing the matrix
 * \param m_pID The packedID for the matrix to be removed
 *
 * If the input n is not NULL, then n will be used; if it is NULL, the
 * function will search the hash table for the matrix.
 *
 * \return NULL if n null, otherwise the pointer to the next node in the hash chain
 */
hash_node_t *remove_matrix_from_scalar_store_by_recordPTR_andor_nodePTR (
  record_ptr_t r_ptr, hash_node_t *n, uint64_t supplied_hash, int64_t m_pID);

#endif // #ifndef SWIG

#ifndef MAR
#ifndef SWIG

mats_ptr_t find_scalar_SPRmode(uint64_t hash, scalarType scalar);

mats_ptr_t insert_scalar_SPRmode(uint64_t hash, scalarType scalar);
#endif // SWIG
#endif // SPR

#ifdef MAR
#ifndef SWIG
/*!
 * \ingroup larc
 * \brief determine whether a tile has a representative stored
 *
 * This function is used in mar.c and therefore cannot be made static.
 *
 * \param hash_chain The index (hash) of the hash chain to be searched
 * \param index_PTR A pointer to the tile index we are looking for
 * \return scalar pointer to record, or SCALAR_PTR_INVALID if not found
 */
mats_ptr_t find_record_from_tile_index(uint64_t hash_chain, MAR_tile_index_t *index_PTR);

/*!
 * \ingroup larc
 * \brief create a new scalar record 
 *
 * This function is used in mar.c and therefore cannot be made static.
 *
 * \param target_tile_index_PTR The tile index for the scalar
 * \param hash_value The hash chain which will receive the new record
 * \param target_scalar The scalar value to be stored
 * \return pointer to new primary record
 */
 mats_ptr_t  insert_primary_scalar_record(const MAR_tile_index_t *target_tile_index_PTR,
				  const uint64_t hash_value, const scalarType target_scalar);

/*!
 * \ingroup larc
 * \brief Create a new record for a neighbor tile claimed by a primary tile
 *
 * This function is used in mar.c and therefore cannot be made static.
 *
 * \param nhbr_hash The hash chain where the neighbor record will be stored
 * \param nhbr_tile_index_PTR Pointer to the tile index for the neighbor tile
 * \param primary_ptr Record pointer for the primary tile that claims this neighbor
 * \param tile_offset Two bits indicating the position of the neighbor relative to the primary tile
 * \return Pointer to the hash node containing the neighbor record
 */
 hash_node_t  *insert_nhbr_record(uint64_t nhbr_hash,
                                  MAR_tile_index_t *nhbr_tile_index_PTR,
                                  mats_ptr_t primary_ptr,
                                  unsigned int tile_offset);

  // void   insert_nhbr_node_ptr_in_primary_record(
  //					       hash_node_t * nhbr_node_ptr,
  //					       matns_ptr_t primary_ptr);

#endif // not SWIG
#endif // MAR

/*!
 * \ingroup larc
 * \brief utility routine that confirms an input packedID has a valid pointer
 *
 * If the packedID is MATRIX_ID_INVALID, or if the pointer associeated to
 * that ID in the ptrs_indexed_by_matrixID array is null, the program exits.
 *
 * \param A_pID A packedID 
 * \param callingRoutine The name of the routine that called the check
 * \param A_Rptr Pointer to the record pointer associated to A_pID
 */
void check_validity_one_input(int64_t A_pID,
     const char *callingRoutine, record_ptr_t *A_Rptr);

/*!
 * \ingroup larc
 * \brief utility routine that confirms input packedIDs have valid pointers
 *
 * If a packedID is MATRIX_ID_INVALID, or if the pointer associeated to
 * an ID in the ptrs_indexed_by_matrixID array is null, the program exits.
 *
 * \param A_pID A packedID 
 * \param B_pID A packedID 
 * \param callingRoutine The name of the routine that called the check
 * \param A_Rptr Pointer to the record pointer associated to A_pID
 * \param B_Rptr Pointer to the record pointer associated to B_pID
 */
void check_validity_two_input(int64_t A_pID, int64_t B_pID,
     const char *callingRoutine, record_ptr_t *A_Rptr, record_ptr_t *B_Rptr);

#if ALERT_ON_SNAP
#ifndef SWIG
void handle_snap(scalarType original_scalar, scalarType new_scalar, const char *tile_layout);
#endif // #ifndef SWIG
#endif // #if ALERT_ON_SNAP

#ifdef __cplusplus
}
#endif // #ifdef __cplusplus

#endif    
