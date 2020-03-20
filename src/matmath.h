//                      matmath.h
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

#ifndef LARC_MATMATH_H
#define LARC_MATMATH_H

#ifdef __cplusplus
extern "C" {
#endif // #ifdef __cplusplus

#include <gmp.h>  

#include "larc.h"
#include "global.h"
#include "op_store.h"
#include "matrix_store.h"
#include "organize.h"

/*!
 * \ingroup larc
 * 
 * \brief Add two matrices
 * 
 * \param ptr_A A handle to the first matrix
 *
 * \param ptr_B A handle to the second matrix
 * 
 * \return A handle to the matrix that is the sum of the matrices
 * referenced by the arguments.
 */
mat_ptr_t matrix_add(mat_ptr_t ptr_A, mat_ptr_t ptr_B);

/*!
 * \ingroup larc
 * 
 * \brief Subtract matrix B from matrix A. 
 * 
 * \param ptr_A A handle to the first matrix (`minuend')
 *
 * \param ptr_B A handle to the second matrix (`subtrahend')
 * 
 * \return A handle to the matrix that is the difference of the matrices 
 * referenced by the arguments.
 */
mat_ptr_t matrix_diff(mat_ptr_t ptr_A, mat_ptr_t ptr_B);

/*!
 * \ingroup larc
 *
 * \brief Multiply two matrices.
 *
 * Empty matrix store and op store after all internal recursive
 * calls to matrix_mult where input matrices have row or column
 * level equal to cleanThresh. Also remove extra matrices produced
 * during recursive calls to matrix_mult where input matrices have row
 * or column level above cleanThresh. For efficiency, cleanThresh should be
 * as large possible such that the calculation can finish. To set to no cleaning, 
 * set cleanThresh to higher than the levels of the input matrices. 
 *
 * \param A_ptr A handle to the first matrix.
 * \param B_ptr A handle to the second matrix.
 * \param cleanThresh A matrix level threshold. 
 *
 * \return A handle to the matrix that is the product of the matrices 
 * referenced by the arguments.
 */
mat_ptr_t matrix_mult_clean(mat_ptr_t A_ptr, mat_ptr_t B_ptr, mat_level_t cleanThresh);
/*!
 * \ingroup larc
 *
 * \brief Multiply two matrices. Wraps matrix_mult_clean to set no cleaning
 * or for backwards compatibility. 
 *
 * \param ptr_A A handle to the first matrix.
 * \param ptr_B A handle to the second matrix.
 *
 * \return A handle to the matrix that is the product of the matrices 
 * referenced by the arguments.
 */
mat_ptr_t matrix_mult(mat_ptr_t ptr_A, mat_ptr_t ptr_B);

/*!
 *  \ingroup larc
 * 
 *  \brief Multiply a matrix by a scalar.
 * 
 *  \param ptr_A A handle to a matrix with a single element, which is considered the scalar.
 *  \param ptr_B A handle to a matrix of any size.
 * 
 *  \return A handle to the scaled matrix.
 */

mat_ptr_t scalar_mult(mat_ptr_t ptr_A, mat_ptr_t ptr_B);

/*!
 *  \ingroup larc
 * 
 *  \brief Divide a matrix by a non-zero scalar.
 * 
 *  \param ptr_A A handle to a matrix of any size.
 *  \param ptr_B A handle to a matrix with a single element, which is considered the scalar,
 *         which must be non-zero.
 * 
 *  \return A handle to the scaled matrix.
 */

mat_ptr_t scalar_divide(mat_ptr_t ptr_A, mat_ptr_t ptr_B);

// static  mat_ptr_t tensor_with_identity_on_left(mat_ptr_t A_ptr, mat_level_t levelI);

/*!
 * \ingroup larc
 *
 * \brief Compute the Kronecker product of two matrices.
 * 
 * \param ptr_A A handle to the first operand, i.e., the operand whose
 * values will be replaced by scaled copies of the the second operand.
 * \param ptr_B A handle to the second operand.
 * 
 * \return A handle to the Kroncker product.
 */
mat_ptr_t kronecker_product(mat_ptr_t ptr_A, mat_ptr_t ptr_B);

/*!
 * \ingroup larc
 *
 * \brief Form a new matrix by joining two matrices of the same size
 * side by side.
 *
 * \param A_ptr: A handle to the first matrix, which will be on the left.
 * \param B_ptr: A handle to the second matrix, which will be on the right.
 *
 * \return A handle to the matrix formed by joining the two matrices.
 */
mat_ptr_t join(mat_ptr_t A_ptr, mat_ptr_t B_ptr);

/*!
 * \ingroup larc
 *
 * \brief Form a new matrix by stacking two matrices of the same size
 * side one on top of the other.
 *
 * \param A_ptr: A handle to the first matrix, which will be on top.
 * \param B_ptr: A handle to the second matrix, which will be on the bottom.
 *
 * \return A handle to the matrix formed by stacking the two matrices.
 */
mat_ptr_t stack(mat_ptr_t A_ptr, mat_ptr_t B_ptr);

/*!
 * \ingroup larc
 * \brief With unitary matrix B, change the basis in which the A matrix is expressed
 * \param B_ptr A matrix pointer (not necessarily unitary)
 * \param A_ptr A matrix pointer
 * \result A matrix pointer to C = B*A*conjugate_transpose(B)
 */
mat_ptr_t matrix_basischange_A_by_B(mat_ptr_t B_ptr, mat_ptr_t A_ptr);

/*!
 * \ingroup larc
 * \brief The standard linear algebra saxpy operation
 * \param A_ptr A matrix pointer
 * \param B_ptr A matrix pointer
 * \param scalar_a_ptr A pointer to a scalar value
 * \param scalar_b_ptr A pointer to a scalar value
 * \result A matrix pointer to C = a*A + b*B
 */
mat_ptr_t matrix_saxpy(mat_ptr_t A_ptr, mat_ptr_t B_ptr,
        mat_ptr_t scalar_a_ptr, mat_ptr_t scalar_b_ptr);

/*!
 * \ingroup larc
 * \brief calculates the dot product of two vectors
 * \param A_ptr A matrix pointer (usually a ROW_VECTOR)
 * \param B_ptr A matrix pointer (usually a COL_VECTOR)
 * \param verbose When above BASIC, much information printed to stderr
 * \result A matrix pointer to the (always) scalar result
 */
mat_ptr_t vector_dot_product(mat_ptr_t A_ptr, mat_ptr_t B_ptr, int verbose);

/*!
 * \ingroup larc
 * \brief A nonlinear operation which squares and scales the norm of every element of a matrix
 * \param ptr_A A matrix pointer
 * \param scale_ptr A matrix pointer to a scalar value
 * \result A pointer to a matrix B such that B_{ij} = scale*|A_{ij}|^2
 */
mat_ptr_t matrix_entrySquared(mat_ptr_t ptr_A, mat_ptr_t scale_ptr);

/*!
 * \ingroup larc
 * \brief An efficient routine for left-multiplying a matrix by an integer Hadamard matrix
 * \param A_ptr A matrix pointer
 * \result A pointer to a matrix B = intH*A
 */
mat_ptr_t iHadamard_times_matrix(mat_ptr_t A_ptr);

/*!
 * \ingroup larc
 * \brief An efficient routine for right-multiplying a matrix by an integer Hadamard matrix
 * \param A_ptr A matrix pointer
 * \result A pointer to a matrix B = A*intH
 */
mat_ptr_t matrix_times_iHadamard(mat_ptr_t A_ptr);

#ifndef SWIG
/*!
 * \ingroup larc
 * \brief Find the max norm of a matrix. 
 * \param mat_ptr A handle to a matrix. 
 * \param max     Parameter to return the max norm of matrix A where 
 *                |A|_{max}=max_{ij}|a_{ij}|. 
 */ 
void matrix_maxnorm(scalarType *max, mat_ptr_t mat_ptr);
/*!
 * \ingroup larc
 * \brief Find the max norm of a matrix. 
 * \param mat_ptr A handle to a matrix. 
 * \param norm A function pointer, allows the user to define their own norm 
 * \param max     Parameter to return the max norm of matrix A where 
 *                |A|_{max}=max_{ij}norm(a_{ij}). 
 */
void matrix_maxnorm_custom_norm(scalarType *max, mat_ptr_t mat_ptr, void (*norm)(scalarType*, const scalarType));

/*!
 * \ingroup larc
 * \brief Find the L2 norm of a matrix.
 * \param mat_ptr A handle to a matrix.
 * \param l2norm  Parameter to return the L2 norm of matrix A where
 *                |A|_{L2}=sqrt(sum_{ij}|a_{ij}|^2).
 */
void matrix_l2norm(scalarType *l2norm, mat_ptr_t mat_ptr);
/*!
 * \ingroup larc
 * \brief Find the L2 norm of a matrix.
 * \param mat_ptr A handle to a matrix.
 * \param norm A function pointer, allows the user to define their own norm
 * \param running_sum     Parameter to return the L2 norm of matrix A where
 *                |A|_{L2}=sqrt(sum_{ij}norm(a_{ij}))
 *                 or
 *                |A|_{L2}=sqrt(sum_{ij}|a_{ij}|^2).
 */
void matrix_sumnorm_custom_norm(scalarType *running_sum, mat_ptr_t mat_ptr, void (*norm)(scalarType*, const scalarType));

/*!
 * \ingroup larc
 * \brief Calculates the scaled tracenorm of a matrix
 * \param ret Returns a scalarType pointer to c = scale_factor*tr(A*conjugate_transpose(A))
 * \param mat_ptr A matrix pointer
 * \param scale_factor A scalarType value
 */
void matrix_tracenorm(scalarType *ret, mat_ptr_t mat_ptr, scalarType scale_factor);

/*!
 * \ingroup larc
 * \brief Find the number of entries with a given value in a matrix.
 *
 * Returns answer through sca_count parameter, which should already be initialized by 
 * mpz_init. Note that locate_entries_larcMatrixFile is probably more efficient for large matrices. 
 *
 * \param count   The count is returned through this parameter. 
 * \param mat_ptr A matrix pointer to a matrix in LARC format. 
 * \param scalar  A scalar value to count. 
 */
void matrix_count_entries(mpz_t count, mat_ptr_t mat_ptr, scalarType scalar);
/*!
 * \ingroup larc
 *
 * \brief Generate a random matrix by randomly distributing a given 
 * number of scalars into quadrants of zeros. 
 *
 * \param row_level row level of matrix
 * \param col_level col level of matrix
 * \param scalarNum   number of entries with value one desired in matrix
 * \param scalarPtr   matrix pointer to scalar value - increases efficiency.
 *
 * \return matrix pointer to matrix
 */
mat_ptr_t random_bool_matrix_from_count(mat_level_t row_level, mat_level_t col_level, mpz_t scalarNum, mat_ptr_t scalarPtr);
#endif


/*!
 * \ingroup larc
 *
 * \brief Allocates and fills an array with the distinct scalar entries in a
 * matrix from its LARCMatrix json file. 
 *
 * \param scalars_ptr Pointer to hold an array of scalars. 
 * \param numScalars  Pointer to hold the count of scalars. 
 * \param path        Path to the LARCMatrix file.
 */
void matrix_get_scalars_larcMatrix(scalarType **scalars_ptr, int64_t *numScalars, char *path);

// user interface (such as python) work with matrixIDs
/*!
 * \ingroup larc
 * \brief Add two matrices
 * \param A_mID A handle to the first matrix
 * \param B_mID A handle to the second matrix
 * \return A handle to the matrix that is the sum of the matrices
 * referenced by the arguments.
 */
int64_t matrix_add_matrixID(int64_t A_mID, int64_t B_mID);
/*!
 * \ingroup larc
 * \brief Subtract matrix B from matrix A
 * \param A_mID A handle to the first matrix ('minuend')
 * \param B_mID A handle to the second matrix ('subtrahend')
 * \return A handle to the matrix that is the difference of the matrices
 * referenced by the arguments.
 */
int64_t matrix_diff_matrixID(int64_t A_mID, int64_t B_mID);

/*!
 * \ingroup larc
 *
 * \brief Multiply two matrices. Wraps matrix_mult_clean to set no cleaning
 * or for backwards compatibility. 
 *
 * \param A_mID A handle to the first matrix.
 * \param B_mID A handle to the second matrix.
 *
 * \return A handle to the matrix that is the product of the matrices 
 * referenced by the arguments.
 */
int64_t matrix_mult_matrixID(int64_t A_mID, int64_t B_mID);

/*!
 * \ingroup larc
 *
 * \brief Multiply two matrices.
 *
 * Empty matrix store and op store after all internal recursive
 * calls to matrix_mult where input matrices have row or column
 * level equal to cleanThresh. Also remove extra matrices produced
 * during recursive calls to matrix_mult where input matrices have row
 * or column level above cleanThresh. For efficiency, cleanThresh should be
 * as large possible such that the calculation can finish. To set to no cleaning, 
 * set cleanThresh to higher than the levels of the input matrices. 
 *
 * \param A_mID A handle to the first matrix.
 * \param B_mID A handle to the second matrix.
 * \param cleanThresh A matrix level threshold. 
 *
 * \return A handle to the matrix that is the product of the matrices 
 * referenced by the arguments.
 */
int64_t matrix_mult_clean_matrixID(int64_t A_mID, int64_t B_mID, mat_level_t cleanThresh);

/*!
 * \ingroup larc
 * 
 * \brief Multiply a matrix by a scalar.
 * 
 * \param A_mID A handle to a matrix with a single element, which is considered the scalar.
 * \param B_mID A handle to a matrix of any size.
 * 
 * \return A handle to the scaled matrix.
 */
int64_t scalar_mult_matrixID(int64_t A_mID, int64_t B_mID);

/*!
 * \ingroup larc
 * 
 * \brief Divide a matrix by a non-zero scalar.
 * 
 * \param A_mID A handle to a matrix of any size.
 * \param B_mID A handle to a matrix with a single element, which is considered the scalar,
 *        which must be non-zero.
 * 
 * \return A handle to the scaled matrix.
 */
int64_t scalar_divide_matrixID(int64_t A_mID, int64_t B_mID);

/*!
 * \ingroup larc
 *
 * \brief Compute the Kronecker product of two matrices.
 * 
 * \param A_mID A handle to the first operand, i.e., the operand whose
 * values will be replaced by scaled copies of the the second operand.
 * \param B_mID A handle to the second operand.
 * 
 * \return A handle to the Kroncker product.
 */
int64_t kronecker_product_matrixID(int64_t A_mID, int64_t B_mID);

/*!
 * \ingroup larc
 *
 * \brief Form a new matrix by joining two matrices of the same size
 * side by side.
 *
 * \param A_mID A handle to the first matrix, which will be on the left.
 * \param B_mID A handle to the second matrix, which will be on the right.
 *
 * \return A handle to the matrix formed by joining the two matrices.
 */
int64_t join_matrixID(int64_t A_mID, int64_t B_mID);

/*!
 * \ingroup larc
 *
 * \brief Form a new matrix by stacking two matrices of the same size
 * side one on top of the other.
 *
 * \param A_mID A handle to the first matrix, which will be on top.
 * \param B_mID A handle to the second matrix, which will be on the bottom.
 *
 * \return A handle to the matrix formed by stacking the two matrices.
 */
int64_t stack_matrixID(int64_t A_mID, int64_t B_mID);

/*!
 * \ingroup larc
 * \brief A nonlinear operation which squares and scales the norm of every element of a matrix
 * \param m_mID A matrixID 
 * \param scale_factor A scalarType value in string format
 * \result A matrixID for a matrix B such that B_{ij} = scale*|A_{ij}|^2
 */
int64_t matrix_entrySquared_matrixID(int64_t m_mID, char *scale_factor);

/*!
 * \ingroup larc
 * \brief An efficient routine for left-multiplying a matrix by an integer Hadamard matrix
 * \param A_mID A matrixID
 * \result A matrixID for a matrix B = intH*A
 */
int64_t iHadamard_times_matrix_matrixID(int64_t A_mID);

/*!
 * \ingroup larc
 * \brief With unitary matrix B, change the basis in which the A matrix is expressed
 * \param B_mID A matrixID (matrix not necessarily unitary)
 * \param A_mID A matrixID
 * \result A matrixID for C = B*A*conjugate_transpose(B)
 */
int64_t matrix_basischange_A_by_B_matrixID(int64_t B_mID, int64_t A_mID);

/*!
 * \ingroup larc
 * \brief The standard linear algebra saxpy operation
 * \param A_mID A matrixID
 * \param B_mID A matrixID
 * \param scalar_a_mID A matrixID for a scalar value
 * \param scalar_b_mID A matrixID for a scalar value
 * \result A matrixID for C = a*A + b*B
 */
int64_t matrix_saxpy_matrixID(int64_t A_mID, int64_t B_mID, 
        int64_t scalar_a_mID, int64_t scalar_b_mID);

/*!
 * \ingroup larc
 * \brief calculates the dot product of two vectors
 * \param A_mID A matrixID (usually for a ROW_VECTOR)
 * \param B_mID A matrixID (usually for a COL_VECTOR)
 * \param verbose When above BASIC, much information printed to stderr
 * \result A matrixID for the (always) scalar result
 */
int64_t vector_dot_product_matrixID(int64_t A_mID, int64_t B_mID, int verbose);

/*!
 * \ingroup larc
 * \brief An efficient routine for right-multiplying a matrix by an integer Hadamard matrix
 * \param A_mID A matrixID
 * \result A matrixID for a matrix B = A*intH
 */
int64_t matrix_times_iHadamard_matrixID(int64_t A_mID);

/*!
 * \ingroup larc
 * \brief Find the max norm of a matrix. 
 * \param mat_mID A handle to a matrix. 
 * \return A scalarType in string format containing the max norm of matrix A where 
 *                |A|_{max}=max_{ij}|a_{ij}|. 
 */
char *matrix_maxnorm_matrixID(int64_t mat_mID);

/*!
 * \ingroup larc
 * \brief Find the L2 norm of a matrix.
 * \param mat_mID A handle to a matrix.
 * \return A scalarType in string format containing the L2 norm of matrix A where
 *                |A|_{L2}=sqrt(sum_{ij}|a_{ij}|^2).
 */
char *matrix_l2norm_matrixID(int64_t mat_mID);

/*!
 * \ingroup larc
 * \brief Calculates the scaled tracenorm of a matrix
 * \param m_mID A matrixID
 * \param scale_factor A scalarType value
 * \return A scalarType in string format containing c = scale_factor*tr(A*conjugate_transpose(A))
 */
char *matrix_tracenorm_matrixID(int64_t m_mID, char *scale_factor);

/*!
 * \ingroup larc
 * \brief Returns the trace of a matrix
 * \param m_mID A matrixID
 * \return The scalarType value of the trace in string format
 */
char *matrix_trace_from_matrixID(int64_t m_mID);

/*!
 * \ingroup larc
 * \brief Returns the scalarType value of a stored scalar
 * \param m_mID A matrixID
 * \return The scalarType value of the scalar in string format
 */
char* scalar_string_from_matrixID(int64_t m_mID);

/*!
 * \ingroup larc
 *
 * \brief Find the number of entries with a given value in a matrix. 
 *
 * \param mat_mID A matrixID for a matrix. 
 * \param scalar_str A scalar value to count. 
 *
 * \return A string of the multi-precision count. 
 */
char *matrix_count_entries_matrixID(int64_t mat_mID, char *scalar_str);

/*!
 * \ingroup larc
 *
 * \brief Returns a printed list of the distinct scalar entries in a LARCMatrix json file. 
 *
 * \param path    The path to the LARCmatrix file. 
 *
 * \return A string of the scalar entries. 
 */
char *matrix_list_scalars_larcMatrix(char *path);

/*!
 * \ingroup larc
 *
 * \brief Generate a random binary matrix by randomly distributing a given 
 * number of 1's into quadrants. 
 *
 * \param row_level row level of matrix
 * \param col_level col level of matrix
 * \param numOnes   number of entries with value one desired in matrix, as string.
 *
 * \return matrixID of the matrix
 */
int64_t random_bool_matrixID_from_count(mat_level_t row_level, mat_level_t col_level, char* numOnes);

#ifdef __cplusplus
}
#endif // #ifdef __cplusplus

#endif
