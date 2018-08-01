//                      matmath.h
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

#ifndef LARC_MATMATH_H
#define LARC_MATMATH_H

#ifdef __cplusplus
extern "C" {
#endif // #ifdef __cplusplus

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
mat_add_t matrix_add(mat_add_t ptr_A, mat_add_t ptr_B);

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
mat_add_t matrix_diff(mat_add_t ptr_A, mat_add_t ptr_B);

/*!
 * \ingroup larc
 *
 * \brief Multiply two matrices.
 *
 * \param ptr_A A handle to the first matrix.
 * \param ptr_B A handle to the second matrix.
 *
 * \return A handle to the matrix that is the product of the matrices 
 * referenced by the arguments.
 */
mat_add_t matrix_mult(mat_add_t ptr_A, mat_add_t ptr_B);

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

mat_add_t scalar_mult(mat_add_t ptr_a, mat_add_t ptr_B);
// static  mat_add_t tensor_with_identity_on_left(mat_add_t A_ptr, mat_level_t levelI);

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
mat_add_t kronecker_product(mat_add_t ptr_A, mat_add_t ptr_B);

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
mat_add_t join(mat_add_t A_ptr, mat_add_t B_ptr);

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
mat_add_t stack(mat_add_t A_ptr, mat_add_t B_ptr);

mat_add_t matrix_entrySquared(mat_add_t ptr_A, double scale_factor,
        mat_add_t scale_ptr);
// user interface (such as python) work with matrixIDs
int64_t matrix_add_matrixID(int64_t A_mID, int64_t B_mID);
int64_t matrix_diff_matrixID(int64_t A_mID, int64_t B_mID);
int64_t matrix_mult_matrixID(int64_t A_mID, int64_t B_mID);
int64_t scalar_mult_matrixID(int64_t A_mID, int64_t B_mID);
int64_t kronecker_product_matrixID(int64_t A_mID, int64_t B_mID);
int64_t join_matrixID(int64_t A_mID, int64_t B_mID);
int64_t stack_matrixID(int64_t A_mID, int64_t B_mID);
int64_t matrix_entrySquared_matrixID(int64_t m_mID, double scale_factor);

mat_add_t iHadamard_times_matrix(mat_add_t A_ptr);
int64_t iHadamard_times_matrix_matrixID(int64_t A_mID);

mat_add_t matrix_times_iHadamard(mat_add_t A_ptr);
int64_t matrix_times_iHadamard_matrixID(int64_t A_mID);

ScalarType matrix_tracenorm(mat_add_t mat_ptr, ScalarType scale);

/*!
 * \ingroup larc
 * 
 * \brief Find the sparcity of a matrix. 
 * 
 * \param mat_ptr A handle to a matrix. 
 *
 * \return The sparcity of matrix \f$A_{n\times m}\f$ where sparcity is 
 * defined to be \f[\frac{|\{(i,j):a_{i,j}=0\}|}{m\times n}.\f] 
 */ 
double matrix_sparsity(mat_add_t mat_ptr);
double matrix_sparsity_matrixID(int64_t mat_mID);

/*!
 * \ingroup larc
 * 
 * \brief Find the max norm of a matrix. 
 * 
 * \param mat_ptr A handle to a matrix. 
 *
 * \return The max norm of matrix \f$A\f$ where 
 * \f[\|A\|_{\max }=\max _{ij}|a_{ij}|.\f] 
 */ 
double matrix_maxnorm_matrixID(int64_t mat_mID);
double matrix_maxnorm(mat_add_t mat_ptr);

#ifdef __cplusplus
}
#endif // #ifdef __cplusplus

#endif
