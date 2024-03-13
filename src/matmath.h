//                      matmath.h
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

#ifndef SWIG

/*!
 * \ingroup larc
 *
 * \brief Allocates and fills an array with the distinct scalar entries in a
 * matrix from its LARCMatrix json file, in order of appearance. 
 *
 * \param scalars_ptr The returned list of scalars found in the file
 * \param numScalars  The number of distinct scalars in our list
 * \param path        Path to the LARCMatrix file.
 */
void get_array_of_scalars_in_larcMatrixFile(scalarType **scalars_ptr, int64_t *numScalars, char *path);
#endif // ifndef SWIG

// user interface (such as python) works with packedIDs
/*!
 * \ingroup larc
 * \brief Add two matrices
 * \param A_pID The packedID of the first matrix
 * \param B_pID The packedID of the second matrix
 * \return The packedID of the matrix that is the sum of the matrices
 * referenced by the arguments.
 */
int64_t matrix_add(int64_t A_pID, int64_t B_pID);

/*!
 * \ingroup larc
 * \brief Subtract matrix B from matrix A
 * \param A_pID The packedID of the first matrix ('minuend')
 * \param B_pID The packedID of the second matrix ('subtrahend')
 * \return The packedID of the matrix that is the difference of the matrices
 * referenced by the arguments.
 */
int64_t matrix_diff(int64_t A_pID, int64_t B_pID);

/*!
 * \ingroup larc
 *
 * \brief Multiply two matrices. Wraps matrix_mult_clean to set no cleaning
 * or for backwards compatibility. 
 *
 * \param A_pID The packedID of the first matrix.
 * \param B_pID The packedID of the second matrix.
 *
 * \return The packedID of the matrix that is the product of the matrices 
 * referenced by the arguments.
 */
int64_t matrix_mult(int64_t A_pID, int64_t B_pID);

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
 * \param A_pID The packedID of the first matrix.
 * \param B_pID The packedID of the second matrix.
 * \param cleanThresh A matrix level threshold. 
 *
 * \return The packedID of the matrix that is the product of the matrices 
 * referenced by the arguments.
 */
int64_t matrix_mult_clean(int64_t A_pID, int64_t B_pID, mat_level_t cleanThresh);

/*!
 * \ingroup larc
 * 
 * \brief Multiply a matrix by a scalar.
 * 
 * \param A_pID The packedID of a matrix with a single element, which is considered the scalar.
 * \param B_pID The packedID of a matrix of any size.
 * 
 * \return The packedID of the scaled matrix.
 */
int64_t scalar_mult(int64_t A_pID, int64_t B_pID);

/*!
 * \ingroup larc
 * 
 * \brief Divide a matrix by a non-zero scalar.
 * 
 * \param A_pID The packedID of a matrix of any size.
 * \param B_pID The packedID of a matrix with a single element, which is considered the scalar,
 *        which must be non-zero.
 * 
 * \return The packedID of the scaled matrix.
 */
int64_t scalar_divide(int64_t A_pID, int64_t B_pID);

/*!
 * \ingroup larc
 *
 * \brief Compute the Kronecker product of two matrices.
 * 
 * \param A_pID The packedID of the first operand, i.e., the operand whose
 * values will be replaced by scaled copies of the the second operand.
 * \param B_pID The packedID of the second operand.
 * 
 * \return The packedID of the Kroncker product.
 */
int64_t kronecker_product(int64_t A_pID, int64_t B_pID);

/*!
 * \ingroup larc
 *
 * \brief Form a new matrix by joining two matrices of the same size
 * side by side.
 *
 * \param A_pID The packedID of the first matrix, which will be on the left.
 * \param B_pID The packedID of the second matrix, which will be on the right.
 *
 * \return The packedID of the matrix formed by joining the two matrices.
 */
int64_t join(int64_t A_pID, int64_t B_pID);

/*!
 * \ingroup larc
 *
 * \brief Form a new matrix by stacking two matrices of the same size
 * side one on top of the other.
 *
 * \param A_pID The packedID of the first matrix, which will be on top.
 * \param B_pID The packedID of the second matrix, which will be on the bottom.
 *
 * \return The packedID of the matrix formed by stacking the two matrices.
 */
int64_t stack(int64_t A_pID, int64_t B_pID);

/*!
 * \ingroup larc
 * \brief A nonlinear operation which squares and scales the norm of every element of a matrix
 * \param m_pID The packedID of the matrix
 * \param scale_factor A scalarType value in string format
 * \return A handle for a matrix B such that B_{ij} = scale*|A_{ij}|^2
 */
int64_t matrix_entrySquared(int64_t m_pID, char *scale_factor);

/*!
 * \ingroup larc
 * \brief An efficient routine for left-multiplying a matrix by an integer Hadamard matrix
 * \param A_pID The packedID of the matrix
 * \return A handle for a matrix B = intH*A
 */
int64_t iHadamard_times_matrix(int64_t A_pID);

/*!
 * \ingroup larc
 * \brief With unitary matrix B, change the basis in which the A matrix is expressed
 * \param B_pID The packedID of a matrix (matrix not necessarily unitary)
 * \param A_pID The packedID of a matrix
 * \return A handle for C = B*A*conjugate_transpose(B)
 */
int64_t matrix_basischange_A_by_B(int64_t B_pID, int64_t A_pID);

/*!
 * \ingroup larc
 * \brief The standard linear algebra saxpy operation
 * \param A_pID A packedID
 * \param B_pID A packedID
 * \param scalar_a_pID A packedID for a scalar value
 * \param scalar_b_pID A packedID for a scalar value
 * \return A packedID for C = a*A + b*B
 */
int64_t matrix_saxpy(int64_t A_pID, int64_t B_pID, 
        int64_t scalar_a_pID, int64_t scalar_b_pID);

/*!
 * \ingroup larc
 * \brief calculates the dot product of two vectors
 * \param A_pID A packedID (usually for a ROW_VECTOR)
 * \param B_pID A packedID (usually for a COL_VECTOR)
 * \param verbose When above BASIC, much information printed to stderr
 * \return A packedID for the (always) scalar result
 */
int64_t vector_dot_product(int64_t A_pID, int64_t B_pID, int verbose);

/*!
 * \ingroup larc
 * \brief An efficient routine for right-multiplying a matrix by an integer Hadamard matrix
 * \param A_pID A packedID
 * \return A packedID for a matrix B = A*intH
 */
int64_t matrix_times_iHadamard(int64_t A_pID);

/*!
 * \ingroup larc
 * \brief Calculates the scaled tracenorm of a matrix
 * \param m_pID A packedID
 * \param scale_factor A scalarType value
 * \return A scalarType in string format containing c = scale_factor*tr(A*conjugate_transpose(A))
 */
char *tracenorm(int64_t m_pID, char *scale_factor);

/*!
 * \ingroup larc
 * \brief Returns the trace of a matrix
 * \param m_pID A packedID
 * \return The scalarType value of the trace in string format
 */
char *traceID(int64_t m_pID);

/*!
 * \ingroup larc
 * \brief Returns the scalarType value of a stored scalar
 * \param m_pID A packedID
 * \return The scalarType value of the scalar in string format
 */
char *get_scalar_value_string(int64_t m_pID);

/*!
 * \ingroup larc
 *
 * \brief Find the number of entries with a given value in a matrix. 
 *
 * \param mat_pID A packedID for a matrix. 
 * \param scalar_str A scalar value to count. 
 *
 * \return A string of the multi-precision count. 
 */
char *matrix_count_entries(int64_t mat_pID, char *scalar_str);

/*!
 * \ingroup larc
 *
 * \brief Returns a printed list of the distinct scalar entries in a LARCMatrix json file. 
 *
 * \param path    The path to the LARCmatrix file. 
 *
 * \return A string of the scalar entries. 
 */
char *get_list_of_scalars_in_larcMatrixFile(char *path);

/*!
 * \ingroup larc
 *
 * \brief Find the (row,col) coordinates of every entry with a given value in a
 * a matrix from its LARCMatrix file.
 *
 * Returns a nx2 array where entry
 *      row(ith matching entry) = array[i][0], and
 *      col(ith matching entry) = array[i][1],
 * where n is the number of entries with the given value. Inner dimension
 * boundaries of the array are marked with -1's and the outer dimension
 * boundary is marked with a NULL so that SWIG can convert to a list of list
 * of ints.
 *
 * \param path  A path to a LARCmatrix json file
 * \param scalar_str A scalar value to count, given as a string.
 * \param maxno The maximum number of locations to be reported.
 *
 * \return The array of coordinates.
 */
int64_t **locate_entries_larcMatrixFile(char *path, char *scalar_str, unsigned maxno);


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
 * \return packedID of the matrix
 */
int64_t random_Boolean_matrix_from_count(mat_level_t row_level, mat_level_t col_level, char* numOnes);

/*!
 * \ingroup larc
 *
 * \brief Return the packedID for the matrix element of largest norm value in a given matrix
 *
 * \param mat_pID The packedID of the input matrix.
 * \return The packedID for the scalar element in the matrix with the largest norm value.
 */
int64_t matrix_element_with_maxNorm(int64_t mat_pID);

/*!
 * \ingroup larc
 *
 * \brief Returns the packedID of the norm of the input matrix
 * \param mat_pID The packedID of the matrix for which we want the norm
 * \param whichNorm Input 0 for L_infty norm, 1 for L_1 norm, 2 for L_2 norm; any other input value is treated as 0.
 * \return The packedID for the norm value for the input matrix
 */
int64_t normID(int64_t mat_pID, int whichNorm);

/*!
 * \ingroup larc
 * \brief Generates a matrix with all elements equal to a constant
 * \param constID The packedID of the constant value
 * \param row_level The log of the row dimension of the desired matrix
 * \param col_level The log of the column dimension of the desired matrix
 * \return The packedID of a level-(row_level,col_level) matrix filled by the constant value
 */
int64_t create_const_matrix(int64_t constID, mat_level_t row_level,
	mat_level_t col_level);

/*!
 * \ingroup larc
 * \brief Applies a user-defined function to every element of a matrix
 *
 * The user-defined function can be memoized, preferably using one of the 
 * operation_type values FUNC_A, FUNC_B or FUNC_C.
 *
 * \param m_pID The packedID of the input matrix
 * \param func The function to apply to the scalarType values
 * \param op_memz A number corresponding to an operation_type enum value
 * \return The packedID of the output matrix
 */
int64_t apply_function_to_matrix_values(int64_t m_pID,
        void (*func)(scalarType*, const scalarType), op_type_t op_memz);


#ifdef __cplusplus
}
#endif // #ifdef __cplusplus

#endif
