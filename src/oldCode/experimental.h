//                       experimental.h 
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

#ifndef EXPERIMENTAL_H
#define EXPERIMENTAL_H

// Standard Libaries
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif


/*!
 * \ingroup larc
 *
 * \brief Find the number of entries with a given value in a matrix from its 
 * LARCmatrix file. 
 *
 * \param larcMatrixFile_path A path to a LARCmatrix json file
 * \param scalar_str A scalar value to count. 
 *
 * \return A string of the multi-precision count. 
 */
char *count_entries_larcMatrixFile(char *larcMatrixFile_path, char *scalar_str);

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
 * \param maxno_str The maximum number of locations to be reported, as a string
 *
 * \return The array of coordinates. 
 */
int64_t **locate_entries_larcMatrixFile(char *path, char *scalar_str, char *maxno_str);

#ifndef SWIG

/*!
 * \ingroup larc
 *
 * \brief Find the number of entries with a given value in a LARCmatrix from its json file.
 *
 * Returns answer through sca_count parameter, which should already 
 * be initialized by mpz_init. 
 *
 * \param sca_count The count is returned through this parameter. 
 * \param larcMatrixFile_path A path to a LARCmatrix json file
 * \param scalar    A scalar value to count. 
 */
void mpz_count_entries_larcMatrixFile(mpz_t sca_count, char *larcMatrixFile_path, scalarType scalar);

/*!
 * \ingroup larc
 *
 * \brief Find the (row,col) coordinates of every entry with a given value in a
 * matrix from its LARCMatrix file.
 *
 * Returns answer through locations parameter, which points to a nx2 array
 * where entry 
 *      row(ith matching entry) = array[i][0], and 
 *      col(ith matching entry) = array[i][1],
 * where n is the number of entries with the given value. Also returns n 
 * through count paramter. 
 *
 * \param locations_ptr The array of coordinates is returned though this pointer.
 * \param count         The count is returned through this parameter. 
 * \param path          A path to a LARCMatrix json file
 * \param scalar        A scalar value to count. 
 * \param maxno         The maximum number of locations to report
 */
void mpz_locate_entries_larcMatrixFile(int64_t ***locations_ptr, mpz_t count, char *path, scalarType scalar, unsigned maxno);

/*!
 * \ingroup larc
 * \brief Reads a LARCMatrix formatted compressed matrix from a json file and
 * applies func to scalar entries as they are read in.
 *
 * The function func() is applied to each scalar in the file before passing
 * it to the routines which put it into the matrix store (or find that the
 * modified value is already there...) This function is only used in the
 * possibly obsolete code matrix_read_larcMatrix_file_flatten_small_scalars_PTR
 * and whatever routines call it.
 *
 * \param path The filename of the LARCMatrix file to be read into LARC.
 * \param func A function that will be applied to each scalar before storing it
 * \return The pointer to the stored matrix.
 */
mat_ptr_t read_and_alter_vals_larcMatrix_file_return_matPTR(char *path, void (*func)(scalarType*, const scalarType));

/*!
 * \ingroup larc
 * \brief Probably-obsolete function which reads in a matrix from a LARCMatrix json file, turning small scalars to zero before storing the matrix
 * \param path The location of the json file to be read
 * \param threshold The smallest scalarType value to not be zeroed
 * \return The pointer to the stored matrix
 */
mat_ptr_t matrix_read_larcMatrix_file_zeroize_small_scalars_PTR(char *path, scalarType threshold);

/*!
 * \ingroup larc
 * \brief Probably-obsolete function which reads in a matrix from a LARCMatrix json file and performs a flattening function to the data before storage
 * \param path The location of the json file to be read
 * \param threshold
 * \param flatValue
 * \return The pointer to the stored matrix
 */
mat_ptr_t matrix_read_larcMatrix_file_flatten_small_scalars_PTR(char *path, scalarType threshold, scalarType flatValue);

/*!
 * \ingroup larc
 * \brief Stores a formatted compressed LARCmatrix to a json file and applies
 * func to scalar entries before they are written.
 *
 * The function *func is applied to each scalar as it is written out to the
 * LARCMatrix json file. The value of the scalar stored in the matrixStore
 * is unchanged.
 *
 * This code has been removed from the main LARC code in io.c, as it was never
 * really used. We are preserving it here in case the idea turns out to be
 * useful at some later date.
 *
 * \param m_ptr The pointer to the matrix to be written
 * \param path The file to which the LARCmatrix will be written
 * \param func A function that will be applied to each scalar before writing
 * \return A value indicating success or the type of error
 */
int writePTR_and_alter_vals_larcMatrix_file(mat_ptr_t m_ptr, char *path,
        void (*func)(scalarType*, const scalarType));

#endif

// TODO: THE FOLLOWING FUNCTIONS NEED BETTER DOXYGEN HEADERS IF WE KEEP THEM
/*!
 * \ingroup larc
 * \brief Wrapper for probably-obsolete function
 * \param path The location of the json file to be read
 * \param threshold
 * \return The matrixID of the stored matrix
 */
int64_t matrix_read_larcMatrix_file_zeroize_small_scalars_matrixID(char *path, char *threshold);

/*!
 * \ingroup larc
 * \brief Wrapper for probably-obsolete function
 * \param path The location of the json file to be read
 * \param threshold
 * \param flatValue
 * \return The matrixID of the stored matrix
 */
int64_t matrix_read_larcMatrix_file_flatten_small_scalars_matrixID(char *path, char *threshold, char *flatValue);


#ifdef __cplusplus
}
#endif

#endif
