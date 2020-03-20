//                       experimental.h 
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

#endif


#ifdef __cplusplus
}
#endif

#endif
