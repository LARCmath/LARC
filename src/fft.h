//                      fft.h
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

#ifndef LARC_FFT_H
#define LARC_FFT_H

#ifdef __cplusplus
extern "C" {
#endif 

#include <stdint.h>       // int64_t 

#include "larc.h"         // mat_level_t
#include "matrix_store.h" // mats_ptr_t, matns_ptr_t
#include "scalars.h"      // conversion function

/*!
 * \ingroup larc
 * \brief Generates or finds the inverse shuffle permutation matrix of the given level
 * \param m_lev The level of the inverse shuffle matrix to generate or find
 * \return The packedID for the inverse shuffle matrix
 */
int64_t create_invShufMat(mat_level_t m_lev);

/*!
 * \ingroup larc
 * \brief Generates or finds the diagonal D matrix of the desired level, for use in a sparse block recursive DFT
 * \param m_lev The level of the D matrix to generate or find
 * \return The packedID for the D matrix
 */
int64_t create_FFT_DMat(mat_level_t m_lev);

/*!
 * \ingroup larc
 * \brief Prints the (2^pow)-th roots of unity 
 * \param pow The power of two to use 
 * \return 0 on success
 */
int print_pow2_roots_unity(uint32_t pow);

/*!
 * \ingroup larc
 * \brief  Returns the principal (2^pow)-th root of unity, e^{i 2 pi/2^pow} and either loads or finds it in the matrix_store 
 * \param pow The log base 2 of the power of unity to generate
 * \return The packedID for the generated root of unity
 */
int64_t principal_pow2_root_unity_pID(uint32_t pow);

/*!
 * \ingroup larc
 * \brief Generates or finds the block C matrix of the desired level, for use in a sparse block recursive DFT
 * \param m_lev The level of the C matrix (to generate or find)
 * \return The packedID for the C matrix
 */
int64_t create_FFT_CMat(mat_level_t m_lev);

/*!
 * \ingroup larc
 * \brief Generates an FFT matrix of the desired size
 * \param m_lev The level of the matrix generated
 * \return The packedID for the generated matrix
 */
int64_t create_FFTMat(mat_level_t m_lev);

/*!
 * \ingroup larc
 * \brief Returns the packedID for larc matrix containing the scalar e^{i 2 pi k/n}, the  k-th power of n-th principal root-of-unity
 * \param k The integer exponent to which to put the n-th root of unity
 * \param n The integer specifying the n-th root-of-unity
 * \param verbose
 * \return The packedID for the larc matrix record for scalar e^{i 2 pi k/n}
 */
int64_t k_th_power_of_n_th_root_of_unity_pID(int k, int n, int verbose);
  

#ifdef __cplusplus
}
#endif // #ifdef __cplusplus


#endif
