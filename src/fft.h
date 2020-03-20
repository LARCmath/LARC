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
#include "matrix_store.h" // mat_ptr_t

/* The create_perm_inv(mat_level_t n) function
   takes the size of the matrix and returns a matrix PTR
   corresponding the inverse permution matrix */
mat_ptr_t create_perm_inv(mat_level_t n);

/* The create_fft_D(mat_level_t n) function
   takes the size of the matrix and returns a matrix PTR
   corresponding the D matrix for the fft recursion*/
mat_ptr_t create_fft_D(mat_level_t n); 

/* return the matrix pointer for the first (2^n)-th root of unity */
mat_ptr_t get_first_pow2_root_unity(mat_level_t n);

/* The create_fft_C(mat_level_t n) function
   takes the size of the matrix and returns a matrix PTR
   corresponding the C matrix for the fft recursion*/
mat_ptr_t create_fft_C(mat_level_t n); 

/* The create_fft_matrix(mat_level_t n) function
   takes the size of the matrix and returns a matrix PTR
   corresponding to the fft matrix */
mat_ptr_t create_fft_matrix(mat_level_t n); 

/* Python Interface for the create_perm_inv(mat_level_t n) function.
   This takes the size of the matrix and returns a matrix ID
   corresponding the inverse permution matrix */
int64_t create_perm_inv_matrixID(mat_level_t n);

/* Python Interface for the create_fft_D(mat_level_t n) function.
   This takes the size of the matrix and returns a matrix ID
   corresponding the D matrix needed for the fft recursion */
int64_t create_fft_D_matrixID(mat_level_t n);

/* prints the (2^n)-th roots of unity */
int print_pow2_roots_unity(mat_level_t n);

/* returns the matrixID of the principal (2^n)-th root of unity */ 
int64_t principal_pow2_root_unity(mat_level_t n);

/* Python Interface for the create_fft_C(mat_level_t n) function.
   This takes the size of the matrix and returns a matrix ID
   corresponding the C matrix needed for the fft recursion */
int64_t create_fft_C_matrixID(mat_level_t n);

/* Python Interface for the create_fft_matrix(mat_level_t n) function.
   This takes the size of the matrix and returns a matrix ID
   corresponding to the fft matrix */
int64_t create_fft_matrix_matrixID(mat_level_t n);

/* prints the n-th roots of unity */
int print_n_th_roots_of_unity(int n, int verbose);
  
/* return the matrix ID for the principal n-th root of unity */
int64_t principal_n_th_root_of_unity_matID(int n, int verbose);

/* return matrix ID for k-th power of the n-th principal root of unity */
int64_t k_th_power_of_n_th_root_of_unity_matID(int k, int n, int verbose);
  

#ifdef __cplusplus
}
#endif // #ifdef __cplusplus


#endif
