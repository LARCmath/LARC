//                      fft.h
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

#ifndef LARC_FFT_H
#define LARC_FFT_H

#include "larc.h"
#include "global.h"
#include "op_store.h"
#include "matrix_store.h"
#include "organize.h"
#include "matmath.h"


/* The create_perm_inv(mat_level_t n) function
   takes the size of the matrix and returns a matrix PTR
   corresponding the inverse permution matrix */
mat_add_t create_perm_inv(mat_level_t n);

/* Python Interface for the create_perm_inv(mat_level_t n) function.
   This takes the size of the matrix and returns a matrix ID
   corresponding the inverse permution matrix */
int64_t create_perm_inv_matrixID(mat_level_t n);


/* Python Interface for the create_fft_D(mat_level_t n) function.
   This takes the size of the matrix and returns a matrix ID
   corresponding the D matrix needed for the fft recursion */
int64_t create_fft_D_matrixID(mat_level_t n);

/* The create_fft_D(mat_level_t n) function
   takes the size of the matrix and returns a matrix PTR
   corresponding the D matrix for the fft recursion*/
mat_add_t create_fft_D(mat_level_t n); 


/* prints the (2^n)-th roots of unity */
int print_pow2_roots_unity(mat_level_t n);

/* returns the matrixID of the principal (2^n)-th root of unity */ 
int64_t principal_pow2_root_unity(mat_level_t n);

/* return the matrix pointer for the first (2^n)-th root of unity */
mat_add_t get_first_pow2_root_unity(mat_level_t n);

/* Python Interface for the create_fft_C(mat_level_t n) function.
   This takes the size of the matrix and returns a matrix ID
   corresponding the C matrix needed for the fft recursion */
int64_t create_fft_C_matrixID(mat_level_t n);

/* The create_fft_C(mat_level_t n) function
   takes the size of the matrix and returns a matrix PTR
   corresponding the C matrix for the fft recursion*/
mat_add_t create_fft_C(mat_level_t n); 


/* Python Interface for the create_fft_matrix(mat_level_t n) function.
   This takes the size of the matrix and returns a matrix ID
   corresponding to the fft matrix */
int64_t create_fft_matrix_matrixID(mat_level_t n);

/* The create_fft_matrix(mat_level_t n) function
   takes the size of the matrix and returns a matrix PTR
   corresponding to the fft matrix */
mat_add_t create_fft_matrix(mat_level_t n); 


#endif
