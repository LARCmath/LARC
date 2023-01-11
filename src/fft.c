//                          fft.c 
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

// From Van Loan:
   /*   F_k = C_k * (I_1 @ F_(k-1) ) * PI_k */
   /* where  */
   /*   PI_k is the 2^k by 2^k inverse shuffle matrix  */
   /*       (its transverse is its inverse and would shuffle a column vector) */
   /*       PI_2 = (1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1) */  //swap gate
   /*       PI_3 = (10000000; 00100000; 00001000; 00000010;  */
   /*              01000000; 00010000; 00000100; 00000001)  */
   /*   C_k is the matrix constructed by having  */
   /*       UL = I_(k-1), UR = D_(k-1);  */
   /*       LL = I_(k-1), LR = -D_(k-1)  */
   /*   D_k is the diagonal matrix which had 1, w, w^2, ....w^(2^k - 1)  */
   /*       on the diagonal, where w is the root of unity cooresponding  */
   /*       to the next largest size fourier transform, that is  */
   /*       the 2^(k+1)th root unity.   e.g. D_1 = (1, 0: 0 , i); */
   /* We end up getting a formula for the FFT in terms of a product of  */
   /* block matrices: */

   /*   F_k =    (product_for j= 0   to k-1)  I_j @ C_(k-j)   */
   /*                                       * I_k @ F_0 */
   /*          * (product_for j= k-1 to 0  )  I_j @ PI_(k-j)  */

   /* Since F_0 = I_0 and PI_1 = I_1  two of these terms are identity  */
   /* matrices (also interesting to note:  C_1 = H), so we have: */

   /*   F_k =    (product_for j= 0   to k-1)  I_j @ C_(k-j)   */
   /*          * (product_for j= k-2 to 0  )  I_j @ PI_(k-j)  */

   /* where @ is the tensor product, and * is matrix multiply. */

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

// TWO NEW INCLUDES FOR THE FFT
#include <complex.h>
#include <math.h>
// for some reason we're not getting this from math.h...
#ifndef M_PIl
# define M_PIl          3.141592653589793238462643383279502884L /* pi */
#endif

#include "fft.h"
#include "io.h"		  // print_naive
#include "larc.h"         // mat_level_t
#include "matrix_store.h" // matns_ptr_t, mats_ptr_t
#include "matmath.h"      // kronecker_product, matrix_mult, etc.

/*!
 * \file fft.c
 * \brief This file contains routines useful for Fast Fourier Transform (FFT)
 * calculations
 */


/* This program returns the packedID of the k-th power of the n-th principal
 * root of unity */
/* This program is often used to create a list of all the roots of unity for
 * fixed n  */
/* see the test code in tests/python/test_roots_unity.py */
int64_t k_th_power_of_n_th_root_of_unity_pID(int k, int n, int verbose){

// must use complex or mpcomplex or mpratcomplex scalarType
#ifndef IS_COMPLEX
  fprintf(stderr,"\nERROR in %s: not using COMPLEX,\n",__func__);
  fprintf(stderr,"MPCOMPLEX, MPRATCOMPLEX, or a complex Clifford type!\n");
  exit(1);
#endif
  
  if (verbose) {
    printf("The routine %s\n",__func__);
    printf("will calculate and print the\n");
    printf("%d-th power of the principal %d-th root of unity.\n",k,n);
    printf("Then after loading this into the matrix store, we\n");
    printf("print the associated value in the matrix store.\n");
  }

  scalarType rootkn;
  sca_init(&rootkn);

#if defined(USE_MPCOMPLEX)
  mpc_rootofunity(rootkn, n, k, MPC_RNDNN);
#else
  mpc_t mpc_rootkn;
  mpc_init2(mpc_rootkn, mpreal_precision);
  mpc_rootofunity(mpc_rootkn, n, k, MPC_RNDNN);
  convert_from_two_mpfr_to_complex_scalarType(&rootkn,
	mpc_realref(mpc_rootkn), mpc_imagref(mpc_rootkn));
  mpc_clear(mpc_rootkn);
#endif  

  mats_ptr_t rootkn_ptr = get_scalarPTR_for_scalarVal(rootkn);
  int64_t rootkn_pID = rootkn_ptr->packedID;

  if (verbose) {
    printf("The calculation of the root gives\n");
    char *s = sca_get_readable_approx_str(rootkn);
    printf("Scalar: %s\n", s);
    printf("\nand its matrixID in the matrix store is %"
	PRId64 "\n\n",MID_FROM_PID(rootkn_pID));
    free(s);
  }

  sca_clear(&rootkn);

  return rootkn_pID;
}

/* This takes the size of the matrix and returns a packedID
   corresponding to the inverse permutation matrix */
int64_t create_invShufMat(mat_level_t m_lev)
{
  int64_t PI_pID;

  // EXCEPTION CHECKING
  if (m_lev < 0) {
    fprintf(stderr,"ERROR in %s: the matrix level=%d is too small\n", __func__,
	m_lev);
    exit(1);
  }
  if (m_lev > max_level_allowed_matrixStore()) {
    fprintf(stderr,"ERROR in %s: the matrix level=%d > max_level\n", __func__,
	m_lev);
    exit(1);
  }

  // SMALL CASES
  if (m_lev == 0) PI_pID = get_identity_pID(0);
  else if (m_lev == 1) PI_pID = get_identity_pID(1);
  // RECURSIVE CASE
  else {
    int64_t PIsmall_pID = create_invShufMat(m_lev-1);
    matns_ptr_t PIsmall_ptr = (matns_ptr_t)get_recordPTR_from_pID(PIsmall_pID,
	"",__func__,0);
    int64_t *small_panel = PIsmall_ptr->subMatList;
    int64_t small_zeroID = get_zero_pID(m_lev-2,m_lev-2);

    int64_t panel[4];
    panel[0] = get_pID_from_four_sub_pIDs(small_panel[0],small_panel[1],
	small_zeroID, small_zeroID, m_lev-1, m_lev-1);
    panel[1] = get_pID_from_four_sub_pIDs(small_zeroID, small_zeroID,
	small_panel[0],small_panel[1], m_lev-1, m_lev-1);
    panel[2] = get_pID_from_four_sub_pIDs(small_panel[2],small_panel[3],
	small_zeroID, small_zeroID, m_lev-1, m_lev-1);
    panel[3] = get_pID_from_four_sub_pIDs(small_zeroID, small_zeroID,
	small_panel[2],small_panel[3], m_lev-1, m_lev-1);

    PI_pID = get_pID_from_array_of_four_sub_pIDs(panel, m_lev, m_lev);
  }
  return (PI_pID);
}


/* Prints the (2^pow)-th roots of unity */
int print_pow2_roots_unity(uint32_t pow){

// must use complex or mpcomplex or mpratcomplex scalarType
#ifndef IS_COMPLEX
  fprintf(stderr,"\nERROR in %s: not using COMPLEX,\n",__func__);
  fprintf(stderr,"MPCOMPLEX, MPRATCOMPLEX, or a complex Clifford type!\n");
  exit(1);
#endif

  if (pow>10)
  {
    fprintf(stderr,"\nERROR in %s, attempting to print out\n",__func__);
    fprintf(stderr,"information about all powers of the (2**%d)th\n",pow);
    fprintf(stderr,"root of unity. This is not a reasonable request.\n");
    return -1;
  }

  /* double complex root; */
  /* long double complex rootlong; */
 
  /* printf("The size of double complex is %zd\n",sizeof(root)); */
  /* printf("The size of long double complex is %zd\n\n",sizeof(rootlong)); */

  // create n = 2^pow
  uint64_t n = (uint64_t)1 << pow;

  printf("The (2^%d)-th roots of unity are:\n",pow);

  for (int k = 0; k< n ; k++) {
    // the following 4 lines are just for comparison purposes
    /* root = cexp(2 * M_PIl * I * k / n); */
    /* printf("k=%d, %lg + %lg i\n",k, creal(root), cimag(root)); */
    /* rootlong = cexpl(2 * M_PIl * I * k / n); */
    /* printf(" LONG:  k=%d, %Lg + %Lg i\n\n",k, creall(rootlong), cimagl(rootlong)); */

    // this calculates the value to full precision
    int64_t root_pID = k_th_power_of_n_th_root_of_unity_pID(k,n,0);

    printf("\tHave loaded the root into the store, it is now:\n\t");
    print_naive(root_pID);
    printf("\n");
  }
  return 0;
}


/* Returns the principal (2^pow)-th root of unity, *
 * e^{i 2 pi/2^pow} and either loads or finds it in the matrix_store */
int64_t principal_pow2_root_unity_pID(uint32_t pow)
{
// must use complex or mpcomplex or mpratcomplex scalarType
#ifndef IS_COMPLEX
  fprintf(stderr,"\nERROR in %s: not using COMPLEX,\n",__func__);
  fprintf(stderr,"MPCOMPLEX, MPRATCOMPLEX, or a complex Clifford type!\n");
  exit(1);
#endif

  // create n = 2^pow
  uint64_t n = (uint64_t)1 << pow;

  printf("The (2^%d)-th principal root of unity is:\n",pow);  

  int64_t root_pID = k_th_power_of_n_th_root_of_unity_pID(1,n,0);

  print_naive(root_pID);

  return root_pID;
  
}


/* Python Interface for the create_fft_D(mat_level_t m_lev) function.
   This takes the size of the matrix and returns a matrix ID
   corresponding the D matrix needed for the fft recursion */
int64_t create_FFT_DMat(mat_level_t m_lev)
{

// must use complex or mpcomplex or mpratcomplex scalarType
#ifndef IS_COMPLEX
  fprintf(stderr,"\nERROR in %s: not using COMPLEX,\n",__func__);
  fprintf(stderr,"MPCOMPLEX, MPRATCOMPLEX, or a complex Clifford type!\n");
  exit(1);
#endif

  int64_t D_pID;

  // EXCEPTION CHECKING
  if (m_lev < 0) {
    fprintf(stderr,"ERROR in %s: the matrix level=%d is too small\n", __func__,
	m_lev);
    exit(1);
  }
  if (m_lev > max_level_allowed_matrixStore()) {
    fprintf(stderr,"ERROR in %s: the matrix level=%d > max_level\n", __func__,
	m_lev);
    exit(1);
  }

  // SMALL CASES
  if (m_lev == 0) {
    D_pID = get_identity_pID(0);
  }
  // RECURSIVE CASE
  else {
    // construct the matrix [1,0;0,(2^(m_lev+1))-th root of unity]
    uint64_t npow2 = (uint64_t)1 << (m_lev+1);
    int64_t panel[4];
    panel[0] = get_identity_pID(0);
    panel[1] = get_zero_pID(0,0);
    panel[2] = get_zero_pID(0,0);
    panel[3] = k_th_power_of_n_th_root_of_unity_pID(1,npow2,0);
    int64_t A_pID = get_pID_from_array_of_four_sub_pIDs(panel, 1, 1);
    int64_t Dminus1_pID = create_FFT_DMat(m_lev-1); 
    D_pID = kronecker_product(Dminus1_pID,A_pID);
  }
  return (D_pID);
}

/* Python Interface for the create_fft_C(mat_level_t m_lev) function.
   This takes the size of the matrix and returns a matrix ID
   corresponding the C matrix needed for the fft recursion */
int64_t create_FFT_CMat(mat_level_t m_lev)
{

// must use complex or mpcomplex or mpratcomplex scalarType
#ifndef IS_COMPLEX
  fprintf(stderr,"\nERROR in %s: not using COMPLEX,\n",__func__);
  fprintf(stderr,"MPCOMPLEX, MPRATCOMPLEX, or a complex Clifford type!\n");
  exit(1);
#endif

  int64_t C_pID;

  // EXCEPTION CHECKING
  if (m_lev < 1) {
    fprintf(stderr,"ERROR in %s: the matrix level=%d is too small\n", __func__,
	m_lev);
    exit(1);
  }
  if (m_lev > max_level_allowed_matrixStore()) {
    fprintf(stderr,"ERROR in %s: the matrix level=%d > max_level\n", __func__,
	m_lev);
    exit(1);
  }

  // SMALL CASE
  if (m_lev == 1) C_pID = get_iHadamard_pID(1);
  // ELSE build C from panels [I(m_lev-1),D(m_lev-1);I(m_lev-1),-D(m_lev-1)]
  else {
    int64_t Dminus1_pID = create_FFT_DMat(m_lev-1); 
    scalarType scalarM1; //'minus 1'
    sca_init(&scalarM1);
    sca_set_str(&scalarM1, "-1");
    mats_ptr_t ptr_scalarM1 = get_scalarPTR_for_scalarVal(scalarM1);
    sca_clear(&scalarM1);
    int64_t matID_scalarM1 = ptr_scalarM1->packedID;
    int64_t minusDminus1_pID = scalar_mult(matID_scalarM1,Dminus1_pID);
    int64_t panel[4];
    panel[0] = get_identity_pID(m_lev-1);
    panel[1] = Dminus1_pID;
    panel[2] = get_identity_pID(m_lev-1);
    panel[3] = minusDminus1_pID;
    C_pID = get_pID_from_array_of_four_sub_pIDs(panel, m_lev, m_lev);
  }
  return (C_pID);
}


/* Python Interface for the create_fft_matrix(mat_level_t m_lev) function.
   This takes the size of the matrix and returns a matrix ID
   corresponding to the fft matrix */
int64_t 
create_FFTMat(mat_level_t m_lev) {

// must use complex or mpcomplex or mpratcomplex scalarType
#ifndef IS_COMPLEX
  fprintf(stderr,"\nERROR in %s: not using COMPLEX,\n",__func__);
  fprintf(stderr,"MPCOMPLEX, MPRATCOMPLEX, or a complex Clifford type!\n");
  exit(1);
#endif

  int64_t F_pID;

  // EXCEPTION CHECKING
  if (m_lev < 0) {
    fprintf(stderr,"ERROR in %s: the matrix level=%d is too small\n", __func__,
	m_lev);
    exit(1);
  }
  if (m_lev > max_level_allowed_matrixStore()) {
    fprintf(stderr,"ERROR in %s: the matrix level=%d > max_level\n", __func__,
	m_lev);
    exit(1);
  }

  // SMALL CASES
  if (m_lev == 0) F_pID = get_identity_pID(0);
  // RECURSIVE CASE   F(n) = C(n) * [I(1) @ F(n-1)] * PI(n)
  else {
    int64_t Cn_pID = create_FFT_CMat(m_lev);
    //printf("C matrix:\n");
    //print_naive(Cn_pID);
    int64_t PIn_pID = create_invShufMat(m_lev);
    //printf("\nPI matrix:\n");
    //print_naive(PIn_pID);
    int64_t Fn_minus1_pID = create_FFTMat(m_lev-1);
    //printf("\nsmaller F matrix:\n");
    //print_naive(Fn_minus1_pID);
    //printf("\n");
    int64_t blockF_pID = kronecker_product(
                       get_identity_pID(1),Fn_minus1_pID);
    F_pID = matrix_mult(Cn_pID, matrix_mult(blockF_pID, PIn_pID));
  }
  return (F_pID);
}
