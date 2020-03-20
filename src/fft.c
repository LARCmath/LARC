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
   /*   F_k = C_k * (I_2 @ F_(k-1) ) * PI_k */
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

#include "fft.h"
#include "io.h"           // print_naive_by_matPTR
#include "larc.h"         // mat_level_t
#include "matrix_store.h" // mat_ptr_t
#include "matmath.h"      // kronecker_product, matrix_mult, etc.

/* Python Interface for the create_perm_inv(mat_level_t n) function.
   This takes the size of the matrix and returns a matrix ID
   corresponding to the inverse permutation matrix */
int64_t 
create_perm_inv_matrixID(mat_level_t n) {

  // get the matrix pointers from the matrixID, and see if still in store 
  mat_ptr_t PI_ptr = create_perm_inv(n);

  // call matrix pointer version of function
  int64_t PI_mID = get_matID_from_matPTR(PI_ptr);
  return PI_mID;
}


/* The create_perm_inv(mat_level_t n) function
   takes the size of the matrix and returns a matrix PTR
   corresponding the inverse permution matrix */
mat_ptr_t
create_perm_inv(mat_level_t n) 
{

  mat_ptr_t PI_ptr;

  // EXCEPTION CHECKING
  if (n < 0) {
    fprintf(stderr,"ERROR in %s: the matrix level n=%d is too small\n", __func__,n);
    exit(1);
  }
  if (n > maximum_level()) {
    fprintf(stderr,"ERROR in %s: the matrix level n=%d > max_level\n", __func__,n);
    exit(1);
  }

  // SMALL CASES
  if (n == 0) {
    PI_ptr = get_identity_matrix_ptr(0);
    return (PI_ptr);
  }

  if (n == 1) {
    PI_ptr = get_identity_matrix_ptr(1);
    return (PI_ptr);
  }

  // RECURSIVE CASE
  else {
    mat_ptr_t PIsmall_ptr = create_perm_inv(n-1);
    mat_ptr_t small_panel[4];
    for (int i=0;i<4;++i) {
      small_panel[i] = matrix_sub(PIsmall_ptr, i);
    }
    mat_ptr_t panel[4];
    mat_ptr_t temp_array[4];
    temp_array[0] = small_panel[0];
    temp_array[1] = small_panel[1];
    temp_array[2] = get_zero_matrix_ptr(n-2,n-2);
    temp_array[3] = get_zero_matrix_ptr(n-2,n-2);
    panel[0] = get_matPTR_from_array_of_four_subMatPTRs(temp_array,n-1,n-1);

    temp_array[0] = get_zero_matrix_ptr(n-2,n-2);
    temp_array[1] = get_zero_matrix_ptr(n-2,n-2);
    temp_array[2] = small_panel[0];
    temp_array[3] = small_panel[1];
    panel[1] = get_matPTR_from_array_of_four_subMatPTRs(temp_array,n-1,n-1);

    temp_array[0] = small_panel[2];
    temp_array[1] = small_panel[3];
    temp_array[2] = get_zero_matrix_ptr(n-2,n-2);
    temp_array[3] = get_zero_matrix_ptr(n-2,n-2);
    panel[2] = get_matPTR_from_array_of_four_subMatPTRs(temp_array,n-1,n-1);

    temp_array[0] = get_zero_matrix_ptr(n-2,n-2);
    temp_array[1] = get_zero_matrix_ptr(n-2,n-2);
    temp_array[2] = small_panel[2];
    temp_array[3] = small_panel[3];
    panel[3] = get_matPTR_from_array_of_four_subMatPTRs(temp_array,n-1,n-1);

    PI_ptr = get_matPTR_from_array_of_four_subMatPTRs(panel,n,n);
    return (PI_ptr);
  }
}


/* prints the (2^n)-th roots of unity */
int print_pow2_roots_unity(mat_level_t n){

// must use complex or mpcomplex or mpratcomplex scalarType
#if !defined(USE_COMPLEX) && !defined(USE_MPCOMPLEX) && !defined(USE_MPRATCOMPLEX)
  fprintf(stderr,"\nERROR in %s not using COMPLEX or MPCOMPLEX or MPRATCOMPLEX type!\n", __func__);
  exit(1);
#endif

  double complex root;
  long double complex rootlong;
 
  printf("The size of double complex is %zd\n",sizeof(root));
  printf("The size of long double complex is %zd\n\n",sizeof(rootlong));

  // create 2^n
  int64_t pow = 1L<<n;

  printf("The (2^%d)-th roots of unity are:\n",n);

  for (int k = 0; k< pow ; k++) {
    root = cexp(2 * M_PI * I * k / pow);
    printf("k=%d, %lg + %lg i\n",k, creal(root), cimag(root));
    rootlong = cexpl(2 * M_PI * I * k / pow);
    printf(" LONG:  k=%d, %Lg + %Lg i\n\n",k, creall(rootlong), cimagl(rootlong));
    // load the values into the store and then print them out again

    scalarType scalarRoot;
    sca_init(&scalarRoot);
    sca_set_2ldoubles(&scalarRoot,creall(rootlong),cimagl(rootlong));
    mat_ptr_t root_ptr = get_valMatPTR_from_val(scalarRoot);
    // int64_t root_mID = get_matID_from_matPTR(root_ptr);
    sca_clear(&scalarRoot);

    printf("\tHave loaded the root into the store, it is now:\n\t");
    print_naive_by_matPTR(root_ptr);
    printf("\n");
  }
  return 0;
  
}


/* return the principal (2^n)-th root of unity */
int64_t principal_pow2_root_unity(mat_level_t n){

// must use complex or mpcomplex or mpratcomplex scalarType
#if !defined(USE_COMPLEX) && !defined(USE_MPCOMPLEX) && !defined(USE_MPRATCOMPLEX)
  fprintf(stderr,"\nERROR in %s not using COMPLEX or MPCOMPLEX or MPRATCOMPLEX type!\n", __func__);
  exit(1);
#endif

  double complex root;
  long double complex rootlong;

  printf("The size of double complex is %zd\n",sizeof(root));
  printf("The size of long double complex is %zd\n\n",sizeof(rootlong));

  // create 2^n
  int64_t pow = 1L<<n;

  printf("The (2^%d)-th principal root of unity is:\n",n);  

  root = cexp(2 * M_PI * I  / pow);
  printf("%lg + %lg i\n", creal(root), cimag(root));
  rootlong = cexpl(2 * M_PI * I / pow);
  printf(" LONG:  %Lg + %Lg i\n\n",creall(rootlong), cimagl(rootlong));
  // load the values into the store and then print them out again
  scalarType scalarRoot;
  sca_init(&scalarRoot);
  sca_set_2ldoubles(&scalarRoot,creall(rootlong),cimagl(rootlong));
  mat_ptr_t root_ptr = get_valMatPTR_from_val(scalarRoot);
  int64_t root_mID = get_matID_from_matPTR(root_ptr);
  sca_clear(&scalarRoot);

  return root_mID;
  
}




/* return the matrix pointer for the first (2^n)-th root of unity */
mat_ptr_t
get_first_pow2_root_unity(mat_level_t n){

// must use complex or mpcomplex or mpratcomplex scalarType
#if !defined(USE_COMPLEX) && !defined(USE_MPCOMPLEX) && !defined(USE_MPRATCOMPLEX)
  fprintf(stderr,"\nERROR in %s not using COMPLEX or MPCOMPLEX or MPRATCOMPLEX type!\n", __func__);
  exit(1);
#endif

  // create 2^n
  int64_t pow = 1L<<n;

  //double complex root;
  //root = cexp(2 * M_PI * I  / pow);

  long double complex rootlong;  
  rootlong = cexpl(2 * M_PI * I  / pow);
  
  scalarType scalarRoot;
  sca_init(&scalarRoot);
  sca_set_2ldoubles(&scalarRoot,creall(rootlong),cimagl(rootlong));
  
  mat_ptr_t root_ptr = get_valMatPTR_from_val(scalarRoot);
  sca_clear(&scalarRoot);

  return root_ptr;
  
}



/* Python Interface for the create_fft_D(mat_level_t n) function.
   This takes the size of the matrix and returns a matrix ID
   corresponding the D matrix needed for the fft recursion */
int64_t 
create_fft_D_matrixID(mat_level_t n) {

// must use complex or mpcomplex or mpratcomplex scalarType
#if !defined(USE_COMPLEX) && !defined(USE_MPCOMPLEX) && !defined(USE_MPRATCOMPLEX)
  fprintf(stderr,"\nERROR in %s not using COMPLEX or MPCOMPLEX or MPRATCOMPLEX type!\n", __func__);
  exit(1);
#endif

  // get the matrix pointers from the matrixID, and see if still in store 
  mat_ptr_t D_ptr = create_fft_D(n);

  // call matrix pointer version of function
  int64_t D_mID = get_matID_from_matPTR(D_ptr);
  return D_mID;
}


/* The create_fft_D(mat_level_t n) function
   takes the size of the matrix and returns a matrix PTR
   corresponding the D matrix for the fft recursion*/
mat_ptr_t create_fft_D(mat_level_t n) 
{

// must use complex or mpcomplex or mpratcomplex scalarType
#if !defined(USE_COMPLEX) && !defined(USE_MPCOMPLEX) && !defined(USE_MPRATCOMPLEX)
  fprintf(stderr,"\nERROR in %s not using COMPLEX or MPCOMPLEX or MPRATCOMPLEX type!\n", __func__);
  exit(1);
#endif

  mat_ptr_t D_ptr;

  // EXCEPTION CHECKING
  if (n < 0) {
    fprintf(stderr,"ERROR in %s: the matrix level n=%d is too small\n", __func__,n);
    exit(1);
  }
  if (n > maximum_level()) {
    fprintf(stderr,"ERROR in %s: the matrix size n=%d > max_level\n", __func__,n);
    exit(1);
  }

  // SMALL CASES
  if (n == 0) {
    D_ptr = get_identity_matrix_ptr(0);
    return (D_ptr);
  }

  // RECURSIVE CASE
  else {
    // construct the matrix [1,0;0,(2^(n+1))-th root of unity]
    mat_ptr_t panel[4];
    panel[0] = get_identity_matrix_ptr(0);
    panel[1] = get_zero_matrix_ptr(0,0);
    panel[2] = get_zero_matrix_ptr(0,0);
    panel[3] = get_first_pow2_root_unity(n+1);
    mat_ptr_t A_ptr = get_matPTR_from_array_of_four_subMatPTRs(panel,1,1);
    mat_ptr_t Dminus1_ptr = create_fft_D(n-1); 
    // D_ptr = kronecker_product(A_ptr, Dminus1_ptr);
    D_ptr = kronecker_product(Dminus1_ptr,A_ptr);

    return (D_ptr);
  }
}



/* Python Interface for the create_fft_C(mat_level_t n) function.
   This takes the size of the matrix and returns a matrix ID
   corresponding the C matrix needed for the fft recursion */
int64_t create_fft_C_matrixID(mat_level_t n) {

// must use complex or mpcomplex or mpratcomplex scalarType
#if !defined(USE_COMPLEX) && !defined(USE_MPCOMPLEX) && !defined(USE_MPRATCOMPLEX)
  fprintf(stderr,"\nERROR in %s not using COMPLEX or MPCOMPLEX or MPRATCOMPLEX type!\n", __func__);
  exit(1);
#endif

  // get the matrix pointers from the matrixID, and see if still in store 
  mat_ptr_t C_ptr = create_fft_C(n);

  // call matrix pointer version of function
  int64_t C_mID = get_matID_from_matPTR(C_ptr);
  return C_mID;
}


/* The create_fft_C(mat_level_t n) function
   takes the size of the matrix and returns a matrix PTR
   corresponding the C matrix for the fft recursion*/
mat_ptr_t create_fft_C(mat_level_t n) 
{

// must use complex or mpcomplex or mpratcomplex scalarType
#if !defined(USE_COMPLEX) && !defined(USE_MPCOMPLEX) && !defined(USE_MPRATCOMPLEX)
  fprintf(stderr,"\nERROR in %s not using COMPLEX or MPCOMPLEX or MPRATCOMPLEX type!\n", __func__);
  exit(1);
#endif

  mat_ptr_t C_ptr;

  // EXCEPTION CHECKING
  if (n < 1) {
    fprintf(stderr,"ERROR in %s: the matrix level n=%d is too small\n", __func__,n);
    exit(1);
  }
  if (n > maximum_level()) {
    fprintf(stderr,"ERROR in %s: the matrix size n=%d > max_level\n", __func__,n);
    exit(1);
  }

  // SMALL CASE
  if (n == 1) {
    C_ptr = get_iHadamard_matrix_ptr(1);
    return (C_ptr);
  }

  // ELSE build C from panels [I(n-1),D(n-1);I(n-1),-D(n-1)]
  else {
    mat_ptr_t Dminus1_ptr = create_fft_D(n-1); 
    scalarType scalarM1; //'minus 1'
    sca_init(&scalarM1);
    sca_set_str(&scalarM1, "-1");
    mat_ptr_t matptr_scalarM1 = get_valMatPTR_from_val(scalarM1);
    mat_ptr_t minusDminus1_ptr = scalar_mult(matptr_scalarM1,Dminus1_ptr); 
    mat_ptr_t panel[4];
    panel[0] = get_identity_matrix_ptr(n-1);
    panel[1] = Dminus1_ptr;
    panel[2] = get_identity_matrix_ptr(n-1);
    panel[3] = minusDminus1_ptr;
    mat_ptr_t C_ptr = get_matPTR_from_array_of_four_subMatPTRs(panel,n,n);

    return (C_ptr);
  }
}


/* Python Interface for the create_fft_matrix(mat_level_t n) function.
   This takes the size of the matrix and returns a matrix ID
   corresponding to the fft matrix */
int64_t 
create_fft_matrix_matrixID(mat_level_t n) {

// must use complex or mpcomplex or mpratcomplex scalarType
#if !defined(USE_COMPLEX) && !defined(USE_MPCOMPLEX) && !defined(USE_MPRATCOMPLEX)
  fprintf(stderr,"\nERROR in %s not using COMPLEX or MPCOMPLEX or MPRATCOMPLEX type!\n", __func__);
  exit(1);
#endif

  // get the matrix pointers from the matrixID, and see if still in store 
  mat_ptr_t F_ptr = create_fft_matrix(n);

  // call matrix pointer version of function
  int64_t F_mID = get_matID_from_matPTR(F_ptr);
  return F_mID;
}


/* The create_fft_matrix(mat_level_t n) function
   takes the size of the matrix and returns a matrix PTR
   corresponding to the fft matrix */
mat_ptr_t create_fft_matrix(mat_level_t n) 
{

// must use complex or mpcomplex or mpratcomplex scalarType
#if !defined(USE_COMPLEX) && !defined(USE_MPCOMPLEX) && !defined(USE_MPRATCOMPLEX)
  fprintf(stderr,"\nERROR in %s not using COMPLEX or MPCOMPLEX or MPRATCOMPLEX type!\n", __func__);
  exit(1);
#endif

  mat_ptr_t F_ptr;

  // EXCEPTION CHECKING
  if (n < 0) {
    fprintf(stderr,"ERROR in %s: the matrix level n=%d is too small\n", __func__,n);
    exit(1);
  }
  if (n > maximum_level()) {
    fprintf(stderr,"ERROR in %s: the matrix size n=%d > max_level\n", __func__,n);
    exit(1);
  }

  // SMALL CASES
  if (n == 0) {
    F_ptr = get_identity_matrix_ptr(0);
    return (F_ptr);
  }

  // RECURSIVE CASE   F(n) = C(n) * [I(1) @ F(n-1)] * PI(n)
  else {
    mat_ptr_t Cn_ptr = create_fft_C(n);
    mat_ptr_t PIn_ptr = create_perm_inv(n);
    mat_ptr_t Fn_minus1_ptr = create_fft_matrix(n-1);
    mat_ptr_t blockF_ptr = kronecker_product(
                       get_identity_matrix_ptr(1),Fn_minus1_ptr);
    F_ptr = matrix_mult(Cn_ptr,matrix_mult(blockF_ptr,PIn_ptr));

    return (F_ptr);
  }
}

////////////////////
// Added Oct 2019 //
////////////////////

/* prints the n-th roots of unity and complex and 
   long double complex calculations */
int print_n_th_roots_of_unity(int n, int verbose){

// must use complex or mpcomplex or mpratcomplex scalarType
#if !defined(USE_COMPLEX) && !defined(USE_MPCOMPLEX) && !defined(USE_MPRATCOMPLEX)
  fprintf(stderr,"\nERROR in %s not using COMPLEX or MPCOMPLEX or MPRATCOMPLEX type!\n", __func__);
  exit(1);
#endif
  
  // compare complex vs long double complex;
  scalarType root1;
  sca_init(&root1);
  scalarType root2;
  sca_init(&root2);

  // these should be double
  complex root_double;
  // long double complex rootlong
  long double complex root_long;
 
  if (verbose) {
    printf("The routine %s will calculate and print a root of unity\n",__func__);
    printf("Then after loading it into the matrix store, print again from the store\n");
    printf("\tThe effectiveness of this may be impacted by the size of scalarType %zd\n",sizeof(root1));
    printf("\tthe size of char is %zd bytes, and\n",sizeof(char));
    printf("\tthe size of complex is %zd bytes, and\n",sizeof(root_double));
    printf("\tthe size of long double complex is %zd bytes.\n\n",sizeof(root_long));
  }

#ifdef DEBUG  
  printf("\n\nPRINTING TEST OF WHAT SHOULD ALREADY BE PRELOADED\n");
  printf("Loading 1,-1, i, -i into store, just in case something was wrong\n");
  scalarType temp;
  sca_init(&temp);
  // sca_set(&temp,1);
  sca_set_2ldoubles(&temp,1.0, 0.0);
  mat_ptr_t temp_ptr = get_valMatPTR_from_val(temp);
  int64_t temp_mID = get_matID_from_matPTR(temp_ptr);
  printf("matrixID for 1 is %zd\n",temp_mID);
  sca_set_2ldoubles(&temp,-1.0, 0.0);
  temp_ptr = get_valMatPTR_from_val(temp);
  temp_mID = get_matID_from_matPTR(temp_ptr);
  printf("matrixID for -1 is %zd\n",temp_mID);
  sca_set_2ldoubles(&temp,0.0, 1.0);
  temp_ptr = get_valMatPTR_from_val(temp);
  temp_mID = get_matID_from_matPTR(temp_ptr);
  printf("matrixID for I is %zd\n",temp_mID);
  sca_set_2ldoubles(&temp,0.0, -1.0);
  temp_ptr = get_valMatPTR_from_val(temp);
  temp_mID = get_matID_from_matPTR(temp_ptr);
  printf("matrixID for -I is %zd\n",temp_mID);
#endif

  printf("\nThe %d-th roots of unity are:\n",n);

  for (int k = 0; k< n ; k++) {
    root_double = cexp(2 * M_PI * I * k / n);
    printf("complex: k=%d, %lg + %lg i\n",k, creal(root_double), cimag(root_double));
    printf("complex hex: k=%d, %la + %la i\n",k, creal(root_double), cimag(root_double));
    root_long = cexpl(2 * M_PI * I * k / n);
    printf("long double complex:  k=%d, %Lg + %Lg i\n\n",k, creall(root_long), cimagl(root_long));
    printf("long double complex hex:  k=%d, %La + %La i\n\n",k, creall(root_long), cimagl(root_long));
    
    // load the values into the store and then print them by matrixPTR
    sca_set_2ldoubles(&root1,creal(root_double), cimag(root_double));
    mat_ptr_t root1_ptr = get_valMatPTR_from_val(root1);
    printf("\tThe root is loaded into the store. Printing value root from store gives:\n\t");
    print_naive_by_matPTR(root1_ptr);
    printf("\n");

    sca_set_2ldoubles(&root2,creall(root_long), cimagl(root_long));
    mat_ptr_t root2_ptr = get_valMatPTR_from_val(root2);
    printf("\tThe root is loaded into the store. Printing value root from store gives:\n\t");
    print_naive_by_matPTR(root2_ptr);
    if (root1_ptr != root2_ptr) {
      printf("the complex calculation gives a different scalar than the long double complex calculation\n");
    }
    else {
      printf("the complex calculation gives the same scalar as the long double complex calculation\n");
    }
    printf("\n");

  }

  return 0;

}



/* return the matrix ID for the principal n-th root of unity */
int64_t principal_n_th_root_of_unity_matID(int n, int verbose){

// must use complex or mpcomplex or mpratcomplex scalarType
#if !defined(USE_COMPLEX) && !defined(USE_MPCOMPLEX) && !defined(USE_MPRATCOMPLEX)
  fprintf(stderr,"\nERROR in %s not using COMPLEX or MPCOMPLEX or MPRATCOMPLEX type!\n", __func__);
  exit(1);
#endif

  // compare complex vs long double complex;
  scalarType root1;
  sca_init(&root1);
  scalarType root2;
  sca_init(&root2);

  // these should be double
  complex root_double;
  // long double complex rootlong
  long double complex root_long;
 
  if (verbose) {
    printf("The routine %s will calculate and print a root of unity\n",__func__);
    printf("Then after loading it into the matrix store, print again from the store\n");
    printf("\tThe effectiveness of this may be impacted by the size of scalarType %zd\n",sizeof(root1));
    printf("\tthe size of complex %zd, and\n",sizeof(root_double));
    printf("\tthe size of long double complex %zd,\n\n",sizeof(root_long));
  }
  
  //for (int k = 0; k< n ; k++) {
  root_double = cexp(2 * M_PI * I / n);
  root_long = cexpl(2 * M_PI * I / n);

  if (verbose) {
    printf("\nThe principal %d-th root of unity is:\n",n);
    printf("complex: %lg + %lg i\n", creal(root_double), cimag(root_double));
    printf("long double complex:  %Lg + %Lg i\n\n", creall(root_long), cimagl(root_long));
  }

  // load the values into the store get the matrixID and return that to the user
  sca_set_2ldoubles(&root1,creal(root_double), cimag(root_double));
  mat_ptr_t root1_ptr = get_valMatPTR_from_val(root1);
  int64_t root1_mID = get_matID_from_matPTR(root1_ptr);

  if (verbose) {
    printf("\nand its matrixID in the matrix store is %zd\n\n",root1_mID);
  }
  
  // check this against what happens if we use long double complex
  sca_set_2ldoubles(&root2,creall(root_long), cimagl(root_long));
  mat_ptr_t root2_ptr = get_valMatPTR_from_val(root2);
  int64_t root2_mID = get_matID_from_matPTR(root2_ptr);
  if (root1_mID != root2_mID) {
    fprintf(stderr,"WARNING: in %s the complex calculation gives a different scalar\n",__func__);
    fprintf(stderr,"\tthan the long double complex calculation!\n");
  }
  printf("\n");
  
  return root1_mID;
  
}

/* This program returns the matrixID of the k-th power of the n-th principal */
/* root of unity */
/* This program is often used to create a list of all the roots of unity for fixed n  */
/* see the test code in tests/python/test_roots_unity.py */
int64_t k_th_power_of_n_th_root_of_unity_matID(int k, int n, int verbose){

// must use complex or mpcomplex or mpratcomplex scalarType
#if !defined(USE_COMPLEX) && !defined(USE_MPCOMPLEX) && !defined(USE_MPRATCOMPLEX)
  fprintf(stderr,"\nERROR in %s not using COMPLEX or MPCOMPLEX or MPRATCOMPLEX type!\n", __func__);
  exit(1);
#endif

  // compare complex vs long double complex;
  scalarType root1;
  sca_init(&root1);
  scalarType root2;
  sca_init(&root2);

  // these should be double
  complex root_double;
  // long double complex rootlong
  long double complex root_long;

  // printf("In %s verbose is %d\n",__func__,verbose);
 
  if (verbose) {
    printf("The routine %s\n",__func__);
    printf("will calculate and print the\n");
    printf("%d-th power of the principal %d-th root of unity.\n",k,n);
    printf("Then after loading this into the matrix store, we\n");
    printf("print the associated value in the matrix store.\n");
    printf("\tThe closeness of the calculated and stored values may be\n");
    printf("\timpacted by the size of scalarType %zd\n",sizeof(root1));
    printf("\tthe size of complex %zd, and\n",sizeof(root_double));
    printf("\tthe size of long double complex %zd,\n\n",sizeof(root_long));
  }

  root_double = cexp(2 * M_PI * I * k / n);
  root_long = cexpl(2 * M_PI * I * k / n);
  if (verbose) {
    printf("The calculation of the root gives\n");
    printf("complex: k=%d, %lg + %lg i\n",k, creal(root_double), cimag(root_double));
    printf("complex hex: k=%d, %la + %la i\n",k, creal(root_double), cimag(root_double));
    printf("long double complex:  k=%d, %Lg + %Lg i\n\n",k, creall(root_long), cimagl(root_long));
    printf("long double complex hex:  k=%d, %La + %La i\n\n",k, creall(root_long), cimagl(root_long));
  }

  // load the values into the store get the matrixID and return that to the user
  sca_set_2ldoubles(&root1,creal(root_double), cimag(root_double));
  mat_ptr_t root1_ptr = get_valMatPTR_from_val(root1);
  int64_t root1_mID = get_matID_from_matPTR(root1_ptr);

  if (verbose) {
    printf("\nand its matrixID in the matrix store is %zd\n\n",root1_mID);
  }
  
  // check this against what happens if we use long double complex
  sca_set_2ldoubles(&root2,creall(root_long), cimagl(root_long));
  mat_ptr_t root2_ptr = get_valMatPTR_from_val(root2);
  int64_t root2_mID = get_matID_from_matPTR(root2_ptr);
  if (root1_mID != root2_mID) {
    fprintf(stderr,"WARNING: in %s the complex calculation gives a different scalar\n",__func__);
    fprintf(stderr,"\tthan the long double complex calculation!\n\n");
  }
  
  return root1_mID;
  
}



