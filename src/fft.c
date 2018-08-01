//                          fft.c 
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

// From Van Loan:
   /*   F_k = C_k * (I_2 @ F_(k-1) ) * PI_k */
   /* where  */
   /*   PI_k is the 2^k by 2^k inverse shuffle matrix  */
   /*       (its transverse is its inverse and would shuffle a column vector) */
   /*       PI_2 = (1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1) is CNOT */
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
#ifdef __cplusplus
#include <complex>
#else
#include <complex.h>
#endif
#include <math.h>

#include "fft.h"
#include "io.h" // print_matrix_naive

/* Python Interface for the create_perm_inv(mat_level_t n) function.
   This takes the size of the matrix and returns a matrix ID
   corresponding to the inverse permutation matrix */
int64_t 
create_perm_inv_matrixID(mat_level_t n) {

#ifndef USE_COMPLEX
  printf("\nERROR in %s not using COMPLEX type!\n", __func__);
  exit(1);
#endif

  // get the matrix pointers from the matrixID, and see if still in store 
  mat_add_t PI_ptr = create_perm_inv(n);

  // call matrix pointer version of function
  int64_t PI_mID = get_matrixID_from_ptr(PI_ptr);
  return PI_mID;
}


/* The create_perm_inv(mat_level_t n) function
   takes the size of the matrix and returns a matrix PTR
   corresponding the inverse permution matrix */
mat_add_t create_perm_inv(mat_level_t n) 
{

#ifndef USE_COMPLEX
  printf("\nERROR in %s not using COMPLEX type!\n", __func__);
  exit(1);
#endif


  mat_add_t PI_ptr;

  // EXCEPTION CHECKING
  if (n < 0) {
    printf("ERROR in %s: the matrix level n=%d is too small\n", __func__,n);
    exit(1);
  }
  if (n > maximum_level()) {
    printf("ERROR in %s: the matrix level n=%d > max_level\n", __func__,n);
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
    mat_add_t PIsmall_ptr = create_perm_inv(n-1);
    mat_add_t small_panel[4];
    for (int i=0;i<4;++i) {
      small_panel[i] = matrix_sub(PIsmall_ptr, i);
    }
    mat_add_t panel[4];
    mat_add_t temp_array[4];
    temp_array[0] = small_panel[0];
    temp_array[1] = small_panel[1];
    temp_array[2] = get_zero_matrix_ptr(n-2,n-2);
    temp_array[3] = get_zero_matrix_ptr(n-2,n-2);
    panel[0] = matrix_get_ptr_panel(temp_array,n-1,n-1);

    temp_array[0] = get_zero_matrix_ptr(n-2,n-2);
    temp_array[1] = get_zero_matrix_ptr(n-2,n-2);
    temp_array[2] = small_panel[0];
    temp_array[3] = small_panel[1];
    panel[1] = matrix_get_ptr_panel(temp_array,n-1,n-1);

    temp_array[0] = small_panel[2];
    temp_array[1] = small_panel[3];
    temp_array[2] = get_zero_matrix_ptr(n-2,n-2);
    temp_array[3] = get_zero_matrix_ptr(n-2,n-2);
    panel[2] = matrix_get_ptr_panel(temp_array,n-1,n-1);

    temp_array[0] = get_zero_matrix_ptr(n-2,n-2);
    temp_array[1] = get_zero_matrix_ptr(n-2,n-2);
    temp_array[2] = small_panel[2];
    temp_array[3] = small_panel[3];
    panel[3] = matrix_get_ptr_panel(temp_array,n-1,n-1);

    PI_ptr = matrix_get_ptr_panel(panel,n,n);
    return (PI_ptr);
  }
}


/* prints the (2^n)-th roots of unity */
int print_pow2_roots_unity(mat_level_t n){

#ifndef USE_COMPLEX
  printf("\nERROR in %s not using COMPLEX type!\n", __func__);
  exit(1);
#endif

  double complex root;
  // long double complex rootlong;
 
  // printf("The size of double complex is %zd\n",sizeof(root));
  // printf("The size of long double complex is %zd\n\n",sizeof(rootlong));

  // create 2^n
  int64_t pow = 1L<<n;

  printf("The (2^%d)-th roots of unity are:\n",n);

  for (int k = 0; k< pow ; k++) {
    root = cexp(2 * M_PI * I * k / pow);
    printf("k=%d, %lg + %lg i\n",k, creal(root), cimag(root));
    // rootlong = cexpl(2 * M_PI * I * k / pow);
    // printf(" LONG:  k=%d, %Lg + %Lg i\n\n",k, creall(rootlong), cimagl(rootlong));
    // load the values into the store and then print them out again
    mat_add_t root_ptr = matrix_get_ptr_scalar(root);
    printf("\tHave loaded the root into the store, it is now:\n\t");
    print_matrix_naive(root_ptr);
    printf("\n");

  }

  return 0;

}


/* return the principal (2^n)-th root of unity */
int64_t principal_pow2_root_unity(mat_level_t n){

#ifndef USE_COMPLEX
  printf("\nERROR in %s not using COMPLEX type!\n", __func__);
  exit(1);
#endif

  double complex root;
  // long double complex rootlong;
 
  // printf("The size of double complex is %zd\n",sizeof(root));
  // printf("The size of long double complex is %zd\n\n",sizeof(rootlong));

  // create 2^n
  int64_t pow = 1L<<n;

  printf("The (2^%d)-th principal root of unity is:\n",n);

  root = cexp(2 * M_PI * I  / pow);
  printf("%lg + %lg i\n", creal(root), cimag(root));
  // rootlong = cexpl(2 * M_PI * I * k / pow);
  // printf(" LONG:  k=%d, %Lg + %Lg i\n\n",k, creall(rootlong), cimagl(rootlong));
  // load the values into the store and then print them out again
  mat_add_t root_ptr = matrix_get_ptr_scalar(root);
  int64_t root_mID = get_matrixID_from_ptr(root_ptr);

  return root_mID;

}




/* return the matrix pointer for the first (2^n)-th root of unity */
mat_add_t get_first_pow2_root_unity(mat_level_t n){

#ifndef USE_COMPLEX
  printf("\nERROR in %s not using COMPLEX type!\n", __func__);
  exit(1);
#endif

  // create 2^n
  int64_t pow = 1L<<n;

  ScalarType root;
  root = cexp(2 * M_PI * I  / pow);
  mat_add_t root_ptr = matrix_get_ptr_scalar(root);

  return root_ptr;
}



/* Python Interface for the create_fft_D(mat_level_t n) function.
   This takes the size of the matrix and returns a matrix ID
   corresponding the D matrix needed for the fft recursion */
int64_t 
create_fft_D_matrixID(mat_level_t n) {

#ifndef USE_COMPLEX
  printf("\nERROR in %s not using COMPLEX type!\n", __func__);
  exit(1);
#endif

  // get the matrix pointers from the matrixID, and see if still in store 
  mat_add_t D_ptr = create_fft_D(n);

  // call matrix pointer version of function
  int64_t D_mID = get_matrixID_from_ptr(D_ptr);
  return D_mID;
}


/* The create_fft_D(mat_level_t n) function
   takes the size of the matrix and returns a matrix PTR
   corresponding the D matrix for the fft recursion*/
mat_add_t create_fft_D(mat_level_t n) 
{

#ifndef USE_COMPLEX
  printf("\nERROR in %s not using COMPLEX type!\n", __func__);
  exit(1);
#endif

  mat_add_t D_ptr;

  // EXCEPTION CHECKING
  if (n < 0) {
    printf("ERROR in %s: the matrix level n=%d is too small\n", __func__,n);
    exit(1);
  }
  if (n > maximum_level()) {
    printf("ERROR in %s: the matrix size n=%d > max_level\n", __func__,n);
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
    mat_add_t panel[4];
    panel[0] = get_identity_matrix_ptr(0);
    panel[1] = get_zero_matrix_ptr(0,0);
    panel[2] = get_zero_matrix_ptr(0,0);
    panel[3] = get_first_pow2_root_unity(n+1);
    mat_add_t A_ptr = matrix_get_ptr_panel(panel,1,1);
    mat_add_t Dminus1_ptr = create_fft_D(n-1); 
    // D_ptr = kronecker_product(A_ptr, Dminus1_ptr);
    D_ptr = kronecker_product(Dminus1_ptr,A_ptr);

    return (D_ptr);
  }
}



/* Python Interface for the create_fft_C(mat_level_t n) function.
   This takes the size of the matrix and returns a matrix ID
   corresponding the C matrix needed for the fft recursion */
int64_t 
create_fft_C_matrixID(mat_level_t n) {

#ifndef USE_COMPLEX
  printf("\nERROR in %s not using COMPLEX type!\n", __func__);
  exit(1);
#endif

  // get the matrix pointers from the matrixID, and see if still in store 
  mat_add_t C_ptr = create_fft_C(n);

  // call matrix pointer version of function
  int64_t C_mID = get_matrixID_from_ptr(C_ptr);
  return C_mID;
}


/* The create_fft_C(mat_level_t n) function
   takes the size of the matrix and returns a matrix PTR
   corresponding the C matrix for the fft recursion*/
mat_add_t create_fft_C(mat_level_t n) 
{

#ifndef USE_COMPLEX
  printf("\nERROR in %s not using COMPLEX type!\n", __func__);
  exit(1);
#endif

  mat_add_t C_ptr;

  // EXCEPTION CHECKING
  if (n < 1) {
    printf("ERROR in %s: the matrix level n=%d is too small\n", __func__,n);
    exit(1);
  }
  if (n > maximum_level()) {
    printf("ERROR in %s: the matrix size n=%d > max_level\n", __func__,n);
    exit(1);
  }

  // SMALL CASE
  if (n == 1) {
    C_ptr = get_iHadamard_matrix_ptr(1);
    return (C_ptr);
  }

  // ELSE build C from panels [I(n-1),D(n-1);I(n-1),-D(n-1)]
  else {
    mat_add_t Dminus1_ptr = create_fft_D(n-1); 
    ScalarType scalarM1 = -1 + 0*I;
    mat_add_t matptr_scalarM1 = matrix_get_ptr_scalar(scalarM1);
    mat_add_t minusDminus1_ptr = scalar_mult(matptr_scalarM1,Dminus1_ptr); 
    mat_add_t panel[4];
    panel[0] = get_identity_matrix_ptr(n-1);
    panel[1] = Dminus1_ptr;
    panel[2] = get_identity_matrix_ptr(n-1);
    panel[3] = minusDminus1_ptr;
    mat_add_t C_ptr = matrix_get_ptr_panel(panel,n,n);

    return (C_ptr);
  }
}





/* Python Interface for the create_fft_matrix(mat_level_t n) function.
   This takes the size of the matrix and returns a matrix ID
   corresponding to the fft matrix */
int64_t 
create_fft_matrix_matrixID(mat_level_t n) {

#ifndef USE_COMPLEX
  printf("\nERROR in %s not using COMPLEX type!\n", __func__);
  exit(1);
#endif

  // get the matrix pointers from the matrixID, and see if still in store 
  mat_add_t F_ptr = create_fft_matrix(n);

  // call matrix pointer version of function
  int64_t F_mID = get_matrixID_from_ptr(F_ptr);
  return F_mID;
}


/* The create_fft_matrix(mat_level_t n) function
   takes the size of the matrix and returns a matrix PTR
   corresponding to the fft matrix */
mat_add_t create_fft_matrix(mat_level_t n) 
{

#ifndef USE_COMPLEX
  printf("\nERROR in %s not using COMPLEX type!\n", __func__);
  exit(1);
#endif

  mat_add_t F_ptr;

  // EXCEPTION CHECKING
  if (n < 0) {
    printf("ERROR in %s: the matrix level n=%d is too small\n", __func__,n);
    exit(1);
  }
  if (n > maximum_level()) {
    printf("ERROR in %s: the matrix size n=%d > max_level\n", __func__,n);
    exit(1);
  }

  // SMALL CASES
  if (n == 0) {
    F_ptr = get_identity_matrix_ptr(0);
    return (F_ptr);
  }

  // RECURSIVE CASE   F(n) = C(n) * [I(1) @ F(n-1)] * PI(n)
  else {
    mat_add_t Cn_ptr = create_fft_C(n);
    mat_add_t PIn_ptr = create_perm_inv(n);
    mat_add_t Fn_minus1_ptr = create_fft_matrix(n-1);
    mat_add_t blockF_ptr = kronecker_product(
                       get_identity_matrix_ptr(1),Fn_minus1_ptr);
    F_ptr = matrix_mult(Cn_ptr,matrix_mult(blockF_ptr,PIn_ptr));

    return (F_ptr);
  }
}




