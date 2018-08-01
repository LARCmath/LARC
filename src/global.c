//                        global.c 
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

#include "larc.h"
#include "global.h"
#include "matrix_store.h"

pthread_t thd;

// this passes information to python about current scalarType
char scalarTypeDef;
#ifdef USE_INTEGER
char scalarTypeDef = 'i';
#endif
#ifdef USE_COMPLEX
char scalarTypeDef = 'c';
#endif
#ifdef USE_REAL
char scalarTypeDef = 'r';
#endif

/***************************************************
 * FREQUENTLY USED MATRICES ARE GIVEN GLOBAL NAMES *
 ***************************************************/

/*********************************************
 *       matrices of SCALAR type             * 
 *********************************************/
int64_t matID_scalar0;
int64_t matID_scalar1;
int64_t matID_scalarM1;
#ifndef USE_INTEGER //REAL OR COMPLEX
int64_t matID_scalar0_5;
int64_t matID_inv_sqrt_2;
int64_t matID_neg_inv_sqrt_2;
int64_t matID_scalarM0_5;
#ifdef USE_COMPLEX
int64_t matID_scalar0i1;
int64_t matID_scalar0iM1;
int64_t matID_plus_root_i;
int64_t matID_minus_root_i;
int64_t matID_i_inv_sqrt_2;
int64_t matID_scalar0i0_5;
int64_t matID_scalar0iM0_5;
#endif
#endif //REAL OR COMPLEX

/*****************************************************
 *   square matrices of MATRIX type level > 0        *
 *****************************************************/
int64_t matID_I1;
int64_t matID_NOT;
int64_t matID_HH1;

#ifndef USE_INTEGER //COMPLEX OR REAL
int64_t matID_H1;
#endif

void init_globals(void) {

  int verbose = 0;
/*********************************************
 *       matrices of SCALAR type             * 
 ********************************************/
  //integer with no imaginary part
  mat_add_t matptr_scalar0;
  mat_add_t matptr_scalar1;
  mat_add_t matptr_scalarM1;

#ifndef USE_INTEGER //REAL OR COMPLEX
  //noninteger with no imaginary part
  mat_add_t matptr_scalar0_5;
  mat_add_t matptr_inv_sqrt_2;
  mat_add_t matptr_neg_inv_sqrt_2;
  mat_add_t matptr_scalarM0_5;

#ifdef USE_COMPLEX
  //values that have imaginary parts
  mat_add_t matptr_scalar0i1;
  mat_add_t matptr_scalar0iM1;
  mat_add_t matptr_plus_root_i;
  mat_add_t matptr_minus_root_i;
  mat_add_t matptr_i_inv_sqrt_2;
  mat_add_t matptr_scalar0i0_5;
  mat_add_t matptr_scalar0iM0_5;
#endif

#endif //REAL OR COMPLEX
/***********************************
 *       2x2 matrices              * 
 ***********************************/
  mat_add_t matptr_I1; // 1-bit identity matrix
  mat_add_t matptr_NOT; // 1-bit NOT / bitflip matrix
  // sqrt(2) * Hadamard matrix 
  // aka iHadamard or integer Hadamard
  mat_add_t matptr_HH1;
#ifndef USE_INTEGER //COMPLEX OR REAL
  mat_add_t matptr_H1; //Hadamard matrix  
#endif

#ifndef USE_COMPLEX //REAL OR INTEGER

/*****************************************************************
 *   matrices of SCALAR type only used when scalars are real     *
 ****************************************************************/

  ScalarType scalar0 = 0;
  matptr_scalar0 = matrix_get_ptr_scalar(scalar0);
  matID_scalar0 = matptr_scalar0->matrixID;
  // we have already locked this value while creating the matrix store

  ScalarType scalar1 = 1;
  matptr_scalar1 = matrix_get_ptr_scalar(scalar1);
  matID_scalar1 = matptr_scalar1->matrixID;
  // we have already locked this value while creating the matrix store

  ScalarType scalarM1 = -1;
  matptr_scalarM1 = matrix_get_ptr_scalar(scalarM1);
  matID_scalarM1 = matptr_scalarM1->matrixID;
  // we have already locked this value while loading iHadamards

#ifdef USE_REAL
  //Start: these values have decimal components, but no imaginary component

  ScalarType scalar0_5 = 0.5;
  matptr_scalar0_5 = matrix_get_ptr_scalar(scalar0_5);
  lock_matrix(matptr_scalar0_5);
  matID_scalar0_5 = matptr_scalar0_5->matrixID;

  ScalarType inv_sqrt_2 = 1/M_SQRT2;
  matptr_inv_sqrt_2 = matrix_get_ptr_scalar(inv_sqrt_2);
  lock_matrix(matptr_inv_sqrt_2);
  matID_inv_sqrt_2 = matptr_inv_sqrt_2->matrixID;

  ScalarType neg_inv_sqrt_2 = -1/M_SQRT2;
  matptr_neg_inv_sqrt_2 = matrix_get_ptr_scalar(neg_inv_sqrt_2);
  lock_matrix(matptr_neg_inv_sqrt_2);
  matID_neg_inv_sqrt_2 = matptr_neg_inv_sqrt_2->matrixID;

  ScalarType scalarM0_5 = -0.5;
  matptr_scalarM0_5 = matrix_get_ptr_scalar(scalarM0_5);
  lock_matrix(matptr_scalarM0_5);
  matID_scalarM0_5 = matptr_scalarM0_5->matrixID;
  //End: these values have no imaginary component

#endif

#else // USE_COMPLEX is defined

/*****************************************************************
 *   matrices of SCALAR type only used when scalars are complex  *
 ****************************************************************/

  //Start: these values have no imaginary component

  ScalarType scalar0 = 0 + 0*I;
  matptr_scalar0 = matrix_get_ptr_scalar(scalar0);
  matID_scalar0 = matptr_scalar0->matrixID;
  // we have already locked this value while creating the matrix store

  ScalarType scalar1 = 1 + 0*I;
  matptr_scalar1 = matrix_get_ptr_scalar(scalar1);
  matID_scalar1 = matptr_scalar1->matrixID;
  // we have already locked this value while creating the matrix store

  ScalarType scalarM1 = -1 + 0*I;
  matptr_scalarM1 = matrix_get_ptr_scalar(scalarM1);
  matID_scalarM1 = matptr_scalarM1->matrixID;
  // we have already locked this value while loading iHadamards

  ScalarType scalar0_5 = 0.5 + 0*I;
  matptr_scalar0_5 = matrix_get_ptr_scalar(scalar0_5);
  lock_matrix(matptr_scalar0_5);
  matID_scalar0_5 = matptr_scalar0_5->matrixID;

  ScalarType inv_sqrt_2 = 1/M_SQRT2 + 0*I;
  matptr_inv_sqrt_2 = matrix_get_ptr_scalar(inv_sqrt_2);
  lock_matrix(matptr_inv_sqrt_2);
  matID_inv_sqrt_2 = matptr_inv_sqrt_2->matrixID;

  ScalarType neg_inv_sqrt_2 = -1/M_SQRT2 + 0*I;
  matptr_neg_inv_sqrt_2 = matrix_get_ptr_scalar(neg_inv_sqrt_2);
  lock_matrix(matptr_neg_inv_sqrt_2);
  matID_neg_inv_sqrt_2 = matptr_neg_inv_sqrt_2->matrixID;

  ScalarType scalarM0_5 = -0.5 + 0*I;
  matptr_scalarM0_5 = matrix_get_ptr_scalar(scalarM0_5);
  lock_matrix(matptr_scalarM0_5);
  matID_scalarM0_5 = matptr_scalarM0_5->matrixID;
  //End: these values have no imaginary component

  ScalarType scalar0i1 = 0 + 1*I;
  matptr_scalar0i1 = matrix_get_ptr_scalar(scalar0i1);
  lock_matrix(matptr_scalar0i1);
  matID_scalar0i1 = matptr_scalar0i1->matrixID;

  ScalarType scalar0iM1 = 0 - 1*I;
  matptr_scalar0iM1 = matrix_get_ptr_scalar(scalar0iM1);
  lock_matrix(matptr_scalar0iM1);
  matID_scalar0iM1 = matptr_scalar0iM1->matrixID;

  ScalarType plus_root_i = 1/M_SQRT2 + I/M_SQRT2; 
  matptr_plus_root_i = matrix_get_ptr_scalar(plus_root_i);
  lock_matrix(matptr_plus_root_i);
  matID_plus_root_i = matptr_plus_root_i->matrixID;

  ScalarType minus_root_i = 1/M_SQRT2 - I/M_SQRT2;
  matptr_minus_root_i = matrix_get_ptr_scalar(minus_root_i);
  lock_matrix(matptr_minus_root_i);
  matID_minus_root_i = matptr_minus_root_i->matrixID;

  ScalarType i_inv_sqrt_2 = 1/M_SQRT2 * I;
  matptr_i_inv_sqrt_2 = matrix_get_ptr_scalar(i_inv_sqrt_2);
  lock_matrix(matptr_i_inv_sqrt_2);
  matID_i_inv_sqrt_2 = matptr_i_inv_sqrt_2->matrixID;

  ScalarType scalar0i0_5 = 0.0 + 0.5*I;
  matptr_scalar0i0_5 = matrix_get_ptr_scalar(scalar0i0_5);
  lock_matrix(matptr_scalar0i0_5);
  matID_scalar0i0_5 = matptr_scalar0i0_5->matrixID;

  ScalarType scalar0iM0_5 = 0.0 - 0.5*I;
  matptr_scalar0iM0_5 = matrix_get_ptr_scalar(scalar0iM0_5);
  lock_matrix(matptr_scalar0iM0_5);
  matID_scalar0iM0_5 = matptr_scalar0iM0_5->matrixID;

#endif //COMPLEX OR REAL
  
  /****************************************************************************
   * Defining and loading content of frequently used matrices of MATRIX type  *
   ****************************************************************************/
  if (verbose) {
    printf("  Defining and loading frequently used matrices.\n");
  }
 
/***********************************
 * 2 by 2 square matrices, level 1 *
 ***********************************/
  if (verbose) {
    printf("  Defining content of frequently used 2 by 2  matrices\n");
  }
  
  // Load into matrix store and set global names for matrix ids
  
// 1-bit identity matrix
  matptr_I1 = get_identity_matrix_ptr(1);
  matID_I1 = get_identity_matrixID(1);

// 1-bit NOT / bitflip matrix
  mat_add_t submat[4];
  submat[0] = submat[3] = matptr_scalar0;
  submat[1] = submat[2] = matptr_scalar1;
  matptr_NOT = matrix_get_ptr_panel(submat,1,1);
  lock_matrix(matptr_NOT);
  matID_NOT = matptr_NOT->matrixID;

#ifndef USE_INTEGER 
//Hadamard matrix  
  submat[0] = submat[1] = submat[2] = matptr_inv_sqrt_2;
  submat[3] = matptr_neg_inv_sqrt_2; 
  matptr_H1 = matrix_get_ptr_panel(submat,1,1);
  lock_matrix(matptr_H1);
  matID_H1 = matptr_H1->matrixID;
#endif

// sqrt(2) * Hadamard matrix 
// aka iHadamard or integer Hadamard
  matptr_HH1 = get_iHadamard_matrix_ptr(1);
  matID_HH1 = matptr_HH1->matrixID;

 if (verbose) {
  printf("finishing routine %s\n",__func__);
 }

}

