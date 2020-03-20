//                        global.c 
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

#include "larc.h"
#include "global.h"
#include "matrix_store.h"
#include "scalars.h"

pthread_t thd;

// this passes information to python about current scalarType
char scalarTypeDef;
char *scalarTypeStr;
#ifdef USE_INTEGER
char scalarTypeDef = 'i';
char *scalarTypeStr = "Integer";
#endif
#ifdef USE_COMPLEX
char scalarTypeDef = 'c';
char *scalarTypeStr = "Complex";
#endif
#ifdef USE_REAL
char scalarTypeDef = 'r';
char *scalarTypeStr = "Real";
#endif
#ifdef USE_MPINTEGER
char scalarTypeDef = 'z';
char *scalarTypeStr = "MPInteger";
#endif
#ifdef USE_MPRATIONAL
char scalarTypeDef = 'q';
char *scalarTypeStr = "MPRational";
#endif
#ifdef USE_MPRATCOMPLEX
char scalarTypeDef = 'v';
char *scalarTypeStr = "MPRatComplex";
#endif
#ifdef USE_MPREAL
char scalarTypeDef = 'm';
char *scalarTypeStr = "MPReal";
#endif
#ifdef USE_MPCOMPLEX
char scalarTypeDef = 'j';
char *scalarTypeStr = "MPComplex";
#endif


scratchVars_t scratchVars = {0};

void scratchVars_init()
{
    // Note that scalar routines need to be loaded before this can be run. 
    if (1 != check_scalarOps()){
        fprintf(stderr,"Scalar ops must be initialized before running matrix math routines.\n");
        exit(1);
    }

    scratchVars.top_level = -1;

    sca_init(&scratchVars.submit_to_store);
    sca_init(&scratchVars.calc_conj);
    sca_init(&scratchVars.quick_use);
    sca_init(&scratchVars.misc);
    sca_init(&scratchVars.approx_value);

    mpz_init(scratchVars.counter);
#if defined(USE_MPRATIONAL) || defined(USE_MPRATCOMPLEX) || defined(USE_MPREAL) || defined(USE_MPCOMPLEX)
    mpz_init(scratchVars.mpinteger);
#endif
#if defined(USE_MPRATIONAL) || defined(USE_MPRATCOMPLEX)
    mpq_init(scratchVars.mprational);
#endif
#if defined(USE_MPINTEGER) || defined(USE_MPREAL) || defined(USE_MPCOMPLEX)
    mpfr_prec_t mpreal_precision = 256;
    mpfr_init2(scratchVars.mpreal, mpreal_precision);
#endif
}

void scratchVars_exitroutine(mat_level_t current_level)
{
    if (scratchVars.top_level == current_level)
        scratchVars.top_level = -1;
}


/***************************************************
 * FREQUENTLY USED SCALARS ARE GIVEN GLOBAL NAMES *
 ***************************************************/

#define M_SQRT1_2L      0.7071067811865475244008443621048490L  /* 1/sqrt(2) */
#define M_SQRT1_2S     "0.7071067811865475244008443621048490392848359376884740365883398689953662392310510"  /* 1/sqrt(2) */


/*********************************************
 *       values of SCALAR type               * 
 *********************************************/
// All SCALAR types have 0 and 1 and can handle -1.
scalarType scalar0;
scalarType scalar1;
scalarType scalarM1; //'minus 1'
#if defined(USE_REAL) || defined(USE_COMPLEX) || defined(USE_MPREAL) || defined(USE_MPCOMPLEX) || defined(USE_MPRATIONAL) || defined(USE_MPRATCOMPLEX)
scalarType scalar0_5;
scalarType scalarM0_5;
scalarType scalar_inv_sqrt_2;
scalarType scalar_neg_inv_sqrt_2;
#endif
#if defined(USE_COMPLEX) || defined(USE_MPCOMPLEX) || defined(USE_MPRATCOMPLEX)
scalarType scalar0i1;
scalarType scalar0iM1;
scalarType scalar0i0_5;
scalarType scalar0iM0_5;
scalarType scalar_i_inv_sqrt_2;
scalarType scalar_plus_root_i;
scalarType scalar_minus_root_i;
#endif 


/***************************************************
 * FREQUENTLY USED MATRICES ARE GIVEN GLOBAL NAMES *
 ***************************************************/

/*********************************************
 *       matrices of SCALAR type             * 
 *********************************************/
// All SCALAR types have 0 and 1 and can handle -1.
int64_t matID_scalar0;
int64_t matID_scalar1;
int64_t matID_scalarM1; //'minus 1'
#if defined(USE_REAL) || defined(USE_COMPLEX) || defined(USE_MPREAL) || defined(USE_MPCOMPLEX) || defined(USE_MPRATIONAL) || defined(USE_MPRATCOMPLEX)
int64_t matID_scalar0_5;
int64_t matID_scalarM0_5;
int64_t matID_inv_sqrt_2;
int64_t matID_neg_inv_sqrt_2;
#endif
#if defined(USE_COMPLEX) || defined(USE_MPCOMPLEX) || defined(USE_MPRATCOMPLEX)
int64_t matID_scalar0i1;
int64_t matID_scalar0iM1;
int64_t matID_scalar0i0_5;
int64_t matID_scalar0iM0_5;
int64_t matID_i_inv_sqrt_2;
int64_t matID_plus_root_i;
int64_t matID_minus_root_i;
#endif 

/*****************************************************
 *   square matrices of MATRIX type level > 0        *
 *****************************************************/
int64_t matID_I1;  // 1-bit identity matrix
int64_t matID_NOT; // 1-bit NOT / bitflip matrix
  // sqrt(2) * Hadamard matrix 
  // aka iHadamard or integer Hadamard
int64_t matID_HH1;
#if defined(USE_REAL) || defined(USE_COMPLEX) || defined(USE_MPREAL) || defined(USE_MPCOMPLEX) || defined(USE_MPRATIONAL) || defined(USE_MPRATCOMPLEX)
int64_t matID_H1; //Hadamard matrix  
#endif

void init_globals(void) {

  int verbose = 0;
/*********************************************
 *       matrices of SCALAR type             * 
 ********************************************/
  mat_ptr_t matptr_scalar0;
  mat_ptr_t matptr_scalar1;
  mat_ptr_t matptr_scalarM1;
#if defined(USE_REAL) || defined(USE_COMPLEX) || defined(USE_MPREAL) || defined(USE_MPCOMPLEX) || defined(USE_MPRATIONAL) || defined(USE_MPRATCOMPLEX)
  mat_ptr_t matptr_scalar0_5;
  mat_ptr_t matptr_scalarM0_5;
  mat_ptr_t matptr_inv_sqrt_2;
  mat_ptr_t matptr_neg_inv_sqrt_2;
#endif
#if defined(USE_COMPLEX) || defined(USE_MPCOMPLEX) || defined(USE_MPRATCOMPLEX)
  mat_ptr_t matptr_scalar0i1;
  mat_ptr_t matptr_scalar0iM1;
  mat_ptr_t matptr_scalar0i0_5;
  mat_ptr_t matptr_scalar0iM0_5;
  mat_ptr_t matptr_i_inv_sqrt_2;
  mat_ptr_t matptr_plus_root_i;
  mat_ptr_t matptr_minus_root_i;
#endif 

/***********************************
 *       2x2 matrices              * 
 ***********************************/
  //mat_ptr_t matptr_I1; 
  mat_ptr_t matptr_NOT;
  mat_ptr_t matptr_HH1;
#if defined(USE_REAL) || defined(USE_COMPLEX) || defined(USE_MPREAL) || defined(USE_MPCOMPLEX) || defined(USE_MPRATIONAL) || defined(USE_MPRATCOMPLEX)
  mat_ptr_t matptr_H1;
#endif

/*****************************************************************
 *   matrices of SCALAR type                                     *
 ****************************************************************/
  sca_init(&scalar0);
  sca_set_str(&scalar0, "0");
  matptr_scalar0 = get_valMatPTR_from_val(scalar0);
  matID_scalar0 = matptr_scalar0->matrixID;
  // we have already locked this value while preloading the matrix store

  sca_init(&scalar1);
  sca_set_str(&scalar1, "1");
  matptr_scalar1 = get_valMatPTR_from_val(scalar1);
  matID_scalar1 = matptr_scalar1->matrixID;
  // we have already locked this value while preloading the matrix store

  sca_init(&scalarM1);
  sca_set_str(&scalarM1, "-1");
  matptr_scalarM1 = get_valMatPTR_from_val(scalarM1);
  matID_scalarM1 = matptr_scalarM1->matrixID;
  // we have already locked this value while preloading iHadamards

#if defined(USE_REAL) || defined(USE_COMPLEX) || defined(USE_MPREAL) || defined(USE_MPCOMPLEX) || defined(USE_MPRATIONAL) || defined(USE_MPRATCOMPLEX)
  sca_init(&scalar0_5);
  sca_set_2ldoubles(&scalar0_5, 0.5L, 0.0L);
  matptr_scalar0_5 = get_valMatPTR_from_val(scalar0_5);
  lock_matrix(matptr_scalar0_5);
  matID_scalar0_5 = matptr_scalar0_5->matrixID;

  sca_init(&scalarM0_5);
  sca_set_2ldoubles(&scalarM0_5, -0.5L, 0.0L);
  matptr_scalarM0_5 = get_valMatPTR_from_val(scalarM0_5);
  lock_matrix(matptr_scalarM0_5);
  matID_scalarM0_5 = matptr_scalarM0_5->matrixID;

  sca_init(&scalar_inv_sqrt_2);
#if defined(USE_MPRATIONAL) || defined(USE_MPRATCOMPLEX)
  sca_set_2ldoubles(&scalar_inv_sqrt_2, M_SQRT1_2L, 0.0L);
#else
  sca_set_str(&scalar_inv_sqrt_2, M_SQRT1_2S);
#endif
  matptr_inv_sqrt_2 = get_valMatPTR_from_val(scalar_inv_sqrt_2);
  lock_matrix(matptr_inv_sqrt_2);
  matID_inv_sqrt_2 = matptr_inv_sqrt_2->matrixID;

  sca_init(&scalar_neg_inv_sqrt_2);
  sca_mult(&scalar_neg_inv_sqrt_2, scalar_inv_sqrt_2, scalarM1);
  matptr_neg_inv_sqrt_2 = get_valMatPTR_from_val(scalar_neg_inv_sqrt_2);
  lock_matrix(matptr_neg_inv_sqrt_2);
  matID_neg_inv_sqrt_2 = matptr_neg_inv_sqrt_2->matrixID;
#endif

#if defined(USE_COMPLEX) || defined(USE_MPCOMPLEX) || defined(USE_MPRATCOMPLEX)
/*****************************************************************
 *   matrices of SCALAR type only used when scalars are complex  *
 ****************************************************************/
  sca_init(&scalar0i1);
  sca_set_2ldoubles(&scalar0i1, 0.0L, 1.0L);
  matptr_scalar0i1 = get_valMatPTR_from_val(scalar0i1);
  lock_matrix(matptr_scalar0i1);
  matID_scalar0i1 = matptr_scalar0i1->matrixID;

  sca_init(&scalar0iM1);
  sca_set_2ldoubles(&scalar0iM1, 0.0L, -1.0L);
  matptr_scalar0iM1 = get_valMatPTR_from_val(scalar0iM1);
  lock_matrix(matptr_scalar0iM1);
  matID_scalar0iM1 = matptr_scalar0iM1->matrixID;

  sca_init(&scalar0i0_5);
  sca_set_2ldoubles(&scalar0i0_5, 0.0L, 0.5L);
  matptr_scalar0i0_5 = get_valMatPTR_from_val(scalar0i0_5);
  lock_matrix(matptr_scalar0i0_5);
  matID_scalar0i0_5 = matptr_scalar0i0_5->matrixID;

  sca_init(&scalar0iM0_5);
  sca_set_2ldoubles(&scalar0iM0_5, 0.0L, -0.5L);
  matptr_scalar0iM0_5 = get_valMatPTR_from_val(scalar0iM0_5);
  lock_matrix(matptr_scalar0iM0_5);
  matID_scalar0iM0_5 = matptr_scalar0iM0_5->matrixID;

  sca_init(&scalar_i_inv_sqrt_2);
  sca_mult(&scalar_i_inv_sqrt_2, scalar_inv_sqrt_2, scalar0i1);
  matptr_i_inv_sqrt_2 = get_valMatPTR_from_val(scalar_i_inv_sqrt_2);
  lock_matrix(matptr_i_inv_sqrt_2);
  matID_i_inv_sqrt_2 = matptr_i_inv_sqrt_2->matrixID;

  sca_init(&scalar_plus_root_i);
  sca_add(&scalar_plus_root_i, scalar_inv_sqrt_2, scalar_i_inv_sqrt_2);
  matptr_plus_root_i = get_valMatPTR_from_val(scalar_plus_root_i);
  lock_matrix(matptr_plus_root_i);
  matID_plus_root_i = matptr_plus_root_i->matrixID;

  sca_init(&scalar_minus_root_i);
  sca_mult(&scalar_minus_root_i, scalar_plus_root_i, scalar0iM1);
  matptr_minus_root_i = get_valMatPTR_from_val(scalar_minus_root_i);
  lock_matrix(matptr_minus_root_i);
  matID_minus_root_i = matptr_minus_root_i->matrixID;
#endif 
  
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
  //matptr_I1 = get_identity_matrix_ptr(1);
  matID_I1 = get_identity_matrixID(1);

  // 1-bit NOT / bitflip matrix
  mat_ptr_t submat[4];
  submat[0] = submat[3] = matptr_scalar0;
  submat[1] = submat[2] = matptr_scalar1;
  matptr_NOT = get_matPTR_from_array_of_four_subMatPTRs(submat, 1, 1);
  lock_matrix(matptr_NOT);
  matID_NOT = matptr_NOT->matrixID;

  // sqrt(2) * Hadamard matrix 
  // aka iHadamard or integer Hadamard
  matptr_HH1 = get_iHadamard_matrix_ptr(1);
  matID_HH1 = matptr_HH1->matrixID;
#if defined(USE_REAL) || defined(USE_COMPLEX) || defined(USE_MPREAL) || defined(USE_MPCOMPLEX) || defined(USE_MPRATIONAL) || defined(USE_MPRATCOMPLEX)
  // Hadamard matrix  
  submat[0] = submat[1] = submat[2] = matptr_inv_sqrt_2;
  submat[3] = matptr_neg_inv_sqrt_2; 
  matptr_H1 = get_matPTR_from_array_of_four_subMatPTRs(submat, 1, 1);
  lock_matrix(matptr_H1);
  matID_H1 = matptr_H1->matrixID;
#endif

 if (verbose) {
  printf("finishing routine %s\n",__func__);
 }

}

