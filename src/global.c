//                        global.c 
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

#include "larc.h"
#include "global.h"
#include "matrix_store.h"
#include "scalars.h"
#include "mar.h"

/*!
 * \file global.c
 * \brief This file contains the global variables used by LARC and also
 * functions which initialize them.
 */

const mpfr_prec_t mpreal_precision = 256;

// this passes information to python about current scalarType
char scalarTypeDef;
char *scalarTypeStr;
#ifdef USE_INTEGER
char scalarTypeDef = 'i';
char *scalarTypeStr = "Integer";
#endif // USE_INTEGER
#ifdef USE_BOOLEAN
char scalarTypeDef = 'b';
char *scalarTypeStr = "Boolean";
#endif // USE_BOOLEAN
#ifdef USE_COMPLEX
char scalarTypeDef = 'c';
char *scalarTypeStr = "Complex";
#endif // USE_COMPLEX
#ifdef USE_REAL
char scalarTypeDef = 'r';
char *scalarTypeStr = "Real";
#endif // USE_REAL
#ifdef USE_MPINTEGER
char scalarTypeDef = 'z';
char *scalarTypeStr = "MPInteger";
#endif // USE_MPINTEGER
#ifdef USE_MPRATIONAL
char scalarTypeDef = 'q';
char *scalarTypeStr = "MPRational";
#endif // USE_MPRATIONAL
#ifdef USE_MPRATCOMPLEX
char scalarTypeDef = 'v';
char *scalarTypeStr = "MPRatComplex";
#endif // USE_MPRATCOMPLEX
#ifdef USE_MPREAL
char scalarTypeDef = 'm';
char *scalarTypeStr = "MPReal";
#endif // USE_MPREAL
#ifdef USE_MPCOMPLEX
char scalarTypeDef = 'j';
char *scalarTypeStr = "MPComplex";
#endif // USE_MPCOMPLEX
#ifdef USE_CLIFFORD
char scalarTypeDef = 'a';
char *scalarTypeStr = "Clifford";
#endif // USE_CLIFFORD
#ifdef USE_LOWER
char scalarTypeDef = 'l';
char *scalarTypeStr = "Lower";
#endif // USE_LOWER
#ifdef USE_UPPER
char scalarTypeDef = 'u';
char *scalarTypeStr = "Upper";
#endif // USE_UPPER

scratchVars_t scratchVars = {0};

void scratchVars_init()
{
    // Note that scalar routines need to be loaded before this can be run. 
    if (sca_init == NULL){
        fprintf(stderr,"Scalar ops must be initialized before running matrix math routines.\n");
        exit(1);
    }

    scratchVars.top_level = -1;

    sca_init(&scratchVars.submit_to_store);
    scratchVars.submit_to_store_in_use = 0;
    sca_init(&scratchVars.calc_conj);
    scratchVars.calc_conj_in_use = 0;
    sca_init(&scratchVars.quick_use);
    scratchVars.quick_use_in_use = 0;
    sca_init(&scratchVars.misc);
    scratchVars.misc_in_use = 0;
    sca_init(&scratchVars.approx_value);
    scratchVars.approx_value_in_use = 0;

    mpz_init(scratchVars.counter);
    scratchVars.counter_in_use = 0;
    mpz_init(scratchVars.mpinteger);
    scratchVars.mpinteger_in_use = 0;
    mpq_init(scratchVars.mprational);
    scratchVars.mprational_in_use = 0;
    mpq_init(scratchVars.mprational2);
    scratchVars.mprational2_in_use = 0;
    mpfr_init2(scratchVars.mpreal, mpreal_precision);
    scratchVars.mpreal_in_use = 0;
#ifdef USE_CLIFFORD
    larc_sca_init_clifford(&(scratchVars.clifford));
    scratchVars.clifford_in_use = 0;
    for (int i=0; i<CLIFFORD_DIMENSION-1; ++i)
    {
        mpfr_init(mpfr_algebraic_approx[i]);
        mpq_init(mpq_algebraic_approx[i]);
    }
#endif // USE_CLIFFORD

#ifdef MAR
    scratchVars_mar_init();
#endif // MARmode
    
}

void scratchVars_exitroutine(mat_level_t current_level)
{
    if (scratchVars.top_level == current_level)
        scratchVars.top_level = -1;
}


/***************************************************
 * FREQUENTLY USED SCALARS ARE GIVEN GLOBAL NAMES *
 ***************************************************/

/*********************************************
 *       values of SCALAR type               * 
 *********************************************/
// All SCALAR types have 0 and 1 and can handle -1.
scalarType scalar0;
scalarType scalar1;
scalarType scalarM1; //'minus 1'
#ifndef IS_INTEGER
scalarType scalar0_5;
scalarType scalarM0_5;
#endif
#ifdef IS_COMPLEX
scalarType scalar0i1;
scalarType scalar0iM1;
scalarType scalar0i0_5;
scalarType scalar0iM0_5;
#endif 


/***************************************************
 * FREQUENTLY USED MATRICES ARE GIVEN GLOBAL NAMES *
 ***************************************************/

/*********************************************
 *       matrices of SCALAR type             * 
 *********************************************/
// All SCALAR types have 0 and 1 and can handle -1.
int64_t packedID_scalar0;
int64_t packedID_scalar1;
int64_t packedID_scalarM1; //'minus 1'
#ifndef IS_INTEGER
int64_t packedID_scalar0_5;
int64_t packedID_scalarM0_5;
#endif
#ifdef IS_COMPLEX
int64_t packedID_scalar0i1;
int64_t packedID_scalar0iM1;
int64_t packedID_scalar0i0_5;
int64_t packedID_scalar0iM0_5;
#endif 

/*****************************************************
 *   square matrices of MATRIX type level > 0        *
 *****************************************************/
int64_t packedID_I1;  // 1-bit identity matrix
int64_t packedID_NOT; // 1-bit NOT / bitflip matrix
  // sqrt(2) * Hadamard matrix 
  // aka iHadamard or integer Hadamard
int64_t packedID_HH1;

void init_globals(void) {

  int verbose = 0;

/*********************************************
 *       matrices of SCALAR type             * 
 ********************************************/
  mats_ptr_t scaptr_scalar0;
  mats_ptr_t scaptr_scalar1;
  mats_ptr_t scaptr_scalarM1;
#ifndef IS_INTEGER
  mats_ptr_t scaptr_scalar0_5;
  mats_ptr_t scaptr_scalarM0_5;
#endif
#ifdef IS_COMPLEX
  mats_ptr_t scaptr_scalar0i1;
  mats_ptr_t scaptr_scalar0iM1;
  mats_ptr_t scaptr_scalar0i0_5;
  mats_ptr_t scaptr_scalar0iM0_5;
#endif 

/*****************************************************************
 *   matrices of SCALAR type                                     *
 ****************************************************************/
  sca_init(&scalar0);
  sca_set_2ldoubles(&scalar0, 0.0L, 0.0L);
  scaptr_scalar0 = get_scalarPTR_for_scalarVal(scalar0);
  packedID_scalar0 = scaptr_scalar0->packedID;
  // we have already locked this value while preloading the matrix store

  sca_init(&scalar1);
  sca_set_2ldoubles(&scalar1, 1.0L, 0.0L);
  scaptr_scalar1 = get_scalarPTR_for_scalarVal(scalar1);
  packedID_scalar1 = scaptr_scalar1->packedID;
  // we have already locked this value while preloading the matrix store

  sca_init(&scalarM1);
  sca_set_2ldoubles(&scalarM1, -1.0L, 0.0L);
  scaptr_scalarM1 = get_scalarPTR_for_scalarVal(scalarM1);
  packedID_scalarM1 = scaptr_scalarM1->packedID;
  // we have already locked this value while preloading iHadamards

#ifndef IS_INTEGER
  sca_init(&scalar0_5);
  sca_set_2ldoubles(&scalar0_5, 0.5L, 0.0L);
  scaptr_scalar0_5 = get_scalarPTR_for_scalarVal(scalar0_5);
  packedID_scalar0_5 = scaptr_scalar0_5->packedID;
  lock_matrix(packedID_scalar0_5);

  sca_init(&scalarM0_5);
  sca_set_2ldoubles(&scalarM0_5, -0.5L, 0.0L);
  scaptr_scalarM0_5 = get_scalarPTR_for_scalarVal(scalarM0_5);
  packedID_scalarM0_5 = scaptr_scalarM0_5->packedID;
  lock_matrix(packedID_scalarM0_5);

#endif

#ifdef IS_COMPLEX
/*****************************************************************
 *   matrices of SCALAR type only used when scalars are complex  *
 ****************************************************************/
  sca_init(&scalar0i1);
  sca_set_2ldoubles(&scalar0i1, 0.0L, 1.0L);
  scaptr_scalar0i1 = get_scalarPTR_for_scalarVal(scalar0i1);
  packedID_scalar0i1 = scaptr_scalar0i1->packedID;
  lock_matrix(packedID_scalar0i1);

  sca_init(&scalar0iM1);
  sca_set_2ldoubles(&scalar0iM1, 0.0L, -1.0L);
  scaptr_scalar0iM1 = get_scalarPTR_for_scalarVal(scalar0iM1);
  packedID_scalar0iM1 = scaptr_scalar0iM1->packedID;
  lock_matrix(packedID_scalar0iM1);

  sca_init(&scalar0i0_5);
  sca_set_2ldoubles(&scalar0i0_5, 0.0L, 0.5L);
  scaptr_scalar0i0_5 = get_scalarPTR_for_scalarVal(scalar0i0_5);
  packedID_scalar0i0_5 = scaptr_scalar0i0_5->packedID;
  lock_matrix(packedID_scalar0i0_5);

  sca_init(&scalar0iM0_5);
  sca_set_2ldoubles(&scalar0iM0_5, 0.0L, -0.5L);
  scaptr_scalar0iM0_5 = get_scalarPTR_for_scalarVal(scalar0iM0_5);
  packedID_scalar0iM0_5 = scaptr_scalar0iM0_5->packedID;
  lock_matrix(packedID_scalar0iM0_5);

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
  packedID_I1 = get_identity_pID(1);

  // 1-bit NOT / bitflip matrix
  int64_t sub_pID[4];
  sub_pID[0] = sub_pID[3] = packedID_scalar0;
  sub_pID[1] = sub_pID[2] = packedID_scalar1;
  packedID_NOT = get_pID_from_array_of_four_sub_pIDs(sub_pID, 1, 1);
  lock_matrix(packedID_NOT);

  // sqrt(2) * Hadamard matrix 
  // aka iHadamard or integer Hadamard
  packedID_HH1 = get_iHadamard_pID(1);


//  fprintf(stderr,"matrixID of last loaded by init_globals()");
//  fprintf(stderr," (2x2 Hadamard) is %" PRId64 "\n",matID_H1);
//  fprintf(stderr,"At this point, num_matrices_in_store returns %" PRId64 "\n",
//		num_matrices_in_store());

  if (verbose) {
   printf("finishing routine %s\n",__func__);
  }

}

void clear_globals(void)
{
    // Free memory allocated to global scalarType variables. (Unless
    // scalarType is multiprecision, sca_clear() is a no-op.) This
    // routine should only be called as part of a LARC shutdown procedure.

    sca_clear(&scalar0);
    sca_clear(&scalar1);
    sca_clear(&scalarM1);

#ifndef IS_INTEGER
    sca_clear(&scalar0_5);
    sca_clear(&scalarM0_5);
#endif

#ifdef IS_COMPLEX
    sca_clear(&scalar0i1);
    sca_clear(&scalar0iM1);
    sca_clear(&scalar0i0_5);
    sca_clear(&scalar0iM0_5);
#endif

}

void scratchVars_clear(void)
{
    // Clear the memory allocated to the global scalarType variables.
    // This routine should only be called as part of a LARC shutdown
    // procedure.

    fprintf(stderr,"freeing memory for scratchVars structure\n");
    sca_clear(&scratchVars.submit_to_store);
    sca_clear(&scratchVars.calc_conj);
    sca_clear(&scratchVars.quick_use);
    sca_clear(&scratchVars.misc);
    sca_clear(&scratchVars.approx_value);

    mpz_clear(scratchVars.counter);
    mpz_clear(scratchVars.mpinteger);
    mpq_clear(scratchVars.mprational);
    mpq_clear(scratchVars.mprational2);
    mpfr_clear(scratchVars.mpreal);
#ifdef USE_CLIFFORD
    larc_sca_clear_clifford(&(scratchVars.clifford));
    for (int i=0; i<CLIFFORD_DIMENSION-1; ++i)
    {
        mpfr_clear(mpfr_algebraic_approx[i]);
        mpq_clear(mpq_algebraic_approx[i]);
    }
#endif // USE_CLIFFORD
#ifdef MAR
    scratchVars_mar_clear();
#endif // MAR
}

