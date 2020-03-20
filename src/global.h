//                      global.h
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

#ifndef GLOBAL_H
#define GLOBAL_H

#include <math.h>
#include <pthread.h> // pthread_t
#include <stdio.h>
#include <gmp.h>  

#include "io.h"

// Unfortunately math.h doesn't contain a max and min
#define MAX(x,y) ((x > y) ? x : y)
#define MIN(x,y) ((x < y) ? x : y)

#ifndef SWIG

// Temporary variables for calculations that would otherwise initialize and
// clear local variables many times. In particular, this is a problem when 
// the scalarType happens to be one of the multi-precision ones.  
//
// In general, we should leave these ready to go (initialized) since routines
// like scalar_mult use without checking. Routines like matrix_mult can re-init
// (for clearing) at top level if desired. 
//
// Use these variables as you like. But don't be surprised if they change after
// calling a routine because that routine might also make use of them!
// These variables are NOT thread safe. Future work could include a locking
// mechanism to make them more safe. 
//
// See dev_notes for more info. 
typedef struct scratchVars_s {
    // for saving the level of the original call of a recursive routine
    mat_level_t top_level;

    // designed [to be calculated and then] submitted to get_valMatPTR_from_val
    // guarantee that submit_to_store is never used in get_valMatPTR_from_val
    scalarType submit_to_store;
    // for use when computing the conjugate
    scalarType calc_conj;
    // for use when you're sure no subroutines in it's scope also use quick_use
    scalarType quick_use;
    // for general use or for some unknown purpose we haven't thought of yet. 
    scalarType misc;

    // for counting some value that accrues over multiple [recursive] calls
    mpz_t counter;     
    // for use in locality approximation rounding 
    scalarType approx_value;
#if defined(USE_MPRATIONAL) || defined(USE_MPRATCOMPLEX) || defined(USE_MPREAL) || defined(USE_MPCOMPLEX)
    mpz_t mpinteger;
#endif
#if defined(USE_MPRATIONAL) || defined(USE_MPRATCOMPLEX)
    mpq_t mprational;
#endif
#if defined(USE_MPINTEGER) || defined(USE_MPREAL) || defined(USE_MPCOMPLEX)
    mpfr_t mpreal;
#endif

} scratchVars_t;

extern scratchVars_t scratchVars;

#endif

/*!
 * \ingroup larc
 *
 * \brief Initialize struct of frequent use variables in matrix math routines.
 * 
 * Note that scalar routines must be loaded before this can run. 
 */
void scratchVars_init();
/*!
 * \ingroup larc
 *
 * \brief Used in cleaning routines to indicate recursion is finished and
 * scratch variables can be used again.
 *
 * Sets the top_level member of a scratchVars_t strucure to -1.
 * 
 * \param current_level The current level of recursion.
 */
void scratchVars_exitroutine(mat_level_t current_level);

extern pthread_t thd;
/*!
 * \ingroup larc
 * \brief Gives names to certain global constants and 2x2 matrices, and when necessary stores and locks them.
 */
void init_globals(void);

// this passes information to python about current scalarType
extern char scalarTypeDef;
extern char* scalarTypeStr;


#ifndef SWIG
/***************************************************
 * FREQUENTLY USED SCALARS ARE GIVEN GLOBAL NAMES *
 ***************************************************/

/*********************************************
 *       values of SCALAR type               * 
 *********************************************/
// All SCALAR types have 0 and 1 and can handle -1.
extern scalarType scalar0;
extern scalarType scalar1;
extern scalarType scalarM1; //'minus 1'
#if defined(USE_REAL) || defined(USE_COMPLEX) || defined(USE_MPREAL) || defined(USE_MPCOMPLEX) || defined(USE_MPRATIONAL) || defined(USE_MPRATCOMPLEX)
extern scalarType scalar0_5;
extern scalarType scalarM0_5;
extern scalarType scalar_inv_sqrt_2;
extern scalarType scalar_neg_inv_sqrt_2;
#endif
#if defined(USE_COMPLEX) || defined(USE_MPCOMPLEX) || defined(USE_MPRATCOMPLEX)
extern scalarType scalar0i1;
extern scalarType scalar0iM1;
extern scalarType scalar0i0_5;
extern scalarType scalar0iM0_5;
extern scalarType scalar_i_inv_sqrt_2;
extern scalarType scalar_plus_root_i;
extern scalarType scalar_minus_root_i;
#endif 
#endif


          /***************************************************
           * FREQUENTLY USED MATRICES ARE GIVEN GLOBAL NAMES *
           ***************************************************/

/*********************************************
 *       matrices of SCALAR type             * 
 ********************************************/
extern int64_t matID_scalar0;
extern int64_t matID_scalar1;
extern int64_t matID_scalarM1;
#if defined(USE_REAL) || defined(USE_COMPLEX) || defined(USE_MPREAL) || defined(USE_MPCOMPLEX) || defined(USE_MPRATIONAL) || defined(USE_MPRATCOMPLEX)
extern int64_t matID_scalar0_5;
extern int64_t matID_scalarM0_5;
extern int64_t matID_inv_sqrt_2;
extern int64_t matID_neg_inv_sqrt_2;
#endif  
#if defined(USE_COMPLEX) || defined(USE_MPCOMPLEX) || defined(USE_MPRATCOMPLEX)
extern int64_t matID_scalar0i1;   // scalar 0.0 + i*1.0
extern int64_t matID_scalar0iM1;  // scalar 0.0 - i*1.0
extern int64_t matID_scalar0i0_5;
extern int64_t matID_scalar0iM0_5;
extern int64_t matID_i_inv_sqrt_2;
extern int64_t matID_plus_root_i;
extern int64_t matID_minus_root_i;
#endif  

/*****************************************************
 *   square matrices of MATRIX type level > 0        *
 *****************************************************/
extern int64_t matID_I1;     // 1-bit identity matrix
extern int64_t matID_NOT;     // 1-bit NOT / bitflip matrix
extern int64_t matID_HH1;    // sqrt(2) * Hadamard matrix
#if defined(USE_REAL) || defined(USE_COMPLEX) || defined(USE_MPREAL) || defined(USE_MPCOMPLEX) || defined(USE_MPRATIONAL) || defined(USE_MPRATCOMPLEX)
extern int64_t matID_H1;     // Hadamard matrix
#endif



#endif

