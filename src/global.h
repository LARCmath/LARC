//                      global.h
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

#ifndef GLOBAL_H
#define GLOBAL_H

#include <math.h>
#include <pthread.h> // pthread_t
#include <stdio.h>

#include "io.h"

// Unfortunately math.h doesn't contain a max and min
#define MAX(x,y) ((x > y) ? x : y)
#define MIN(x,y) ((x < y) ? x : y)

extern pthread_t thd;

void init_globals(void);

// this passes information to python about current scalarType
extern char scalarTypeDef;


          /***************************************************
           * FREQUENTLY USED MATRICES ARE GIVEN GLOBAL NAMES *
           ***************************************************/

/*********************************************
 *       matrices of SCALAR type             * 
 ********************************************/
extern int64_t matID_scalar0;
extern int64_t matID_scalar1;
extern int64_t matID_scalarM1;

#ifndef USE_INTEGER //REAL OR COMPLEX
extern int64_t matID_scalar0_5;
extern int64_t matID_scalarM0_5;
extern int64_t matID_inv_sqrt_2;
extern int64_t matID_neg_inv_sqrt_2;
#endif

#ifdef USE_COMPLEX

extern int64_t matID_scalar0i1;   // scalar 0.0 + i*1.0
extern int64_t matID_scalar0iM1;  // scalar 0.0 - i*1.0
extern int64_t matID_scalar0i0_5;
extern int64_t matID_scalar0iM0_5;
extern int64_t matID_plus_root_i;
extern int64_t matID_minus_root_i;
extern int64_t matID_i_inv_sqrt_2;

#endif   // USE_COMPLEX

/*****************************************************
 *   square matrices of MATRIX type level > 0        *
 *****************************************************/
extern int64_t matID_I1;     // 1-bit identity matrix
extern int64_t matID_NOT;     // 1-bit NOT / bitflip matrix
extern int64_t matID_HH1;    // sqrt(2) * Hadamard matrix

#ifndef USE_INTEGER //COMPLEX OR REAL
extern int64_t matID_H1;     // Hadamard matrix
#endif

#endif

