//                         spr.h
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

#ifndef LARC_SPR_H
#define LARC_SPR_H

#include <stdio.h> // FILE
#include <stddef.h> // size_t
#include <stdint.h> // int8_t, int32_t, uint32_t, uint64_t, uintptr_t
#include <float.h> // LDBL_MANT_DIG
#include <complex.h> // complex numbers

#include <gmp.h>
#include <mpfr.h>
#include <mpc.h>
#include "type.h"
#include "larc.h"
#include "hash.h"
#include "global.h"
#include "matrix_store.h"

/*!
 * \file spr.c
 * \brief contains routines for the Single-tile Probabilistic Retrieval technique
 */

#ifndef MAR

#ifndef SWIG

/*!
 * \ingroup larc
 * \brief The function that finds the SPR region center
 * \param output The SPR region center associated to the input value
 * \param input The input value
 */
void return_SPR_region_center(scalarType *output, const scalarType input);

/*!
 * \brief Returns a PTR to a matrixRecord for a scalar that either is equal to the target_scalar, or is a previously stored scalar that LARC deems to be sufficiently close.
 *
 * This routine first looks for the SPR tile associated with a scalar in the
 * matrix_store. If such a tile is found, the pointer to the representative
 * scalar is returned. If no such tile is found, then a new scalar record is
 * formed, and a pointer to the new record is returned.
 * 
 * \param target_scalar The input scalar.
 * \return A record pointer for a scalar
 */
mats_ptr_t get_PTR_scalar_record(scalarType scalar);

/*!
 * \ingroup larc
 * \brief Clean scalar store for SPR mode.
 *
 *  Rremoves all eligible matrices from the scalar store.
 *  WARNING: LARC does not toggle counters to show that a scalar is
 *  needed for norms or traces, so these scalars will be cleaned
 *  unless locked, held, or used by some 2 by 2 matrix 
 *
 * \param matrix_store_ptr A pointer to the static matrix store object.
 * \result Returns 1.
 */
int clean_scalar_matrices_SPR(struct matrix_store_t * matrix_store_ptr);

#endif // #ifndef SWIG

#endif // MAR not defined 

#endif  // LARC_SPR_H



