//                         mar.h
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

#ifndef LARC_MAR_H
#define LARC_MAR_H

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
 * \file mar.c
 * \brief contains routines for the Multi-tile Assured Retrieval technique
 */

#ifdef MAR

#ifndef SWIG
extern MAR_tile_index_t *GLOBAL_primary_tile_index_PTR;
extern char Gpti_in_use;
#endif // #ifndef SWIG

void scratchVars_mar_init(void);
void scratchVars_mar_clear(void);

/*!
 * \ingroup larc
 * \brief compares two tile indices to see if they are equal
 * \param tile1_index_PTR Pointer to the first tile index
 * \param tile2_index_PTR Pointer to the second tile index
 * \return 1 if the tiles contain the same index values
 */
int tile_indices_equal(MAR_tile_index_t *tile1_index_PTR, MAR_tile_index_t *tile2_index_PTR);

/*!
 * \ingroup larc
 * \brief checks tile indices to see if primary + offset = target
 * \param primary_index_PTR Pointer to a primary tile index
 * \param target_index_PTR Pointer to a possible secondary tile index
 * \param tile_offset_flag A specified offset for the target tile
 * \return 1 if primary tile index + tile_offset is the target tile index, 0 otherwise
 */
int tile_indices_equal_after_offset(MAR_tile_index_t *primary_index_PTR,
                                    MAR_tile_index_t *target_index_PTR,
                                    unsigned int tile_offset_flag);

/*!
 * \ingroup larc
 * \brief Prints out the tile index
 *
 * The tile index is one multiprecision integer for non-complex types, and
 * two multiprecision integers for complex types.
 *
 * \param tile_index_PTR The pointer to the index to be printed
 */
void print_tile_index(MAR_tile_index_t *tile_index_PTR);

/*!
 * \ingroup larc
 * \brief Uses two (single-precision) integers to set a tile index
 * \param tile_index_PTR The pointer to the tile index to be set
 * \param i The real_index (or for non-complex types the only index)
 * \param j The imag_index (ignored for non-complex types)
 */
void set_tile_index(MAR_tile_index_t *tile_index_PTR, int i, int j);

/*!
 * \ingroup larc
 * \brief Returns a PTR to a MAR_tile_index structure for the tile that contains this scalar.
 * \param tile_index_PTR The pointer to the structure containing the tile index.
 * \param tile_exp The integer determining tile widths of 1/2^tile_exp
 * \param scalar The scalarType input to the function.
 */
void get_tile_index(MAR_tile_index_t *tile_index_PTR,
		    int tile_exp, const scalarType scalar);

#ifndef SWIG


/*!
 * \ingroup larc
 * \brief Returns a PTR to a matrixRecord for a scalar that either is equal to the target_scalar, or is a previously stored scalar that LARC deems to be sufficiently close.
 *
 * This routine first looks for the MAR tile associated with a scalar in the matrix_store.
 * If such a tile is found, it is either a primary record (containing a scalar in the same tile)
 * or a neighbor record (which was claimed by an adjacent primary tile, to create a
 * multi-tile MAR region), and matrix PTR to the primary record is returned.
 * If no such tile is found, then a new primary matrix record is formed (whose matrix
 * record PTR will be returned).   Before returning the record PTR an attempt is made
 * to expand the MAR region to include previously unclaimed neighbor tiles which
 *  are adjacent to the quarter of the primary tile that the target scalar sits in 
 * (or in the half of the primary tile for real types).
 * 
 * \param target_scalar The input scalar.
 * \return A record pointer for a scalar
 */
mats_ptr_t retrieve_PTR_scalar_record(scalarType target_scalar);

/*!
 * \ingroup larc
 * \brief Clean scalar store for MAR mode.
 *
 *  Removes all eligible matrices from the scalar store.
 *  WARNING: LARC does not toggle counters to show that a scalar is
 *  needed for norms or traces, so these scalars will be cleaned
 *  unless locked, held, or used by some 2 by 2 matrix 
 *
 * \param matrix_store_ptr A pointer to the static matrix_store object.
 * \result Returns 1.
 */
int clean_scalar_matrices_MAR(struct matrix_store_t * matrix_store_ptr);


#endif // not SWIG


#endif // MAR defined 


#endif  // LARC_MAR_H
