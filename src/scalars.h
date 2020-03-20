//                       scalars.h
/******************************************************************
 *                                                                *
 * Copyright (C) 2014, Institute for Defense Analyses             *
 * 4850 Mark Center Drive, Alexandria, VA; 703-845-2500           *
 * This material may be reproduced by or for the US Government    *
 * pursuant to the copyright license under the clauses at DFARS   *
 * 252.227-7013 and 252.227-7014.                                 *
 *                                                                *
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

#ifndef SCALARS_H
#define SCALARS_H

#ifdef __cplusplus
extern "C" {
#endif

// Standard Libaries
#include <stdint.h>

#include "larc.h"

// See larc.h for prototypes of scalar operations

// for the user to set their own functions
/*!
 * \ingroup larc
 * \brief override the default function for allocating a new scalarType
 * \param func The new function
 */
void define_sca_init(void func(scalarType *));
/*!
 * \ingroup larc
 * \brief override the default function for deallocating a scalarType
 * \param func The new function
 */
void define_sca_clear(void func(scalarType *));
/*!
 * \ingroup larc
 * \brief override the default function for setting a scalarType
 * \param func The new function
 */
void define_sca_set(void func(scalarType *, const scalarType));
/*!
 * \ingroup larc
 * \brief override the default function for setting a scalarType using an input string
 * \param func The new function
 */
void define_sca_set_str(void func(scalarType *, const char *));
/*!
 * \ingroup larc
 * \brief override the default function for setting a scalarType using two input double precision numbers
 * \param func The new function
 */
void define_sca_set_2ldoubles(void func(scalarType *, long double, long double));
/*!
 * \ingroup larc
 * \brief override the default function for expressing a scalarType in string format
 * \param func The new function
 */
void define_sca_get_str(char* func(const scalarType));
/*!
 * \ingroup larc
 * \brief override the default function for hashing a scalarType
 * \param func The new function
 */
void define_sca_hash(uint64_t func(const scalarType, uint64_t));
/*!
 * \ingroup larc
 * \brief override the default function for adding two scalarType variables
 * \param func The new function
 */
void define_sca_add(void func(scalarType *, const scalarType, const scalarType));
/*!
 * \ingroup larc
 * \brief override the default function for multiplying two scalarType variables
 * \param func The new function
 */
void define_sca_mult(void func(scalarType *, const scalarType, const scalarType));
/*!
 * \ingroup larc
 * \brief override the default function for dividing two scalarType variables
 * \param func The new function
 */
void define_sca_divide(void func(scalarType *, const scalarType, const scalarType));
/*!
 * \ingroup larc
 * \brief override the default function for comparing two scalarType variables
 * \param func The new function
 */
void define_sca_cmp(int func(const scalarType, const scalarType));
/*!
 * \ingroup larc
 * \brief override the default function for determining whether two scalarType variables are equal
 * \param func The new function
 */
void define_sca_eq(int func(const scalarType, const scalarType));
/*!
 * \ingroup larc
 * \brief override the default function for determining whether two scalarType variables have the same neighborhood approximation
 * \param func The new function
 */
void define_sca_eq_approx(int func(const scalarType, const scalarType));
/*!
 * \ingroup larc
 * \brief override the default function for the sqrt of a scalarType 
 * \param func The new function
 */
void define_sca_sqrt(void func(scalarType *, const scalarType));
/*!
 * \ingroup larc
 * \brief override the default function for the norm of a scalarType 
 * \param func The new function
 */
void define_sca_norm(void func(scalarType *, const scalarType));
/*!
 * \ingroup larc
 * \brief override the default function for conjugating a scalarType 
 * \param func The new function
 */
void define_sca_conj(void func(scalarType *, const scalarType));
/*!
 * \ingroup larc
 * \brief override the default function for determining whether a scalarType has no imaginary part
 * \param func The new function
 */
void define_sca_is_real(int func(const scalarType));


//NOTE: while it is currently possible to switch between these if (for a fixed
//scalarType, it is not recommended. In particular, it is difficult to
//predict how memoizing could unexpectedly cause errors. 
//E.g. if scalarType = INTEGER and we go from arithmetic to boolean. Then 
//reading in a new matrix might NOT produce a boolean matrix. 

/*!
 * \ingroup larc
 * \brief Tells LARC to use Boolean versions of the usual arithmetic operations
 * \param verbose an enum value in {SILENT, BASIC, CHATTY, DEBUG}; determines which messsages are printed
 */
void init_boolean_scalarOps(int verbose);
/*!
 * \ingroup larc
 * \brief Tells LARC to use the usual arithmetic operations
 * \param verbose an enum value in {SILENT, BASIC, CHATTY, DEBUG}; determines which messsages are printed
 */
void init_arithmetic_scalarOps(int verbose);

/*!
 * \ingroup larc
 *
 * \brief Check to see if scalar operations have already been set.
 *
 * \return  0 if no scalar ops have been set. 
 *          1 if all scalar ops have been set. 
 *         -1 if some, but not all, scalar ops have been set. 
 */
int check_scalarOps();

#ifdef __cplusplus
}
#endif

#endif
