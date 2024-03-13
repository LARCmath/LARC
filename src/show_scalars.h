//                       show_scalars.h
/******************************************************************
 *                                                                *
 * Copyright (C) 2014-2024, Institute for Defense Analyses        *
 * 4850 Mark Center Drive, Alexandria, VA; 703-845-2500           *
 * This material may be reproduced by or for the US Government    *
 * pursuant to the copyright license under the clauses at DFARS   *
 * 252.227-7013 and 252.227-7014.                                 *
 *                                                                *
 * LARC (Linear Algebra via Recursive Compression)                *
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

#ifndef SHOW_SCALARS_H
#define SHOW_SCALARS_H

#include <inttypes.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpc.h>
#include "larc.h"


/*!
 * \ingroup larc
 *
 * \brief Return the versions of GMP, MPFR, and MPC.
 * \return A string giving the current versions of the GMP software packages.
 */
const char* get_string_gmp_versions (void);

/*!
 * \ingroup larc
 *
 * \brief Print the versions of GMP, MPFR, and MPC to standard output.
 */
void print_gmp_versions (void);

/*!
 * \ingroup larc
 *
 * \brief Print the full internal details of a GNU multi-precision
 * integer structure to standard output.
 *
 * This prints out the address of the base structure, the address and
 * amount of any heap memory allocated, the binary data stored, and
 * the value of the represented multi-precision integer.
 * 
 * \param v     The multi-precision integer structure to debug.
 * \param title A title string to add to the debugging output.
 */
void debug_out_mpz (const mpz_t v, const char * title);

/*!
 * \ingroup larc
 *
 * \brief Print the full internal details of a GNU multi-precision
 * rational structure to standard output.
 *
 * This prints out the address of the base structure, the address and
 * amount of any heap memory allocated, the binary data stored, and
 * the value of the represented multi-precision rational number.
 * 
 * \param v     The multi-precision rational structure to debug.
 * \param title A title string to add to the debugging output.
 */
void debug_out_mpq (const mpq_t v, const char * title);

#ifdef USE_MPRATCOMPLEX
/*!
 * \ingroup larc
 *
 * \brief Print the full internal details of a LARC multi-precision
 * rational complex structure to standard output.
 *
 * This prints out the address of the base structure, the address and
 * amount of any heap memory allocated, the binary data stored, and
 * the value of the represented multi-precision rational complex number.
 * 
 * \param v     The multi-precision rational complex structure to debug.
 * \param title A title string to add to the debugging output.
 */
void debug_out_RatComplex (const larc_mpratcomplex_t v, const char * title);
#endif // #ifdef USE_MPRATCOMPLEX

#ifdef USE_CLIFFORD
/*!
 * \ingroup larc
 *
 * \brief Print the full internal details of a LARC Clifford Algebra
 * structure to standard output.
 *
 * This prints out the address of the base structure, the address and
 * amount of any heap memory allocated, the binary data stored, and
 * the value of the represented multi-precision Clifford Algebra value.
 * 
 * \param v     The multi-precision Clifford Algebra structure to debug.
 * \param title A title string to add to the debugging output.
 */
void debug_out_Clifford (const clifford_t v, const char * title);
#endif // #ifdef USE_CLIFFORD

/*!
 * \ingroup larc
 *
 * \brief Print the full internal details of a GNU multi-precision
 * real structure to standard output.
 *
 * This prints out the address of the base structure, the address and
 * amount of any heap memory allocated, the binary data stored, and
 * the value of the represented multi-precision real number.
 * 
 * \param v     The multi-precision real structure to debug.
 * \param title A title string to add to the debugging output.
 */
void debug_out_mpfr (const mpfr_t v, const char * title);

/*!
 * \ingroup larc
 *
 * \brief Print the full internal details of a GNU multi-precision
 * complex structure to standard output.
 *
 * This prints out the address of the base structure, the address and
 * amount of any heap memory allocated, the binary data stored, and
 * the value of the represented multi-precision complex number.
 * 
 * \param v     The multi-precision complex structure to debug.
 * \param title A title string to add to the debugging output.
 */
void debug_out_mpc (const mpc_t v, const char * title);

/*!
 * \ingroup larc
 *
 * \brief Print the internal representation of the 64-bit integer
 * structure to standard output.
 *
 * \param v     The 64-bit integer structure to debug.
 * \param title A title string to add to the debugging output.
 */
void debug_out_integer (const int64_t v, const char * title);

/*!
 * \ingroup larc
 *
 * \brief Print the internal representation of the C long double
 * structure to standard output.
 *
 * \param v     The C long double structure to debug.
 * \param title A title string to add to the debugging output.
 */
void debug_out_real (const long double v, const char * title);

/*!
 * \ingroup larc
 *
 * \brief Print the internal representation of the C long double
 * complex structure to standard output.
 *
 * \param v     The C long double complex structure to debug.
 * \param title A title string to add to the debugging output.
 */
void debug_out_complex (const long double complex v, const char * title);

/*!
 * \ingroup larc
 *
 * \brief Print the full internal details of whichever ScalarType
 * LARC is currently using to standard output.
 *
 * \param v     The LARC scalar structure to debug.
 * \param title A title string to add to the debugging output.
 */
void debug_out_scalar (const scalarType v, const char * title);


#endif  // #ifndef SHOW_SCALARS_H

