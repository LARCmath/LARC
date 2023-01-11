//                       clifford.h
/******************************************************************
 *                                                                *
 * Copyright (C) 2014, Institute for Defense Analyses             *
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

#ifndef CLIFFORD_H
#define CLIFFORD_H

#ifdef USE_CLIFFORD

#ifdef __cplusplus
extern "C" {
#endif // #ifdef __cplusplus

// Standard Libaries
#include <stdint.h>
#include "larc.h"
#include "gmp.h"

/*!
 * \file clifford.h
 * \brief This is the header file for Clifford algebra scalarTypes.
 *
 * See the clifford.c file for more information about Clifford algebras.
 */


/*!
 * \brief This is a short description of the particular
 * Clifford algebra that has been selected.
 */
  extern char *clifford_description;


/*!
 * \brief This specifies the names for each of the
 * adjoined constant values for the Clifford algebra.
 */
  extern int io_string_base_dup;
  extern char *clifford_consts[CLIFFORD_DIMENSION];

  typedef struct clifford_term_s
  {
     int numerator;
     int denominator;
     int const_index;
  } clifford_term_t;
  extern clifford_term_t clifford_mult[CLIFFORD_DIMENSION][CLIFFORD_DIMENSION];

  extern mpfr_t mpfr_algebraic_approx[CLIFFORD_DIMENSION-1];
  extern mpq_t mpq_algebraic_approx[CLIFFORD_DIMENSION-1];

#ifndef SWIG

/*!
 * \brief This defines larc scalar type for the Clifford algebra.
 *
 * For example, if there are two adjoined constant values S1 and S2,
 * then a Clifford algebra scalar value is given by four rational
 * coefficients, two for the real part and two for the imaginary part.
 * If the real coefficients are "a" and "b", and the imaginary
 * coefficients are "c" and "d", then the value represented is:
 *
 *    (a*S1 + b*S2) + I*(c*S1 + d*S2)
 */
  typedef struct clifford_s
  {
      mpq_t real_coeffs[CLIFFORD_DIMENSION];
#ifdef IS_COMPLEX
      mpq_t imag_coeffs[CLIFFORD_DIMENSION];
#endif // IS_COMPLEX
  } clifford_t[1];

/* Operation declarations */

/*!
 * \ingroup LARC
 * \brief checks a Clifford value to see if it could be represented as a rational number a/b (or for complex numbers, is in the field Q[\sqrt{-1}], e.g. a/b + \sqrt{-1}c/d)
 * \param val the Clifford value to be tested
 * \return 0 if the number is irrational, 1 if it is rational
 */
int clifford_test_pure_rational(const clifford_t val);

int lookup_clifford_name(const char *s);

/*!
 * \ingroup LARC
 * \brief Gets a 'good' rational approximation for a Clifford number.
 * \param approx An approximation for the number in Clifford format, but with coefficients for all algebraic terms set to zero
 * \param scalar The number in Clifford format
 *
 */
void convert_clifford_to_mprational_approx(clifford_t approx, const clifford_t scalar);

#endif // #ifndef SWIG

#ifdef __cplusplus
}
#endif // #ifdef __cplusplus

#endif // #ifdef USE_CLIFFORD

#endif   // #define CLIFFORD_H
