//                        clifford.c
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


// Standard Libaries
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <pthread.h>
#include <inttypes.h> //for printing
#include "gmp.h"
#include "mpfr.h"

// Our header files structures and functions
#include "larc.h"
#include "matrix_store.h"
#include "op_store.h"
#include "global.h"
#include "organize.h"
#include "version.h"
#include "info_store.h"
#include "scalars.h"
#include "hash.h"
#include "mar.h"
#ifdef USE_CLIFFORD
#include "clifford.h"
#endif // #ifdef USE_CLIFFORD

/*!
 * \file clifford.c
 * \brief This file contains the routines that handle Clifford algebra
 * scalar types.
 *
 * In particular, it contains routines which translate string values into
 * scalarType and vice versa, and functions which wrap basic arithmetic
 * operations and value comparisons for scalarType so that the underlying
 * Clifford algebra scalar type is handled correctly.
 *
 */

#ifdef USE_CLIFFORD

int io_string_base_dup = 10;

/*
 * For each implemented Clifford algebra scalar type, we need to specify
 * the corresponding data.
 */

#if defined(MPC_CLIFFORD_S2) || defined(MPR_CLIFFORD_S2)

/*
 * This section is only needed for the Clifford algebra MPC_CLIFFORD_S2 = Q(i, sqrt(2)).
 * It is not needed for the other Clifford algebras.
 */

#if defined(MPC_CLIFFORD_S2)
char *clifford_description = "Q(i, sqrt(2))";
#else
char *clifford_description = "Q(sqrt(2))";
#endif

/*!
 * 1  = square root of 1 (sometimes S1 in comments)
 * S2 = square root of 2
 */
char *clifford_consts[CLIFFORD_DIMENSION] = {"1", "S2"};


#elif defined(MPC_CLIFFORD_S2_S3) || defined(MPR_CLIFFORD_S2_S3)

/*
 * This section is only needed for the Clifford algebra MPC_CLIFFORD_S2_S3 = Q(i, sqrt(2), sqrt(3)).
 * It is not needed for the other Clifford algebras.
 */

#if defined(MPC_CLIFFORD_S2_S3)
char *clifford_description = "Q(i, sqrt(2), sqrt(3))";
#else
char *clifford_description = "Q(sqrt(2), sqrt(3))";
#endif

/*!
 * 1  = square root of 1 (sometimes S1 in comments)
 * S2 = square root of 2
 * S3 = square root of 3
 * S6 = square root of 6
 */
char *clifford_consts[CLIFFORD_DIMENSION] = {"1", "S2", "S3", "S6"};


#elif defined(MPC_CLIFFORD_C2) || defined(MPR_CLIFFORD_C2)

/*
 * This section is only needed for the Clifford algebra MPC_CLIFFORD_C2 = Q(i, 2^(1/3)).
 * It is not needed for the other Clifford algebras.
 */

#if defined(MPC_CLIFFORD_C2)
char *clifford_description = "Q(i, 2^(1/3))";
#else
char *clifford_description = "Q(2^(1/3))";
#endif

/*!
 * 1  = cube root of 1 (sometimes C1 in comments)
 * C2 = cube root of 2
 * C4 = cube root of 4
 */
char *clifford_consts[CLIFFORD_DIMENSION] = {"1", "C2", "C4"};


#endif // #if defined(MPC_CLIFFORD_S2)


/*
 * For each Clifford algebra scalar type, we also need to
 * specify how to multiply the adjoined constants and how
 * to take the inverse of a Clifford scalar value.
 */

/*
 * For example, in MPC_CLIFFORD_S2, but just the real part,
 * a general multiplication looks like:
 *
 *   (a*sqrt(1) + b*sqrt(2)) * (c*sqrt(1) + d*sqrt(2))
 *     = (a*c + 2*b*d)*sqrt(1) + (a*d + b*c)*sqrt(2)
 *
 * This equation can be derived from the distributive law once the
 * following atomic facts are known:
 *
 *   sqrt(1) * sqrt(1) = sqrt(1)
 *   sqrt(1) * sqrt(2) = sqrt(2)
 *   sqrt(2) * sqrt(1) = sqrt(2)
 *   sqrt(2) * sqrt(2) = 2*sqrt(1)
 *
 * Note the above rules are the rules for MPC_CLIFFORD_S2: Q(i, sqrt(2)).
 * They are not needed for the other Clifford algebras.
 */

/*
 * These are the rules for MPC_CLIFFORD_S2_S3: Q(i, sqrt(2), sqrt(3)).
 * This is not needed for the other Clifford algebras.
 *
 * For MPC_CLIFFORD_S2_S3, there are 16 atomic facts:
 *
 *   sqrt(1) * sqrt(1) = sqrt(1)
 *   sqrt(1) * sqrt(2) = sqrt(2)
 *   sqrt(1) * sqrt(3) = sqrt(3)
 *   sqrt(1) * sqrt(6) = sqrt(6)
 *
 *   sqrt(2) * sqrt(1) = sqrt(2)
 *   sqrt(2) * sqrt(2) = 2*sqrt(1)
 *   sqrt(2) * sqrt(3) = sqrt(6)
 *   sqrt(2) * sqrt(6) = 2*sqrt(3)
 *
 *   sqrt(3) * sqrt(1) = sqrt(3)
 *   sqrt(3) * sqrt(2) = sqrt(6)
 *   sqrt(3) * sqrt(3) = 3*sqrt(1)
 *   sqrt(3) * sqrt(6) = 3*sqrt(2)
 *
 *   sqrt(6) * sqrt(1) = sqrt(6)
 *   sqrt(6) * sqrt(2) = 2*sqrt(3)
 *   sqrt(6) * sqrt(3) = 3*sqrt(2)
 *   sqrt(6) * sqrt(6) = 6*sqrt(1)
 *
 * These atomic facts yield the following multiplication rule:
 *
 *     (a*sqrt(1) + b*sqrt(2) + c*sqrt(3) + d*sqrt(6))
 *   * (e*sqrt(1) + f*sqrt(2) + g*sqrt(3) + h*sqrt(6))
 *
 *     = (a*e + 2*b*f + 3*c*g + 6*d*h) * sqrt(1)
 *     + (a*f +   b*e + 3*c*h + 3*d*g) * sqrt(2)
 *     + (a*g + 2*b*h +   c*e + 2*d*f) * sqrt(3)
 *     + (a*h +   b*g +   c*f +   d*e) * sqrt(6)
 */

/*
 * These are the rules for MPC_CLIFFORD_C2: Q(i, 2^(1/3)).
 * This is not needed for the other Clifford algebras.
 *
 * In MPC_CLIFFORD_C2, there are 9 atomic facts:
 *
 *   C1 * C1 = C1
 *   C1 * C2 = C2
 *   C1 * C4 = C4
 *
 *   C2 * C1 = C2
 *   C2 * C2 = C4
 *   C2 * C4 = 2*C1
 *
 *   C4 * C1 = C4
 *   C4 * C2 = 2*C1
 *   C4 * C4 = 2*C2
 */

#if defined(MPC_CLIFFORD_S2) || defined(MPR_CLIFFORD_S2)

clifford_term_t clifford_mult[CLIFFORD_DIMENSION][CLIFFORD_DIMENSION] =
  {
    {
      {1, 1, 0},     // clifford_mult[S1][S1] = S1 * S1 = (1/1)*S1
      {1, 1, 1}      // clifford_mult[S1][S2] = S1 * S2 = (1/1)*S2
    },
    {
      {1, 1, 1},     // clifford_mult[S2][S1] = S2 * S1 = (1/1)*S2
      {2, 1, 0}      // clifford_mult[S2][S2] = S2 * S2 = (2/1)*S1
    }
  };


#elif defined(MPC_CLIFFORD_S2_S3) || defined(MPR_CLIFFORD_S2_S3)

clifford_term_t clifford_mult[CLIFFORD_DIMENSION][CLIFFORD_DIMENSION] =
  {
    {
      {1, 1, 0},     // clifford_mult[S1][S1] = S1 * S1 = (1/1)*S1
      {1, 1, 1},     // clifford_mult[S1][S2] = S1 * S2 = (1/1)*S2
      {1, 1, 2},     // clifford_mult[S1][S3] = S1 * S3 = (1/1)*S3
      {1, 1, 3}      // clifford_mult[S1][S6] = S1 * S6 = (1/1)*S6
    },
    {
      {1, 1, 1},     // clifford_mult[S2][S1] = S2 * S1 = (1/1)*S2
      {2, 1, 0},     // clifford_mult[S2][S2] = S2 * S2 = (2/1)*S1
      {1, 1, 3},     // clifford_mult[S2][S3] = S2 * S3 = (1/1)*S6
      {2, 1, 2}      // clifford_mult[S2][S6] = S2 * S6 = (2/1)*S3
    },
    {
      {1, 1, 2},     // clifford_mult[S3][S1] = S3 * S1 = (1/1)*S3
      {1, 1, 3},     // clifford_mult[S3][S2] = S3 * S2 = (1/1)*S6
      {3, 1, 0},     // clifford_mult[S3][S3] = S3 * S3 = (3/1)*S1
      {3, 1, 1}      // clifford_mult[S3][S6] = S3 * S6 = (3/1)*S2
    },
    {
      {1, 1, 3},     // clifford_mult[S6][S1] = S6 * S1 = (1/1)*S6
      {2, 1, 2},     // clifford_mult[S6][S2] = S6 * S2 = (2/1)*S3
      {3, 1, 1},     // clifford_mult[S6][S3] = S6 * S3 = (3/1)*S2
      {6, 1, 0}      // clifford_mult[S6][S6] = S6 * S6 = (6/1)*S1
    }
  };


#elif defined(MPC_CLIFFORD_C2) || defined(MPR_CLIFFORD_C2)

clifford_term_t clifford_mult[CLIFFORD_DIMENSION][CLIFFORD_DIMENSION] =
  {
    {
      {1, 1, 0},     // clifford_mult[C1][C1] = C1 * C1 = (1/1)*C1
      {1, 1, 1},     // clifford_mult[C1][C2] = C1 * C2 = (1/1)*C2
      {1, 1, 2}      // clifford_mult[C1][C4] = C1 * C4 = (1/1)*C4
    },
    {
      {1, 1, 1},     // clifford_mult[C2][C1] = C2 * C1 = (1/1)*C2
      {1, 1, 2},     // clifford_mult[C2][C2] = C2 * C2 = (1/1)*C4
      {2, 1, 0}      // clifford_mult[C2][C4] = C2 * C4 = (2/1)*C1
    },
    {
      {1, 1, 2},     // clifford_mult[C4][C1] = C4 * C1 = (1/1)*C4
      {2, 1, 0},     // clifford_mult[C4][C2] = C4 * C2 = (2/1)*C1
      {2, 1, 1}      // clifford_mult[C4][C4] = C4 * C4 = (2/1)*C2
    }
  };


#endif // #if defined(MPC_CLIFFORD_S2) || defined(MPR_CLIFFORD_S2)




/*
 * Division is accomplished by taking the inverse and then multiplying.
 *
 *     1 / (a*sqrt(1) + b*sqrt(2))
 *   =
 *     ((a/(a*a-2*b*b))*sqrt(1) - (b/(a*a-2*b*b))*sqrt(2))
 */




/*
 * Other Clifford Operations Follow:
 */


/*!
 * \brief This contains the multiprecision real values
 * of each of the adjoined Clifford constants.
 */
  mpfr_t mpfr_algebraic_approx[CLIFFORD_DIMENSION-1];

/*!
 * \brief This contains the multiprecision rational values
 * of each of the adjoined Clifford constants.
 */
  mpq_t mpq_algebraic_approx[CLIFFORD_DIMENSION-1];

// The SPR region or MAR tile that contains a point expressed in terms of a
// Clifford algebra is most easily found by first calculating an MPRATIONAL
// approximation to the original point. This is done by obtaining an MPREAL
// representation of the adjoining algebraic numbers (eg sqrt2) to whatever
// precision is needed, converting these values to MPRATIONAL type, then
// performing MPRATIONAL multiplications and additions to get a value very 
// close to the desired point in MPRATIONAL form. This approximation is then
// fed into the MPRATIONAL routine that returns the proper tile or region
// label. As all the rational operations are exact, the approximation accuracy
// is limited only by the precision of the MPREAL approximation of the
// adjoining algebraic numbers.

/*!
 * 
 * \ingroup LARC
 * \brief Estimate the number of bits of precision needed for a good floating-point approximation to a number in Clifford form
 *
 * \param scalar The Clifford number to be approximated
 * \result The number of bits needed for the approximation
 */
static size_t get_clifford_algebraic_precision(const clifford_t scalar)
{
    // Given an input Clifford scalar C(a,b,c,...) representing a value
    // a*1 + b*A1 + c*A2 + ..., this function estimates the number of bits
    // needed for floating-point approximations A1f, A2f, ... to the algebraic
    // numbers A1, A2, ...  such that C(a+b*A1f+c*A2f+...,0,0,...) is a good
    // approximation to the input scalar.
    //
    // Ideally we would not use any more precision than necessary. In practice,
    // we estimate a safe (more than likely too large) number and use that.
    size_t num_bits = 0;
    for (int i=0; i<CLIFFORD_DIMENSION; ++i)
    {
        num_bits += mpz_sizeinbase(mpq_numref(scalar->real_coeffs[i]),2);
        num_bits += mpz_sizeinbase(mpq_denref(scalar->real_coeffs[i]),2);
    }
#ifdef IS_COMPLEX
    size_t num_bits_i = 0;
    for (int i=0; i<CLIFFORD_DIMENSION; ++i)
    {
        num_bits_i += mpz_sizeinbase(mpq_numref(scalar->imag_coeffs[i]),2);
        num_bits_i += mpz_sizeinbase(mpq_denref(scalar->imag_coeffs[i]),2);
    }
    num_bits = MAX(num_bits,num_bits_i);
#endif
    // this should ensure several digits of precision in our final answers
    num_bits += 20;
    // always use as a minimum the LARC standard-sized mpfr_t
    if (num_bits<256) num_bits = 256;
    return num_bits;
}

// not sure where the following text should go...
//
    // The algebraic approximation global arrays have dimension one less than
    // CLIFFORD_DIMENSION (since we don't need an approximation for 1 or i).
    // Their names are mpfr_algebraic_approx[] and mpq_algebraic_approx[].
    //
    // As we find new and larger worst cases, we may need to allocate more
    // memory for the MPFR variables. We use mpfr_set_prec( ) to clear
    // old values and only allocate more memory when necessary (it does not
    // reduce the size allocated even when possible). Using the smallest
    // precision known to give correct results is what we are striving for.
    // We do not need to do anything special with the MPQ variables, which 
    // automatically adjust their memory.
    //
    // Since we only have a few of these variables, it is not a problem that
    // more memory may be allocated than needed at a particular time. Our tests
    // used ~2million bit numbers to add two numbers of O(10^{401 370}) and get
    // a result of O(10^{-401 370}), and it's unlikely that we'll be dealing
    // with a case as bad as this in real applications.
    //
    // We will need to recompute the values to the desired precision, but that
    // doesn't seem to be too much of a time sink. This should be tested.
    // An alternative is to round the high-precision value to lower precisions,
    // requiring a seperate low-precision variable and whatever time that
    // takes. I'm guessing recomputation is better.
//

/*!
 *
 * \ingroup LARC
 * \brief Gets a 'good' rational approximation to the algebraic terms in a Clifford field
 *
 * The approximation is obtained by first getting multiprecision floating point
 * numbers with nbits_scalar bits of precision, using these to approximate the
 * algebraic terms, then exactly converting these approximations to
 * multiprecision rational numbers
 *
 * \param nbits_scalar The number of bits needed for the floating-point conversions
 *
 */
static void set_values_for_clifford_rational_approx(size_t nbits_scalar)
{
    // This function puts approximate values for the algebraic numbers in
    // our Clifford algebra into multiprecision floating point format, then
    // converts these floats into exact rational numbers. The accuracy of the
    // approximation is thus determined by the number of bits of precision
    // used for the floats. We use global variables to hold these values to
    // minimize repeated allocating/deallocating.

    // Ensure there is enough memory for the desired floating point precision
    // variables (also sets value to NaN). There is no reallocation performed
    // unless the new precision requires more memory than currently allocated.
    for (int i=0;i<CLIFFORD_DIMENSION-1;++i)
        mpfr_set_prec(mpfr_algebraic_approx[i],nbits_scalar);

    // Calculate floating-point approximations for the algebraic numbers which
    // are used in the field. Calculations are done to nbits_scalar precision.
#if defined(MPC_CLIFFORD_S2) || defined(MPR_CLIFFORD_S2)
    // the algebraic number is sqrt(2)
    mpfr_sqrt_ui(mpfr_algebraic_approx[0],2,MPFR_RNDN);
#elif defined(MPC_CLIFFORD_S2_S3) || defined(MPR_CLIFFORD_S2_S3)
    // the algebraic numbers are sqrt(2) and sqrt(3) 
    // we also need their product sqrt(6)
    mpfr_sqrt_ui(mpfr_algebraic_approx[0],2,MPFR_RNDN);
    mpfr_sqrt_ui(mpfr_algebraic_approx[1],3,MPFR_RNDN);
    mpfr_sqrt_ui(mpfr_algebraic_approx[2],6,MPFR_RNDN);
#elif defined(MPC_CLIFFORD_C2) || defined(MPR_CLIFFORD_C2)
    // the algebraic number is cube_root(2)
    // we also need its square, cube_root(4)
    if (scratchVars.mpreal_in_use)
        fprintf(stderr,"%s reusing scratchVars.mpreal!\n",__func__);
    scratchVars.mpreal_in_use = 1;

    mpfr_t *x = &scratchVars.mpreal; // this precision does not need to change
    mpfr_set_ui(*x,2,MPFR_RNDN);
    mpfr_rootn_ui(mpfr_algebraic_approx[0],*x,3,MPFR_RNDN);
    mpfr_set_ui(*x,4,MPFR_RNDN);
    mpfr_rootn_ui(mpfr_algebraic_approx[1],*x,3,MPFR_RNDN);
#endif

    // Convert these floats into rationals (this conversion is exact).
    // GMP automatically handles any memory reallocations
    for (int i=0;i<CLIFFORD_DIMENSION-1;++i)
        mpfr_get_q(mpq_algebraic_approx[i],mpfr_algebraic_approx[i]);

    scratchVars.mpreal_in_use = 0;
    return; 
}

void convert_clifford_to_mprational_approx(clifford_t approx, const clifford_t scalar)
{
    // This function converts a Clifford number C(a,b,c,...) to an approximate
    // value C(d,0,0,...); in other words, the first rational number in the
    // approximation is nearly equal to the exact value C(a,b,c,...). For 
    // complex numbers, our representation C(a,b,c,...)+sqrt(-1)*C(d,e,f,...)
    // is approximated as C(g,0,0,...)+sqrt(-1)*C(h,0,0,....).

    // the mprational approximation is stored in the first positions of the
    // clifford structure (real and imaginary calculated separately)
    if (scratchVars.mprational_in_use)
        fprintf(stderr,"%s reusing scratchVars.mprational!\n",__func__);
    scratchVars.mprational_in_use = 1;

    mpq_t *temp_q = &scratchVars.mprational;

    // In the case where we have C(a,b) == a + b*A1 with a, b large in
    // magnitude and C(a,b) small in magnitude, Bad Things can happen when the
    // algebraic number A1 is approximated to insufficient precision. The
    // problem generalizes to more complex Clifford algebras. We use the
    // following function to determine an adequate precision for the MPFR
    // approximations to the algebraic numbers.
    size_t nbits_scalar = get_clifford_algebraic_precision(scalar);

    // call another routine which knows what the algebraic numbers are and
    // how to calculate them; the routine sets global variables used below
    set_values_for_clifford_rational_approx(nbits_scalar);

    // set the coefficient of the unit position of the approximate number
    // equal to the coefficient of the unit position of the input.
    mpq_set(approx->real_coeffs[0],scalar->real_coeffs[0]);

    // the product of the approximation to an algebraic number and its
    // coefficient are summed into the unit position coefficient. To be safe,
    // the algebraic position coefficients of the approximate value are zeroed.
    for (int i=1;i<CLIFFORD_DIMENSION;++i)
    {
        mpq_mul(*temp_q,scalar->real_coeffs[i],mpq_algebraic_approx[i-1]);
        mpq_add(approx->real_coeffs[0],approx->real_coeffs[0],*temp_q);
        mpq_set_ui(approx->real_coeffs[i],0,1);
    }
#ifdef IS_COMPLEX
    mpq_set(approx->imag_coeffs[0],scalar->imag_coeffs[0]);
    for (int i=1;i<CLIFFORD_DIMENSION;++i)
    {
        mpq_mul(*temp_q,scalar->imag_coeffs[i],mpq_algebraic_approx[i-1]);
        mpq_add(approx->imag_coeffs[0],approx->imag_coeffs[0],*temp_q);
        mpq_set_ui(approx->imag_coeffs[i],0,1);
    }
#endif // IS_COMPLEX
    scratchVars.mprational_in_use = 0;
    return;
}

int clifford_test_pure_rational(const clifford_t val){
    for (int i=1; i<CLIFFORD_DIMENSION; ++i) {
        if (0 != mpq_cmp_ui(val->real_coeffs[i], 0, 1)) return 0;
#ifdef IS_COMPLEX
        if (0 != mpq_cmp_ui(val->imag_coeffs[i], 0, 1)) return 0;
#endif // IS_COMPLEX
    }
    return 1;
}


int lookup_clifford_name(const char *s) {
    if (s[0] == '\0') {
        return 0;
    }

    for (int i = 0; i < CLIFFORD_DIMENSION; ++i) {
        if (0 == strcmp(clifford_consts[i], s)) {
            return i;
        }
    }

    return -1;
}

#endif  // #ifdef USE_CLIFFORD
