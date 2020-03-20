//                        scalars.c
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



// See larc.h for prototypes of scalar operations.

// To add a new scalar type, modify the following files:
//    * Makefile
//        - Add new scalar type to build and unittests.
//    * src/typeNEWSCALARTYPENAME.h
//        - Create this file by copy & paste & edit of typeREAL.h.
//    * src/larc.h
//        - Add the type definition of the new scalar type.
//    * src/organize.c
//        - Add the new scalar type to the list of scalar types.
//    * src/scalars.c
//        - Implement the operations for the new scalar type.
//    * src/global.h
//        - Specify which global constants exist in the new scalar type.
//    * src/global.c
//        - Initialize the global constants for the new scalar type.
//    * src/io.c
//        - Add the new scalar type to the list of scalar types that
//          don't have a legacy format.
//    * src/exampleLARC.c
//        - Add example code for the new scalar type.


// NOTE: The code for the various scalar types is scattered below,
// and contains a lot of copy & paste of similar code fragments.
// We have tried to re-organize it several times, but every attempt
// ends up not compiling, usually due to an issue with SWIG.



/*****************************************
 * Required Constants                    *
 ****************************************/
#if defined(USE_REAL) || defined(USE_COMPLEX)
#define LDBL_MANT_DEC_DIG ((int) lround(log10(pow(2, (double) LDBL_MANT_DIG))))
#endif

#if defined(USE_MPINTEGER) || defined(USE_MPRATIONAL) || defined(USE_MPRATCOMPLEX) || defined(USE_MPREAL) || defined(USE_MPCOMPLEX)
int io_string_base = 10;
#endif

#if defined(USE_MPREAL) || defined(USE_MPCOMPLEX)
mpfr_prec_t mpreal_precision = 256;
#endif



/*****************************************
 *  Locality Approximation Functions     *
 ****************************************/
// Each scalarType needs its own version of the two locality functions
// collapse_near_zero_<scalarType> and round_sig_fig_<scalarType>. In
// addition, the function larc_nbhd_approx() is defined to use
// the correct version of the two functions. This is all handled by
// using #ifdef USE_<scalarType> so that the preprocessor sees only the
// needed functions.


#if defined(USE_MPRATIONAL) || defined(USE_MPRATCOMPLEX)
/*!
 * \ingroup larc
 * \brief Compares a multiprecision rational input to a threshold value; returns 0 if it is less than the threshold, returns the input if it is not
 * \param rsc The multiprecision rational result of the function
 * \param sc The multiprecision rational input
 */
void collapse_near_zero_mprational(mpq_t *rsc, const mpq_t sc)
{
  long double zerorealthresh = get_zerorealthresh();
  mpq_t zeromprationalthresh;
  mpq_init(zeromprationalthresh);
  mpq_set_d(zeromprationalthresh, zerorealthresh);

  mpq_t abs_sc;
  mpq_init(abs_sc);
  mpq_abs(abs_sc, sc);

  if (mpq_cmp(abs_sc, zeromprationalthresh) < 0) {
    mpq_set_ui(*rsc, 0, 1);
  } else {
    mpq_set(*rsc, sc);
  }
  mpq_clear(abs_sc);
  mpq_clear(zeromprationalthresh);
}

/*!
 * \ingroup larc
 * \brief Returns the multiprecision rational input rounded to the SIGHASH most significant bits
 * \param rsc The multiprecision rational result of the function
 * \param sc The multiprecision rational input
 */
void round_sig_fig_mprational(mpq_t *rsc, const mpq_t sc)
{
    int sighash = get_sighash();

    // Multiply by 2^sighash:
    mpq_t *temp_q = &scratchVars.mprational;
    mpq_mul_2exp(*temp_q, sc, sighash);  // shifts rational left by sighash bits
    mpq_set(*rsc,*temp_q);

    // Change from Mark's original code: add \pm0.5, to change subsequent
    // truncation into round towards 0
    mpq_set_ui(*temp_q, 1, 2);
    if (mpq_cmp_ui(*rsc, 0, 1) >= 0) {
        mpq_add(*rsc, *rsc, *temp_q);
    } else {
        mpq_sub(*rsc, *rsc, *temp_q);
    }

    // Truncate *rsc:
    mpz_set_q(scratchVars.mpinteger, *rsc);   // takes rational, sets integer
    mpq_set_z(*rsc, scratchVars.mpinteger);   // takes integer, sets rational

    // Divide by 2^sighash:
    mpq_div_2exp(*temp_q, *rsc, sighash);  // shifts rational right by sighash bits
    mpq_set(*rsc,*temp_q);
}
#endif

#if defined(USE_MPREAL) || defined(USE_MPCOMPLEX)
/*!
 * \ingroup larc
 * \brief Compares the multiprecision real input to a threshold value; returns 0 if it is less than the threshold, returns the input if it is not
 * \param rsc The multiprecision real result of the function
 * \param sc The multiprecision real input
 */
void collapse_near_zero_mpreal(mpfr_t *rsc, const mpfr_t sc)
{
  if (VERBOSE>DEBUG) printf("in %s\n",__func__);
  long double zerorealthresh = get_zerorealthresh();
  mpfr_t zeromprealthresh;
  mpfr_init2(zeromprealthresh, mpreal_precision);
  mpfr_set_ld(zeromprealthresh, zerorealthresh, MPFR_RNDN);

  mpfr_t abs_sc;
  mpfr_init2(abs_sc, mpreal_precision);
  mpfr_abs(abs_sc, sc, MPFR_RNDN);

#ifdef DEBUG
  char *larc_sca_get_str_mpreal(const mpfr_t);
  char *temp1 = larc_sca_get_str_mpreal(zeromprealthresh);
  char *temp2 = larc_sca_get_str_mpreal(abs_sc);
  printf("threshold is %s\n",temp1);
  printf("abs(scalar) is %s\n",temp2);
  free(temp1);
  free(temp2);
#endif

  if (mpfr_cmp(abs_sc, zeromprealthresh) < 0) {
    if (VERBOSE>DEBUG) printf("\tsetting scalar to zero\n\n");
    mpfr_set_ui(*rsc, 0, MPFR_RNDN);
  } else {
    if (VERBOSE>DEBUG) printf("\tnot setting scalar to zero\n\n");
    mpfr_set(*rsc, sc, MPFR_RNDN);
  }
  mpfr_clear(abs_sc);
  mpfr_clear(zeromprealthresh);
  if (VERBOSE>DEBUG) printf("exiting %s\n",__func__);
#undef DEBUG
}

/*!
 * \ingroup larc
 * \brief Returns the multiprecision input rounded to the SIGHASH most significant bits
 * \param rsc The multiprecision real result of the function
 * \param sc The multiprecision real input
 */
void round_sig_fig_mpreal(mpfr_t *rsc, const mpfr_t sc)
{
    int sighash = get_sighash();

    // Multiply by 2^sighash:
    mpfr_t *temp_q = &scratchVars.mpreal;
    mpfr_mul_2exp(*temp_q, sc, sighash, MPFR_RNDN);  // shifts real left by sighash bits
    mpfr_set(*rsc, *temp_q, MPFR_RNDN);

    // Change from Mark's original code: add \pm0.5, to change subsequent
    // truncation into round towards 0
    mpfr_set_d(*temp_q, 0.5, MPFR_RNDN);
    if (mpfr_cmp_ui(*rsc, 0) >= 0) {
        mpfr_add(*rsc, *rsc, *temp_q, MPFR_RNDN);
    } else {
        mpfr_sub(*rsc, *rsc, *temp_q, MPFR_RNDN);
    }

    // Truncate *rsc:
    mpfr_get_z(scratchVars.mpinteger, *rsc, MPFR_RNDN);   // takes real, sets integer
    mpfr_set_z(*rsc, scratchVars.mpinteger, MPFR_RNDN);   // takes integer, sets real

    // Divide by 2^sighash:
    mpfr_div_2exp(*temp_q, *rsc, sighash, MPFR_RNDN);  // shifts real right by sighash bits
    mpfr_set(*rsc, *temp_q, MPFR_RNDN);
}
#endif

#if defined(USE_REAL) || defined(USE_COMPLEX)
/*!
 * \ingroup larc
 * \brief Returns the double precision input rounded to the SIGHASH most significant bits
 * \param output The double precision floating point result of the function
 * \param input The double precision floating point input
 */
void round_sig_fig_real(long double *output, const long double input)
{
// This uses standard C functions (frexp, modf, and ldexp) to perform rounding
// in way that is independent of the specific floating point representation
// that is used by the underlying machine architecture.
  int sighash = get_sighash();

  if (isnanl(input) ) {
    fprintf(stderr,"Inside %s a NaN showed up\n", __func__);
    exit (1);
  }

  // Shift the input left by the number of significant bits to keep
  //   using the ldexp() function.
  long double shifted_real = ldexpl(input, sighash);
  // Post-condition: shifted_real == input * (2 ^ sighash)

  // Account for the next bit after the number of significant bits
  //   by adding or subtracting 1/2.
  // After doing this, truncating to the number of significant bits
  //   will actually be doing a rounding operation.
  if (input >= 0) {
      shifted_real += 0.5;
  } else {
      shifted_real -= 0.5;
  }

  // Truncate to the number of significant bits using the modf() function.
  long double shifted_real_int;
  //double shifted_real_frac = modf(shifted_real, &shifted_real_int);
  (void) modfl(shifted_real, &shifted_real_int);
  // Post-condition: shifted_real == shifted_real_int + shifted_real_frac
  // Post-condition: [shifted_real >= 0.0] --> [0.0 <= shifted_real_frac < 1.0]
  // Post-condition: [shifted_real < 0.0] --> [-1.0 < shifted_real_frac <= 0.0]
  // Post-condition: The value of shifted_real_int is an integer.

  // Shift the rounded real back to where it belongs using the ldexp() function.
  long double rounded_answer = ldexpl(shifted_real_int, - sighash);
  // Post-condition: rounded_answer == shifted_real_int * (2 ^ (- sighash))

  // printf("Hash debug:  input = %25.16la   output = %25.16la   sighash = %d\n", input, rounded_answer, sighash);

  *output = rounded_answer;
}

/*!
 * \ingroup larc
 * \brief Compares the double precision input to a threshold value; returns 0 if it is less than the threshold, returns the input if it is not
 * \param output The double precision floating point result of the function
 * \param input The double precision floating point input
 */
void collapse_near_zero_real(long double *output, const long double input)
{
  if (VERBOSE>DEBUG) printf("in %s\n",__func__);
  long double threshold = get_zerorealthresh();
  long double abs_input = fabsl(input);

#ifdef DEBUG
  char *larc_sca_get_str_real(const long double);
  char *temp1 = larc_sca_get_str_real(threshold);
  char *temp2 = larc_sca_get_str_real(abs_input);
  printf("threshold is %s\n",temp1);
  printf("abs(scalar) is %s\n",temp2);
  free(temp1);
  free(temp2);
#endif

  if (abs_input < threshold) {
     if (VERBOSE>DEBUG) printf("\tsetting scalar to zero\n\n");
     *output = 0.0;
  }
  else
  {
     if (VERBOSE>DEBUG) printf("\tnot setting scalar to zero\n\n");
     *output = input;
  }
  if (VERBOSE>DEBUG) printf("exiting %s\n",__func__);
  return;
}
#endif

/*!
 * \ingroup larc
 * \brief Compares the input to a threshold value; returns 0 if the input is less than the threshold, or the input value if it is not. Integer types do not collapse.
 * \param output The scalarType result of the function
 * \param input The scalarType input to the function
 */
void collapse_near_zero(scalarType *output, const scalarType input)
{
#if defined(USE_INTEGER)
   *output = input;
#elif defined(USE_MPINTEGER)
   mpz_set(*output,input);
#elif defined(USE_REAL)
   collapse_near_zero_real(output,input);
#elif defined(USE_COMPLEX)
   long double realval, imagval;
   collapse_near_zero_real(&realval,creall(input));
   collapse_near_zero_real(&imagval,cimagl(input));
   *output = realval + I*imagval;
#elif defined(USE_MPRATIONAL)
   collapse_near_zero_mprational(output,input);
#elif defined(USE_MPRATCOMPLEX)
   collapse_near_zero_mprational(&((*output)->real),input->real);
   collapse_near_zero_mprational(&((*output)->imag),input->imag);
#elif defined(USE_MPREAL)
   collapse_near_zero_mpreal(output,input);
#elif defined(USE_MPCOMPLEX)
   collapse_near_zero_mpreal(&(mpc_realref(*output)),mpc_realref(input));
   collapse_near_zero_mpreal(&(mpc_imagref(*output)),mpc_imagref(input));
#else
   fprintf(stderr,"no valid scalarType for collapse_near_zero!\n");
   exit(-1);
#endif
   return;
}

/*!
 * \ingroup larc
 * \brief Returns the input rounded to the SIGHASH most significant bits. For integer types, returns the input without rounding.
 * \param output The scalarType output of the function
 * \param input The scalarType input to the function
 */
void round_sig_fig(scalarType *output, const scalarType input)
{
#if defined(USE_INTEGER)
   *output = input;
#elif defined(USE_REAL)
   round_sig_fig_real(output,input);
#elif defined(USE_COMPLEX)
   long double realval, imagval;
   round_sig_fig_real(&realval,creall(input));
   round_sig_fig_real(&imagval,cimagl(input));
   *output = realval + I*imagval;
#elif defined(USE_MPINTEGER)
   mpz_set(*output,input);
#elif defined(USE_MPRATIONAL)
   round_sig_fig_mprational(output,input);
#elif defined(USE_MPRATCOMPLEX)
   round_sig_fig_mprational(&((*output)->real),input->real);
   round_sig_fig_mprational(&((*output)->imag),input->imag);
#elif defined(USE_MPREAL)
   round_sig_fig_mpreal(output,input);
#elif defined(USE_MPCOMPLEX)
   round_sig_fig_mpreal(&(mpc_realref(*output)),mpc_realref(input));
   round_sig_fig_mpreal(&(mpc_imagref(*output)),mpc_imagref(input));
#else
   fprintf(stderr,"no valid scalarType for round_sig_fig!\n");
   exit(-1);
#endif
   return;
}

/*!
 * \ingroup larc
 * \brief The neighborhood hash approximation function
 * \param output The neighborhood representative for the input value
 * \param input The input value
 */
void larc_nbhd_approx(scalarType *output, const scalarType input)
{
   collapse_near_zero(output,input);
   round_sig_fig(output,*output);
   return;
}


#if defined(USE_INTEGER) || defined(USE_REAL) || defined(USE_COMPLEX)
/******************************************************************************
 * empty wrappers
 ******************************************************************************/
/*!
 * \ingroup larc
 * \brief This is a placeholder for any sca_* function we do not want to define
 * \param scalar a pointer to a scalarType
 */
void empty1(scalarType *scalar){;}
#endif


// Hashes for multi-precision types.

#if defined(USE_MPINTEGER) || defined(USE_MPRATIONAL) || defined(USE_MPRATCOMPLEX)
/*!
 * \ingroup larc
 * \brief Hashes a multiprecision integer to an integer in range [0,2**exponent)
 * \param sc The multiprecision integer
 * \param exponent Determines the size of the hash table (2**exponent)
 */
uint64_t larc_mpz_hash(const mpz_t sc, uint64_t exponent)
{
/*
 * This function, modeled off Jenny's recursive hash from integer list, takes
 * a GMP MP Integer and a hash table size exponent e (2^e = size of table).
 * We take 64-bit portions of the integer and apply a 64-bit mult-golden_hash
 * (which multiplies by the fractional part of the golden mean and returns
 * the full 64 bits of the result). Then we bit-wise mod 2 add in the next
 * 64-bit portion. We alternate these steps until the integer is exhausted and
 * then end with an application of mult_golden_hash and only return exponent
 * bits.
 * Author: Andy Wills
 */
    uint64_t ping, pong;
    // GMP integers are stored in (32 or) 64 bit words called `limbs'.
    size_t limbNum = mpz_size(sc);
    // a limb is "currently a long, but on some systems it's an int for
    // efficiency, and on some systems it will be long long in the future."
    // -GMP 6.1.2 docs.
    //if (64 != mp_bits_per_limb){printf("WARNING: hash function for GMP integers not optimized for implementations with limb size %ld instead of %ld.\n", mp_bits_per_limb, 64);}
    mp_limb_t limb;
    int i;

    // initializing ping with the sign of sc: 0 or 1.
    ping = (mpz_sgn(sc)<0);

    // Go from most significant limb to least significant limb. Reversing the
    // order allows coalescence between a, b if b = a * 2^(k64). In particular
    // h(0 + h(0 + ...))) = 0 for any number of interations: two zero inputs
    // that are hashed different numbers of times are still zero. Reversing the
    // direction guarantees that the first value hashed is nonzero. Two inputs
    // that are hashed a different number of times are different (when nonzero).
    // Note that this type of coalescence would be impossible if the number of
    // limbs was fixed.
    // Verified by experiment that mpz_t allocated for n limbs but only holding
    // m limbs (where m < n) has limbNum = m. (Theoretically limbNum could
    // be the allocation size and unused limbs could be set to zero. If that was
    // the case, the two value could be the same with different hashes, which
    // would be bad.)
    for (i = limbNum-1; i >= 0; i--){
        pong = mult_golden_hash(ping, 64);
        limb = mpz_getlimbn(sc, i);
        ping = pong ^ limb;
    }
    pong = mult_golden_hash(ping, exponent);

    return pong;
}
#endif

#if defined(USE_MPREAL) || defined(USE_MPCOMPLEX)
/*!
 * \ingroup larc
 * \brief Hashes a multiprecision real to an integer in range [0,2**exponent)
 * \param sc The multiprecision real
 * \param exponent Determines the size of the hash table (2**exponent)
 */
uint64_t larc_mpfr_hash(const mpfr_t sc, uint64_t exponent)
{
    if (VERBOSE > DEBUG) printf("in %s:\n",__func__);
    uint64_t ping, pong;
    size_t limbNum = (mpreal_precision + mp_bits_per_limb - 1) / mp_bits_per_limb;
    mp_limb_t limb;
    int i;

    if (VERBOSE>DEBUG) printf("\tlimbNum = %d\n",(int)limbNum);
    // initializing ping with the sign of sc: 0 or 1.
    ping = (mpfr_sgn(sc)<0);
    if (VERBOSE>DEBUG) printf("\tping initialized to %lu\n",ping);

    for (i = limbNum-1; i >= 0; i--){
    //for (i = 0; i < limbNum; ++i){
        pong = mult_golden_hash(ping, 64);
        limb = sc->_mpfr_d[i];
        ping = pong ^ limb;
#ifdef DEBUG
        printf("in loop with i = %d:",i);
        printf(">> ping = %lu, pong = %lu\n",ping,pong);
        printf("limb = %lu\n",(long unsigned)limb);
#endif
    }
    pong = mult_golden_hash(ping, exponent);
    if (VERBOSE>DEBUG) {
        printf("\tfinal value of pong is %lu\n",pong);
        printf("exiting %s\n",__func__);
    }
    return pong;
#undef DEBUG
}
#endif


/*************************************************************
 *  These functions define the various scalarType functions  *
 *  used in LARC, so that all definitions depending on a     *
 *  particular type (eg, multiprecision integers) are hidden *
 *  behind a scalarType sca_* function. This is done to keep *
 *  all the #defined functions in scalars.c.                 *
 *************************************************************/

#if defined(USE_MPINTEGER)
/******************************************************************************
 * mpz_t (GMP) scalar operations:
 * these wrap appropriate mpz routines and convert return pointers to regular
 * mpz_t structs
 ******************************************************************************/
void larc_sca_init_mpinteger(mpz_t *sc){mpz_init(*sc);}
void larc_sca_clear_mpinteger(mpz_t *sc){mpz_clear(*sc);}
void larc_sca_set_mpinteger(mpz_t *rsc, const mpz_t sc){mpz_set(*rsc, sc);}
void larc_sca_set_str_mpinteger(mpz_t *rsc, const char *input_str){mpz_set_str(*rsc, input_str, io_string_base);}
void larc_sca_set_2ldoubles_mpinteger(mpz_t *rsc, long double real_val, long double imag_val){
    if (0.0 != imag_val) {
        fprintf(stderr,"ERROR in %s: MPINTEGER cannot have imaginary component.\n", __func__);
        exit(1);
    }
    mpfr_set_ld(scratchVars.mpreal, real_val, MPFR_RNDN);
    mpfr_get_z(*rsc, scratchVars.mpreal, MPFR_RNDN);
}
char *larc_sca_get_str_mpinteger(const mpz_t sc){return mpz_get_str(NULL, io_string_base, sc);}
void larc_sca_add_mpinteger (mpz_t *rsc, const mpz_t sc1, const mpz_t sc2){mpz_add(*rsc, sc1, sc2);}
void larc_sca_mult_mpinteger(mpz_t *rsc, const mpz_t sc1, const mpz_t sc2){mpz_mul(*rsc, sc1, sc2);}
// NOTE: mpz_tdiv_q returns the quotient of the division (not remainder)
void larc_sca_divide_mpinteger(mpz_t *rsc, const mpz_t sc1, const mpz_t sc2){mpz_tdiv_q(*rsc, sc1, sc2);}
void larc_sca_sqrt_mpinteger(mpz_t *sqroot, const mpz_t sc){mpz_sqrt(*sqroot, sc);}
void larc_sca_norm_mpinteger(mpz_t *norm, const mpz_t sc){mpz_abs(*norm, sc);}
int  larc_sca_eq_mpinteger(const mpz_t sc1, const mpz_t sc2){return 0 == mpz_cmp(sc1, sc2);}
int  larc_sca_eq_approx_mpinteger(const mpz_t sc1, const mpz_t sc2){return 0 == mpz_cmp(sc1, sc2);}
#endif

#if defined(USE_MPRATIONAL) || defined(USE_MPRATCOMPLEX)
/******************************************************************************
 * mpq_t (GMP) scalar operations:
 * these wrap appropriate mpz routines and convert return pointers to regular
 * mpq_t structs
 ******************************************************************************/
void larc_sca_init_mprational(mpq_t *sc){mpq_init(*sc);}
void larc_sca_clear_mprational(mpq_t *sc){mpq_clear(*sc);}
void larc_sca_set_mprational(mpq_t *rsc, const mpq_t sc){mpq_set(*rsc, sc);}

int larc_sca_eq_approx_mprational(const mpq_t val1, const mpq_t val2) {
  if (VERBOSE>DEBUG) printf("in %s\n",__func__);
  mpq_t val1_approx;
  mpq_t val2_approx;
  mpq_t *temp = &scratchVars.mprational;

  larc_sca_init_mprational(&val1_approx);
  larc_sca_init_mprational(&val2_approx);
  collapse_near_zero_mprational(temp,val1);
  round_sig_fig_mprational(&val1_approx,*temp);
  collapse_near_zero_mprational(temp,val2);
  round_sig_fig_mprational(&val2_approx,*temp);
  // compare calculated values
  int approx_same = mpq_equal(val1_approx,val2_approx);
  larc_sca_clear_mprational(&val1_approx);
  larc_sca_clear_mprational(&val2_approx);
  if (VERBOSE>DEBUG) printf("exiting %s\n",__func__);

  return approx_same;
}
void larc_sca_set_str_mprational(mpq_t *rsc, const char *input_str){
    mpq_set_str(*rsc, input_str, io_string_base);
    mpq_canonicalize(*rsc);
}
char *larc_sca_get_str_mprational(const mpq_t sc){return mpq_get_str(NULL, io_string_base, sc);}
void larc_sca_set_2ldoubles_mprational(mpq_t *rsc, long double real_val, long double imag_val){
    if (0.0 != imag_val) {
        fprintf(stderr,"ERROR in %s: MPRATIONAL cannot have imaginary component.\n", __func__);
        exit(1);
    }
    mpq_set_d(*rsc, real_val);
}
void larc_sca_add_mprational(mpq_t *rsc, const mpq_t sc1, const mpq_t sc2){mpq_add(*rsc, sc1, sc2);}
void larc_sca_mult_mprational(mpq_t *rsc, const mpq_t sc1, const mpq_t sc2){mpq_mul(*rsc, sc1, sc2);}
void larc_sca_divide_mprational(mpq_t *rsc, const mpq_t sc1, const mpq_t sc2){mpq_div(*rsc, sc1, sc2);}
uint64_t larc_sca_hash_mprational(const mpq_t sc, uint64_t exponent)
{
    uint64_t num_hash = larc_mpz_hash(mpq_numref(sc), 64);
    uint64_t den_hash = larc_mpz_hash(mpq_denref(sc), 64);

    return recursive_hash_from_two_integers(num_hash, den_hash, exponent);
}
void larc_sca_sqrt_mprational(mpq_t *sqroot, const mpq_t sc){
    long double sc_dbl, sqroot_dbl;
    sc_dbl = mpq_get_d(sc);
    sqroot_dbl = sqrtl(sc_dbl);
    mpq_set_d(*sqroot, sqroot_dbl);
}
void larc_sca_norm_mprational(mpq_t *norm, const mpq_t sc){mpq_abs(*norm, sc);}
#endif

#if defined(USE_MPRATCOMPLEX)
/******************************************************************************
 * Multi-precision complex numbers
 ******************************************************************************/
void larc_sca_init_mpratcomplex(larc_mpratcomplex_t *val){
    mpq_init((*val)->real);
    mpq_init((*val)->imag);
}

void larc_sca_clear_mpratcomplex(larc_mpratcomplex_t *val){
    mpq_clear((*val)->real);
    mpq_clear((*val)->imag);
}

void larc_sca_set_mpratcomplex(larc_mpratcomplex_t *rsc, const larc_mpratcomplex_t sc){
    larc_sca_set_mprational(&(*rsc)->real, sc->real);
    larc_sca_set_mprational(&(*rsc)->imag, sc->imag);
}

void larc_sca_set_2ldoubles_mpratcomplex(larc_mpratcomplex_t *rsc, long double real_val, long double imag_val){
    mpq_set_d((*rsc)->real, real_val);
    mpq_set_d((*rsc)->imag, imag_val);
}

int larc_sca_eq_mpratcomplex(const larc_mpratcomplex_t val1, const larc_mpratcomplex_t val2){
    return (   mpq_equal(val1->real, val2->real)
            && mpq_equal(val1->imag, val2->imag) );
}

int larc_sca_eq_approx_mpratcomplex(const larc_mpratcomplex_t val1, const larc_mpratcomplex_t val2){
    return (   larc_sca_eq_approx_mprational(val1->real, val2->real)
            && larc_sca_eq_approx_mprational(val1->imag, val2->imag) );
}

int larc_sca_cmp_mpratcomplex(const larc_mpratcomplex_t val1, const larc_mpratcomplex_t val2){
    int cmpval = mpq_cmp(val1->real, val2->real);
    if (0 != cmpval) {
        return cmpval;
    } else {
        return mpq_cmp(val1->imag, val2->imag);
    }
}

void larc_sca_add_mpratcomplex(larc_mpratcomplex_t *rsc, const larc_mpratcomplex_t sc1, const larc_mpratcomplex_t sc2){
    mpq_add((*rsc)->real, sc1->real, sc2->real);
    mpq_add((*rsc)->imag, sc1->imag, sc2->imag);
}

void larc_sca_mult_mpratcomplex(larc_mpratcomplex_t *rsc, const larc_mpratcomplex_t sc1, const larc_mpratcomplex_t sc2){
    mpq_t temp_r1_r2, temp_r1_i2, temp_i1_r2, temp_i1_i2;
    mpq_init(temp_r1_r2);
    mpq_init(temp_r1_i2);
    mpq_init(temp_i1_r2);
    mpq_init(temp_i1_i2);
    mpq_mul(temp_r1_r2, sc1->real, sc2->real);
    mpq_mul(temp_r1_i2, sc1->real, sc2->imag);
    mpq_mul(temp_i1_r2, sc1->imag, sc2->real);
    mpq_mul(temp_i1_i2, sc1->imag, sc2->imag);
    mpq_sub((*rsc)->real, temp_r1_r2, temp_i1_i2);
    mpq_add((*rsc)->imag, temp_r1_i2, temp_i1_r2);
    mpq_clear(temp_r1_r2);
    mpq_clear(temp_r1_i2);
    mpq_clear(temp_i1_r2);
    mpq_clear(temp_i1_i2);
}

void larc_sca_sqrt_mpratcomplex(larc_mpratcomplex_t *sqroot, const larc_mpratcomplex_t sc){
    long double complex sc_dbl, sqroot_dbl;
    sc_dbl = mpq_get_d(sc->real) + I*mpq_get_d(sc->imag);
    sqroot_dbl = csqrtl(sc_dbl);
    mpq_set_d((*sqroot)->real, creall(sqroot_dbl));
    mpq_set_d((*sqroot)->imag, cimagl(sqroot_dbl));
}

void larc_sca_norm_mpratcomplex(larc_mpratcomplex_t *norm, const larc_mpratcomplex_t sc){
    mpq_t temp_real, temp_imag;
    mpq_init(temp_real);
    mpq_init(temp_imag);
    mpq_mul(temp_real, sc->real, sc->real);
    mpq_mul(temp_imag, sc->imag, sc->imag);
    mpq_add((*norm)->real, temp_real, temp_imag);
    larc_sca_sqrt_mprational(&(*norm)->real, (*norm)->real);
    mpq_set_ui((*norm)->imag, 0, 1);
    mpq_clear(temp_real);
    mpq_clear(temp_imag);
}

void larc_sca_conj_mpratcomplex(larc_mpratcomplex_t *conj_out, const larc_mpratcomplex_t sc){
    mpq_set((*conj_out)->real, sc->real);
    mpq_neg((*conj_out)->imag, sc->imag);
}

int larc_sca_is_real_mpratcomplex(const larc_mpratcomplex_t sc){
    // A complex number is real iff the imaginary part is equal to zero.
    return (0 == mpq_cmp_ui(sc->imag, 0, 1));
}

int larc_sca_is_imag_mpratcomplex(const larc_mpratcomplex_t sc){
    // A complex number is imaginary iff the real part is equal to zero.
    return (0 == mpq_cmp_ui(sc->real, 0, 1));
}

void larc_sca_divide_mpratcomplex(larc_mpratcomplex_t *rsc, const larc_mpratcomplex_t sc1, const larc_mpratcomplex_t sc2){
    if (larc_sca_is_real_mpratcomplex(sc2)) {
        mpq_div((*rsc)->real, sc1->real, sc2->real);
        mpq_div((*rsc)->imag, sc1->imag, sc2->real);
    } else if (larc_sca_is_imag_mpratcomplex(sc2)) {
        mpq_div((*rsc)->real, sc1->imag, sc2->imag);
        mpq_div((*rsc)->imag, sc1->real, sc2->imag);
        mpq_neg((*rsc)->imag, (*rsc)->imag);
    } else {
        larc_sca_conj_mpratcomplex(rsc, sc2);
        larc_sca_mult_mpratcomplex(rsc, sc1, *rsc);
        mpq_t temp_r2_r2, temp_i2_i2, temp_norm_squared;
        mpq_init(temp_r2_r2);
        mpq_init(temp_i2_i2);
        mpq_init(temp_norm_squared);
        mpq_mul(temp_r2_r2, sc2->real, sc2->real);
        mpq_mul(temp_i2_i2, sc2->imag, sc2->imag);
        mpq_add(temp_norm_squared, temp_r2_r2, temp_i2_i2);
        mpq_div((*rsc)->real, (*rsc)->real, temp_norm_squared);
        mpq_div((*rsc)->imag, (*rsc)->imag, temp_norm_squared);
        mpq_clear(temp_r2_r2);
        mpq_clear(temp_i2_i2);
        mpq_clear(temp_norm_squared);
    }
}

uint64_t larc_sca_hash_mpratcomplex(const larc_mpratcomplex_t sc, uint64_t exponent)
{
    uint64_t real_hash = larc_sca_hash_mprational(sc->real, 64);
    uint64_t imag_hash = larc_sca_hash_mprational(sc->imag, 64);

    return recursive_hash_from_two_integers(real_hash, imag_hash, exponent);
}

char *larc_sca_get_str_mpratcomplex(const larc_mpratcomplex_t val){
// In mpq_get_str, "if str is NULL, the result string is allocated using the
// current allocation function (see Custom Allocation). The block will be
// strlen(str)+1 bytes, that being exactly enough for the string and
// null-terminator." -GMP6.1.2
    char *real = mpq_get_str(NULL, io_string_base, val->real);
    char *imag = mpq_get_str(NULL, io_string_base, val->imag);
    char *num = calloc(strlen(real) + strlen(imag) + strlen("+I*")+1, sizeof(char));
    strcat(num, real);
    char *imag_start = imag;
    if (imag_start[0] == '+') {
        ++imag_start;
        strcat(num, "+I*");
    } else if (imag_start[0] == '-') {
        ++imag_start;
        strcat(num, "-I*");
    } else {
        strcat(num, "+I*");
    }
    strcat(num, imag_start);
    free(real);
    free(imag);
    return num;
}

void larc_sca_set_str_mpratcomplex(larc_mpratcomplex_t *rsc, const char *input_str){
    char *ipos = strchr(input_str, 'I');
    if (ipos == NULL) {
        mpq_set_str((*rsc)->real, input_str, io_string_base);
        mpq_canonicalize((*rsc)->real);
        mpq_set_ui((*rsc)->imag, 0, 1);
    } else if ((ipos[1] == '*') && ((ipos == input_str) || (ipos[-1] == '+') || (ipos[-1] == '-'))) {
        mpq_set_ui((*rsc)->real, 0, 1);
        mpq_set_str((*rsc)->imag, ipos + 2, io_string_base);
        mpq_canonicalize((*rsc)->imag);
        if ((ipos != input_str) && (ipos[-1] == '-')) {
            mpq_neg((*rsc)->imag, (*rsc)->imag);
        }
        if (ipos > input_str + 1) {
            char *realpart = calloc(ipos - input_str, sizeof(char));
            strncpy(realpart, input_str, ipos - input_str - 1);
            realpart[ipos - input_str - 1] = '\0';
            mpq_set_str((*rsc)->real, realpart, io_string_base);
            mpq_canonicalize((*rsc)->real);
            free(realpart);
        }
    } else {
        fprintf(stderr,"ERROR in %s: could not convert '%s' to rational complex number.\n", __func__, input_str);
        exit(1);
    }
}
#endif

#if defined(USE_MPREAL) || defined(USE_MPCOMPLEX)
/******************************************************************************
 * mpfr_t (MPFR) scalar operations:
 * these wrap appropriate mpfr routines and convert return pointers to regular
 * mpfr_t structs
 ******************************************************************************/
void larc_sca_init_mpreal(mpfr_t *sc){mpfr_init2(*sc, mpreal_precision);}
void larc_sca_clear_mpreal(mpfr_t *sc){mpfr_clear(*sc);}
void larc_sca_set_mpreal(mpfr_t *rsc, const mpfr_t sc){mpfr_set(*rsc, sc, MPFR_RNDN);}

int larc_sca_eq_approx_mpreal(const mpfr_t val1, const mpfr_t val2) {
  if (VERBOSE>DEBUG) printf("in %s\n",__func__);
  mpfr_t val1_approx;
  mpfr_t val2_approx;
  mpfr_t *temp = &scratchVars.mpreal;

  larc_sca_init_mpreal(&val1_approx);
  larc_sca_init_mpreal(&val2_approx);
  collapse_near_zero_mpreal(temp,val1); 
  round_sig_fig_mpreal(&val1_approx,*temp);
  collapse_near_zero_mpreal(temp,val2); 
  round_sig_fig_mpreal(&val2_approx,*temp);
  // compare calculated values
  int approx_same = mpfr_equal_p(val1_approx,val2_approx);
  larc_sca_clear_mpreal(&val1_approx);
  larc_sca_clear_mpreal(&val2_approx);
  if (VERBOSE>DEBUG) printf("exiting %s\n",__func__);

  return approx_same;
}

char *larc_sca_get_str_mpreal(const mpfr_t sc){
    mpfr_exp_t scexp;
    char *raw = mpfr_get_str(NULL, &scexp, io_string_base, 0, sc, MPFR_RNDN);
    char exponent_str[25];
    sprintf(exponent_str, "e%ld", scexp);
    char *s = calloc(strlen(raw) + strlen(exponent_str) + 4, sizeof(char));
    s[0] = '\0';
    int is_negative = 0;
    int is_numeric = 0;
    if (raw[0] == '-') {
        is_negative = 1;
    }
    if (('0' <= raw[is_negative]) && (raw[is_negative] <= '9')) {
        is_numeric = 1;
    }
    if (is_negative) {
        strcat(s, "-");
    }
    if (is_numeric) {
        strcat(s, "0.");
    }
    strcat(s, &raw[is_negative]);
    if (is_numeric) {
        while ((s[strlen(s)-1] == '0') && ('0' <= s[strlen(s)-2]) && (s[strlen(s)-2] <= '9')) {
            s[strlen(s)-1] = '\0';
        }
    }
    strcat(s, exponent_str);
    mpfr_free_str(raw);
    return s;
}
void larc_sca_set_str_mpreal(mpfr_t *rsc, const char *input_str){
    mpfr_set_str(*rsc, input_str, 0, MPFR_RNDN);
}
void larc_sca_set_2ldoubles_mpreal(mpfr_t *rsc, long double real_val, long double imag_val){
    if (0.0 != imag_val) {
        fprintf(stderr,"ERROR in %s: MPREAL cannot have imaginary component.\n", __func__);
        exit(1);
    }
    mpfr_set_ld(*rsc, real_val, MPFR_RNDN);
}
void larc_sca_add_mpreal(mpfr_t *rsc, const mpfr_t sc1, const mpfr_t sc2){mpfr_add(*rsc, sc1, sc2, MPFR_RNDN);}
void larc_sca_mult_mpreal(mpfr_t *rsc, const mpfr_t sc1, const mpfr_t sc2){mpfr_mul(*rsc, sc1, sc2, MPFR_RNDN);}
void larc_sca_divide_mpreal(mpfr_t *rsc, const mpfr_t sc1, const mpfr_t sc2){mpfr_div(*rsc, sc1, sc2, MPFR_RNDN);}
uint64_t larc_sca_hash_mpreal(const mpfr_t sc, uint64_t exponent)
{
    return larc_mpfr_hash(sc, exponent);
}
void larc_sca_sqrt_mpreal(mpfr_t *sqroot, const mpfr_t sc){mpfr_sqrt(*sqroot, sc, MPFR_RNDN);}
void larc_sca_norm_mpreal(mpfr_t *norm, const mpfr_t sc){mpfr_abs(*norm, sc, MPFR_RNDN);}

int larc_sca_is_real_mpreal(const mpfr_t sc){return 1;}
void larc_sca_conj_mpreal(mpfr_t *conj_out, const mpfr_t sc){mpfr_set(*conj_out, sc, MPFR_RNDN);}

int larc_sca_eq_mpreal(const mpfr_t val1, const mpfr_t val2){
    return mpfr_equal_p(val1, val2);
}

int larc_sca_cmp_mpreal(const mpfr_t val1, const mpfr_t val2){
    return mpfr_cmp(val1, val2);
}
#endif

#if defined(USE_MPCOMPLEX)
/******************************************************************************
 * mpc_t (MPC) scalar operations:
 * these wrap appropriate mpfr routines and convert return pointers to regular
 * mpc_t structs
 ******************************************************************************/
void larc_sca_init_mpcomplex(mpc_t *sc){mpc_init2(*sc, mpreal_precision);}
void larc_sca_clear_mpcomplex(mpc_t *sc){mpc_clear(*sc);}
void larc_sca_set_mpcomplex(mpc_t *rsc, const mpc_t sc){mpc_set(*rsc, sc, MPC_RNDNN);}

int larc_sca_eq_approx_mpcomplex(const mpc_t val1, const mpc_t val2) {
  return larc_sca_eq_approx_mpreal(mpc_realref(val1), mpc_realref(val2))
      && larc_sca_eq_approx_mpreal(mpc_imagref(val1), mpc_imagref(val2));
}
char *larc_sca_get_str_mpcomplex(const mpc_t sc){
    // mpc_set_str returns complex numbers in the following format:
    //   "(" ++ real_part ++ " " ++ imag_part ++ ")"
    char *raw = mpc_get_str(io_string_base, 0, sc, MPC_RNDNN);

    // The following code converts this format to real_part+I*imag_part
    // or real_part-I*imag_part as appropriate.
    char *s = calloc(strlen(raw) + 6, sizeof(char));
    int rindex = 0;  // Index into raw[].
    int sindex = 0;  // Index into s[].
    if (raw[rindex] == '(') {
        ++rindex;
    }
    while ((rindex < strlen(raw)) && (raw[rindex] != ' ') && (raw[rindex] != ')')) {
        s[sindex] = raw[rindex];
        ++rindex;
        ++sindex;
    }
    s[sindex] = '+';
    int imag_sign_index = sindex;  // Save location of sign in case we need to change it.
    ++sindex;
    s[sindex] = 'I';
    ++sindex;
    s[sindex] = '*';
    ++sindex;
    ++rindex;
    if ((rindex < strlen(raw) && ((raw[rindex] == '+') || (raw[rindex] == '-')))) {
        s[imag_sign_index] = raw[rindex];
        ++rindex;
    }
    while ((rindex < strlen(raw)) && (raw[rindex] != ' ') && (raw[rindex] != ')')) {
        s[sindex] = raw[rindex];
        ++rindex;
        ++sindex;
    }
    s[sindex] = '\0';
    mpfr_free_str(raw);
    return s;
}
void larc_sca_set_str_mpcomplex(mpc_t *rsc, const char *input_str){
    char *ipos = strchr(input_str, 'I');
    if (ipos == NULL) {
        mpfr_set_str(mpc_realref(*rsc), input_str, 0, MPFR_RNDN);
        mpfr_set_ui(mpc_imagref(*rsc), 0, MPFR_RNDN);
    } else if ((ipos[1] == '*') && ((ipos == input_str) || (ipos[-1] == '+') || (ipos[-1] == '-'))) {
        mpfr_set_ui(mpc_realref(*rsc), 0, MPFR_RNDN);
        mpfr_set_str(mpc_imagref(*rsc), ipos + 2, 0, MPFR_RNDN);
        if ((ipos != input_str) && (ipos[-1] == '-')) {
            mpfr_neg(mpc_imagref(*rsc), mpc_imagref(*rsc), MPFR_RNDN);
        }
        if (ipos > input_str + 1) {
            char *realpart = calloc(ipos - input_str, sizeof(char));
            strncpy(realpart, input_str, ipos - input_str - 1);
            realpart[ipos - input_str - 1] = '\0';
            mpfr_set_str(mpc_realref(*rsc), realpart, 0, MPFR_RNDN);
            free(realpart);
        }
    } else {
        fprintf(stderr,"ERROR in %s: could not convert '%s' to mpcomplex number.\n", __func__, input_str);
        exit(1);
    }
}
void larc_sca_set_2ldoubles_mpcomplex(mpc_t *rsc, long double real_val, long double imag_val){
    mpc_set_ld_ld(*rsc, real_val, imag_val, MPC_RNDNN);
}
void larc_sca_add_mpcomplex(mpc_t *rsc, const mpc_t sc1, const mpc_t sc2){mpc_add(*rsc, sc1, sc2, MPC_RNDNN);}
void larc_sca_mult_mpcomplex(mpc_t *rsc, const mpc_t sc1, const mpc_t sc2){mpc_mul(*rsc, sc1, sc2, MPC_RNDNN);}
void larc_sca_divide_mpcomplex(mpc_t *rsc, const mpc_t sc1, const mpc_t sc2){mpc_div(*rsc, sc1, sc2, MPC_RNDNN);}
uint64_t larc_sca_hash_mpcomplex(const mpc_t sc, uint64_t exponent)
{
    uint64_t realhash = larc_mpfr_hash(mpc_realref(sc), 64);
    uint64_t imaghash = larc_mpfr_hash(mpc_imagref(sc), 64);
#ifdef DEBUG
    printf("in %s, realhash = %lu, imaghash = %lu\n",
        __func__, realhash, imaghash);
#endif
    uint64_t hash = recursive_hash_from_two_integers(realhash, imaghash, exponent);
#ifdef DEBUG
    printf("result of recursive_hash_from_two_integers is %lu\n",hash);
#endif
    return hash;
#undef DEBUG
}
void larc_sca_sqrt_mpcomplex(mpc_t *sqroot, const mpc_t sc){
    mpc_sqrt(*sqroot, sc, MPC_RNDNN);
 }
void larc_sca_norm_mpcomplex(mpc_t *norm, const mpc_t sc){
    mpc_abs(mpc_realref(*norm), sc, MPC_RNDNN);
    mpfr_set_ui(mpc_imagref(*norm), 0, MPFR_RNDN);
}

int larc_sca_is_real_mpcomplex(const mpc_t sc){
    return mpfr_zero_p(mpc_imagref(sc));
}
void larc_sca_conj_mpcomplex(mpc_t *conj_out, const mpc_t sc){mpc_conj(*conj_out, sc, MPC_RNDNN);}

int larc_sca_eq_mpcomplex(const mpc_t val1, const mpc_t val2){
    return mpfr_equal_p(mpc_realref(val1), mpc_realref(val2))
        && mpfr_equal_p(mpc_imagref(val1), mpc_imagref(val2));
}

int larc_sca_cmp_mpcomplex(const mpc_t val1, const mpc_t val2){
    int cmpval = mpfr_cmp(mpc_realref(val1), mpc_realref(val2));
    if (0 != cmpval) {
        return cmpval;
    } else {
        return mpfr_cmp(mpc_imagref(val1), mpc_imagref(val2));
    }
}
#endif

#if defined(USE_INTEGER) || defined(USE_REAL) || defined(USE_COMPLEX)
/******************************************************************************
 * Integer, Complex, Real:
 * these wrap the usual operations
 * TODO: for real/complex - do we need to worry about FPCLASSIFY here?
 ******************************************************************************/
void larc_sca_set_arith (scalarType *ret, const scalarType a){*ret = a;}
void larc_sca_add_arith (scalarType *ret, const scalarType a, const scalarType b){*ret = a + b;}
void larc_sca_mult_arith(scalarType *ret, const scalarType a, const scalarType b){*ret = a * b;}
void larc_sca_divide_arith(scalarType *ret, const scalarType a, const scalarType b){*ret = a / b;}
#endif

#if defined(USE_REAL) || defined(USE_COMPLEX)
// Now using strtold() instead of sscanf() for best coverage of numeric format possibilities.
void larc_sca_set_str_real(long double *ret, const char *input_str){
    char *endptr = NULL;
    errno = 0;  // Zero out errno so that we can tell if strtold returns an error.
    *ret = strtold(input_str, &endptr);
    int rc = 0;  // Assume success to start with.
    if ((errno != 0) || (endptr == input_str)) {
        rc = 1;  // Failure
    }
    if (rc) {
        fprintf(stderr,"ERROR in %s: could not convert '%s' to real number.\n", __func__, input_str);
        exit(1);
    }
    return;
}
char *larc_sca_get_str_real(const long double n){
    char *out = calloc(LDBL_MANT_DEC_DIG + 21, sizeof(char));
    snprintf(out, LDBL_MANT_DEC_DIG + 20, "%.20Lg", n);
    return out;
}
void larc_sca_2ldoubles_real(long double *ret, long double real_val, long double imag_val){
    if (0.0L != imag_val) {
        fprintf(stderr,"ERROR in %s: REAL cannot have imaginary component.\n", __func__);
        exit(1);
    }
    *ret = real_val;
}
void larc_sca_sqrt_real(long double *sqroot, const long double a){*sqroot = sqrtl(a);}
void larc_sca_norm_real(long double *norm, const long double a){*norm = fabsl(a);}
#endif

#if defined(USE_COMPLEX)
void larc_sca_set_str_complex(long double complex *ret, const char *input_str){
    long double real, imag;
    char *ipos = strchr(input_str, 'I');
    if (ipos == NULL) {
        larc_sca_set_str_real(&real, input_str);
        imag = 0.0L;
    } else if ((ipos[1] == '*') && ((ipos == input_str) || (ipos[-1] == '+') || (ipos[-1] == '-'))) {
        real = 0.0L;
        larc_sca_set_str_real(&imag, ipos + 2);
        if ((ipos != input_str) && (ipos[-1] == '-')) {
            imag = -imag;
        }
        if (ipos > input_str + 1) {
            char *realpart = calloc(ipos - input_str, sizeof(char));
            strncpy(realpart, input_str, ipos - input_str - 1);
            realpart[ipos - input_str - 1] = '\0';
            larc_sca_set_str_real(&real, realpart);
            free(realpart);
        }
    } else {
        fprintf(stderr,"ERROR in %s: could not convert '%s' to complex number.\n", __func__, input_str);
        exit(1);
    }
    *ret = real+imag*I;
}
char *larc_sca_get_str_complex(const long double complex n){
    char *out = calloc(2*LDBL_MANT_DEC_DIG + 44, sizeof(char));
    snprintf(out, 2*LDBL_MANT_DEC_DIG + 43,
	       "%.20Lg+I*%.20Lg", creall(n), cimagl(n));
    return out;
}
void larc_sca_2ldoubles_complex(long double complex *ret, long double real_val, long double imag_val){
    *ret = real_val+imag_val*I;
}
uint64_t larc_sca_hash_complex(long double complex sc, uint64_t exponent){return hash_from_two_longdoubles(creall(sc), cimagl(sc), exponent);}
void larc_sca_sqrt_complex(long double complex *sqroot, const long double complex a){*sqroot= csqrtl(a);}
void larc_sca_norm_complex(long double complex *norm, const long double complex a){*norm= cabsl(a) + I*0.0L;}
void larc_sca_conj_complex(long double complex *conj_out, const long double complex sc){*conj_out = conjl(sc);}
int larc_sca_is_real_complex(const long double complex sc){return (cimagl(sc) == 0.0L);}
#endif

#if defined(USE_REAL) || defined(USE_COMPLEX)
int larc_sca_cmp_real(const long double a, const long double b){return (a > b) - (a < b);}
#endif

#if defined(USE_COMPLEX)
int larc_sca_cmp_complex(const long double complex a, const long double complex b){
    long double complex a_norm, b_norm;
    larc_sca_norm_complex(&a_norm, a);
    larc_sca_norm_complex(&b_norm, b);
    return larc_sca_cmp_real(creall(a_norm), creall(b_norm));
}
#endif

#if defined(USE_REAL) || defined(USE_COMPLEX)
int rc_eq(const scalarType a, const scalarType b){
  // Check to see if two scalars are equal
  return a == b;
}

int rc_eq_approx(const scalarType a, const scalarType b){
  // Two scalars have equal approximation functions
  // This approximation is called neighborhood or locality-approximations
  scalarType marker_of_a, marker_of_b;
  // these will go away, to be replaced by nbhd_approx
  sca_nbhd_approx(&marker_of_a, a);
  sca_nbhd_approx(&marker_of_b, b);
  return marker_of_a == marker_of_b;
}
#endif


#if defined(USE_INTEGER)
void int_set_str(int64_t *ret, const char *input_str){*ret = atol(input_str);}
//NOTE: consider strtol() for input with base
char *int_get_str(const int64_t n){
    char *out = calloc(31, sizeof(char));
    snprintf(out, 30, "%ld", n);
    return out;
}
void int_set_2ldoubles(int64_t *ret, long double real_val, long double imag_val){
    if (0.0L != imag_val) {
        fprintf(stderr,"ERROR in %s: INTEGER cannot have imaginary component.\n", __func__);
        exit(1);
    }
    *ret = (int64_t) real_val;
    if (real_val != (long double) *ret) {
        fprintf(stderr,"ERROR in %s: Can't convert %.20Lg to INTEGER.\n", __func__, real_val);
        exit(1);
    }
}
int int_cmp(const scalarType a, const scalarType b){return (a > b) - (a < b);}
int int_eq(const scalarType a, const scalarType b){return a == b;}
int int_eq_approx(const scalarType a, const scalarType b){return a == b;}
uint64_t int_hash(int64_t sc, uint64_t exponent){return mult_golden_hash(sc, exponent);}
void larc_sca_sqrt_int(int64_t *sqroot, const int64_t sc){*sqroot = sqrt(sc);}
void larc_sca_norm_int(int64_t *norm, const int64_t sc){*norm = abs(sc);}


/******************************************************************************
 * Boolean arithmetic for integers:
 * As of 5-2019 we have no Boolean scalarType yet, however, once we do
 * create one, there are moe efficient ways of defining Boolean versions
 * of the functions sca_cmp, sca_add, sca_mult, sca_set, and sca_set_str
 ******************************************************************************/
// we'll say these only apply to integer type for the sake of examples
int  bool_cmp (const int64_t a, const int64_t b){return a && b;}
void bool_add (int64_t *ret, const int64_t a, const int64_t b){*ret = a || b;}
void bool_mult(int64_t *ret, const int64_t a, const int64_t b){*ret = a && b;}
void bool_divide(int64_t *ret, const int64_t a, const int64_t b){*ret = a;}
void bool_set (int64_t *ret, const int64_t a){*ret = (a!=0);}
// automatically convert all nonzero integers to 1.
void bool_set_str(int64_t *ret, const char *input_str){*ret = (atol(input_str) != 0);}
void bool_set_2ldoubles(int64_t *ret, long double real_val, long double imag_val){
    if (0.0 != imag_val) {
        fprintf(stderr,"ERROR in %s: INTEGER/boolean cannot have imaginary component.\n", __func__);
        exit(1);
    }
    *ret = (0.0 != real_val);
}
#endif

int always_true_is_real (const scalarType a){return 1;}


/******************************************************************************
 * Initialize scalar routines
 ******************************************************************************/
// to do standard arithmetic, use this.
void init_arithmetic_scalarOps(int verbose){
    int ops_set = check_scalarOps();
    if (ops_set == -1 && verbose > SILENT)
        fprintf(stderr,"WARNING: user was previously operating without all scalar ops defined.\n");

    if (verbose>BASIC){
        if (ops_set == 1)
            printf("Switching ");
        else
            printf("Setting ");
        printf("scalar ops to standard arithmetic for all scalarTypes.\n");
    }

/*****************************************************************/
// We need to define the sca_ functions (ops) for each scalarType*
/*****************************************************************/
    sca_nbhd_approx = larc_nbhd_approx;
#if defined(USE_MPINTEGER)
    sca_init     = larc_sca_init_mpinteger;
    sca_clear    = larc_sca_clear_mpinteger;
    sca_cmp      = mpz_cmp;
    sca_eq       = larc_sca_eq_mpinteger;
    sca_eq_approx = larc_sca_eq_approx_mpinteger; // currently no approximation implemented
    sca_set      = larc_sca_set_mpinteger;
    sca_set_str  = larc_sca_set_str_mpinteger;
    sca_set_2ldoubles = larc_sca_set_2ldoubles_mpinteger;
    sca_get_str  = larc_sca_get_str_mpinteger;
    if (64 != mp_bits_per_limb && verbose > SILENT)
    {
      fprintf(stderr,"WARNING: hash function for GMP integers not optimized");
      fprintf(stderr," for implementations with limb size %d instead of %d.\n",
              mp_bits_per_limb, 64);
    }
    sca_hash     = larc_mpz_hash;
    sca_add      = larc_sca_add_mpinteger;
    sca_mult     = larc_sca_mult_mpinteger;
    sca_divide   = larc_sca_divide_mpinteger;
    sca_sqrt     = larc_sca_sqrt_mpinteger;
    sca_norm     = larc_sca_norm_mpinteger;
    sca_conj     = larc_sca_set_mpinteger;
    sca_is_real  = always_true_is_real;
#endif
#if defined(USE_MPRATIONAL)
    sca_init     = larc_sca_init_mprational;
    sca_clear    = larc_sca_clear_mprational;
    sca_cmp      = mpq_cmp;
    sca_eq       = mpq_equal;
    sca_eq_approx    = larc_sca_eq_approx_mprational;
    sca_set      = larc_sca_set_mprational;
    sca_set_str  = larc_sca_set_str_mprational;
    sca_set_2ldoubles = larc_sca_set_2ldoubles_mprational;
    sca_get_str  = larc_sca_get_str_mprational;
    if (64 != mp_bits_per_limb && verbose > SILENT)
    {
      fprintf(stderr,"WARNING: hash function for GMP rationals not optimized");
      fprintf(stderr," for implementations with limb size %d instead of %d.\n",
              mp_bits_per_limb, 64);
    }
    sca_hash     = larc_sca_hash_mprational;
    sca_add      = larc_sca_add_mprational;
    sca_mult     = larc_sca_mult_mprational;
    sca_divide   = larc_sca_divide_mprational;
    sca_sqrt     = larc_sca_sqrt_mprational;
    sca_norm     = larc_sca_norm_mprational;
    sca_conj     = larc_sca_set_mprational;
    sca_is_real  = always_true_is_real;
#endif
#if defined(USE_MPRATCOMPLEX)
    sca_init     = larc_sca_init_mpratcomplex;
    sca_clear    = larc_sca_clear_mpratcomplex;
    sca_cmp      = larc_sca_cmp_mpratcomplex;
    sca_eq       = larc_sca_eq_mpratcomplex;
    sca_eq_approx    = larc_sca_eq_approx_mpratcomplex;
    sca_set      = larc_sca_set_mpratcomplex;
    sca_set_str  = larc_sca_set_str_mpratcomplex;
    sca_set_2ldoubles = larc_sca_set_2ldoubles_mpratcomplex;
    sca_get_str  = larc_sca_get_str_mpratcomplex;
    if (64 != mp_bits_per_limb && verbose > SILENT)
    {
      fprintf(stderr,"WARNING: hash function for GMP rationals not optimized");
      fprintf(stderr," for implementations with limb size %d instead of %d.\n",
              mp_bits_per_limb, 64);
    }
    sca_hash     = larc_sca_hash_mpratcomplex;
    sca_add      = larc_sca_add_mpratcomplex;
    sca_mult     = larc_sca_mult_mpratcomplex;
    sca_divide   = larc_sca_divide_mpratcomplex;
    sca_sqrt     = larc_sca_sqrt_mpratcomplex;
    sca_norm     = larc_sca_norm_mpratcomplex;
    sca_conj     = larc_sca_conj_mpratcomplex;
    sca_is_real  = larc_sca_is_real_mpratcomplex;
#endif
#if defined(USE_MPREAL)
    sca_init     = larc_sca_init_mpreal;
    sca_clear    = larc_sca_clear_mpreal;
    sca_cmp      = larc_sca_cmp_mpreal;
    sca_eq       = larc_sca_eq_mpreal;
    sca_eq_approx    = larc_sca_eq_approx_mpreal;
    sca_set      = larc_sca_set_mpreal;
    sca_set_str  = larc_sca_set_str_mpreal;
    sca_set_2ldoubles = larc_sca_set_2ldoubles_mpreal;
    sca_get_str  = larc_sca_get_str_mpreal;
    sca_hash     = larc_sca_hash_mpreal;
    sca_add      = larc_sca_add_mpreal;
    sca_mult     = larc_sca_mult_mpreal;
    sca_divide   = larc_sca_divide_mpreal;
    sca_sqrt     = larc_sca_sqrt_mpreal;
    sca_norm     = larc_sca_norm_mpreal;
    sca_conj     = larc_sca_conj_mpreal;
    sca_is_real  = larc_sca_is_real_mpreal;
#endif
#if defined(USE_MPCOMPLEX)
    sca_init     = larc_sca_init_mpcomplex;
    sca_clear    = larc_sca_clear_mpcomplex;
    sca_cmp      = larc_sca_cmp_mpcomplex;
    sca_eq       = larc_sca_eq_mpcomplex;
    sca_eq_approx    = larc_sca_eq_approx_mpcomplex;
    sca_set      = larc_sca_set_mpcomplex;
    sca_set_str  = larc_sca_set_str_mpcomplex;
    sca_set_2ldoubles = larc_sca_set_2ldoubles_mpcomplex;
    sca_get_str  = larc_sca_get_str_mpcomplex;
    sca_hash     = larc_sca_hash_mpcomplex;
    sca_add      = larc_sca_add_mpcomplex;
    sca_mult     = larc_sca_mult_mpcomplex;
    sca_divide   = larc_sca_divide_mpcomplex;
    sca_sqrt     = larc_sca_sqrt_mpcomplex;
    sca_norm     = larc_sca_norm_mpcomplex;
    sca_conj     = larc_sca_conj_mpcomplex;
    sca_is_real  = larc_sca_is_real_mpcomplex;
#endif
#if defined(USE_INTEGER) || defined(USE_REAL) || defined(USE_COMPLEX)
    // these five definitions are the same for COMPLEX, REAL, and INTEGER
    sca_init     = empty1;
    sca_clear    = empty1;
    sca_set      = larc_sca_set_arith;
    sca_add      = larc_sca_add_arith;
    sca_mult     = larc_sca_mult_arith;
    sca_divide   = larc_sca_divide_arith;
#endif
#if defined(USE_INTEGER)
    sca_cmp      = int_cmp;
    sca_eq    = int_eq;
    sca_eq_approx  = int_eq_approx; // currently no approximation implemented
    sca_set_str  = int_set_str;
    sca_set_2ldoubles = int_set_2ldoubles;
    sca_get_str  = int_get_str;
    sca_hash     = int_hash;
    sca_sqrt     = larc_sca_sqrt_int;
    sca_norm     = larc_sca_norm_int;
    sca_conj     = larc_sca_set_arith;
    sca_is_real  = always_true_is_real;
#endif
#if defined(USE_REAL)
    sca_cmp      = larc_sca_cmp_real;
    sca_eq    = rc_eq;
    sca_eq_approx    = rc_eq_approx;
    sca_set_str  = larc_sca_set_str_real;
    sca_set_2ldoubles = larc_sca_2ldoubles_real;
    sca_get_str  = larc_sca_get_str_real;
    sca_hash     = hash_from_one_longdouble;
    sca_sqrt     = larc_sca_sqrt_real;
    sca_norm     = larc_sca_norm_real;
    sca_conj     = larc_sca_set_arith;
    sca_is_real  = always_true_is_real;
#endif
#if defined(USE_COMPLEX)
    sca_cmp      = larc_sca_cmp_complex;
    sca_eq    = rc_eq;
    sca_eq_approx    = rc_eq_approx;
    sca_set_str  = larc_sca_set_str_complex;
    sca_set_2ldoubles = larc_sca_2ldoubles_complex;
    sca_get_str  = larc_sca_get_str_complex;
    sca_hash     = larc_sca_hash_complex;
    sca_sqrt     = larc_sca_sqrt_complex;
    sca_norm     = larc_sca_norm_complex;
    sca_conj     = larc_sca_conj_complex;
    sca_is_real  = larc_sca_is_real_complex;
#endif
}

// to do Boolean arithmetic, use this.
void init_boolean_scalarOps(int verbose){
    int ops_set = check_scalarOps();
    if (ops_set == -1 && verbose > SILENT)
        fprintf(stderr,"WARNING: user was previously operating without all scalar ops defined.\n");
    if (verbose>BASIC){
        if (ops_set == 1)
            fprintf(stderr,"Switching ");
        if (ops_set == 0)
            fprintf(stderr,"Setting ");
        fprintf(stderr,"scalar ops to Boolean Arithmetic for Integers.\n");
    }
    sca_nbhd_approx = larc_nbhd_approx;
#if defined(USE_INTEGER)
    sca_init     = empty1;
    sca_clear    = empty1;
    sca_set      = larc_sca_set_arith;
    sca_set_str  = bool_set_str;
    sca_set_2ldoubles = bool_set_2ldoubles;
    sca_get_str  = int_get_str;
    sca_hash     = int_hash;
    sca_add      = bool_add;
    sca_mult     = bool_mult;
    sca_divide   = bool_divide;
    // consider if setting to Boolean functions after LARC is initialized with
    // arithmetic functions. bool_cmp would consider e.g. 1 and -1 equal,
    // breaking the rule that a matrix is in LARC only once.
    // Proposition: if narrowing domain, sca_cmp should be aware of wider domain.
    // Also, using num_cmp preserves matrix memoizing, I think.
    sca_cmp      = int_cmp;
    sca_eq    = int_eq;
    sca_eq_approx    = int_eq_approx; // currently no approximation implemented
    //sca_cmp      = bool_cmp;
    sca_norm     = bool_set;
    sca_conj     = bool_set;
    sca_is_real  = always_true_is_real;
#else
    fprintf(stderr,"ERROR in %s: matrix store must be INTEGER for Boolean arithmetic.\n", __func__);
    exit(1);
#endif
}


/**********************************************************
 * Andy added this check_scalarOps to make sure that either     *
 * all the sca_ values were undefined (returns 0),        *
 * or that all the sca_ values  are defined (returns 1).  *
 * A minus one is returned if neither of these cases hold.*
 **********************************************************/
int check_scalarOps(){
  int totalOps = 18;
  // The prototypes for these are declared in larc.h.
  int not_defined = (sca_init == NULL)
                    + (sca_clear == NULL)
                    + (sca_set == NULL)
                    + (sca_set_str == NULL)
                    + (sca_set_2ldoubles == NULL)
                    + (sca_nbhd_approx == NULL)
                    + (sca_get_str == NULL)
                    + (sca_hash == NULL)
                    + (sca_add == NULL)
                    + (sca_mult == NULL)
                    + (sca_divide == NULL)
                    + (sca_cmp == NULL)
                    + (sca_eq == NULL)
                    + (sca_eq_approx == NULL)
                    + (sca_sqrt == NULL)
                    + (sca_norm == NULL)
                    + (sca_conj == NULL)
                    + (sca_is_real == NULL);
  if (not_defined == 0)
      return 1;
  else if (not_defined == totalOps)
      return 0;
  else 
  {
    if (VERBOSE>BASIC) {
      fprintf(stderr,"not_defined =  %d, totalOps = %d\n",not_defined,totalOps);
      if (sca_init == NULL) fprintf(stderr, "sca_init is NULL\n");
      if (sca_clear == NULL) fprintf(stderr, "sca_clear is NULL\n");
      if (sca_set == NULL) fprintf(stderr, "sca_set is NULL\n");
      if (sca_set_str == NULL) fprintf(stderr, "sca_set_str is NULL\n");
      if (sca_set_2ldoubles == NULL) fprintf(stderr, "sca_set_2ldoubles is NULL\n");
      if (sca_nbhd_approx == NULL) fprintf(stderr, "sca_nbhd_approx is NULL\n");
      if (sca_get_str == NULL) fprintf(stderr, "sca_get_str is NULL\n");
      if (sca_hash == NULL) fprintf(stderr, "sca_hash is NULL\n");
      if (sca_add == NULL) fprintf(stderr, "sca_add is NULL\n");
      if (sca_mult == NULL) fprintf(stderr, "sca_mult is NULL\n");
      if (sca_divide == NULL) fprintf(stderr, "sca_divide is NULL\n");
      if (sca_cmp == NULL) fprintf(stderr, "sca_cmp is NULL\n");
      if (sca_eq == NULL) fprintf(stderr, "sca_eq is NULL\n");
      if (sca_eq_approx == NULL) fprintf(stderr, "sca_eq_approx is NULL\n");
      if (sca_sqrt == NULL) fprintf(stderr, "sca_sqrt is NULL\n");
      if (sca_norm == NULL) fprintf(stderr, "sca_norm is NULL\n");
      if (sca_conj == NULL) fprintf(stderr, "sca_conj is NULL\n");
      if (sca_is_real == NULL) fprintf(stderr, "sca_is_real is NULL\n");
    }
    return -1;
  }
}



// for the user to set their own function for one of these
void define_sca_init(void func(scalarType *)){sca_init = func;}
void define_sca_clear(void func(scalarType *)){sca_clear = func;}
void define_sca_set(void func(scalarType *, const scalarType)){sca_set = func;}
void define_sca_set_str(void func(scalarType *, const char *)){sca_set_str = func;}
void define_sca_set_2ldoubles(void func(scalarType *, long double real_val, long double imag_val)){sca_set_2ldoubles = func;}
void define_sca_get_str(char *func(const scalarType)){sca_get_str = func;}
void define_sca_hash(uint64_t func(const scalarType, uint64_t)){sca_hash = func;}
void define_sca_add(void func(scalarType *, const scalarType, const scalarType)){sca_add = func;}
void define_sca_mult(void func(scalarType *, const scalarType, const scalarType)){sca_mult = func;}
void define_sca_divide(void func(scalarType *, const scalarType, const scalarType)){sca_divide = func;}
void define_sca_cmp(int func(const scalarType, const scalarType)){sca_cmp = func;}
void define_sca_eq(int func(const scalarType, const scalarType)){sca_eq = func;}
void define_sca_eq_approx(int func(const scalarType, const scalarType)){sca_eq_approx = func;}
void define_sca_sqrt(void func(scalarType *, const scalarType)){sca_sqrt = func;}
void define_sca_norm(void func(scalarType *, const scalarType)){sca_norm = func;}
void define_sca_conj(void func(scalarType *, const scalarType)){sca_conj = func;}
void define_sca_is_real(int func(const scalarType)){sca_is_real = func;}


