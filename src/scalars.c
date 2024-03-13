//                        scalars.c
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
#include <signal.h>

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
#include "spr.h"
#include "clifford.h"

/*!
 * \file scalars.c
 * \brief This file contains the routines that handle scalarType data.
 *
 * In particular, it contains routines which translate string values into
 * scalarType and vice versa, and functions which wrap basic arithmetic
 * and value comparisons for scalarType so that the underlying native type
 * is handled correctly.
 *
 */

// this is used to pass information to python about whether
// we are in MAR (1) or SPR (0) mode for locality hashing
#ifdef MAR
int MARmode = 1;
#else
int MARmode = 0;
#endif

// See scalars.h for doxygen comments on the prototypes of scalar operations.

void     (*sca_init)    (scalarType* s_ptr);
void     (*sca_clear)   (scalarType* s_ptr);
void     (*sca_set)     (scalarType* d_ptr, const scalarType s_ptr);
void     (*sca_set_str) (scalarType* s_ptr, const char *input_str);
void     (*sca_set_2ldoubles) (scalarType*, long double real_val, long double imag_val);
void     (*sca_set_enum) (scalarType*, const enum sca_constant_spec enum_val);
char    *(*sca_get_readable_approx_str) (const scalarType s);
char    *(*sca_get_exact_str) (const scalarType s);
uint64_t (*sca_hash)    (const scalarType s, uint64_t exp);
void     (*sca_add)   (scalarType* s, const scalarType s1, const scalarType s2);
void     (*sca_mult)  (scalarType* p, const scalarType s1, const scalarType s2);
void     (*sca_divide)  (scalarType* s, const scalarType n, const scalarType d);
int      (*sca_cmp)     (const scalarType a, const scalarType b); //
int      (*sca_eq)   (const scalarType a, const scalarType b); //
void     (*sca_sqrt)    (scalarType* a, const scalarType b);
void     (*sca_norm)    (scalarType* a, const scalarType b);
void     (*sca_conj)    (scalarType* a, const scalarType b);
int      (*sca_is_pure_real)  (const scalarType a);
int      (*sca_is_pure_imag)  (const scalarType a);

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
#define LDBL_MANT_DEC_DIG ((int) lroundl(log10l(powl(2, LDBL_MANT_DIG))))
#endif

#ifdef IS_MP
int io_string_base = 10;
#endif

/*****************************************
 *  Locality Approximation Functions     *
 ****************************************/
// LARC has two locality hashes, and the functions which are needed for these
// hashes to work differ. Both MAR and SPR use round_sig_fig_<scalarType>.
// SPR also uses collapse_near_zero_<scalarType> whenever zeroregionbitparam
// is less than regionbitparam-1. 
//
// SPR labels its regions by the geometric center of the region. The function
// return_SPR_region_center() is used to determine these values. (Equality
// testing uses sca_eq().)
//
// MAR labels its tiles by indexing them. The indices are calculated by
// round_sig_fig_<scalarType>, and are multiprecision integers (for real
// scalarTypes, a single index suffices, but complex numbers require two
// indices). 

void return_SPRregion_or_MARtile_label(void* ptr,
	int __attribute__((unused)) tile_exp, const scalarType input)
{
#ifdef MAR // MARmode
   get_tile_index((MAR_tile_index_t *)ptr, tile_exp, input);
#else // SPRmode
   return_SPR_region_center((scalarType *)ptr, input);
#endif // #ifdef MAR
}


/******************************************************************************
 * empty wrappers
 ******************************************************************************/
//#ifndef IS_MP
/*!
 * \ingroup larc
 * \brief This is a placeholder for any sca_* function we do not want to define
 * \param scalar a pointer to a scalarType
 */
void empty1(scalarType *scalar){;}
// only used if USE_CLIFFORD
/*!
 * \ingroup larc
 * \brief This is a placeholder for any sca_* function we do not want to define
 * \param scalar a scalarType
 * \param scalar2 another scalarType
 * \result 0
 */
int empty2(const scalarType scalar, const scalarType scalar2){return 0;}
/*!
 * \ingroup larc
 * \brief This is a placeholder for any sca_* function we do not want to define
 * \param scalar a pointer to a scalarType
 * \param x type enum sca_constant_spec, whatever that is
 */
void empty3(scalarType *scalar, const enum sca_constant_spec x){;}
/*!
 * \ingroup larc
 * \brief This is a placeholder for any sca_* function we do not want to define
 * \param scalar a pointer to a scalarType
 */
char * empty4(const scalarType scalar){return "";}
//#endif

// Hashes for multi-precision types.


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
 * Authors: Andy Wills and Jenny Zito
 */
    uint64_t ping, pong;
    // GMP integers are stored in (32 or) 64 bit words called `limbs'.
    size_t limbNum = mpz_size(sc);
    // a limb is "currently a long, but on some systems it's an int for
    // efficiency, and on some systems it will be long long in the future."
    // -GMP 6.1.2 docs.
    if (64 != mp_bits_per_limb){printf("WARNING: hash function for GMP integers not optimized for implementations with limb size %d instead of %d.\n", mp_bits_per_limb, 64);}
    mp_limb_t limb;
    int i;

    // in case there is something strange about GMP integer zero, we
    // test for it and deal with it in a standard way
    if (mpz_sgn(sc)==0) return mult_golden_hash((uint64_t)0,exponent);

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


#if defined(USE_MPREAL) || defined(USE_MPCOMPLEX)
/*!
 * \ingroup larc
 * \brief Hashes a multiprecision real to an integer in range [0,2**exponent)
 *
 *     Multiprecision floating point zero is weird, as it doesn't seem to
 *     care what's in its limbs. There is also the possibility that we have
 *     stored NaN or infinity, which may behave likewise. To avoid trouble,
 *     we use mpfr_regular_p() to determine if we are in one of these cases,
 *     and if so we just hash the 64-bit integer zero. Otherwise, we hash
 *     the sign bit, XOR the result with the first limb and hash that, and
 *     repeat XORing and hashing until all limbs are included.
 *
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

    // Multiprecision floating point zero is weird, as it doesn't seem to
    // care what's in its limbs. To avoid trouble, we hash a standard
    // integer zero value instead.
    if (mpfr_zero_p(sc)) return mult_golden_hash((uint64_t)0, exponent);

    // "regular" numbers (not NaN, infinity, or zero) are hashed normally
    if (VERBOSE>DEBUG) printf("\tlimbNum = %d\n",(int)limbNum);
    // initializing ping with the sign of sc: 0 or 1.
    ping = (mpfr_sgn(sc)<0);
    if (VERBOSE>DEBUG) printf("\tping initialized to %lu\n",ping);

    for (i = limbNum-1; i >= 0; i--){
    //for (i = 0; i < limbNum; ++i){
        pong = mult_golden_hash(ping, 64);
        limb = sc->_mpfr_d[i];
        ping = pong ^ limb;
        if (VERBOSE>DEBUG) {
            printf("in loop with i = %d:",i);
            printf(">> ping = %lu, pong = %lu\n",ping,pong);
            printf("limb = %lu\n",(long unsigned)limb);
        }
    }
    pong = mult_golden_hash(ping, exponent);
    if (VERBOSE>DEBUG) {
        printf("\tfinal value of pong is %lu\n",pong);
        printf("exiting %s\n",__func__);
    }
    return pong;
}
#endif // USE_MPREAL or USE_MPCOMPLEX


/*************************************************************
 *  These functions define the various scalarType functions  *
 *  used in LARC, so that all definitions depending on a     *
 *  particular type (eg, multiprecision integers) are hidden *
 *  behind a scalarType sca_* function. This is done to keep *
 *  all the #defined functions in scalars.c.                 *
 *************************************************************/

#ifdef USE_MPINTEGER
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
        raise(SIGSEGV);
        exit(1);
    }
    if (scratchVars.mpreal_in_use)
        fprintf(stderr,"%s reusing scratchVars.mpreal!\n",__func__);
    scratchVars.mpreal_in_use = 1;

    mpfr_set_ld(scratchVars.mpreal, real_val, MPFR_RNDN);
    mpfr_get_z(*rsc, scratchVars.mpreal, MPFR_RNDN);
    scratchVars.mpreal_in_use = 0;
}
char *larc_sca_get_str_mpinteger(const mpz_t sc){return mpz_get_str(NULL, io_string_base, sc);}
void larc_sca_add_mpinteger (mpz_t *rsc, const mpz_t sc1, const mpz_t sc2){mpz_add(*rsc, sc1, sc2);}
void larc_sca_mult_mpinteger(mpz_t *rsc, const mpz_t sc1, const mpz_t sc2){mpz_mul(*rsc, sc1, sc2);}
// NOTE: mpz_tdiv_q returns the quotient of the division (not remainder)
void larc_sca_divide_mpinteger(mpz_t *rsc, const mpz_t sc1, const mpz_t sc2){mpz_tdiv_q(*rsc, sc1, sc2);}
void larc_sca_sqrt_mpinteger(mpz_t *sqroot, const mpz_t sc){mpz_sqrt(*sqroot, sc);}
void larc_sca_norm_mpinteger(mpz_t *norm, const mpz_t sc){mpz_abs(*norm, sc);}
int  larc_sca_eq_mpinteger(const mpz_t sc1, const mpz_t sc2){return 0 == mpz_cmp(sc1, sc2);}
#endif // USE_MPINTEGER

#ifdef IS_RATIONAL
/******************************************************************************
 * mpq_t (GMP) scalar operations:
 * these wrap appropriate mpz routines and convert return pointers to regular
 * mpq_t structs
 ******************************************************************************/
void larc_sca_init_mprational(mpq_t *sc){mpq_init(*sc);}
void larc_sca_clear_mprational(mpq_t *sc){mpq_clear(*sc);}
void larc_sca_set_mprational(mpq_t *rsc, const mpq_t sc){mpq_set(*rsc, sc);}

void larc_sca_set_str_mprational(mpq_t *rsc, const char *input_str){
    mpq_set_str(*rsc, input_str, io_string_base);
    mpq_canonicalize(*rsc);
}
char *larc_sca_get_str_mprational(const mpq_t sc){return mpq_get_str(NULL, io_string_base, sc);}
void larc_sca_set_2ldoubles_mprational(mpq_t *rsc, long double real_val, long double imag_val){
    if (0.0 != imag_val) {
        fprintf(stderr,"ERROR in %s: MPRATIONAL cannot have imaginary component.\n", __func__);
        raise(SIGSEGV);
        exit(1);
    }
    // go through mpfr_t to keep full long double precision in conversion
    if (scratchVars.mpreal_in_use)
        fprintf(stderr,"%s reusing scratchVars.mpreal!\n",__func__);
    scratchVars.mpreal_in_use = 1;

    mpfr_set_ld(scratchVars.mpreal,real_val,MPFR_RNDN);
    mpfr_get_q(*rsc,scratchVars.mpreal);
    scratchVars.mpreal_in_use = 0;
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
    // go through mpfr_t to keep full long double precision in conversion
    mpfr_t screal,sqrtreal;
    mpfr_inits2(mpreal_precision,screal,sqrtreal,NULL);
    // convert mpq to mpfr
    mpfr_set_q(screal,sc,MPFR_RNDN);
    // take square root
    mpfr_sqrt(sqrtreal,screal,MPFR_RNDN);
    // convert back
    mpfr_get_q(*sqroot,sqrtreal);
    mpfr_clears(screal,sqrtreal,NULL);
}

void larc_sca_norm_mprational(mpq_t *norm, const mpq_t sc){mpq_abs(*norm, sc);}

#ifdef USE_CLIFFORD

void larc_sca_init_clifford(clifford_t *scalar)
{
   // CLIFFORD_DIMENSION is the dimension of the Clifford Algebra as
   // a vector space over the complex rationals (or real rationals).
   // For example, when we append the two roots sqrt(2) and sqrt(3)
   // to the rationals, then the dimension is 4 since 
   // 1, sqrt(2), sqrt(3), sqrt(6) is a basis.
   int i;
   for (i=0;i<CLIFFORD_DIMENSION;++i)
   {
        mpq_init((*scalar)->real_coeffs[i]);
#ifdef IS_COMPLEX
        // if we change to always initializing imaginary part (to zero for
        // real types, remove the #ifdef
        mpq_init((*scalar)->imag_coeffs[i]);
#endif // IS_COMPLEX
   }
   return;
}

void larc_sca_clear_clifford(clifford_t *scalarClif)
{
    int i;
    for (i=0;i<CLIFFORD_DIMENSION;++i)
    {
        mpq_clear((*scalarClif)->real_coeffs[i]);
#ifdef IS_COMPLEX
        // if we change the code to initialize imaginary parts to
        // zero when type is not complex, remove the #ifdef
        mpq_clear((*scalarClif)->imag_coeffs[i]);
#endif
    }
    return;
}

/*!
 * \ingroup LARC
 * \brief This function adds two Clifford type variables
 *
 * \param sum A pointer to the result of the addition
 * \param a The first addend
 * \param b The second addend
 */
static void larc_sca_add_clifford(clifford_t *sum, const clifford_t a, const clifford_t b)
{
    // Clifford algebra addition is simple - just add the ith term of
    // the two inputs and put it into the ith term of the sum
    for (int i=0;i<CLIFFORD_DIMENSION;++i) 
    {
        mpq_add((*sum)->real_coeffs[i],a->real_coeffs[i],b->real_coeffs[i]);
#ifdef IS_COMPLEX
	// do addition of imaginary parts if type is complex
        mpq_add((*sum)->imag_coeffs[i],a->imag_coeffs[i],b->imag_coeffs[i]);
#endif // IS_COMPLEX
    }
    return;
}

static void larc_sca_set_clifford_to_zero(clifford_t *rsc)
{
    for (int i=0;i<CLIFFORD_DIMENSION;++i)
    {
        mpq_set_ui((*rsc)->real_coeffs[i], 0, 1);
#ifdef IS_COMPLEX
        mpq_set_ui((*rsc)->imag_coeffs[i], 0, 1);
#endif // IS_COMPLEX
    }
    return;
}

/*!
 * \ingroup LARC
 * \brief This function multiplies two Clifford type variables
 *
 * \param product A pointer to the result of the multiplication
 * \param a The multiplicand
 * \param b The multiplier
 */
static void larc_sca_mult_clifford(clifford_t *product, const clifford_t a, const clifford_t b)
{
    // Clifford algebra multiplication is not simple. The rules for individual
    // term multiplication are contained in the array
    // clifford_term_t clifford_mult[CLIFFORD_DIMENSION][CLIFFORD_DIMENSION]
    //
    // The type clifford_term_t consists of 3 ints, the first two the numerator
    // and denominator of any constant factor resulting from multiplying two
    // terms, and the third the index of the resulting term.
    //
    // Since this function could conceivably be called with 'product' being
    // the same as either 'a' or 'b', we protect ourselves by accumulating the
    // result of the multiplication in a scratch variable.
    //
    // This function uses scratchVars.quick_use, so calling functions cannot
    // store data in this variable.
  
    // initialize and zero the variable holding the accumulating product
    if (scratchVars.quick_use_in_use)
        fprintf(stderr,"%s reusing scratchVars.quick_use!\n",__func__);
    scratchVars.quick_use_in_use = 1;

    clifford_t *temp_prod = &scratchVars.quick_use;
    larc_sca_set_clifford_to_zero(temp_prod);
    
    if (scratchVars.mprational_in_use)
        fprintf(stderr,"%s reusing scratchVars.mprational!\n",__func__);
    scratchVars.mprational_in_use = 1;
    if (scratchVars.mprational2_in_use)
        fprintf(stderr,"%s reusing scratchVars.mprational2!\n",__func__);
    scratchVars.mprational2_in_use = 1;

    mpq_t *temp_q = &scratchVars.mprational;
    mpq_t *scale_q = &scratchVars.mprational2;

    for (int i=0;i<CLIFFORD_DIMENSION;++i)
    for (int j=0;j<CLIFFORD_DIMENSION;++j)
    {
        int numer = clifford_mult[i][j].numerator;
        int denom = clifford_mult[i][j].denominator;
        int index = clifford_mult[i][j].const_index;
        mpq_set_si(*scale_q,numer,denom);

        // real-real part is always done
        mpq_mul(*temp_q,a->real_coeffs[i],b->real_coeffs[j]);
        mpq_mul(*temp_q,*temp_q,*scale_q);
        mpq_add((*temp_prod)->real_coeffs[index],
                (*temp_prod)->real_coeffs[index],*temp_q);
#ifdef IS_COMPLEX
        // if one element of the Clifford algebra is sqrt(-1), also do
        // real-imag, imag-real, imag-imag part

        //real*imag = imag
        mpq_mul(*temp_q,a->real_coeffs[i],b->imag_coeffs[j]);
        mpq_mul(*temp_q,*temp_q,*scale_q);
        mpq_add((*temp_prod)->imag_coeffs[index],
                (*temp_prod)->imag_coeffs[index],*temp_q);
        //imag*real = imag (scale_q is the same)
        mpq_mul(*temp_q,a->imag_coeffs[i],b->real_coeffs[j]);
        mpq_mul(*temp_q,*temp_q,*scale_q);
        mpq_add((*temp_prod)->imag_coeffs[index],
                (*temp_prod)->imag_coeffs[index],*temp_q);
        //imag*imag = -1*real (sign of scale_q changes)
        mpq_neg(*scale_q,*scale_q);
        mpq_mul(*temp_q,a->imag_coeffs[i],b->imag_coeffs[j]);
        mpq_mul(*temp_q,*temp_q,*scale_q);
        mpq_add((*temp_prod)->real_coeffs[index],
                (*temp_prod)->real_coeffs[index],*temp_q);
#endif // IS_COMPLEX
    } // end loop over j
    scratchVars.mprational_in_use = 0;
    scratchVars.mprational2_in_use = 0;

    // copy accumulated value into argument 'product' and return
    sca_set(product,*temp_prod);
    scratchVars.quick_use_in_use = 0;
    return;
}

static void larc_sca_set_clifford(clifford_t *rsc, const clifford_t sc)
{
    for (int i=0;i<CLIFFORD_DIMENSION;++i) 
    {
        mpq_set((*rsc)->real_coeffs[i],sc->real_coeffs[i]);
#ifdef IS_COMPLEX
        mpq_set((*rsc)->imag_coeffs[i],sc->imag_coeffs[i]);
#endif // IS_COMPLEX
    }
    return;
}

#ifdef IS_COMPLEX
static char *larc_sca_get_str_clifford_rational_complex(const clifford_t val)
{
// In mpq_get_str, "if str is NULL, the result string is allocated using the
// current allocation function (see Custom Allocation). The block will be
// strlen(str)+1 bytes, that being exactly enough for the string and
// null-terminator." -GMP6.1.2
    char *real = mpq_get_str(NULL, io_string_base_dup, val->real_coeffs[0]);
    char *imag = mpq_get_str(NULL, io_string_base_dup, val->imag_coeffs[0]);
    char *num = calloc(strlen(real) + strlen(imag) + strlen("+I*")+1, sizeof(char));
    if (num==NULL) { ALLOCFAIL(); }
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
#else
static char *larc_sca_get_str_clifford_rational_real(const clifford_t val)
{
// In mpq_get_str, "if str is NULL, the result string is allocated using the
// current allocation function (see Custom Allocation). The block will be
// strlen(str)+1 bytes, that being exactly enough for the string and
// null-terminator." -GMP6.1.2
    return mpq_get_str(NULL, io_string_base_dup, val->real_coeffs[0]);
}
#endif // #ifdef IS_COMPLEX

static char *larc_sca_get_str_clifford_full(const clifford_t sc)
{
    char *real_strs[CLIFFORD_DIMENSION];
#ifdef IS_COMPLEX
    char *imag_strs[CLIFFORD_DIMENSION];
#endif // IS_COMPLEX
    int tot_len = 0;
    for (int i = 0; i < CLIFFORD_DIMENSION; ++i) {
        real_strs[i] = mpq_get_str(NULL, io_string_base_dup, sc->real_coeffs[i]);
        tot_len += strlen(real_strs[i]) + strlen(clifford_consts[i]);
#ifdef IS_COMPLEX
        imag_strs[i] = mpq_get_str(NULL, io_string_base_dup, sc->imag_coeffs[i]);
        tot_len += strlen(imag_strs[i]) + strlen(clifford_consts[i]);
#endif // IS_COMPLEX
    }

    char *num = calloc(tot_len + 12*CLIFFORD_DIMENSION + 64, sizeof(char));
    if (num==NULL) { ALLOCFAIL(); }
    int is_first = 1;
    for (int i = 0; i < CLIFFORD_DIMENSION; ++i) {
        if (0 == mpq_cmp_ui(sc->real_coeffs[i], 0, 1)) continue;
        if (is_first) {
            is_first = 0;
        } else {
            strcat(num, "+");
        }
        strcat(num, "(");
        strcat(num, real_strs[i]);
        strcat(num, ")*{");
        strcat(num, clifford_consts[i]);
        strcat(num, "}");
    }
#ifdef IS_COMPLEX
    for (int i = 0; i < CLIFFORD_DIMENSION; ++i) {
        if (0 == mpq_cmp_ui(sc->imag_coeffs[i], 0, 1)) continue;
        if (is_first) {
            is_first = 0;
        } else {
            strcat(num, "+");
        }
        strcat(num, "(");
        strcat(num, imag_strs[i]);
        strcat(num, ")*{");
        strcat(num, clifford_consts[i]);
        strcat(num, "}*I");
    }
#endif // IS_COMPLEX
    return num;
}

static char *larc_sca_get_str_clifford(const clifford_t val)
{
    if (clifford_test_pure_rational(val)) {
#ifdef IS_COMPLEX
        return larc_sca_get_str_clifford_rational_complex(val);
#else
        return larc_sca_get_str_clifford_rational_real(val);
#endif
    } else {
        return larc_sca_get_str_clifford_full(val);
    }
}

#ifdef IS_COMPLEX
static void larc_sca_set_str_clifford_rational_complex(clifford_t *rsc, const char *input_str)
{
    char *ipos = strchr(input_str, 'I');
    if (ipos == NULL) {
        mpq_set_str((*rsc)->real_coeffs[0], input_str, io_string_base_dup);
        mpq_canonicalize((*rsc)->real_coeffs[0]);
        mpq_set_ui((*rsc)->imag_coeffs[0], 0, 1);
    } else if ((ipos[1] == '*') && ((ipos == input_str) || (ipos[-1] == '+') || (ipos[-1] == '-'))) {
        mpq_set_ui((*rsc)->real_coeffs[0], 0, 1);
        mpq_set_str((*rsc)->imag_coeffs[0], ipos + 2, io_string_base_dup);
        mpq_canonicalize((*rsc)->imag_coeffs[0]);
        if ((ipos != input_str) && (ipos[-1] == '-')) {
            mpq_neg((*rsc)->imag_coeffs[0], (*rsc)->imag_coeffs[0]);
        }
        if (ipos > input_str + 1) {
            char *realpart = calloc(ipos - input_str, sizeof(char));
            if (realpart==NULL) { ALLOCFAIL(); }
            strncpy(realpart, input_str, ipos - input_str - 1);
            realpart[ipos - input_str - 1] = '\0';
            mpq_set_str((*rsc)->real_coeffs[0], realpart, io_string_base_dup);
            mpq_canonicalize((*rsc)->real_coeffs[0]);
            free(realpart);
        }
    } else {
        fprintf(stderr,"ERROR in %s: could not convert '%s' to Clifford rational complex number.\n",
                __func__, input_str);
        exit(1);
    }
}
#else
static void larc_sca_set_str_clifford_rational_real(clifford_t *rsc, const char *input_str)
{
    mpq_set_str((*rsc)->real_coeffs[0], input_str, io_string_base_dup);
    mpq_canonicalize((*rsc)->real_coeffs[0]);
}
#endif // #ifdef IS_COMPLEX

static void larc_sca_set_str_clifford_full(clifford_t *rsc, const char *input_str)
{
/*
Format for CLIFFORD numbers:

(1/3)*{1}-(2/5)*{S2}+3*{S3}*I-(4/5)*{S6}*I
*/
    char *aux = calloc(strlen(input_str) + 1, sizeof(char));
    if (aux==NULL) { ALLOCFAIL(); }
    int term_pos = 0;
    int is_negative = 0;
    int total_length = strlen(input_str);
    while (term_pos < total_length) {

        // Determine Clifford constant.
        char *start_name_ptr = strchr(input_str + term_pos, '{');
        char *end_name_ptr = strchr(input_str + term_pos, '}');
        if ((start_name_ptr == NULL) || (end_name_ptr == NULL)) {
            fprintf(stderr,"ERROR in %s: could not convert '%s' to Clifford number.\n",
                    __func__, input_str);
            exit(1);
        }
        int start_name_pos = start_name_ptr - input_str;
        int end_name_pos = end_name_ptr - input_str;
        strncpy(aux, input_str + start_name_pos + 1, end_name_pos - start_name_pos - 1);
        aux[end_name_pos - start_name_pos - 1] = '\0';
        int coeff_index = lookup_clifford_name(aux);
        if (coeff_index == -1) {
            fprintf(stderr,"ERROR in %s: could not convert '%s' to Clifford number.\n",
                    __func__, input_str);
            exit(1);
        }

        // Determine real or imaginary
        int is_imaginary = 0;
#ifdef IS_COMPLEX
        if (   ((end_name_pos + 2) < total_length)
            && (input_str[end_name_pos+1] == '*')
            && (input_str[end_name_pos+2] == 'I') ) {
            is_imaginary = 1;
        }
#endif

        // Determine coefficient of Clifford constant.
        if (input_str[start_name_pos-1] != '*') {
            fprintf(stderr,"ERROR in %s: could not convert '%s' to Clifford number.\n",
                    __func__, input_str);
            exit(1);
        }
        int start_coeff_pos = term_pos;
        int end_coeff_pos = start_name_pos - 2;
        if (   (input_str[start_coeff_pos] == '(')
            && (input_str[end_coeff_pos]   == ')') ) {
            ++start_coeff_pos;
            --end_coeff_pos;
        }
        strncpy(aux, input_str + start_coeff_pos, end_coeff_pos - start_coeff_pos + 1);
        aux[end_coeff_pos - start_coeff_pos + 1] = '\0';
#ifdef IS_COMPLEX
        if (is_imaginary) {
            mpq_set_str((*rsc)->imag_coeffs[coeff_index], aux, io_string_base_dup);
            if (is_negative) {
                mpq_neg((*rsc)->imag_coeffs[coeff_index], (*rsc)->imag_coeffs[coeff_index]);
            }
        } else {
            mpq_set_str((*rsc)->real_coeffs[coeff_index], aux, io_string_base_dup);
            if (is_negative) {
                mpq_neg((*rsc)->real_coeffs[coeff_index], (*rsc)->real_coeffs[coeff_index]);
            }
        }
#else
        mpq_set_str((*rsc)->real_coeffs[coeff_index], aux, io_string_base_dup);
        if (is_negative) {
            mpq_neg((*rsc)->real_coeffs[coeff_index], (*rsc)->real_coeffs[coeff_index]);
        }
#endif // #ifdef IS_COMPLEX

        // Advance to next term.
        term_pos = end_name_pos + (2 * is_imaginary) + 1;
        if (term_pos < total_length) {
            if (input_str[term_pos] == '+') {
                is_negative = 0;
                ++term_pos;
            } else if (input_str[term_pos] == '-') {
                is_negative = 1;
                ++term_pos;
            } else {
                fprintf(stderr,"ERROR in %s: could not convert '%s' to Clifford number.\n",
                        __func__, input_str);
                exit(1);
            }
        }
    }
    free(aux);
}

static void larc_sca_set_str_clifford(clifford_t *rsc, const char *input_str)
{
    larc_sca_set_clifford_to_zero(rsc);

    if (NULL == strchr(input_str, '{')) {
#ifdef IS_COMPLEX
        return larc_sca_set_str_clifford_rational_complex(rsc, input_str);
#else
        return larc_sca_set_str_clifford_rational_real(rsc, input_str);
#endif
    } else {
        return larc_sca_set_str_clifford_full(rsc, input_str);
    }
}

void sca_set_enum_via_mpfr (scalarType *output, const enum sca_constant_spec enum_val);

static void sca_set_enum_clifford(clifford_t *scalar, const enum sca_constant_spec enum_val)
{
    larc_sca_set_clifford_to_zero(scalar);

    int cliff_index;
    if (enum_val == SCALAR_ENUM_SQRT2) {
        cliff_index = lookup_clifford_name("S2");
        if (cliff_index == -1) {
            sca_set_enum_via_mpfr(scalar, enum_val);
        } else {
            mpq_set_ui((*scalar)->real_coeffs[cliff_index], 1, 1);
        }
    } else if (enum_val == SCALAR_ENUM_INV_SQRT2) {
        cliff_index = lookup_clifford_name("S2");
        if (cliff_index == -1) {
            sca_set_enum_via_mpfr(scalar, enum_val);
        } else {
            mpq_set_ui((*scalar)->real_coeffs[cliff_index], 1, 2);
        }
    } else if (enum_val == SCALAR_ENUM_SQRT3) {
        cliff_index = lookup_clifford_name("S3");
        if (cliff_index == -1) {
            sca_set_enum_via_mpfr(scalar, enum_val);
        } else {
            mpq_set_ui((*scalar)->real_coeffs[cliff_index], 1, 1);
        }
    } else if (enum_val == SCALAR_ENUM_INV_SQRT3) {
        cliff_index = lookup_clifford_name("S3");
        if (cliff_index == -1) {
            sca_set_enum_via_mpfr(scalar, enum_val);
        } else {
            mpq_set_ui((*scalar)->real_coeffs[cliff_index], 1, 3);
        }
    } else if (enum_val == SCALAR_ENUM_SQRT6) {
        cliff_index = lookup_clifford_name("S6");
        if (cliff_index == -1) {
            sca_set_enum_via_mpfr(scalar, enum_val);
        } else {
            mpq_set_ui((*scalar)->real_coeffs[cliff_index], 1, 1);
        }
    } else if (enum_val == SCALAR_ENUM_INV_SQRT6) {
        cliff_index = lookup_clifford_name("S6");
        if (cliff_index == -1) {
            sca_set_enum_via_mpfr(scalar, enum_val);
        } else {
            mpq_set_ui((*scalar)->real_coeffs[cliff_index], 1, 6);
        }
    } else if (enum_val == SCALAR_ENUM_CUBERT2) {
        cliff_index = lookup_clifford_name("C2");
        if (cliff_index == -1) {
            sca_set_enum_via_mpfr(scalar, enum_val);
        } else {
            mpq_set_ui((*scalar)->real_coeffs[cliff_index], 1, 1);
        }
    } else if (enum_val == SCALAR_ENUM_INV_CUBERT2) {
        cliff_index = lookup_clifford_name("C4");
        if (cliff_index == -1) {
            sca_set_enum_via_mpfr(scalar, enum_val);
        } else {
            mpq_set_ui((*scalar)->real_coeffs[cliff_index], 1, 2);
        }
    } else if (enum_val == SCALAR_ENUM_CUBERT4) {
        cliff_index = lookup_clifford_name("C4");
        if (cliff_index == -1) {
            sca_set_enum_via_mpfr(scalar, enum_val);
        } else {
            mpq_set_ui((*scalar)->real_coeffs[cliff_index], 1, 1);
        }
    } else if (enum_val == SCALAR_ENUM_INV_CUBERT4) {
        cliff_index = lookup_clifford_name("C2");
        if (cliff_index == -1) {
            sca_set_enum_via_mpfr(scalar, enum_val);
        } else {
            mpq_set_ui((*scalar)->real_coeffs[cliff_index], 1, 2);
        }
    } else {
        fprintf(stderr, "ERROR: Unrecognized scalar enum %d.\n", (int) enum_val);
    }
}

static int larc_sca_cmp_clifford(const clifford_t val1, const clifford_t val2)
{
    // In general, will need to convert both Clifford values to approximate
    // all-nonalgebraic form and compare those values.
    //
    // When the Clifford algebra is not complex, the result is a
    // straightforward comparison of two rational numbers. In the case of
    // a complex Clifford algebra, the returned value is determined by code
    // copied from the function complex_comparison_logic(real_cmp,imag_cmp).

    clifford_t temp1, temp2;
    larc_sca_init_clifford(&temp1);
    larc_sca_init_clifford(&temp2);
    //
    convert_clifford_to_mprational_approx(temp1, val1);
    convert_clifford_to_mprational_approx(temp2, val2);
    //
    int real_cmp = mpq_cmp(temp1->real_coeffs[0],temp2->real_coeffs[0]);
#ifndef IS_COMPLEX
    larc_sca_clear_clifford(&temp1);
    larc_sca_clear_clifford(&temp2);
    return real_cmp;
#else
    static int warned = 0;
    int retval = 777;
    int imag_cmp = mpq_cmp(temp1->imag_coeffs[0],temp2->imag_coeffs[0]);
    // equality is testable
    if ( !real_cmp && !imag_cmp ) retval=0;
    // if both values real, comparable
    if ( !imag_cmp && (mpq_cmp_ui(temp1->imag_coeffs[0],0,1)==0) )
        retval=real_cmp;
    // if both values imaginary, comparable
    if ( !real_cmp && (mpq_cmp_ui(temp1->real_coeffs[0],0,1)==0) )
        retval=imag_cmp;
    // otherwise not comparable; default to cartesian heuristic.
    if (retval==777 && !warned) {
        warned = 1;
        fprintf(stderr,"in %s: No natural way to compare complex\n",__func__);
        fprintf(stderr,"rational numbers ");
        mpq_out_str(stderr,10,temp1->real_coeffs[0]);
        fprintf(stderr," + I*");
        mpq_out_str(stderr,10,temp1->imag_coeffs[0]);
        fprintf(stderr,"\nand ");
        mpq_out_str(stderr,10,temp2->real_coeffs[0]);
        fprintf(stderr," + I*");
        mpq_out_str(stderr,10,temp2->imag_coeffs[0]);
        fprintf(stderr,"\nWe will use an arbitrary heuristic where the real");
        fprintf(stderr,"\npart comparison takes precedence.\n");
    }
    larc_sca_clear_clifford(&temp1);
    larc_sca_clear_clifford(&temp2);

    if (retval != 777) return retval;
    else if (real_cmp) return real_cmp;
    else return imag_cmp;
#endif
}

static int larc_sca_eq_clifford(const clifford_t val1, const clifford_t val2)
{
    for (int i=0; i<CLIFFORD_DIMENSION; ++i)
    {
        if (!mpq_equal(val1->real_coeffs[i],val2->real_coeffs[i])) return 0;
    }
#ifdef IS_COMPLEX
    for (int i=0; i<CLIFFORD_DIMENSION; ++i)
    {
        if (!mpq_equal(val1->imag_coeffs[i],val2->imag_coeffs[i])) return 0;
    }
#endif
    return 1;
}

static void larc_sca_2ldoubles_clifford(clifford_t *ret, long double real_val,
                                 long double imag_val)
{
    // This function only sets the non-algebraic (i.e., rational) part of
    // the Clifford type. 
#ifndef IS_COMPLEX
    if (fpclassify(imag_val) != FP_ZERO) {
        fprintf(stderr,"ERROR in %s: non-complex CLIFFORD cannot have imaginary component.\n", __func__);
        exit(1);
    }
#endif // not IS_COMPLEX
    larc_sca_set_clifford_to_zero(ret);
    // go through mpfr_t to keep full long double precision in conversion
    if (scratchVars.mpreal_in_use)
        fprintf(stderr,"%s reusing scratchVars.mpreal!\n",__func__);
    scratchVars.mpreal_in_use = 1;

    mpfr_set_ld(scratchVars.mpreal,real_val,MPFR_RNDN);
    mpfr_get_q((*ret)->real_coeffs[0],scratchVars.mpreal);
#ifdef IS_COMPLEX
    mpfr_set_ld(scratchVars.mpreal,imag_val,MPFR_RNDN);
    mpfr_get_q((*ret)->imag_coeffs[0],scratchVars.mpreal);
#endif
    scratchVars.mpreal_in_use = 0;
}

static uint64_t larc_sca_hash_clifford(const clifford_t sc, uint64_t exponent)
{
    // we must hash CLIFFORD_DIMENSION mpq_t values (or twice that for
    // complex Clifford types). We hash the numerators and denominators
    // separately, adding another factor of two to the values to be hashed.
    int num_coeffs = 2*CLIFFORD_DIMENSION;
#ifdef IS_COMPLEX
    num_coeffs *= 2;
#endif
    uint64_t *hash_list = (uint64_t *)malloc(num_coeffs*sizeof(uint64_t));
    if (hash_list == NULL) { ALLOCFAIL(); }

    int j = 0;
    for (int i=0; i<CLIFFORD_DIMENSION; ++i)
    {
        hash_list[j++] = larc_mpz_hash(mpq_numref(sc->real_coeffs[i]),exponent);
        hash_list[j++] = larc_mpz_hash(mpq_denref(sc->real_coeffs[i]),exponent);
#ifdef IS_COMPLEX
        hash_list[j++] = larc_mpz_hash(mpq_numref(sc->imag_coeffs[i]),exponent);
        hash_list[j++] = larc_mpz_hash(mpq_denref(sc->imag_coeffs[i]),exponent);
#endif
    }

    uint64_t result = recursive_hash_from_int_list(hash_list,
                                                   num_coeffs,exponent);
    free(hash_list);
    return result;
}

// in scalars.c: set sca_add = clifford_add_scalar
// in scalars.c: set sca_mult = clifford_mult_scalar

static void larc_sca_divide_clifford(clifford_t *rsc, const clifford_t sc1,
                              const clifford_t sc2)
{
    // This function puts the value sc1/sc2 into rsc. Since this function could
    // conceivably be called with 'rsc' being the same as sc1 or sc2, we do
    // not overwrite rsc until the last step of the algorithm.
    //
    // To divide in Clifford algebras, we implicitly express the denominator
    // sc2 in a fractional form sc2 == A/B such that the expression for B is
    // just a rational number (i.e., the coefficients of the algebraic terms
    // are all zero). Since a rational number is just a fraction, we have no
    // problem finding its inverse, so we set rsc = sc1*A*(1/B).
    // 
    // For any x \in {\mathbb Q}[...], there is a 'conjugate' term y which has
    // the property that x*y \in {\mathbb Q}. More precisely, if we express x
    // in terms of a basis of its coefficients (a,b,...,\omega), there is a 
    // function Y(a,b,...,\omega) which produces the desired value y. Each
    // different Clifford algebra has its own such function. In terms of the
    // A and B we used above, we find that B = x*y and A = y.
    //
    // This function can use up to four intermediate values at the same time.
    // To avoid confusion we always assign as follows:
    //        clifford_t *termA = &(scratchVars.clifford);
    //        clifford_t *termB = &(scratchVars.misc);
    //        clifford_t *y0 = &(scratchVars.calc_conj); // temporary, so OK
    //        clifford_t *y1 = &(scratchVars.approx_value)
    //        clifford_t *inv_sc2 = &(scratchVars.calc_conj);
    //
    // We note that the Clifford multiply uses scratchVars.quick_use, so it is
    // unavailable for us here.
    //

    int debug = 0;

    if (debug) 
    {
        printf("input sc1 is %s\n",larc_sca_get_str_clifford(sc1));
        printf("input sc2 is %s\n",larc_sca_get_str_clifford(sc2));
    }

    if (scratchVars.clifford_in_use)
        fprintf(stderr,"%s reusing scratchVars.clifford!\n",__func__);
    scratchVars.clifford_in_use = 1;
    if (scratchVars.misc_in_use)
        fprintf(stderr,"%s reusing scratchVars.misc!\n",__func__);
    scratchVars.misc_in_use = 1;

    clifford_t *termA = &(scratchVars.clifford);
    clifford_t *termB = &(scratchVars.misc);

    if (testForZero(sc2))
    {
        printf("in %s, divisor is zero!\n",__func__);
        exit(0);
    }
    // calculate x*y for the denominator x(==sc2)
    // this often (but not always) involves calculating y first
#ifdef MPR_CLIFFORD_S2
    // for the algebra Q[sqrt2]:
    //   x = a*1 + b*s2
    //   -> y = a*1 - b*s2, yielding x*y = (a^2 - 2*b^2)*1
    //
    // calculate y, put into termA

    mpq_set((*termA)->real_coeffs[0],sc2->real_coeffs[0]);
    mpq_neg((*termA)->real_coeffs[1],sc2->real_coeffs[1]);
    if (debug) printf("termA is %s\n",larc_sca_get_str_clifford(*termA));
    // calculate x*y, put into termB
    sca_mult(termB,sc2,*termA);
    if (debug) printf("termB is %s\n",larc_sca_get_str_clifford(*termB));

#elif defined(MPC_CLIFFORD_S2)
    // for the algebra Q[i,sqrt2]
    // x = (a*1 + b*s2) + i*(c*1 + d*s2) =  a*1+b*s2+c*i+d*i*s2
    //
    // -> y = (a*1 - b*s2) + i*(c*1 - d*s2)
    // 		* (a*1 + b*s2) - i*(c*1 + d*s2)
    // 		* (a*1 - b*s2) - i*(c*1 - d*s2)
    // 
    // an alternative is to eliminate the (odd powers of) s2 and i in turn
    // in which case we get
    //
    // -> y1 = (a*1+b*s2) - i*(c*1+d*s2)
    //      yielding x*y1 = (a*1+b*s2)^2 + (c*1 + d*s2)^2
    //                   = (a^2+2ab*s2+2b^2) + (c^2+2cd*s2+2d^2)
    //                   = (a^2+2b^2+c^2+2d^2)*1 + 2(ab+cd)*s2
    //                   = A + B*s2
    // -> y2 = A - B*s2, yielding x*y1*y2 = (A^2 - 2*B^2)*1
    //
    // In this form, termB should end up with x*y1*y2 and termA with y1*y2
    //
    // This requires fewer variable manipulations.

    if (scratchVars.approx_value_in_use)
        fprintf(stderr,"%s reusing scratchVars.approx_value!\n",__func__);
    scratchVars.approx_value_in_use = 1;

    clifford_t *y1 = &(scratchVars.approx_value);
    // calculate y1
    sca_conj(y1,sc2);
    // calculate x*y1, using termB as temporary location
    sca_mult(termB,*y1,sc2);
    // calculate y2, using termA as temporary location
    sca_set(termA,*termB);    
    mpq_neg((*termA)->real_coeffs[1],(*termA)->real_coeffs[1]);
    // calculate x*y1*y2 = x*y, writing into termB
    sca_mult(termB,*termA,*termB);
    if (debug) printf("termB is %s\n",larc_sca_get_str_clifford(*termB));
    // calculate y1*y2, writing into termA
    sca_mult(termA,*termA,*y1);
    if (debug) printf("termA is %s\n",larc_sca_get_str_clifford(*termA));
    // y1 = (clifford_t *)0; // do not use again!
    scratchVars.approx_value_in_use = 0;

#elif defined(MPR_CLIFFORD_S2_S3)
    // for the algebra Q[sqrt2,sqrt3]:
    //   x = a*1 + b*s2 + c*s3 + d*s6
    // -> y = (a*1 + b*s2 - c*s3 - d*s6)
    //        * (a*1 - b*s2 + c*s3 - d*s6)
    //        * (a*1 - b*s2 - c*s3 + d*s6)
    // 
    // an alternative is to eliminate the (odd powers of) s2 and s3 in turn,
    // in which case we get
    //
    // -> y1 = (a*1+b*s2-c*s3-d*s6)
    //      yielding x*y1 = (a*1+b*s2+c*s3+d*s6)(a*1+b*s2-c*s3-d*s6)
    //      	= (a^2+2b^2-3c^2-6d^2)*1 + (2ab-6cd)*s2
    //                   = A + B*s2
    // -> y2 = A - B*s2, yielding x*y1*y2 = (A^2 - 2*B^2)*1
    //
    // In this form, termB should end up with x*y1*y2 and termA with y1*y2
    //
    // This requires fewer variable manipulations.

    if (scratchVars.approx_value_in_use)
        fprintf(stderr,"%s reusing scratchVars.approx_value!\n",__func__);
    scratchVars.approx_value_in_use = 1;

    clifford_t *y1 = &(scratchVars.approx_value);
    // calculate y1
    sca_set(y1,sc2);
    mpq_neg((*y1)->real_coeffs[2],(*y1)->real_coeffs[2]);
    mpq_neg((*y1)->real_coeffs[3],(*y1)->real_coeffs[3]);
    // calculate x*y1, using termB as temporary location
    sca_mult(termB,*y1,sc2); // x*y1
    // calculate y2, using termA as temporary location
    sca_set(termA,*termB);
    mpq_neg((*termA)->real_coeffs[1],(*termA)->real_coeffs[1]);
    // calculate x*y1*y2 = x*y, writing into termB
    sca_mult(termB,*termA,*termB);
    if (debug) printf("termB is %s\n",larc_sca_get_str_clifford(*termB));
    // calculate y1*y2 = y, writing into termA
    sca_mult(termA,*termA,*y1);
    if (debug) printf("termA is %s\n",larc_sca_get_str_clifford(*termA));
    // y1 = (clifford_t *)0; // do not use again!
    scratchVars.approx_value_in_use = 0;

#elif defined(MPC_CLIFFORD_S2_S3)
    // for the algebra Q[i,sqrt2,sqrt3]:
    //   x = a*1 + b*s2 + c*s3 + d*s6 + e*i + f*i*s2 + g*i*s3 + h*i*s6
    // -> y = (a*1 + b*s2 + c*s3 + d*s6 - e*i - f*i*s2 - g*i*s3 - h*i*s6) *
    //   (a*1 + b*s2 - c*s3 - d*s6 + e*i + f*i*s2 - g*i*s3 - h*i*s6) *
    //   (a*1 + b*s2 - c*s3 - d*s6 - e*i - f*i*s2 + g*i*s3 + h*i*s6) *
    //   (a*1 - b*s2 + c*s3 - d*s6 + e*i - f*i*s2 + g*i*s3 - h*i*s6) *
    //   (a*1 - b*s2 + c*s3 - d*s6 - e*i + f*i*s2 - g*i*s3 + h*i*s6) *
    //   (a*1 - b*s2 - c*s3 + d*s6 + e*i - f*i*s2 - g*i*s3 + h*i*s6) *
    //   (a*1 - b*s2 - c*s3 + d*s6 - e*i + f*i*s2 + g*i*s3 - h*i*s6)
    //
    // We first eliminate the imaginary part, then continue as with
    // MPR_CLIFFORD_S2_S3 for the new value. We therefore have a y0, y1,
    // and y2 to calculate.
    //
    // termB should end up with x*y0*y1*y2 and termA with y0*y1*y2
    if (scratchVars.approx_value_in_use)
        fprintf(stderr,"%s reusing scratchVars.approx_value!\n",__func__);
    scratchVars.approx_value_in_use = 1;
    if (scratchVars.calc_conj_in_use)
        fprintf(stderr,"%s reusing scratchVars.calc_conj!\n",__func__);
    scratchVars.calc_conj_in_use = 1;

    clifford_t *y0 = &(scratchVars.calc_conj);
    clifford_t *y1 = &(scratchVars.approx_value);
    // calcluate y0
    sca_conj(y0,sc2);
    if (debug) printf("y0 is %s\n",larc_sca_get_str_clifford(*y0));
    // calculate x*y0, using termA as temporary location
    sca_mult(termA,*y0,sc2);
    if (debug) printf("x*y0 is %s\n",larc_sca_get_str_clifford(*termA));
    // calculate y1
    sca_set(y1,*termA);
    mpq_neg((*y1)->real_coeffs[2],(*y1)->real_coeffs[2]);
    mpq_neg((*y1)->real_coeffs[3],(*y1)->real_coeffs[3]);
    if (debug) printf("y1 is %s\n",larc_sca_get_str_clifford(*y1));
    // calculate x*y0*y1, using termB as temporary location
    sca_mult(termB,*termA,*y1);
    if (debug) printf("x*y0*y1 is %s\n",larc_sca_get_str_clifford(*termB));
    // calculate y2, using termA as temporary location
    sca_set(termA,*termB);
    mpq_neg((*termA)->real_coeffs[1],(*termA)->real_coeffs[1]);
    if (debug) printf("y2 is %s\n",larc_sca_get_str_clifford(*termA));
    // calculate x*y0*y1*y2 = x*y, writing into termB
    sca_mult(termB,*termA,*termB);
    if (debug) printf("termB is %s\n",larc_sca_get_str_clifford(*termB));
    // calculate y0*y1*y2 = y, writing into termA
    sca_mult(y1,*y0,*y1);
    sca_mult(termA,*termA,*y1);
    if (debug) printf("termA is %s\n",larc_sca_get_str_clifford(*termA));
    // y0 = (clifford_t *)0; // do not use again!
    // y1 = (clifford_t *)0; // do not use again!
    scratchVars.approx_value_in_use = 0;
    scratchVars.calc_conj_in_use = 0;

#elif defined(MPR_CLIFFORD_C2)
    // for the algebra Q[cuberoot2]:
    // x = a*1 + b*c2 + c*c4
    // -> y = (a*1+w*b*c2+w^2*c*c4)*(a*1+w^2*b*c2 + w*c*c4)
    // where w = exp(2*i*pi/3)
    //
    // To avoid imaginary terms, we work out what y is symbolically
    //    using the identity 1 + w + w^2 = 0:
    // y = (a^2 + 2*b*c*w^2 + 2*b*c*w)
    //   + (a*b*w^2 + b*a*w + 2*c^2)*C2
    //   + (a*c*w + b^2 + c*a*w^2)*C4
    //   = (a^2 - 2*b*c) + (2*c^2 - a*b)*C2 + (b^2 - a*c)*C4
    
    if (scratchVars.mprational_in_use)
        fprintf(stderr,"%s reusing scratchVars.mprational!\n",__func__);
    scratchVars.mprational_in_use = 1;
    if (scratchVars.mprational2_in_use)
        fprintf(stderr,"%s reusing scratchVars.mprational2!\n",__func__);
    scratchVars.mprational2_in_use = 1;

    mpq_t *fac = &(scratchVars.mprational);
    mpq_t *temp = &(scratchVars.mprational2);

    // calculate each coefficient of y into termA
    mpq_mul((*termA)->real_coeffs[0],sc2->real_coeffs[0],sc2->real_coeffs[0]);
    mpq_set_si(*fac,-2,1);
    mpq_mul(*temp,sc2->real_coeffs[1],sc2->real_coeffs[2]);
    mpq_mul(*temp,*fac,*temp);
    mpq_add((*termA)->real_coeffs[0],(*termA)->real_coeffs[0],*temp);

    mpq_set_si(*fac,2,1);
    mpq_mul(*temp,sc2->real_coeffs[2],sc2->real_coeffs[2]);
    mpq_mul((*termA)->real_coeffs[1],*fac,*temp);
    mpq_mul(*temp,sc2->real_coeffs[0],sc2->real_coeffs[1]);
    mpq_neg(*temp,*temp);
    mpq_add((*termA)->real_coeffs[1],(*termA)->real_coeffs[1],*temp);

    mpq_mul((*termA)->real_coeffs[2],sc2->real_coeffs[1],sc2->real_coeffs[1]);
    mpq_mul(*temp,sc2->real_coeffs[0],sc2->real_coeffs[2]);
    mpq_neg(*temp,*temp);
    mpq_add((*termA)->real_coeffs[2],(*termA)->real_coeffs[2],*temp);
    if (debug) printf("termA is %s\n",larc_sca_get_str_clifford(*termA));

    // put x*y into termB
    sca_mult(termB,*termA,sc1);
    if (debug) printf("termB is %s\n",larc_sca_get_str_clifford(*termB));

    fac = (mpq_t *)0; // do not use again!
    temp = (mpq_t *)0; // do not use again!
    scratchVars.mprational_in_use = 0;
    scratchVars.mprational2_in_use = 0;

#elif defined(MPC_CLIFFORD_C2)
    // for the algebra Q[i,cuberoot2]:
    // x = a*1 * b*c2 + c*c4 + d*i + e*i*c2 + f*i*c4
    //
    // We first eliminate the imaginary part with a y1, then continue as with
    // MPR_CLIFFORD_C2 for y2
    //
    // termB should end up with x*y1*y2 and termA with y1*y2

    if (scratchVars.mprational_in_use)
        fprintf(stderr,"%s reusing scratchVars.mprational!\n",__func__);
    scratchVars.mprational_in_use = 1;
    if (scratchVars.mprational2_in_use)
        fprintf(stderr,"%s reusing scratchVars.mprational2!\n",__func__);
    scratchVars.mprational2_in_use = 1;
    if (scratchVars.approx_value_in_use)
        fprintf(stderr,"%s reusing scratchVars.approx_value!\n",__func__);
    scratchVars.approx_value_in_use = 1;

    mpq_t *fac = &(scratchVars.mprational);
    mpq_t *temp = &(scratchVars.mprational2);

    clifford_t *y1 = &(scratchVars.approx_value);
    // calculate y1
    sca_conj(y1,sc2);
    // calculate x*y1, using termB as temporary location
    sca_mult(termB,*y1,sc2);

    // calculate each coefficient of y2 into termA (temporary location)
    mpq_mul((*termA)->real_coeffs[0],(*termB)->real_coeffs[0],
             (*termB)->real_coeffs[0]);
    mpq_set_si(*fac,-2,1);
    mpq_mul(*temp,(*termB)->real_coeffs[1],(*termB)->real_coeffs[2]);
    mpq_mul(*temp,*fac,*temp);
    mpq_add((*termA)->real_coeffs[0],(*termA)->real_coeffs[0],*temp);

    mpq_set_si(*fac,2,1);
    mpq_mul(*temp,(*termB)->real_coeffs[2],(*termB)->real_coeffs[2]);
    mpq_mul((*termA)->real_coeffs[1],*fac,*temp);
    mpq_mul(*temp,(*termB)->real_coeffs[0],(*termB)->real_coeffs[1]);
    mpq_neg(*temp,*temp);
    mpq_add((*termA)->real_coeffs[1],(*termA)->real_coeffs[1],*temp);

    mpq_mul((*termA)->real_coeffs[2],(*termB)->real_coeffs[1],
             (*termB)->real_coeffs[1]);
    mpq_mul(*temp,(*termB)->real_coeffs[0],(*termB)->real_coeffs[2]);
    mpq_neg(*temp,*temp);
    mpq_add((*termA)->real_coeffs[2],(*termA)->real_coeffs[2],*temp);

    // .... put y1*y2 into termA
    sca_mult(termA,*termA,*y1);
    if (debug) printf("termA is %s\n",larc_sca_get_str_clifford(*termA));
    // .... put x*y1*y2 into termB
    sca_mult(termB,*termA,sc2);
    if (debug) printf("termB is %s\n",larc_sca_get_str_clifford(*termB));

    y1 = (clifford_t *)0; // do not use again!
    fac = (mpq_t *)0; // do not use again!
    temp = (mpq_t *)0; // do not use again!
    scratchVars.approx_value_in_use = 0;
    scratchVars.mprational_in_use = 0;
    scratchVars.mprational2_in_use = 0;

#endif // which clifford

    // we now have what we need to calculate 1/sc2 = termB/termA

    // calculate 1/termA knowing that real_coeffs[i] is zero for i>0
    // (and imag_coeffs[i], if it exists, is zero for all i)
    mpq_inv((*termB)->real_coeffs[0],(*termB)->real_coeffs[0]);

    if (debug) printf("termB inverted is %s\n",larc_sca_get_str_clifford(*termB));
    // multiply termB by (1/termA) to get 1/sc2
    if (scratchVars.calc_conj_in_use)
        fprintf(stderr,"%s reusing scratchVars.calc_conj!\n",__func__);
    scratchVars.calc_conj_in_use = 1;

    clifford_t *inv_sc2 = &(scratchVars.calc_conj);
    sca_mult(inv_sc2,*termA,*termB);
    if (debug) printf("the inverse of sc2 is %s\n",larc_sca_get_str_clifford(*inv_sc2));
    scratchVars.clifford_in_use = 0;
    scratchVars.misc_in_use = 0;

    // get answer = sc1/sc2 = sc1*[1/sc2]
    sca_mult(rsc,sc1,*inv_sc2);
    if (debug) printf("sc1/sc2 is %s\n",larc_sca_get_str_clifford(*rsc));
    scratchVars.calc_conj_in_use = 0;
}

static void larc_sca_sqrt_clifford(clifford_t *sqroot, const clifford_t sc)
{
//    printf("clifford square root is not yet implemented\n");
//    exit(0);

    convert_clifford_to_mprational_approx(*sqroot,sc);
#ifndef IS_COMPLEX
    larc_sca_sqrt_mprational(&((*sqroot)->real_coeffs[0]),
                               (*sqroot)->real_coeffs[0]);
#else
    // modified from larc_sca_sqrt_mpratcomplex
    // should consider having global mpc_t scratchVar for complex types
    mpc_t sc_mpc, root_mpc;
    mpc_init2(sc_mpc,mpreal_precision);
    mpc_init2(root_mpc,mpreal_precision);
    mpc_set_q_q(sc_mpc,sc->real_coeffs[0],sc->imag_coeffs[0],MPC_RNDNN);
    mpc_sqrt(root_mpc,sc_mpc,MPC_RNDNN);
    mpfr_get_q((*sqroot)->real_coeffs[0],mpc_realref(root_mpc));
    mpfr_get_q((*sqroot)->imag_coeffs[0],mpc_imagref(root_mpc));
    mpc_clear(sc_mpc);
    mpc_clear(root_mpc);
#endif
    return;
}

static void larc_sca_norm_clifford(clifford_t *norm, const clifford_t sc)
{
#ifdef IS_COMPLEX
    // When finding the norm of complex Clifford scalars, e.g.
    //             | a+c*sqrt(2) + i*(b+d*sqrt(2)) |,
    // there's not an obvious way to preserve the algebraic values, so we do
    // not try. The input value is multiplied by its complex conjugate (so that
    // the result has no imaginary part). This real Clifford value is passed
    // to our routine larc_sca_sqrt_clifford. which first approximates the
    // Clifford as a multiprecision rational, then converts it to a floating
    // point number with 'sufficient' precision, takes the square root of this
    // float, then converts back to a rational value. The end result is a
    // Clifford variable with all algebraic terms zero (this includes the
    // imaginary part as well). The norm is necessarily approximate.

    // If input is zero, so is the norm (avoid divide-by-zero in algorithm)
    if (testForZero(sc))
    {
        sca_set(norm,sc);
        return;
    }
    clifford_t *temp = norm;
//    printf("input: %s\n",sca_get_readable_approx_str(sc));
    sca_conj(temp,sc);
//    printf("conjugate: %s\n",sca_get_readable_approx_str(*temp));
    sca_mult(temp,*temp,sc);
//    printf("product: %s\n",sca_get_readable_approx_str(*temp));
    sca_sqrt(norm,*temp);
//    printf("norm: %s\n",sca_get_readable_approx_str(*temp));
#else
    // Calculation of the "absolute value" of a non-complex Clifford scalar is
    // not trivial; for example, the number x + y*sqrt(2) should be nearly
    // zero when (the rational value) x is sqrt(2)+\epsilon and y is -1, but
    // the sign of \epsilon determines whether the sum is positive or negative.
    // We calculate a good approximation for the scalar that has zero
    // coefficients for the algebraic terms, then use the sign of that
    // approximate value to see whether to leave the scalar unchanged or to
    // change the signs of all its terms. Assuming that the approximation's
    // sign is correct, the returned scalar norm is exact.
    if (scratchVars.clifford_in_use)
        fprintf(stderr,"%s reusing scratchVars.clifford!\n",__func__);
    scratchVars.clifford_in_use = 1;

    clifford_t *temp = &scratchVars.clifford;
    convert_clifford_to_mprational_approx(*temp, sc);
    int neg_flag = mpq_cmp_ui((*temp)->real_coeffs[0],0,1);
    scratchVars.clifford_in_use = 0;

    if (neg_flag >= 0)
        for (int i=0; i<CLIFFORD_DIMENSION; ++i)
            mpq_set((*norm)->real_coeffs[i],sc->real_coeffs[i]);
    else
        for (int i=0; i<CLIFFORD_DIMENSION; ++i)
            mpq_neg((*norm)->real_coeffs[i],sc->real_coeffs[i]);
#endif
}

static void larc_sca_conj_clifford(clifford_t *conj_out, const clifford_t sc)
{
#ifdef IS_COMPLEX
    for (int i=0; i<CLIFFORD_DIMENSION; ++i)
    {
        mpq_set((*conj_out)->real_coeffs[i],sc->real_coeffs[i]);
        mpq_neg((*conj_out)->imag_coeffs[i],sc->imag_coeffs[i]);
    }
#else
    sca_set(conj_out,sc);
#endif
    return;
}

static int larc_sca_is_real_clifford(const clifford_t sc)
{
#ifdef IS_COMPLEX
    for (int i=0; i<CLIFFORD_DIMENSION; ++i)
        if (mpq_cmp_ui(sc->imag_coeffs[i],0,1)) return 0;
#endif // #ifdef IS_COMPLEX
    return 1;
}

static int larc_sca_is_imag_clifford(const clifford_t sc)
{
    for (int i=0; i<CLIFFORD_DIMENSION; ++i)
        if (mpq_cmp_ui(sc->real_coeffs[i],0,1)) return 0;
    return 1;
}

#endif // USE_CLIFFORD
#endif // IS_RATIONAL

#ifdef USE_MPRATCOMPLEX
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
    // go through mpfr_t to keep full long double precision in conversion
    if (scratchVars.mpreal_in_use)
        fprintf(stderr,"%s reusing scratchVars.mpreal!\n",__func__);
    scratchVars.mpreal_in_use = 1;

    mpfr_set_ld(scratchVars.mpreal,real_val,MPFR_RNDN);
    mpfr_get_q((*rsc)->real,scratchVars.mpreal);
    mpfr_set_ld(scratchVars.mpreal,imag_val,MPFR_RNDN);
    mpfr_get_q((*rsc)->imag,scratchVars.mpreal);
    scratchVars.mpreal_in_use = 0;
}

int larc_sca_eq_mpratcomplex(const larc_mpratcomplex_t val1, const larc_mpratcomplex_t val2){
    return (   mpq_equal(val1->real, val2->real)
            && mpq_equal(val1->imag, val2->imag) );
}

int larc_sca_cmp_mpratcomplex(const larc_mpratcomplex_t val1,
                              const larc_mpratcomplex_t val2)
{
    // In general, two complex numbers have no natural ordering in the way
    // that two integers or real numbers have. We produce an ordering for the
    // limited cases where ordering makes sense (e.g., two real numbers which
    // are stored in a variable of complex type), and also return 0 when the
    // two complex numbers are equal. Currently, we warn the user when they
    // try to compare numbers that are not naturally orderable, and use an
    // arbitary method to make the comparison.
    static int warned = 0;
    int real_cmp = mpq_cmp(val1->real, val2->real);
    int imag_cmp = mpq_cmp(val1->imag, val2->imag);

    // equality is testable
    if ( !real_cmp && !imag_cmp ) return 0;
    // if both values real, comparable
    if ( !imag_cmp && (mpq_cmp_ui(val1->imag,0,1)==0) )
        return real_cmp;
    // if both values imaginary, comparable
    if ( !real_cmp && (mpq_cmp_ui(val1->real,0,1)==0) )
        return imag_cmp;
    // otherwise not comparable; default to cartesian heuristic.
    if (!warned) {
        warned = 1;
        fprintf(stderr,"in %s: No natural way to compare complex\n",__func__);
        fprintf(stderr,"rational numbers ");
        mpq_out_str(stderr,10,val1->real);
        fprintf(stderr," + I*");
        mpq_out_str(stderr,10,val1->imag);
        fprintf(stderr,"\nand ");
        mpq_out_str(stderr,10,val2->real);
        fprintf(stderr," + I*");
        mpq_out_str(stderr,10,val2->imag);
        fprintf(stderr,"\nWe will use an arbitrary heuristic where the real");
        fprintf(stderr,"\npart comparison takes precedence.\n");
    }
    if (real_cmp) return real_cmp;
    else return imag_cmp;

}

void larc_sca_add_mpratcomplex(larc_mpratcomplex_t *rsc, const larc_mpratcomplex_t sc1, const larc_mpratcomplex_t sc2){
    mpq_add((*rsc)->real, sc1->real, sc2->real);
    mpq_add((*rsc)->imag, sc1->imag, sc2->imag);
}

void larc_sca_mult_mpratcomplex(larc_mpratcomplex_t *rsc, const larc_mpratcomplex_t sc1, const larc_mpratcomplex_t sc2){
    mpq_t temp_r1_r2, temp_r1_i2, temp_i1_r2, temp_i1_i2;
    mpq_inits(temp_r1_r2, temp_r1_i2, temp_i1_r2, temp_i1_i2, NULL);
    mpq_mul(temp_r1_r2, sc1->real, sc2->real);
    mpq_mul(temp_r1_i2, sc1->real, sc2->imag);
    mpq_mul(temp_i1_r2, sc1->imag, sc2->real);
    mpq_mul(temp_i1_i2, sc1->imag, sc2->imag);
    mpq_sub((*rsc)->real, temp_r1_r2, temp_i1_i2);
    mpq_add((*rsc)->imag, temp_r1_i2, temp_i1_r2);
    mpq_clears(temp_r1_r2, temp_r1_i2, temp_i1_r2, temp_i1_i2, NULL);
}

void larc_sca_sqrt_mpratcomplex(larc_mpratcomplex_t *sqroot, const larc_mpratcomplex_t sc){
    // go through mpc_t, mpfr_t to keep full precision in conversion
    mpc_t sc_mpc, root_mpc;
    mpc_init2(sc_mpc,mpreal_precision);
    mpc_init2(root_mpc,mpreal_precision);
    mpc_set_q_q(sc_mpc,sc->real,sc->imag,MPC_RNDNN);
    mpc_sqrt(root_mpc,sc_mpc,MPC_RNDNN);
    mpfr_get_q((*sqroot)->real,mpc_realref(root_mpc));
    mpfr_get_q((*sqroot)->imag,mpc_imagref(root_mpc));
    mpc_clear(sc_mpc);
    mpc_clear(root_mpc);
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
    if (num == NULL) { ALLOCFAIL(); }
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
            if (realpart == NULL) { ALLOCFAIL(); }
            strncpy(realpart, input_str, ipos - input_str - 1);
            realpart[ipos - input_str - 1] = '\0';
            mpq_set_str((*rsc)->real, realpart, io_string_base);
            mpq_canonicalize((*rsc)->real);
            free(realpart);
        }
    } else {
        fprintf(stderr,"ERROR in %s: could not convert '%s' to rational complex number.\n", __func__, input_str);
        raise(SIGSEGV);
        exit(1);
    }
}
#endif // USE_MPRATCOMPLEX

#if defined(USE_MPREAL) || defined(USE_MPCOMPLEX)
/******************************************************************************
 * mpfr_t (MPFR) scalar operations:
 * these wrap appropriate mpfr routines and convert return pointers to regular
 * mpfr_t structs
 ******************************************************************************/
void larc_sca_init_mpreal(mpfr_t *sc){mpfr_init2(*sc, mpreal_precision);}
void larc_sca_clear_mpreal(mpfr_t *sc){mpfr_clear(*sc);}
void larc_sca_set_mpreal(mpfr_t *rsc, const mpfr_t sc){mpfr_set(*rsc, sc, MPFR_RNDN);}

char *larc_sca_get_str_mpreal(const mpfr_t sc){
    mpfr_exp_t scexp;
    char *raw = mpfr_get_str(NULL, &scexp, io_string_base, 0, sc, MPFR_RNDN);
    char exponent_str[25];
    snprintf(exponent_str, 25,"e%ld", scexp);
    char *s = calloc(strlen(raw) + strlen(exponent_str) + 4, sizeof(char));
    if (s == NULL) { ALLOCFAIL(); }
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
    if (is_numeric) {
        strcat(s, exponent_str);
    }
    mpfr_free_str(raw);
    return s;
}
char *larc_sca_get_exact_str_mpreal(const mpfr_t sc){
    mpfr_exp_t scexp;
    char *raw = mpfr_get_str(NULL, &scexp, 16, 0, sc, MPFR_RNDN);
    char exponent_str[25];
    snprintf(exponent_str, 25, "p%ld", scexp);
    char *s = calloc(strlen(raw) + strlen(exponent_str) + 6, sizeof(char));
    if (s == NULL) { ALLOCFAIL(); }
    s[0] = '\0';
    int is_negative = 0;
    int is_numeric = 0;
    if (raw[0] == '-') {
        is_negative = 1;
    }
    if (('0' <= raw[is_negative]) && (raw[is_negative] <= '9')) {
        is_numeric = 1;
    } else if (('a' <= raw[is_negative]) && (raw[is_negative] <= 'f')) {
        is_numeric = 1;
    } else if (('A' <= raw[is_negative]) && (raw[is_negative] <= 'F')) {
        is_numeric = 1;
    }
    if (is_negative) {
        strcat(s, "-");
    }
    if (is_numeric) {
        strcat(s, "0x0.");
    }
    strcat(s, &raw[is_negative]);
    if (is_numeric) {
        strcat(s, exponent_str);
    }
    mpfr_free_str(raw);
    return s;
}
void larc_sca_set_str_mpreal(mpfr_t *rsc, const char *input_str){
    mpfr_set_str(*rsc, input_str, 0, MPFR_RNDN);
    // Currently LARC does not support NaN or Infinity values
    if (mpfr_number_p(*rsc)==0)
    {
        fprintf(stderr,"ERROR: %s given string which produced NaN\n",__func__);
        fprintf(stderr,"or Infinity value. This is not currently\n");
        fprintf(stderr,"supported by LARC. Exiting...\n");
        raise(SIGSEGV);
        exit(1);
    }
}
void larc_sca_set_2ldoubles_mpreal(mpfr_t *rsc, long double real_val, long double imag_val){
    if (0.0 != imag_val) {
        fprintf(stderr,"ERROR in %s: MPREAL cannot have imaginary component.\n", __func__);
        raise(SIGSEGV);
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

void larc_sca_conj_mpreal(mpfr_t *conj_out, const mpfr_t sc){mpfr_set(*conj_out, sc, MPFR_RNDN);}

int larc_sca_eq_mpreal(const mpfr_t val1, const mpfr_t val2){
    return mpfr_equal_p(val1, val2);
}

int larc_sca_cmp_mpreal(const mpfr_t val1, const mpfr_t val2){
    return mpfr_cmp(val1, val2);
}
#endif // #if defined(USE_MPREAL) || defined(USE_MPCOMPLEX)

#if defined(USE_MPCOMPLEX)
/******************************************************************************
 * mpc_t (MPC) scalar operations:
 * these wrap appropriate mpfr routines and convert return pointers to regular
 * mpc_t structs
 ******************************************************************************/
void larc_sca_init_mpcomplex(mpc_t *sc){mpc_init2(*sc, mpreal_precision);}
void larc_sca_clear_mpcomplex(mpc_t *sc){mpc_clear(*sc);}
void larc_sca_set_mpcomplex(mpc_t *rsc, const mpc_t sc){mpc_set(*rsc, sc, MPC_RNDNN);}

char *larc_sca_get_str_mpcomplex(const mpc_t sc){
    // mpc_set_str returns complex numbers in the following format:
    //   "(" ++ real_part ++ " " ++ imag_part ++ ")"
    char *raw = mpc_get_str(io_string_base, 0, sc, MPC_RNDNN);

    // The following code converts this format to real_part+I*imag_part
    // or real_part-I*imag_part as appropriate.
    char *s = calloc(strlen(raw) + 6, sizeof(char));
    if (s == NULL) { ALLOCFAIL(); }
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
char *larc_sca_get_exact_str_mpcomplex(const mpc_t sc){
    char *rpart = larc_sca_get_exact_str_mpreal(mpc_realref(sc));
    char *ipart = larc_sca_get_exact_str_mpreal(mpc_imagref(sc));
    char *s = calloc(strlen(rpart) + strlen(ipart) + 5, sizeof(char));
    if (s == NULL) { ALLOCFAIL(); }
    s[0] = '\0';
    strcat(s, rpart);
    strcat(s, "+I*");
    strcat(s, ipart);
    free(rpart);
    free(ipart);
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
            if (realpart == NULL) { ALLOCFAIL(); }
            strncpy(realpart, input_str, ipos - input_str - 1);
            realpart[ipos - input_str - 1] = '\0';
            mpfr_set_str(mpc_realref(*rsc), realpart, 0, MPFR_RNDN);
            free(realpart);
        }
    } else {
        fprintf(stderr,"ERROR in %s: could not convert '%s' to mpcomplex number.\n", __func__, input_str);
        raise(SIGSEGV);
        exit(1);
    }
    // Currently LARC does not support NaN or Infinity values
    if ( (mpfr_number_p(mpc_realref(*rsc))==0) ||
		(mpfr_number_p(mpc_imagref(*rsc))==0) )
    {
        fprintf(stderr,"ERROR: %s given string which produced NaN\n",__func__);
        fprintf(stderr,"or Infinity value. This is not currently\n");
        fprintf(stderr,"supported by LARC. Exiting...\n");
        raise(SIGSEGV);
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
    uint64_t hash = recursive_hash_from_two_integers(realhash, imaghash,
                      exponent);
#ifdef DEBUG_SCALARS_C
    printf("in %s, realhash = %lu, imaghash = %lu\n",
        __func__, realhash, imaghash);
    printf("result of recursive_hash_from_two_integers is %lu\n",hash);
#endif // DEBUG_SCALARS_C
    return hash;
#undef DEBUG_SCALARS_C
}
void larc_sca_sqrt_mpcomplex(mpc_t *sqroot, const mpc_t sc){
    mpc_sqrt(*sqroot, sc, MPC_RNDNN);
 }

void larc_sca_norm_mpcomplex(mpc_t *norm, const mpc_t sc){
    // the mpc_abs function takes an mpfr_t as its result,
    // and an mpfr_rnd_t to determine how it rounds
    mpc_abs(mpc_realref(*norm), sc, MPFR_RNDN);
    mpfr_set_ui(mpc_imagref(*norm), 0, MPFR_RNDN);
}

int larc_sca_is_real_mpcomplex(const mpc_t sc){
    return mpfr_zero_p(mpc_imagref(sc));
}

int larc_sca_is_imag_mpcomplex(const mpc_t sc){
    return mpfr_zero_p(mpc_realref(sc));
}

void larc_sca_conj_mpcomplex(mpc_t *conj_out, const mpc_t sc){mpc_conj(*conj_out, sc, MPC_RNDNN);}

int larc_sca_eq_mpcomplex(const mpc_t val1, const mpc_t val2)
{
    // this is taken from the MPC documentation
    return (mpc_cmp(val1, val2) == 0);
}

int larc_sca_cmp_mpcomplex(const mpc_t val1, const mpc_t val2)
{
    // In general, two complex numbers have no natural ordering in the way that
    // two integers or real numbers have. We produce an ordering for the
    // limited cases where ordering makes sense (e.g., two real numbers
    // which are stored in a variable of complex type), and also return 0 when
    // the two complex numbers are equal.
    // Currently, we warn the user when they try to compare numbers
    // that are not naturally orderable, and use an arbitary method to
    // make the comparison.
    static int warned = 0;
    int real_cmp = mpfr_cmp(mpc_realref(val1), mpc_realref(val2));
    int imag_cmp = mpfr_cmp(mpc_imagref(val1), mpc_imagref(val2));
    // equality is testable
    if (!real_cmp && !imag_cmp) return 0;
    // if both values real, comparable
    if ( !imag_cmp && (mpfr_cmp_ui(mpc_imagref(val1),0)==0) )
        return real_cmp;
    // if both values imaginary, comparable
    if ( !real_cmp && (mpfr_cmp_ui(mpc_realref(val1),0)==0) )
        return imag_cmp;
    // otherwise not comparable; default to cartesian heuristic.
    if (!warned) {
        warned = 1;
        fprintf(stderr,"in %s: Cannot compare input complex\n",__func__);
        fprintf(stderr,"numbers ");
        mpfr_out_str(stderr,10,0,mpc_realref(val1),MPFR_RNDN);
        fprintf(stderr," + I*");
        mpfr_out_str(stderr,10,0,mpc_imagref(val1),MPFR_RNDN);
        fprintf(stderr,"\nand ");
        mpfr_out_str(stderr,10,0,mpc_realref(val2),MPFR_RNDN);
        fprintf(stderr," + I*");
        mpfr_out_str(stderr,10,0,mpc_imagref(val2),MPFR_RNDN);
        fprintf(stderr,"\nWe will use an arbitrary heuristic where the real");
        fprintf(stderr,"\npart comparison takes precedence.\n");
    }
    if (real_cmp) return real_cmp;
    else return imag_cmp;
}
#endif // USE_MPCOMPLEX


#ifdef IS_BOUNDING
/******************************************************************************
 * BOUNDING scalar operations:
 ******************************************************************************/
void larc_sca_init_bounding(larc_exponent_scalar_t *sc){
    for (int i = 0; i < NUM_EXPONENTS; ++i) {
        (*sc)->explist[i] = 0;
    }
    (*sc)->is_zero = 1;
    (*sc)->is_one = 0;
    (*sc)->is_nan = 0;
    (*sc)->is_exact = 1;
}

void larc_sca_clear_bounding(larc_exponent_scalar_t *sc){
    for (int i = 0; i < NUM_EXPONENTS; ++i) {
        (*sc)->explist[i] = 0;
    }
    (*sc)->is_zero = 1;
    (*sc)->is_one = 0;
    (*sc)->is_nan = 0;
    (*sc)->is_exact = 1;
}

void larc_sca_set_bounding(larc_exponent_scalar_t *rsc, const larc_exponent_scalar_t sc){
    for (int i = 0; i < NUM_EXPONENTS; ++i) {
        (*rsc)->explist[i] = sc->explist[i];
    }
    (*rsc)->is_zero = sc->is_zero;
    (*rsc)->is_one = sc->is_one;
    (*rsc)->is_nan = sc->is_nan;
    (*rsc)->is_exact = sc->is_exact;
}

// This function takes a list of exponents, in any order, and normalizes them.
// The possible outputs are:
//      n > 0, which means temp_exps[0..(n-1)] are now the normalized non-zero exponents.
//      n == 0, which means the input was all zeros.
//      n == -1, which means the exponents add to exactly one.
//      n == -2, which means the exponents add to more than one.
// In all cases, temp_exps[0..(size_of_temp-1)] contains a normalized
// (except for two or more 1's) list of exponents.
int normalize_exponent_list(bexp_type temp_exps[], int size_of_temp) {
    int i, swap_temp;
    int done = 0;
    while (!done) {
        done = 1;
        for (i = 0; i < (size_of_temp-1); ++i) {
            // Scan for something not normalized, except two 1's next to each other.
            // If you find something, fix it and set done to false.
            if ((temp_exps[i] == 0) && (temp_exps[i+1] != 0)) {
                temp_exps[i] = temp_exps[i+1];
                temp_exps[i+1] = 0;
                done = 0;
            } else if ((temp_exps[i] > temp_exps[i+1]) && (temp_exps[i+1] != 0)) {
                swap_temp = temp_exps[i];
                temp_exps[i] = temp_exps[i+1];
                temp_exps[i+1] = swap_temp;
                done = 0;
            } else if ((temp_exps[i] == temp_exps[i+1]) && (temp_exps[i+1] > 1)) {
                --temp_exps[i];
                temp_exps[i+1] = 0;
                done = 0;
            }
        }
    }

    int count_of_ones = 0;
    int count_of_non_zeros = 0;
    for (int i = 0; i < size_of_temp; ++i) {
        if (temp_exps[i] == 1) {
            ++count_of_ones;
            ++count_of_non_zeros;
        } else if (temp_exps[i] != 0) {
            ++count_of_non_zeros;
        }
    }

    int answer = count_of_non_zeros;

    if (count_of_ones >= 2) {
       if (count_of_non_zeros > 2) {
           answer = -2;
       } else {
           answer = -1;
       }
    }
    return answer;
}

// Set rsc to final answer based on combined exponent list.
// This function assumes that the (*rsc)->is_exact flag has been set correctly,
// that is, the exact flag is 0 if the value represented by the combined exponent list
// is already known to be inexact, and the exact flag is 1 if the value represented by
// the combined exponent list is exact, even if possibly too long.  In the event the
// combined exponent list is too long, it will be shortened and the new exact flag will
// be set to 0.
void translate_explist_to_bounding(larc_exponent_scalar_t *rsc, bexp_type temp_exps[], int size_of_temp) {
    // Zero out everything except the exact flag.
    for (int i = 0; i < NUM_EXPONENTS; ++i) {
        (*rsc)->explist[i] = 0;
    }
    (*rsc)->is_zero = 0;
    (*rsc)->is_one = 0;
    (*rsc)->is_nan = 0;

    // Normalize the combined exponent list.
    int retval = normalize_exponent_list(temp_exps, size_of_temp);

    if (retval == -2) {
        // Result is larger than 1.
        if ((*rsc)->is_exact) {
            (*rsc)->is_nan = 1;
        } else {
#ifdef USE_LOWER
            (*rsc)->is_nan = 1;
#else
            (*rsc)->is_one = 1;
#endif // #ifdef USE_LOWER
        }
    } else if (retval == -1) {
        // Result is equal to 1.
        if ((*rsc)->is_exact) {
            (*rsc)->is_one = 1;
        } else {
#ifdef USE_LOWER
            (*rsc)->is_one = 1;
            (*rsc)->is_exact = 1;
#else
            (*rsc)->is_one = 1;
#endif // #ifdef USE_LOWER
        }
    } else if (retval == 0) {
        // Result is equal to 0.
        (*rsc)->is_zero = 1;
    } else if (retval <= NUM_EXPONENTS) {
        // Put final exponents in answer.
        for (int index = 0; index < retval; ++index) {
            (*rsc)->explist[index] = temp_exps[index];
        }
    } else {
        // Too many exponents to put in answer, need to round.
        (*rsc)->is_exact = 0;
#ifdef USE_LOWER
        for (int index = 0; index < NUM_EXPONENTS; ++index) {
            (*rsc)->explist[index] = temp_exps[index];
        }
#else
        // Duplicate last exponent and renormalize.
        temp_exps[NUM_EXPONENTS] = temp_exps[NUM_EXPONENTS-1];
        retval = normalize_exponent_list(temp_exps, NUM_EXPONENTS+1);
        if (retval < 0) {
            (*rsc)->is_one = 1;
        } else {
            for (int index = 0; index < NUM_EXPONENTS; ++index) {
                (*rsc)->explist[index] = temp_exps[index];
            }
        }
#endif // #ifdef USE_LOWER
    }
}

// This function takes a long double value v, and creates a list of exponents
// to represent v as closely as possible.
// The possible outputs (rc = return code) are:
//      rc == 0, which means v == 0
//      0 < rc <= size_of_temp, which means 0 < v < 1, and the returned exponent
//          list represents v exactly using rc exponents.
//      rc == size_of_temp+1, which means 0 < v < 1, and the returned exponent list
//          represents the best lower bound of v possible using size_of_temp exponents.
//      rc == -1, which means v == 1
//      rc == -2, which means v > 1
//      rc == -3, which means v < 0
int convert_ldouble_to_explist(long double v, bexp_type temp_exps[], int size_of_temp){
    for (int i = 0; i < size_of_temp; ++i) {
        temp_exps[i] = 0;
    }

    if (v == 0.0) {
        return 0;
    } else if (v < 0.0) {
        return -3;
    } else if (v == 1.0) {
        return -1;
    } else if (v > 1.0) {
        return -2;
    } else {
        int found_so_far = 0;
        int shifts_so_far = 0;
        long double vprime = v;
        while ((vprime > 0.0) && (found_so_far < size_of_temp)) {
// New algorithm:
            int new_shift;
            vprime = frexpl(vprime, &new_shift);
            shifts_so_far -= new_shift;
// Old algorithm:
//          while (vprime < 0.5) {
//              vprime *= 2.0;
//              shifts_so_far += 1;
//          }
            temp_exps[found_so_far] = shifts_so_far + 1;
            found_so_far += 1;
            vprime -= 0.5;
        }
        return found_so_far + ((vprime != 0.0) ? 1 : 0);
    }
}

// This function takes a mprational value v, and creates a list of exponents
// to represent v as closely as possible.
// The possible outputs (rc = return code) are:
//      rc == 0, which means v == 0
//      0 < rc <= size_of_temp, which means 0 < v < 1, and the returned exponent
//          list represents v exactly using rc exponents.
//      rc == size_of_temp+1, which means 0 < v < 1, and the returned exponent list
//          represents the best lower bound of v possible using size_of_temp exponents.
//      rc == -1, which means v == 1
//      rc == -2, which means v > 1
//      rc == -3, which means v < 0
int convert_mprational_to_explist(mpq_t v, bexp_type temp_exps[], int size_of_temp){
    for (int i = 0; i < size_of_temp; ++i) {
        temp_exps[i] = 0;
    }

    int rc0 = mpq_cmp_ui(v, 0, 1);
    int rc1 = mpq_cmp_ui(v, 1, 1);
    if (rc0 == 0) {
        return 0;
    } else if (rc0 < 0) {
        return -3;
    } else if (rc1 == 0) {
        return -1;
    } else if (rc1 > 0) {
        return -2;
    } else {
        int found_so_far = 0;
        int shifts_so_far = 0;
        mpq_t vprime, one;
        mpq_init(one);
        mpq_set_ui(one, 1, 1);
        mpq_init(vprime);
        mpq_set(vprime, v);
        while ((rc0 > 0) && (found_so_far < size_of_temp)) {
            rc1 = mpq_cmp_ui(vprime, 1, 1);
            while (rc1 < 0) {
                mpq_mul_2exp(vprime, vprime, 1);
                shifts_so_far += 1;
                rc1 = mpq_cmp_ui(vprime, 1, 1);
            }
            temp_exps[found_so_far] = shifts_so_far;
            found_so_far += 1;
            mpq_sub(vprime, vprime, one);
            rc0 = mpq_cmp_ui(vprime, 0, 1);
        }
        mpq_clear(vprime);
        mpq_clear(one);
        return found_so_far + ((rc0 != 0) ? 1 : 0);
    }
}

// In string representation, zero is "0", one is "1", and all other values are
// "explist[E1,E2,...]" where Ei is the i-th exponent in order E1 < E2 < ...
char *larc_sca_get_str_bounding(const larc_exponent_scalar_t sc){
    char local[4096];
    char oneexp[1024];

#ifdef USE_LOWER
    char *approx = "v";
#else
    char *approx = "^";
#endif // #ifdef USE_LOWER

    local[0] = '\0';
    if (!sc->is_exact) {
        strcat(local, approx);
    }
    if (sc->is_nan) {
        strcat(local, "NaN");
    } else if (sc->is_zero) {
        strcat(local, "0");
    } else if (sc->is_one) {
        strcat(local, "1");
    } else {
        strcat(local, "explist[");
        for (int i = 0; i < NUM_EXPONENTS; ++i) {
            if (sc->explist[i] != 0) {
                if (i > 0) {
                    strcat(local, ",");
                }
                sprintf(oneexp, "%lu", (unsigned long int) sc->explist[i]);
                strcat(local, oneexp);
            }
        }
        strcat(local, "]");
    }

    char *answer = (char *) malloc(strlen(local) + 1);
    if (answer == NULL) { ALLOCFAIL(); }
    strcpy(answer, local);
    return answer;
}

char *larc_sca_get_exact_str_bounding(const larc_exponent_scalar_t sc){
    return larc_sca_get_str_bounding(sc);
}

void store_retval_and_explist_bounding(larc_exponent_scalar_t *rsc, int retval, bexp_type temp_exps[], int size_of_temp){
    if (retval == 0) {
        (*rsc)->is_zero = 1;
    } else if (retval == -1) {
        // Switch flags; leave exponent list all zero.
        (*rsc)->is_zero = 0;
        (*rsc)->is_one = 1;
    } else if (retval < -1) {
        // Invalid value.
        (*rsc)->is_zero = 0;
        (*rsc)->is_nan = 1;
    } else {
        // Value is between zero and one.
        (*rsc)->is_zero = 0;
        if (retval > NUM_EXPONENTS) {
            (*rsc)->is_exact = 0;
        }
        if (retval > size_of_temp) {
            retval = size_of_temp;
        }
        translate_explist_to_bounding(rsc, temp_exps, retval);
    }
}

void larc_sca_set_2ldoubles_bounding(larc_exponent_scalar_t *rsc, long double real_val, long double imag_val){
    if (0.0 != imag_val) {
        fprintf(stderr,"ERROR in %s: BOUNDING cannot have imaginary component.\n", __func__);
        raise(SIGSEGV);
        exit(1);
    }
    larc_sca_clear_bounding(rsc);
    bexp_type temp_exponents[NUM_EXPONENTS+1];
    int retval = convert_ldouble_to_explist(real_val, temp_exponents, NUM_EXPONENTS+1);
    store_retval_and_explist_bounding(rsc, retval, temp_exponents, NUM_EXPONENTS+1);
}

void larc_sca_set_str_real(long double *ret, const char *input_str);

void larc_sca_set_str_bounding(larc_exponent_scalar_t *rsc, const char *input_str){
#ifdef USE_LOWER
    char approx = 'v';
#else
    char approx = '^';
#endif // #ifdef USE_LOWER
    larc_sca_clear_bounding(rsc);
    (*rsc)->is_zero = 0;

    int approximate_found = 0;
    if (input_str[0] == approx) {
        approximate_found = 1;
        input_str += 1;
    }

    if (strncmp(input_str, "NaN", strlen("NaN")) == 0) {
        // input_str begins with "NaN"
        (*rsc)->is_nan = 1;
    } else if (strncmp(input_str, "explist[", strlen("explist[")) == 0) {
        // input_str begins with "explist["
        int ilen = strlen(input_str);
        if (input_str[ilen-1] != ']') {
            (*rsc)->is_nan = 1;
        } else {
            bexp_type temp_exponents[NUM_EXPONENTS+1];
            int index = 0;
            int pos = strlen("explist[");
            char *next;
            while ((input_str[pos] != '\0') && (index <= NUM_EXPONENTS)) {
                long int ev = strtol(input_str+pos, &next, 10);
                temp_exponents[index] = ev;
                ++index;
                pos = (next - input_str) + 1;
            }
            translate_explist_to_bounding(rsc, temp_exponents, index);
            if (approximate_found) {
                (*rsc)->is_exact = 0;
            }
        }
    } else if (strchr(input_str, '/') != NULL) {
        // input_str contains '/', try MPRATIONAL format
        mpq_t rat_input;
        mpq_init(rat_input);
        mpq_set_str(rat_input, input_str, 10);
        mpq_canonicalize(rat_input);
        bexp_type temp_exponents[NUM_EXPONENTS+1];
        int retval = convert_mprational_to_explist(rat_input, temp_exponents, NUM_EXPONENTS+1);
        store_retval_and_explist_bounding(rsc, retval, temp_exponents, NUM_EXPONENTS+1);
        mpq_clear(rat_input);
        if (approximate_found) {
            (*rsc)->is_exact = 0;
        }
    } else {
        // try floating point format
        long double val;
        larc_sca_set_str_real(&val, input_str);
        larc_sca_set_2ldoubles_bounding(rsc, val, 0.0L);
        if (approximate_found) {
            (*rsc)->is_exact = 0;
        }
    }
}

void larc_sca_add_bounding(larc_exponent_scalar_t *rsc, const larc_exponent_scalar_t sc1, const larc_exponent_scalar_t sc2){
    if (sc1->is_nan || sc2->is_nan) {
        // If either of the inputs is a NaN, then so is the output.
        larc_sca_clear_bounding(rsc);
        (*rsc)->is_zero = 0;
        (*rsc)->is_nan = 1;
    } else if (sc1->is_zero) {
        // 0 + sc2 = sc2
        // Note: No such thing as an inexact zero.
        larc_sca_set_bounding(rsc, sc2);
    } else if (sc2->is_zero) {
        // sc1 + 0 = sc1
        // Note: No such thing as an inexact zero.
        larc_sca_set_bounding(rsc, sc1);
    } else if (sc1->is_one && sc1->is_exact) {
        // exact-1 + non-zero-sc2 = NaN
        larc_sca_clear_bounding(rsc);
        (*rsc)->is_zero = 0;
        (*rsc)->is_nan = 1;
    } else if (sc2->is_one && sc2->is_exact) {
        // non-zero-sc1 + exact-1 = NaN
        larc_sca_clear_bounding(rsc);
        (*rsc)->is_zero = 0;
        (*rsc)->is_nan = 1;
    } else if (sc1->is_one || sc2->is_one) {
        // approx-1 + anything-but-exact-1 = approx-1
        larc_sca_clear_bounding(rsc);
        (*rsc)->is_zero = 0;
        (*rsc)->is_one = 1;
        (*rsc)->is_exact = 0;
    } else {
        // Clear is_zero for now; rest of answer will be filled in later.
        larc_sca_clear_bounding(rsc);
        (*rsc)->is_zero = 0;
        // Determine if answer is exact.
        (*rsc)->is_exact = sc1->is_exact && sc2->is_exact;

        // Calclate result temporary exponent list.
        bexp_type temp_exponents[2*NUM_EXPONENTS];
        int index = 0;
        for (int index1 = 0; index1 < NUM_EXPONENTS; ++index1) {
            if (sc1->explist[index1] != 0)  {
                temp_exponents[index++] = sc1->explist[index1];
            }
        }
        for (int index2 = 0; index2 < NUM_EXPONENTS; ++index2) {
            if (sc2->explist[index2] != 0)  {
                temp_exponents[index++] = sc2->explist[index2];
            }
        }

        // Fill in answer.
        translate_explist_to_bounding(rsc, temp_exponents, index);
    }
}

void larc_sca_mult_bounding(larc_exponent_scalar_t *rsc, const larc_exponent_scalar_t sc1, const larc_exponent_scalar_t sc2){
    if (sc1->is_nan || sc2->is_nan) {
        // If either of the inputs is a NaN, then so is the output.
        larc_sca_clear_bounding(rsc);
        (*rsc)->is_zero = 0;
        (*rsc)->is_nan = 1;
    } else if (sc1->is_zero || sc2->is_zero) {
        // 0 * anything = 0
        // Note: No such thing as an inexact zero.
        larc_sca_clear_bounding(rsc);
    } else if (sc1->is_one) {
        // 1 * sc2 = sc2
        larc_sca_set_bounding(rsc, sc2);
        (*rsc)->is_exact = sc1->is_exact && sc2->is_exact;
    } else if (sc2->is_one) {
        // sc1 * 1 = sc1
        larc_sca_set_bounding(rsc, sc1);
        (*rsc)->is_exact = sc1->is_exact && sc2->is_exact;
    } else {
        // Clear is_zero for now; rest of answer will be filled in later.
        larc_sca_clear_bounding(rsc);
        (*rsc)->is_zero = 0;
        // Determine if answer is exact.
        (*rsc)->is_exact = sc1->is_exact && sc2->is_exact;

        // Calclate result temporary exponent list.
        bexp_type temp_exponents[NUM_EXPONENTS*NUM_EXPONENTS];
        int index = 0;
        for (int index1 = 0; index1 < NUM_EXPONENTS; ++index1) {
            for (int index2 = 0; index2 < NUM_EXPONENTS; ++index2) {
                if ((sc1->explist[index1] != 0) && (sc2->explist[index2] != 0)) {
                    temp_exponents[index++] = sc1->explist[index1] + sc2->explist[index2];
                }
            }
        }

        // Fill in answer.
        translate_explist_to_bounding(rsc, temp_exponents, index);
    }
}

void larc_sca_divide_bounding(larc_exponent_scalar_t *rsc, const larc_exponent_scalar_t sc1, const larc_exponent_scalar_t sc2){
    fprintf(stderr,"ERROR in %s: Functionality is not implemented yet.\n", __func__);
    raise(SIGSEGV);
    exit(1);
}

uint64_t larc_sca_hash_bounding(const larc_exponent_scalar_t sc, uint64_t exponent)
{
    uint64_t ping, pong, limb;

    // initializing ping with flag information.
    ping = 16 + sc->is_zero + 2*sc->is_one + 4*sc->is_nan + 8*sc->is_exact;

    // Loop though exponents, starting with likely non-zero's first.
    for (int i = 0; i < NUM_EXPONENTS; ++i) {
        pong = mult_golden_hash(ping, 64);
        limb = sc->explist[i];
        ping = pong ^ limb;
    }
    pong = mult_golden_hash(ping, exponent);

    return pong;
}

void larc_sca_sqrt_bounding(larc_exponent_scalar_t *sqroot, const larc_exponent_scalar_t sc){
    fprintf(stderr,"ERROR in %s: Functionality is not implemented yet.\n", __func__);
    raise(SIGSEGV);
    exit(1);
}

void larc_sca_norm_bounding(larc_exponent_scalar_t *norm, const larc_exponent_scalar_t sc){
    fprintf(stderr,"ERROR in %s: Functionality is not implemented yet.\n", __func__);
    raise(SIGSEGV);
    exit(1);
}

int larc_sca_eq_bounding(const larc_exponent_scalar_t val1, const larc_exponent_scalar_t val2){
    // Currently, equality comparison ignores the isexact flag.
    if (val1->is_zero != val2->is_zero) { return 0; }
    if (val1->is_one != val2->is_one) { return 0; }
    if (val1->is_nan != val2->is_nan) { return 0; }
    for (int i = 0; i < NUM_EXPONENTS; ++i) {
        if (val1->explist[i] != val2->explist[i]) { return 0; }
    }
    return 1;
}

int larc_sca_cmp_bounding(const larc_exponent_scalar_t val1, const larc_exponent_scalar_t val2){
    // Currently, 3-way comparison ignores the isexact flag.
    if (val1->is_nan) {
        if (val2->is_nan) {
            return 0;
        } else {
            return 1;
        }
    } else if (val2->is_nan) {
        return -1;
    } else if (val1->is_zero) {
        if (val2->is_zero) {
            return 0;
        } else {
            return -1;
        }
    } else if (val2->is_zero) {
        return 1;
    } else if (val1->is_one) {
        if (val2->is_one) {
            return 0;
        } else {
            return 1;
        }
    } else if (val2->is_one) {
        return -1;
    } else {
        for (int i = 0; i < NUM_EXPONENTS; ++i) {
            if ((val1->explist[i] == 0) && (val2->explist[i] != 0)) {
                return -1;
            } else if ((val1->explist[i] != 0) && (val2->explist[i] == 0)) {
                return 1;
            } else if (val1->explist[i] > val2->explist[i]) {
                return -1;
            } else if (val1->explist[i] < val2->explist[i]) {
                return 1;
            }
        }
    }
    return 0;
}

void sca_set_enum_bounding(larc_exponent_scalar_t *scalar, const enum sca_constant_spec enum_val){
    sca_set_enum_via_mpfr(scalar, enum_val);
}

#endif // #ifdef IS_BOUNDING


void sca_set_enum_via_mpfr (scalarType *output, const enum sca_constant_spec enum_val) {
    // here are two example of how to store a highly precise irrational number
    if (enum_val == SCALAR_ENUM_SQRT2) {
        // get square root of 2 to mpreal_precision-bit precision, then
        // convert to selected scalarType, rounding if necessary
        mpfr_t mp_sqrt2;
        mpfr_init2(mp_sqrt2, mpreal_precision);
        mpfr_set_d(mp_sqrt2, 2.0, MPFR_RNDN);
        mpfr_sqrt(mp_sqrt2, mp_sqrt2, MPFR_RNDN);
        if (VERBOSE>=DEBUG) {
            printf("calculated value of sqrt(2):\n");
            mpfr_out_str(stdout, 10, 0, mp_sqrt2, MPFR_RNDN);
            printf("\n");
        }
        convert_from_mpfr_to_scalarType(output, mp_sqrt2);
        mpfr_clear(mp_sqrt2);
    } else if (enum_val == SCALAR_ENUM_INV_SQRT2) {
        // get inverse square root of 2 to mpreal_precision-bit precision, then
        // convert to selected scalarType, rounding if necessary
        mpfr_t mp_invsqrt2;
        mpfr_init2(mp_invsqrt2, mpreal_precision);
        mpfr_set_d(mp_invsqrt2, 0.5, MPFR_RNDN);
        mpfr_sqrt(mp_invsqrt2, mp_invsqrt2, MPFR_RNDN);
        if (VERBOSE>=DEBUG) {
            printf("calculated value of sqrt(0.5):\n");
            mpfr_out_str(stdout, 10, 0, mp_invsqrt2, MPFR_RNDN);
            printf("\n");
        }
        convert_from_mpfr_to_scalarType(output, mp_invsqrt2);
        mpfr_clear(mp_invsqrt2);
    } else if (enum_val == SCALAR_ENUM_SQRT3) {
        // get square root of 3 to mpreal_precision-bit precision, then
        // convert to selected scalarType, rounding if necessary
        mpfr_t mp_sqrt3;
        mpfr_init2(mp_sqrt3, mpreal_precision);
        mpfr_set_d(mp_sqrt3, 3.0, MPFR_RNDN);
        mpfr_sqrt(mp_sqrt3, mp_sqrt3, MPFR_RNDN);
        if (VERBOSE>=DEBUG) {
            printf("calculated value of sqrt(3):\n");
            mpfr_out_str(stdout, 10, 0, mp_sqrt3, MPFR_RNDN);
            printf("\n");
        }
        convert_from_mpfr_to_scalarType(output, mp_sqrt3);
        mpfr_clear(mp_sqrt3);
    } else if (enum_val == SCALAR_ENUM_INV_SQRT3) {
        // get inverse square root of 3 to mpreal_precision-bit precision, then
        // convert to selected scalarType, rounding if necessary
        mpfr_t mp_invsqrt3;
        mpfr_init2(mp_invsqrt3, mpreal_precision);
        mpfr_set_d(mp_invsqrt3, 3.0, MPFR_RNDN);
        mpfr_ui_div(mp_invsqrt3, 1, mp_invsqrt3, MPFR_RNDN);
        mpfr_sqrt(mp_invsqrt3, mp_invsqrt3, MPFR_RNDN);
        if (VERBOSE>=DEBUG) {
            printf("calculated value of sqrt(1/3.0):\n");
            mpfr_out_str(stdout, 10, 0, mp_invsqrt3, MPFR_RNDN);
            printf("\n");
        }
        convert_from_mpfr_to_scalarType(output, mp_invsqrt3);
        mpfr_clear(mp_invsqrt3);
    } else if (enum_val == SCALAR_ENUM_SQRT6) {
        // get square root of 6 to mpreal_precision-bit precision, then
        // convert to selected scalarType, rounding if necessary
        mpfr_t mp_sqrt6;
        mpfr_init2(mp_sqrt6, mpreal_precision);
        mpfr_set_d(mp_sqrt6, 6.0, MPFR_RNDN);
        mpfr_sqrt(mp_sqrt6, mp_sqrt6, MPFR_RNDN);
        if (VERBOSE>=DEBUG) {
            printf("calculated value of sqrt(6):\n");
            mpfr_out_str(stdout, 10, 0, mp_sqrt6, MPFR_RNDN);
            printf("\n");
        }
        convert_from_mpfr_to_scalarType(output, mp_sqrt6);
        mpfr_clear(mp_sqrt6);
    } else if (enum_val == SCALAR_ENUM_INV_SQRT6) {
        // get inverse square root of 6 to mpreal_precision-bit precision, then
        // convert to selected scalarType, rounding if necessary
        mpfr_t mp_invsqrt6;
        mpfr_init2(mp_invsqrt6, mpreal_precision);
        mpfr_set_d(mp_invsqrt6, 6.0, MPFR_RNDN);
        mpfr_ui_div(mp_invsqrt6, 1, mp_invsqrt6, MPFR_RNDN);
        mpfr_sqrt(mp_invsqrt6, mp_invsqrt6, MPFR_RNDN);
        if (VERBOSE>=DEBUG) {
            printf("calculated value of sqrt(1/6.0):\n");
            mpfr_out_str(stdout, 10, 0, mp_invsqrt6, MPFR_RNDN);
            printf("\n");
        }
        convert_from_mpfr_to_scalarType(output, mp_invsqrt6);
        mpfr_clear(mp_invsqrt6);
    } else if (enum_val == SCALAR_ENUM_CUBERT2) {
        // get cube root of 2 to mpreal_precision-bit precision, then
        // convert to selected scalarType, rounding if necessary
        mpfr_t mp_cbrt2;
        mpfr_init2(mp_cbrt2, mpreal_precision);
        mpfr_set_d(mp_cbrt2, 2.0, MPFR_RNDN);
        mpfr_cbrt(mp_cbrt2, mp_cbrt2, MPFR_RNDN);
        if (VERBOSE>=DEBUG) {
            printf("calculated value of cbrt(2):\n");
            mpfr_out_str(stdout, 10, 0, mp_cbrt2, MPFR_RNDN);
            printf("\n");
        }
        convert_from_mpfr_to_scalarType(output, mp_cbrt2);
        mpfr_clear(mp_cbrt2);
    } else if (enum_val == SCALAR_ENUM_INV_CUBERT2) {
        // get inverse cube root of 2 to mpreal_precision-bit precision, then
        // convert to selected scalarType, rounding if necessary
        mpfr_t mp_invcbrt2;
        mpfr_init2(mp_invcbrt2, mpreal_precision);
        mpfr_set_d(mp_invcbrt2, 0.5, MPFR_RNDN);
        mpfr_cbrt(mp_invcbrt2, mp_invcbrt2, MPFR_RNDN);
        if (VERBOSE>=DEBUG) {
            printf("calculated value of cbrt(0.5):\n");
            mpfr_out_str(stdout, 10, 0, mp_invcbrt2, MPFR_RNDN);
            printf("\n");
        }
        convert_from_mpfr_to_scalarType(output, mp_invcbrt2);
        mpfr_clear(mp_invcbrt2);
    } else if (enum_val == SCALAR_ENUM_CUBERT4) {
        // get cube root of 4 to mpreal_precision-bit precision, then
        // convert to selected scalarType, rounding if necessary
        mpfr_t mp_cbrt4;
        mpfr_init2(mp_cbrt4, mpreal_precision);
        mpfr_set_d(mp_cbrt4, 4.0, MPFR_RNDN);
        mpfr_cbrt(mp_cbrt4, mp_cbrt4, MPFR_RNDN);
        if (VERBOSE>=DEBUG) {
            printf("calculated value of cbrt(4):\n");
            mpfr_out_str(stdout, 10, 0, mp_cbrt4, MPFR_RNDN);
            printf("\n");
        }
        convert_from_mpfr_to_scalarType(output, mp_cbrt4);
        mpfr_clear(mp_cbrt4);
    } else if (enum_val == SCALAR_ENUM_INV_CUBERT4) {
        // get inverse cube root of 4 to mpreal_precision-bit precision, then
        // convert to selected scalarType, rounding if necessary
        mpfr_t mp_invcbrt4;
        mpfr_init2(mp_invcbrt4, mpreal_precision);
        mpfr_set_d(mp_invcbrt4, 0.25, MPFR_RNDN);
        mpfr_cbrt(mp_invcbrt4, mp_invcbrt4, MPFR_RNDN);
        if (VERBOSE>=DEBUG) {
            printf("calculated value of cbrt(0.25):\n");
            mpfr_out_str(stdout, 10, 0, mp_invcbrt4, MPFR_RNDN);
            printf("\n");
        }
        convert_from_mpfr_to_scalarType(output, mp_invcbrt4);
        mpfr_clear(mp_invcbrt4);
    } else {
        fprintf(stderr, "ERROR: Unrecognized scalar enum %d.\n", (int) enum_val);
    }
}


#if !defined(IS_MP) && !defined(IS_BOUNDING)
/******************************************************************************
 * Integer, Complex, Real:
 * these wrap the usual operations
 * TODO: for real/complex - do we need to worry about FPCLASSIFY here?
 ******************************************************************************/
void larc_sca_set_arith (scalarType *ret, const scalarType a){*ret = a;}
void larc_sca_add_arith (scalarType *ret, const scalarType a, const scalarType b){*ret = a + b;}
void larc_sca_mult_arith(scalarType *ret, const scalarType a, const scalarType b){*ret = a * b;}
void larc_sca_divide_arith(scalarType *ret, const scalarType a, const scalarType b){*ret = a / b;}
#endif // not IS_MP

#if defined(USE_REAL) || defined(USE_COMPLEX) || defined(IS_BOUNDING)
// Now using strtold() instead of sscanf() for best coverage of numeric format
// possibilities, since it will disregard excess precision in the string input.
// We will still get an error on underflow or overflow.
void larc_sca_set_str_real(long double *ret, const char *input_str){
    char *endptr = NULL;
    errno = 0;  // Zero errno so that we can tell if strtold returns an error.
    *ret = strtold(input_str, &endptr);
    int rc = 0;  // Assume success to start with.
    if ((errno != 0) || (endptr == input_str)) {
        rc = 1;  // Failure
    }
    if (rc) {
        fprintf(stderr,"ERROR in %s: could not convert '%s' to real number.\n", __func__, input_str);
        raise(SIGSEGV);
        exit(1);
    }
    // Currently LARC does not support NaN or Infinity values
    if ( (fpclassify(*ret)==FP_NAN) || (fpclassify(*ret)==FP_INFINITE) )
    {
        fprintf(stderr,"ERROR: %s given string which produced NaN\n",__func__);
        fprintf(stderr,"or Infinity value. This is not currently\n");
        fprintf(stderr,"supported by LARC. Exiting...\n");
        raise(SIGSEGV);
	exit(1);
    }
    return;
}
#endif // USE_REAL or USE_COMPLEX or IS_BOUNDING

#if defined(USE_REAL) || defined(USE_COMPLEX)
char *larc_sca_get_str_real(const long double n){
    char *out = calloc(LDBL_MANT_DEC_DIG + 21, sizeof(char));
    if (out == NULL) { ALLOCFAIL(); }
    snprintf(out, LDBL_MANT_DEC_DIG + 20, "%.20Lg", n);
    return out;
}
char *larc_sca_get_exact_str_real(const long double n){
    char *out = calloc(LDBL_MANT_DEC_DIG + 21, sizeof(char));
    if (out == NULL) { ALLOCFAIL(); }
    snprintf(out, LDBL_MANT_DEC_DIG + 20, "%La", n);
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
int larc_sca_cmp_real(const long double a, const long double b){return (a > b) - (a < b);}
int rc_eq(const scalarType a, const scalarType b){
  // Check to see if two scalars are equal
  // NOTE: this is probably not a good thing to do!
  return a == b;
}
#endif // USE_REAL or USE_COMPLEX

#ifdef USE_COMPLEX
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
            if (realpart == NULL) { ALLOCFAIL(); }
            strncpy(realpart, input_str, ipos - input_str - 1);
            realpart[ipos - input_str - 1] = '\0';
            larc_sca_set_str_real(&real, realpart);
            free(realpart);
        }
    } else {
        fprintf(stderr,"ERROR in %s: could not convert '%s' to complex number.\n", __func__, input_str);
        raise(SIGSEGV);
        exit(1);
    }
    *ret = real+imag*I;
}
char *larc_sca_get_str_complex(const long double complex n){
    char *out = calloc(2*LDBL_MANT_DEC_DIG + 44, sizeof(char));
    if (out == NULL) { ALLOCFAIL(); }
    if (cimagl(n) < 0.0) {
        snprintf(out, 2*LDBL_MANT_DEC_DIG + 43,
                 "%.20Lg-I*%.20Lg", creall(n), -cimagl(n));
    } else {
        snprintf(out, 2*LDBL_MANT_DEC_DIG + 43,
                 "%.20Lg+I*%.20Lg", creall(n), cimagl(n));
    }
    return out;
}
char *larc_sca_get_exact_str_complex(const long double complex n){
    char *out = calloc(2*LDBL_MANT_DEC_DIG + 44, sizeof(char));
    if (out == NULL) { ALLOCFAIL(); }
    if (cimagl(n) < 0.0) {
        snprintf(out, 2*LDBL_MANT_DEC_DIG + 43,
                 "%La-I*%La", creall(n), -cimagl(n));
    } else {
        snprintf(out, 2*LDBL_MANT_DEC_DIG + 43,
                 "%La+I*%La", creall(n), cimagl(n));
    }
    return out;
}
void larc_sca_2ldoubles_complex(long double complex *ret, long double real_val, long double imag_val){
    *ret = real_val+imag_val*I;
}
uint64_t larc_sca_hash_complex(long double complex sc, uint64_t exponent){return hash_from_two_longdoubles(creall(sc), cimagl(sc), exponent);}
void larc_sca_sqrt_complex(long double complex *sqroot, const long double complex a){*sqroot= csqrtl(a);}
void larc_sca_norm_complex(long double complex *norm, const long double complex a){*norm= cabsl(a) + I*0.0L;}
void larc_sca_conj_complex(long double complex *conj_out, const long double complex sc){*conj_out = conjl(sc);}
int larc_sca_is_real_complex(const long double complex sc)
	{return fpclassify(cimagl(sc)) == FP_ZERO; }
int larc_sca_is_imag_complex(const long double complex sc)
	{return fpclassify(creall(sc)) == FP_ZERO; }
int larc_sca_cmp_complex(const long double complex a, const long double complex b){
    // In general, two complex numbers have no natural ordering in the way that
    // two integers or real numbers have. We produce an ordering for the
    // limited cases where ordering makes sense (e.g., two real numbers
    // which are stored in a variable of complex type), and also return 0 when
    // the two complex numbers are equal.
    // Currently, we warn the user when they try to compare numbers
    // that are not naturally orderable, and use an arbitary method to
    // make the comparison.
    static int warned = 0;
    int real_cmp = larc_sca_cmp_real(creall(a),creall(b));
    int imag_cmp = larc_sca_cmp_real(cimagl(a),cimagl(b));
    // equality is testable
    if (!real_cmp && !imag_cmp) return 0;
    // if both values real, comparable
    if ( !imag_cmp && (fpclassify(cimagl(a))==FP_ZERO) )
        return real_cmp;
    // if both values imaginary, comparable
    if ( !real_cmp && (fpclassify(creall(a))==FP_ZERO) )
        return imag_cmp;
    // otherwise not comparable; default to cartesian heuristic
    if (!warned) {
        warned = 1;
        fprintf(stderr,"in %s: Cannot compare input complex\n",__func__);
        fprintf(stderr,"numbers %Lg + I*%Lg", creall(a),cimagl(a));
        fprintf(stderr,"and %Lg + I*%Lg\n", creall(b),cimagl(b));
        fprintf(stderr,"\nWe will use an arbitrary heuristic where the real");
        fprintf(stderr,"\npart comparison takes precedence.\n");
    }
    if (real_cmp) return real_cmp;
    else return imag_cmp;
}
#endif // USE_COMPLEX

#if defined(USE_INTEGER) || defined(USE_BOOLEAN)
void int_set_str(int64_t *ret, const char *input_str){*ret = atol(input_str);}
//NOTE: consider strtol() for input with base
char *int_get_str(const int64_t n){
    char *out = calloc(31, sizeof(char));
    if (out == NULL) { ALLOCFAIL(); }
    snprintf(out, 30, "%ld", n);
    return out;
}
void int_set_2ldoubles(int64_t *ret, long double real_val, long double imag_val){
    if (0.0L != imag_val) {
        fprintf(stderr,"ERROR in %s: INTEGER cannot have imaginary component.\n", __func__);
        raise(SIGSEGV);
        exit(1);
    }
    *ret = (int64_t) real_val;
    if (real_val != (long double) *ret) {
        fprintf(stderr,"ERROR in %s: Can't convert %.20Lg to INTEGER.\n", __func__, real_val);
        raise(SIGSEGV);
        exit(1);
    }
}
int int_cmp(const scalarType a, const scalarType b){return (a > b) - (a < b);}
int int_eq(const scalarType a, const scalarType b){return a == b;}
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
void bool_sqrt(int64_t *ret, const int64_t a){*ret = a;}
void bool_set (int64_t *ret, const int64_t a){*ret = (a!=0);}
// automatically convert all nonzero integers to 1.
void bool_set_str(int64_t *ret, const char *input_str){*ret = (atol(input_str) != 0);}
void bool_set_2ldoubles(int64_t *ret, long double real_val, long double imag_val){
    if (0.0 != imag_val) {
        fprintf(stderr,"ERROR in %s: INTEGER/boolean cannot have imaginary component.\n", __func__);
        raise(SIGSEGV);
        exit(1);
    }
    *ret = (0.0 != real_val);
}
#endif // #if defined(USE_INTEGER) || defined(USE_BOOLEAN)

int always_true_is_real (const scalarType a){return 1;}
int always_false_is_imag (const scalarType a){return 0;}

/******************************************************************************
 * Initialize scalar routines
 ******************************************************************************/
// to do standard arithmetic, use this.
void init_scalarOps(int verbose){
    if (verbose>BASIC){
        printf("Setting scalar ops to standard arithmetic for all scalarTypes.\n");
    }

/*****************************************************************/
// We need to define the sca_ functions (ops) for each scalarType*
/*****************************************************************/
#ifdef USE_MPINTEGER
    sca_init     = larc_sca_init_mpinteger;
    sca_clear    = larc_sca_clear_mpinteger;
    sca_cmp      = mpz_cmp;
    sca_eq       = larc_sca_eq_mpinteger;
    sca_set      = larc_sca_set_mpinteger;
    sca_set_str  = larc_sca_set_str_mpinteger;
    sca_set_2ldoubles = larc_sca_set_2ldoubles_mpinteger;
    sca_set_enum = sca_set_enum_via_mpfr;
    sca_get_readable_approx_str  = larc_sca_get_str_mpinteger;
    sca_get_exact_str  = larc_sca_get_str_mpinteger;
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
    sca_is_pure_real  = always_true_is_real;
    sca_is_pure_imag  = always_false_is_imag;
#endif // USE_MPINTEGER

#ifdef USE_MPRATIONAL
    sca_init     = larc_sca_init_mprational;
    sca_clear    = larc_sca_clear_mprational;
    sca_cmp      = mpq_cmp;
    sca_eq       = mpq_equal;
    sca_set      = larc_sca_set_mprational;
    sca_set_str  = larc_sca_set_str_mprational;
    sca_set_2ldoubles = larc_sca_set_2ldoubles_mprational;
    sca_set_enum = sca_set_enum_via_mpfr;
    sca_get_readable_approx_str  = larc_sca_get_str_mprational;
    sca_get_exact_str  = larc_sca_get_str_mprational;
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
    sca_is_pure_real  = always_true_is_real;
    sca_is_pure_imag  = always_false_is_imag;
#endif // USE_MPRATIONAL

#ifdef USE_MPRATCOMPLEX
    sca_init     = larc_sca_init_mpratcomplex;
    sca_clear    = larc_sca_clear_mpratcomplex;
    sca_cmp      = larc_sca_cmp_mpratcomplex;
    sca_eq       = larc_sca_eq_mpratcomplex;
    sca_set      = larc_sca_set_mpratcomplex;
    sca_set_str  = larc_sca_set_str_mpratcomplex;
    sca_set_2ldoubles = larc_sca_set_2ldoubles_mpratcomplex;
    sca_set_enum = sca_set_enum_via_mpfr;
    sca_get_readable_approx_str  = larc_sca_get_str_mpratcomplex;
    sca_get_exact_str  = larc_sca_get_str_mpratcomplex;
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
    sca_is_pure_real  = larc_sca_is_real_mpratcomplex;
    sca_is_pure_imag  = larc_sca_is_imag_mpratcomplex;
#endif // USE_MPRATCOMPLEX

#ifdef USE_MPREAL
    sca_init     = larc_sca_init_mpreal;
    sca_clear    = larc_sca_clear_mpreal;
    sca_cmp      = larc_sca_cmp_mpreal;
    sca_eq       = larc_sca_eq_mpreal;
    sca_set      = larc_sca_set_mpreal;
    sca_set_str  = larc_sca_set_str_mpreal;
    sca_set_2ldoubles = larc_sca_set_2ldoubles_mpreal;
    sca_set_enum = sca_set_enum_via_mpfr;
    sca_get_readable_approx_str  = larc_sca_get_str_mpreal;
    sca_get_exact_str  = larc_sca_get_exact_str_mpreal;
    sca_hash     = larc_sca_hash_mpreal;
    sca_add      = larc_sca_add_mpreal;
    sca_mult     = larc_sca_mult_mpreal;
    sca_divide   = larc_sca_divide_mpreal;
    sca_sqrt     = larc_sca_sqrt_mpreal;
    sca_norm     = larc_sca_norm_mpreal;
    sca_conj     = larc_sca_conj_mpreal;
    sca_is_pure_real  = always_true_is_real;
    sca_is_pure_imag  = always_false_is_imag;
#endif // USE_MPREAL

#ifdef USE_MPCOMPLEX
    sca_init     = larc_sca_init_mpcomplex;
    sca_clear    = larc_sca_clear_mpcomplex;
    sca_cmp      = larc_sca_cmp_mpcomplex;
    sca_eq       = larc_sca_eq_mpcomplex;
    sca_set      = larc_sca_set_mpcomplex;
    sca_set_str  = larc_sca_set_str_mpcomplex;
    sca_set_2ldoubles = larc_sca_set_2ldoubles_mpcomplex;
    sca_set_enum = sca_set_enum_via_mpfr;
    sca_get_readable_approx_str  = larc_sca_get_str_mpcomplex;
    sca_get_exact_str  = larc_sca_get_exact_str_mpcomplex;
    sca_hash     = larc_sca_hash_mpcomplex;
    sca_add      = larc_sca_add_mpcomplex;
    sca_mult     = larc_sca_mult_mpcomplex;
    sca_divide   = larc_sca_divide_mpcomplex;
    sca_sqrt     = larc_sca_sqrt_mpcomplex;
    sca_norm     = larc_sca_norm_mpcomplex;
    sca_conj     = larc_sca_conj_mpcomplex;
    sca_is_pure_real  = larc_sca_is_real_mpcomplex;
    sca_is_pure_imag  = larc_sca_is_imag_mpcomplex;
#endif // USE_MPCOMPLEX

#if !defined(IS_MP) && !defined(IS_BOOLEAN) && !defined(IS_BOUNDING)
    // these seven definitions are the same for COMPLEX, REAL, and INTEGER
    sca_init     = empty1;
    sca_clear    = empty1;
    sca_set      = larc_sca_set_arith;
    sca_set_enum = sca_set_enum_via_mpfr;
    sca_add      = larc_sca_add_arith;
    sca_mult     = larc_sca_mult_arith;
    sca_divide   = larc_sca_divide_arith;
#endif // !defined(IS_MP) && !defined(IS_BOOLEAN) && !defined(IS_BOUNDING)

#ifdef USE_INTEGER
    sca_cmp      = int_cmp;
    sca_eq    = int_eq;
    sca_set_str  = int_set_str;
    sca_set_2ldoubles = int_set_2ldoubles;
    sca_get_readable_approx_str  = int_get_str;
    sca_get_exact_str  = int_get_str;
    sca_hash     = int_hash;
    sca_sqrt     = larc_sca_sqrt_int;
    sca_norm     = larc_sca_norm_int;
    sca_conj     = larc_sca_set_arith;
    sca_is_pure_real  = always_true_is_real;
    sca_is_pure_imag  = always_false_is_imag;
#endif // USE_INTEGER

#ifdef USE_REAL
    sca_cmp      = larc_sca_cmp_real;
    sca_eq    = rc_eq;
    sca_set_str  = larc_sca_set_str_real;
    sca_set_2ldoubles = larc_sca_2ldoubles_real;
    sca_get_readable_approx_str  = larc_sca_get_str_real;
    sca_get_exact_str  = larc_sca_get_exact_str_real;
    sca_hash     = hash_from_one_longdouble;
    sca_sqrt     = larc_sca_sqrt_real;
    sca_norm     = larc_sca_norm_real;
    sca_conj     = larc_sca_set_arith;
    sca_is_pure_real  = always_true_is_real;
    sca_is_pure_imag  = always_false_is_imag;
#endif // USE_REAL

#ifdef USE_COMPLEX
    sca_cmp      = larc_sca_cmp_complex;
    sca_eq    = rc_eq;
    sca_set_str  = larc_sca_set_str_complex;
    sca_set_2ldoubles = larc_sca_2ldoubles_complex;
    sca_get_readable_approx_str  = larc_sca_get_str_complex;
    sca_get_exact_str  = larc_sca_get_exact_str_complex;
    sca_hash     = larc_sca_hash_complex;
    sca_sqrt     = larc_sca_sqrt_complex;
    sca_norm     = larc_sca_norm_complex;
    sca_conj     = larc_sca_conj_complex;
    sca_is_pure_real  = larc_sca_is_real_complex;
    sca_is_pure_imag  = larc_sca_is_imag_complex;
#endif // USE_COMPLEX

#ifdef USE_BOOLEAN
    sca_init     = empty1;
    sca_clear    = empty1;
    sca_set      = larc_sca_set_arith;
    sca_set_str  = bool_set_str;
    sca_set_2ldoubles = bool_set_2ldoubles;
    sca_set_enum = sca_set_enum_via_mpfr;
    sca_get_readable_approx_str  = int_get_str;
    sca_get_exact_str  = int_get_str;
    sca_hash     = int_hash;
    sca_add      = bool_add;
    sca_mult     = bool_mult;
    sca_divide   = bool_divide;
    sca_sqrt     = bool_sqrt;
    //sca_cmp      = bool_cmp;
    sca_cmp      = int_cmp;
    sca_eq    = int_eq;
    sca_norm     = bool_set;
    sca_conj     = bool_set;
    sca_is_pure_real  = always_true_is_real;
    sca_is_pure_imag  = always_false_is_imag;
#endif // USE_BOOLEAN

#ifdef USE_CLIFFORD
    sca_init     = larc_sca_init_clifford;
    sca_clear    = larc_sca_clear_clifford;
    sca_cmp      = larc_sca_cmp_clifford;
    sca_eq       = larc_sca_eq_clifford;
    sca_set      = larc_sca_set_clifford;
    sca_set_str  = larc_sca_set_str_clifford;
    sca_set_2ldoubles = larc_sca_2ldoubles_clifford;
    sca_set_enum = sca_set_enum_clifford;
    sca_get_readable_approx_str  = larc_sca_get_str_clifford;
    sca_get_exact_str  = larc_sca_get_str_clifford;
    sca_hash     = larc_sca_hash_clifford;
    sca_add      = larc_sca_add_clifford;
    sca_mult     = larc_sca_mult_clifford;
    sca_divide   = larc_sca_divide_clifford;
    sca_sqrt     = larc_sca_sqrt_clifford;
    sca_norm     = larc_sca_norm_clifford;
    sca_conj     = larc_sca_conj_clifford;
    sca_is_pure_real  = larc_sca_is_real_clifford;
    sca_is_pure_imag  = larc_sca_is_imag_clifford;
#endif // USE_CLIFFORD

#ifdef IS_BOUNDING
    sca_init     = larc_sca_init_bounding;
    sca_clear    = larc_sca_clear_bounding;
    sca_cmp      = larc_sca_cmp_bounding;
    sca_eq       = larc_sca_eq_bounding;
    sca_set      = larc_sca_set_bounding;
    sca_set_str  = larc_sca_set_str_bounding;
    sca_set_2ldoubles = larc_sca_set_2ldoubles_bounding;
    sca_set_enum = sca_set_enum_bounding;
    sca_get_readable_approx_str  = larc_sca_get_str_bounding;
    sca_get_exact_str  = larc_sca_get_str_bounding;
    sca_hash     = larc_sca_hash_bounding;
    sca_add      = larc_sca_add_bounding;
    sca_mult     = larc_sca_mult_bounding;
    sca_divide   = larc_sca_divide_bounding;
    sca_sqrt     = larc_sca_sqrt_bounding;
    sca_norm     = larc_sca_norm_bounding;
    sca_conj     = larc_sca_set_bounding;
    sca_is_pure_real  = always_true_is_real;
    sca_is_pure_imag  = always_false_is_imag;
#endif // IS_BOUNDING
}



void convert_from_mpfr_to_scalarType(scalarType *output, mpfr_t input)
{

#if defined(IS_INTEGER) || defined(IS_RATIONAL)
  if (0==mpfr_number_p(input))
  {
     fprintf(stderr,"ERROR in %s:\n",__func__);
     fprintf(stderr,"cannon convert NaN or Infinity into non-floating-point\n");
     fprintf(stderr,"scalarType. Exiting...\n");
     raise(SIGSEGV);
     exit(1);
  }
#endif // IS_INTEGER or IS_RATIONAL

#ifdef USE_MPRATIONAL
  mpfr_get_q(*output,input);
#elif defined(USE_MPRATCOMPLEX)
  mpfr_get_q((*output)->real,input);
  mpq_set_d((*output)->imag,0.0);
#elif defined(USE_MPREAL)
  mpfr_set(*output,input,MPFR_RNDN);
#elif defined(USE_MPCOMPLEX)
  mpc_set_fr(*output,input,MPC_RNDNN);
#elif defined(USE_MPINTEGER)
  mpfr_get_z(*output,input,MPFR_RNDN);
#elif defined(USE_COMPLEX)
  *output = mpfr_get_ld(input,MPFR_RNDN) + I*0.0;
#elif defined(USE_REAL)
  *output = mpfr_get_ld(input,MPFR_RNDN);
#elif defined(USE_INTEGER)
  *output = (int64_t)mpfr_get_si(input,MPFR_RNDN);
#elif defined(USE_BOOLEAN)
  *output = (int64_t)mpfr_get_si(input,MPFR_RNDN);
#elif defined(USE_CLIFFORD)
  larc_sca_set_clifford_to_zero(output);
  mpfr_get_q((*output)->real_coeffs[0],input);
#else
  fprintf(stderr,"in %s, unknown scalarType\n",__func__);
  raise(SIGSEGV);
  exit(1);
#endif
}

void convert_from_two_mpfr_to_complex_scalarType(scalarType *output,
	mpfr_t realpart, mpfr_t imagpart)
{
#ifdef USE_MPRATCOMPLEX
  mpfr_get_q((*output)->real,realpart);
  mpfr_get_q((*output)->imag,imagpart);
#elif defined(USE_MPCOMPLEX)
  mpc_set_fr_fr(*output,realpart,imagpart,MPC_RNDNN);
#elif defined(USE_COMPLEX)
  *output = mpfr_get_ld(realpart,MPFR_RNDN)
	+ I*mpfr_get_ld(imagpart,MPFR_RNDN);
#elif defined(USE_CLIFFORD) && defined(IS_COMPLEX)
  larc_sca_set_clifford_to_zero(output);
  mpfr_get_q((*output)->real_coeffs[0],realpart);
  mpfr_get_q((*output)->imag_coeffs[0],imagpart);
#else
#ifdef USE_MPRATIONAL
  mpfr_get_q(*output,realpart);
#elif defined(USE_MPREAL)
  mpfr_set(*output,realpart,MPFR_RNDN);
#elif defined(USE_MPINTEGER)
  mpfr_get_z(*output,realpart,MPFR_RNDN);
#elif defined(USE_REAL)
  *output = mpfr_get_ld(realpart,MPFR_RNDN);
#elif defined(USE_INTEGER)
  *output = (int64_t)mpfr_get_si(realpart,MPFR_RNDN);
#elif defined(USE_BOOLEAN)
  *output = (int64_t)mpfr_get_si(realpart,MPFR_RNDN);
#elif defined(USE_CLIFFORD) && !defined(IS_COMPLEX)
  larc_sca_set_clifford_to_zero(output);
  mpfr_get_q((*output)->real_coeffs[0],realpart);
#else
  fprintf(stderr,"in %s, unknown scalarType\n",__func__);
  raise(SIGSEGV);
  exit(1);
#endif
  if (0==mpfr_zero_p(imagpart))
  {
    fprintf(stderr,"WARNING from %s\n",__func__);
    fprintf(stderr,"tried to assign imaginary part to non-complex number\n");
    fprintf(stderr,"input ignored\n");
  }
#endif
}

int64_t get_pID_for_enum_const(const enum sca_constant_spec enum_val)
{
    if (scratchVars.submit_to_store_in_use)
        fprintf(stderr,"%s reusing scratchVars.submit_to_store!\n",__func__);
    scratchVars.submit_to_store_in_use = 1;

    sca_set_enum(&(scratchVars.submit_to_store),enum_val);
    mats_ptr_t temp = get_scalarPTR_for_scalarVal(scratchVars.submit_to_store);
    scratchVars.submit_to_store_in_use = 0;
    return temp->packedID;
}

int testForNaN(const scalarType variable)
{
   // this routine does not need to be called in the case of INTEGER, BOOLEAN,
   // MPINTEGER, MPRATIONAL, MPRATCOMPLEX, or CLIFFORD, but if it is it should
   // always return 0.
   int notanumber = 0;
#ifdef USE_REAL
   notanumber = (fpclassify(variable)==FP_NAN);
#elif defined(USE_COMPLEX)
   notanumber = (fpclassify(creall(variable))==FP_NAN);
   notanumber |= (fpclassify(cimagl(variable))==FP_NAN);
#elif defined(USE_MPREAL)
   notanumber = mpfr_nan_p(variable);
#elif defined(USE_MPCOMPLEX)
   notanumber = mpfr_nan_p(mpc_realref(variable));
   notanumber |= mpfr_nan_p(mpc_imagref(variable));
#endif

   return notanumber;
}

int testForInfinity(const scalarType variable)
{
   // this routine does not need to be called in the case of INTEGER, BOOLEAN,
   // MPINTEGER, MPRATIONAL, MPRATCOMPLEX, or CLIFFORD, but if it is it should
   // always return 0.

   int isinfinite = 0;
#ifdef USE_REAL
   isinfinite = (fpclassify(variable)==FP_INFINITE);
#elif defined(USE_COMPLEX)
   isinfinite = (fpclassify(creall(variable))==FP_INFINITE);
   isinfinite |= (fpclassify(cimagl(variable))==FP_INFINITE);
#elif defined(USE_MPREAL)
   isinfinite = mpfr_inf_p(variable);
#elif defined(USE_MPCOMPLEX)
   isinfinite = mpfr_inf_p(mpc_realref(variable));
   isinfinite |= mpfr_inf_p(mpc_imagref(variable));
#endif

   return isinfinite;
}

int testForZero(const scalarType variable)
{
   // this routine is a wrapper which contains the code needed for each
   // scalarType to test for zero; it returns 1 if variable == 0.
   int iszero = 0;
#ifdef USE_INTEGER
   iszero = (variable==0);
#elif defined(USE_BOOLEAN)
   iszero = (variable==0);
#elif defined(USE_MPINTEGER)
   // this is not the most efficient function, but does save the trouble
   // of having a multiprecision zero hanging around to compare with
   iszero = (mpz_sgn(variable)==0);
#elif defined(USE_REAL)
   iszero = (fpclassify(variable)==FP_ZERO);
#elif defined(USE_MPREAL)
   iszero = mpfr_zero_p(variable);
#elif defined(USE_COMPLEX)
   iszero = (fpclassify(creall(variable))==FP_ZERO);
   iszero &= (fpclassify(cimagl(variable))==FP_ZERO);
#elif defined(USE_MPCOMPLEX)
   iszero = mpfr_zero_p(mpc_realref(variable));
   iszero &= mpfr_zero_p(mpc_imagref(variable));
#elif defined(USE_MPRATIONAL)
   // this is not the most efficient function, but does save the trouble
   // of having a multiprecision zero hanging around to compare with
   iszero = (mpq_sgn(variable)==0);
#elif defined(USE_MPRATCOMPLEX)
   // this is not the most efficient function, but does save the trouble
   // of having a multiprecision zero hanging around to compare with
   iszero = (mpq_sgn(variable->real)==0);
   iszero &= (mpq_sgn(variable->imag)==0);
#elif defined(USE_CLIFFORD)
   iszero = larc_sca_is_real_clifford(variable); // 0 if nonzero imag
   iszero &= larc_sca_is_imag_clifford(variable); // 0 if nonzero real
#endif
   return iszero;
}

int testForRegular(const scalarType variable)
{
   // this routine is a wrapper which contains the code needed for each
   // scalarType to test that the given variable is not zero, NaN or
   // infinity; it returns 1 if this is true
   // "regular" is a standard definition in MPC/MPFR, but for C floats
   // it is either defined as below or by (!NaN & !infinity & !zero)
   int isregular = 0;
#ifdef USE_INTEGER
   isregular = (variable!=0);
#elif defined(USE_BOOLEAN)
   isregular = (variable!=0);
#elif defined(USE_MPINTEGER)
   // this is not the most efficient function, but does save the trouble
   // of having a multiprecision zero hanging around to compare with
   isregular = (mpz_sgn(variable)!=0);
#elif defined(USE_REAL)
   isregular = (fpclassify(variable)==FP_NORMAL) ||
	(fpclassify(variable)==FP_SUBNORMAL);
#elif defined(USE_MPREAL)
   isregular = mpfr_regular_p(variable);
#elif defined(USE_COMPLEX)
   isregular = (fpclassify(creall(variable))==FP_NORMAL) ||
	(fpclassify(creall(variable))==FP_SUBNORMAL);
   isregular &= (fpclassify(cimagl(variable))==FP_NORMAL) ||
	(fpclassify(cimagl(variable))==FP_SUBNORMAL);
#elif defined(USE_MPCOMPLEX)
   isregular = mpfr_regular_p(mpc_realref(variable));
   isregular &= mpfr_regular_p(mpc_imagref(variable));
#elif defined(USE_MPRATIONAL)
   // this is not the most efficient function, but does save the trouble
   // of having a multiprecision regular hanging around to compare with
   isregular = (mpq_sgn(variable)!=0);
#elif defined(USE_MPRATCOMPLEX)
   // this is not the most efficient function, but does save the trouble
   // of having a multiprecision regular hanging around to compare with
   isregular = (mpq_sgn(variable->real)!=0);
   isregular &= (mpq_sgn(variable->imag)!=0);
#elif defined(USE_CLIFFORD)
   isregular = !larc_sca_is_real_clifford(variable); // !0 if nonzero imag
   isregular |= !larc_sca_is_imag_clifford(variable); // !0 if nonzero real
#endif
   return isregular;
}

int larc_MPC_cmp(const scalarType a, const scalarType b)
{
/*  This routine uses the GNU MPC convention where:
 *	1) the returned int is real_flag + 4*imag_flag
 *	2) real_flag is 2 when Re(a)<Re(b), 0 when Re(a)==Re(b), 1 when
 *		Re(a)>Re(b) (that is, the normal -1 return of a comparison
 *		is converted to 2)
 *	3) imag_flag is 2 when Im(a)<Im(b), 0 when Im(a)==Im(b), 1 when
 *		Im(a)>Im(b) (that is, the normal -1 return of a comparison
 *		is converted to 2)
 *	For non-complex types, imag_flag == 0 and the returned value is
 *	real_flag.
 */

#ifdef USE_MPCOMPLEX
    return mpc_cmp(a,b); // use the MPC library function
#elif defined(USE_CLIFFORD)
  // deal with Clifford type differently, since it can be either COMPLEX or not
    clifford_t approx1, approx2;
    larc_sca_init_clifford(&approx1);
    larc_sca_init_clifford(&approx2);
    convert_clifford_to_mprational_approx(approx1, a);
    convert_clifford_to_mprational_approx(approx2, b);
    int real_flag, imag_flag;
    real_flag = mpq_cmp(approx1->real_coeffs[0],approx2->real_coeffs[0]);
  #ifdef IS_COMPLEX
    imag_flag = mpq_cmp(approx1->imag_coeffs[0],approx2->imag_coeffs[0]);
  #else
    imag_flag = 0;
  #endif
    // emulate the MPC library function
    if (real_flag < 0) real_flag = 2;
    if (imag_flag < 0) imag_flag = 2;
    return(real_flag + 4*imag_flag);
#elif defined(IS_COMPLEX)
  #if defined(USE_COMPLEX)
    int real_flag = larc_sca_cmp_real(creall(a),creall(b));
    int imag_flag = larc_sca_cmp_real(cimagl(a),cimagl(b));
  #elif defined(USE_MPRATCOMPLEX)
    int real_flag = mpq_cmp(a->real, b->real);
    int imag_flag = mpq_cmp(a->imag, b->imag);
  #endif // type-dependent calculation of real_flag and imag_flag
    // emulate the MPC library function
    if (real_flag < 0) real_flag = 2;
    if (imag_flag < 0) imag_flag = 2;
    return(real_flag + 4*imag_flag);
#else // non-complex comparisons straightforward, but we still emulate MPC
    int real_flag = sca_cmp(a,b);
    if (real_flag < 0) return 2;
    else return real_flag;
#endif // #ifdef USE_MPCOMPLEX
}
