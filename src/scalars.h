//                       scalars.h
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

#ifndef SCALARS_H
#define SCALARS_H

#ifdef __cplusplus
extern "C" {
#endif

// Standard Libaries
#include <stdint.h>
#include "larc.h"

// See larc.h for prototypes of scalar operations

// this is used to pass information to python about whether
// we are in MAR (1) or SPR (0) mode for locality hashing
extern int MARmode;

#ifndef SWIG

/*!
 * \ingroup larc
 * \brief Returns a GMP MPC-style comparison of two numbers
 *
 * MPC's comparison function returns real_flag + 4*imag_flag, in which the real
 * and imaginary parts of a and b are compared separately. Where the usual
 * comparison function would return -1, these comparisons return 2. Thus the
 * possible results are:
 *      0 (0000) : a == b
 *      1 (0001) : Re(a)>Re(b), Im(a)==Im(b)
 *      2 (0010) : Re(a)<Re(b), Im(a)==Im(b)
 *      4 (0100) : Re(a)==Re(b), Im(a)>Im(b)
 *      5 (0101) : Re(a)>Re(b), Im(a)>Im(b)
 *      6 (0110) : Re(a)<Re(b), Im(a)>Im(b)
 *      8 (1000) : Re(a)==Re(b), Im(a)<Im(b)
 *      9 (1001) : Re(a)>Re(b), Im(a)<Im(b)
 *     10 (1010) : Re(a)<Re(b), Im(a)<Im(b)
 *
 * LARC's version of this function operates on any defined scalarType and
 * returns results equivalent to what MPC would return if a and b were
 * compared after being converted to mpc_t variables. If LARC is compiled with
 * a non-complex scalarType, the returned value is consistent with
 * Im(a)==Im(b)==0.
 *
 * \param a The first scalarType to be compared
 * \param b The second scalarType to be compared
 * \result One of nine integers encoding the relations between the real parts
 *         of a and b and the imaginary parts of a and b
 */
int larc_MPC_cmp(const scalarType a, const scalarType b);


/*******************************************************************
 *  Prototypes for scalar operations that are parameter based.     *
 *   must be initialized before preloading matrix store            *
 *   If adding additional scalar operations, update scalars.c      *
 *   especially check_scalarOps().                                 *
 ******************************************************************/

/****************************************************************************
 * These functions ensure that multiprecision types (which usually allocate *
 * memory for structures) are properly created, initialized and freed. When *
 * using standard C scalar types, these functions are usually just wrappers *
 * for C operations.                                                        * 
 ***************************************************************************/

/*!
 * \ingroup larc
 * \brief Allocates memory for a new variable of type scalarType (if needed for that scalarType)
 * \param s_ptr A pointer to scalarType (to hold the new variable)
 */
extern void     (*sca_init)    (scalarType* s_ptr);

/*!
 * \ingroup larc
 * \brief Deallocates memory for a variable of type scalarType (if needed for that scalarType)
 * \param s_ptr A pointer to scalarType to be freed
 */
extern void     (*sca_clear)   (scalarType* s_ptr);

/*!
 * \ingroup larc
 * \brief Copies the value of one scalarType to another
 * \param d_ptr A pointer to a scalarType value
 * \param s_ptr A pointer to a scalarType value
 */
extern void     (*sca_set)     (scalarType* d_ptr, const scalarType s_ptr);

/*!
 * \ingroup larc
 * \brief Converts a string into the appropriate scalarType and stores that value
 * \param s_ptr The scalarType value to be written to the store
 * \param input_string The string-form representation to the scalarType value
 */
extern void     (*sca_set_str) (scalarType* s_ptr, const char *input_str);

/*!
 * \ingroup larc
 * \brief Converts two double precision floats into the appropriate scalarType and stores that value
 * \param real_val The real part of the value
 * \param imag_val The imaginary part of the value (if the scalarType is not a complex type, this must be zero)
 */
extern void     (*sca_set_2ldoubles) (scalarType*, long double real_val, long double imag_val);

/*!
 * \ingroup larc
 * \brief Sets the scalar to the preset value specified by the enum
 * \param enum_val Which preset value to use
 */
extern void     (*sca_set_enum) (scalarType*, const enum sca_constant_spec enum_val);

/*!
 * \ingroup larc
 * \brief Converts a scalarType into string format and returns that string; this routine uses malloc()
 * \param s A scalarType variable
 * \return The string form of the value of s
 */
extern char    *(*sca_get_readable_approx_str) (const scalarType s);

/*!
 * \ingroup larc
 * \brief Converts a scalarType into an exact string format and returns that string; this routine uses malloc()
 * \param s A scalarType variable
 * \return The string form of the exact value of s
 */
extern char    *(*sca_get_exact_str) (const scalarType s);

/*!
 * \ingroup larc
 * \brief Hashes a scalarType to a 64-bit integer according to a predefined function
 * \param s A scalarType variable
 * \param exp The log base 2 of the size of the hash table
 * \return The hash value of the input scalar
 */
extern uint64_t (*sca_hash)    (const scalarType s, uint64_t exp);

/*!
 * \ingroup larc
 * \brief Adds two scalarType variables and returns the result
 * \param s A pointer to the sum of the two variables
 * \param s1 The first of two scalarType variables to be added
 * \param s2 The second of two scalarType variables to be added
 */
extern void     (*sca_add)   (scalarType* s, const scalarType s1, const scalarType s2);

/*!
 * \ingroup larc
 * \brief Multiplies two scalarType variables and returns the result
 * \param s A pointer to the product of the two variables
 * \param s1 The first of two scalarType variables to be multiplied
 * \param s2 The second of two scalarType variables to be multiplied
 */
extern void     (*sca_mult)  (scalarType* p, const scalarType s1, const scalarType s2);

/*!
 * \ingroup larc
 * \brief Divides one scalarType variable by another and returns the result
 * \param s A pointer to the result of division 
 * \param n The scalarType variable to be the numerator
 * \param d The scalarType variables to be the denominator
 */
extern void     (*sca_divide)  (scalarType* s, const scalarType n, const scalarType d);

#ifndef SWIG
/*!
 * \ingroup larc
 * \brief Determines which SPRregion or MARtile a given scalar lies within
 *
 * When LARC is compiled in SPR mode, this function finds the region center for
 * the SPR region that the scalar belongs to, and returns a pointer to that 
 * scalarType value. When in MAR mode, this function finds the MARtile index
 * for the MAR region that the scalar belongs to, and returns a pointer to
 * that MAR_tile_index_t value. The returned value must be cast appropriately.
 *
 * \param ptr A void * pointer to the SPRregion or MARtile information
 * \param tile_exp In MARmode, determines the size of the tile (unused in SPRmode)
 * \param input The input scalarType
 */
void     return_SPRregion_or_MARtile_label(void* ptr,
                 int tile_exp, const scalarType input);
#endif // not SWIG

/**********************************************************************
 * sca_cmp: The function's return value is equal to zero, if the two  *
 * scalars are equal; and nonzero otherwise. Negative vs positive     *
 * nonzero values aren't required for LARC yet.                       *
 * If they were, return a positive value if first is greater than     *
 * second, and a negative value if first is less than second.         *
***********************************************************************/
/*!
 * \ingroup larc
 * \brief Compares two scalarType values
 * \param a The first scalarType variable
 * \param b The second scalarType variable
 * \return -1/0/1 if the first is less than/equal to/greater than the second
 */
extern int      (*sca_cmp)     (const scalarType a, const scalarType b); //

/*!
 * \ingroup larc
 * \brief Checks for exact equality of two scalarType variables
 * \param a The first scalarType variable
 * \param b The second scalarType variable
 * \return 1 if the the two scalarType values are equal, 0 otherwise
 * 
 */
extern int      (*sca_eq)   (const scalarType a, const scalarType b); //


/*!
 * \ingroup larc
 * \brief Returns the sqrt of the scalarType variable
 * \param a A pointer to the sqrt of the input scalarType variable
 * \param b The input scalarType variable
 */
extern void     (*sca_sqrt)    (scalarType* a, const scalarType b);
/*!
 * \ingroup larc
 * \brief Returns the norm of the scalarType variable
 * \param a A pointer to the norm of the input scalarType variable
 * \param b The input scalarType variable
 */
extern void     (*sca_norm)    (scalarType* a, const scalarType b);
/*!
 * \ingroup larc
 * \brief Returns the conjugate of the scalarType variable
 * \param a A pointer to the conjugate of the input scalarType variable
 * \param b The input scalarType variable
 */
extern void     (*sca_conj)    (scalarType* a, const scalarType b);
/*!
 * \ingroup larc
 * \brief Tests whether the (complex) input has an imaginary part
 * \param a The input scalarType variable
 * \return 1 if the scalarType variable has no imaginary part, 0 otherwise
 */
extern int      (*sca_is_pure_real)  (const scalarType a);

/*!
 * \ingroup larc
 * \brief Tests whether the (complex) input has a real part
 * \param a The input scalarType variable
 * \return 1 if the scalarType variable has no real part, 0 otherwise
 */
extern int      (*sca_is_pure_imag)  (const scalarType a);


#endif // ifndef SWIG


#ifndef SWIG
/*!
 * \ingroup larc
 * \brief Tells LARC to use the usual scalar operations
 * \param verbose an enum value in {SILENT, BASIC, CHATTY, DEBUG}; determines which messsages are printed
 */
void init_scalarOps(int verbose);

/*!
 * \ingroup larc
 * \brief Hashes a multiprecision integer to an integer in range [0,2**exponent)
 * \param sc The multiprecision integer
 * \param exponent Determines the size of the hash table (2**exponent)
 * \return The hash of the input multiprecision integer value
 */
uint64_t larc_mpz_hash(const mpz_t sc, uint64_t exponent);

  
/*!
 * \brief converts a single multiprecision floating point number to the compiled scalarType, rounding to nearest supported value when necessary
 * \param output A pointer to the scalarType value to be set
 * \param input A mpfr_t multiprecision floating point value
 */
void convert_from_mpfr_to_scalarType(scalarType *output, mpfr_t input);

/*!
 * \brief converts two multiprecision floating point numbers to the compiled complex scalarType, rounding to nearest supported value when necessary
 *
 * \param output A pointer to the scalarType value to be set
 * \param realpart A mpfr_t multiprecision floating point value
 * \param imagpart A mpfr_t multiprecision floating point value
 */
void convert_from_two_mpfr_to_complex_scalarType(scalarType *output,
        mpfr_t realpart, mpfr_t imagpart);


// These functions are used in global.c so must be in scalars.h header
#ifdef USE_CLIFFORD
/*!
 * \ingroup LARC
 * \brief This allocates and initializes a Clifford scalarType variable
 * \param sc A pointer to the scalar to be initialized
 *
 */
void larc_sca_init_clifford(clifford_t *sc);

/*!
 * \ingroup LARC
 * \brief This frees the memory allocated to a Clifford scalarType variable
 * \param sc A pointer to the scalar to be freed
 *
 */
void larc_sca_clear_clifford(clifford_t *sc);
#endif

/*!
 * \brief Tests whether a scalarType variable contains NaN
 * \param variable The scalarType variable to be tested
 * \return 1 if the argument is NaN, 0 otherwise
 */
int testForNaN(const scalarType variable);

/*!
 * \brief Tests whether a scalarType variable contains infinity
 * \param variable The scalarType variable to be tested
 * \return nonzero if the argument is infinity, 0 otherwise
 */
int testForInfinity(const scalarType variable);

/*!
 * \brief Tests whether a scalarType variable contains zero
 * \param variable The scalarType variable to be tested
 * \return 1 if the value pointed at is zero, 0 otherwise
 */
int testForZero(const scalarType variable);

/*!
 * \brief Tests whether a scalarType variable contains a regular value
 * \param variable The scalarType variable to be tested
 * \return 0 if the value pointed at is zero, NaN or infinity, 1 otherwise
 */
int testForRegular(const scalarType variable);

#endif  // ifndef SWIG

/*!
 * \brief Finds/stores the packedID for one of the values specified in the sca_constant_spec enumeration
 *
 * The values in this enumeration are essential for the predefined Clifford
 * types, and will properly set the coefficient for any algebraic value which
 * extends the rational field for that Clifford algebra, but use of this
 * function is not restricted to Clifford types. For * example, if the user
 * desires a good approximation to 1/sqrt(2) in any scalarType, this function
 * will return it. (For integer types, the value 1 is the best approximation
 * possible for 1/sqrt(2), and that is what is returned.)
 *
 * \param A sca_constant_spec enum value e.g. SCALAR_ENUM_INV_SQRT2
 * \return The packedID of the stored value
 */
int64_t get_pID_for_enum_const(const enum sca_constant_spec enum_val);

#ifdef __cplusplus
}
#endif


#endif   // #define SCALARS_H
