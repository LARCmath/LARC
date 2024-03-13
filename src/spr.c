 //                         spr.c
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

#include "spr.h"
#include <signal.h>
#include "clifford.h"

#ifndef MAR // SPRmode

#ifndef SWIG

#ifdef IS_RATIONAL


/*!
 * \ingroup larc
 * \brief Compares a multiprecision rational input to a threshold value; returns 0 if it is less than the threshold, returns the input if it is not
 * \param rsc The multiprecision rational result of the function
 * \param sc The multiprecision rational input
 */
static void collapse_near_zero_mprational(mpq_t *rsc, const mpq_t sc)
{
  const mpq_t* const zeromprationalthresh = get_zerorealthresh();

  if (scratchVars.mprational_in_use)
      fprintf(stderr,"%s reusing scratchVars.mprational!\n",__func__);
  scratchVars.mprational_in_use = 1;

  mpq_t *abs_sc = &scratchVars.mprational;
  mpq_abs(*abs_sc, sc);

  if (mpq_cmp(*abs_sc, *zeromprationalthresh) < 0) {
    mpq_set_ui(*rsc, 0, 1);
  } else {
    mpq_set(*rsc, sc);
  }
  scratchVars.mprational_in_use = 0;
}

/*!
 * \ingroup larc
 * \brief Returns a scalar value that defines the SPR region or MAR tile containing the multiprecision rational input (modified later if special treatment near zero).
 * \param rsc The multiprecision rational result of the function
 * \param sc The multiprecision rational input
 */
static void round_sig_fig_mprational(mpq_t *rsc, const mpq_t sc)
{
    int regionbitparam = get_regionbitparam();

    // Multiply by 2^regionbitparam:
    if (scratchVars.mprational_in_use)
        fprintf(stderr,"%s reusing scratchVars.mprational!\n",__func__);
    scratchVars.mprational_in_use = 1;

    mpq_t *temp_q = &scratchVars.mprational;

    mpq_mul_2exp(*temp_q, sc, regionbitparam);  // shifts rational left by regionbitparam bits
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
    if (scratchVars.mpinteger_in_use)
        fprintf(stderr,"%s reusing scratchVars.mpinteger!\n",__func__);
    scratchVars.mpinteger_in_use = 1;

    mpz_set_q(scratchVars.mpinteger, *rsc);   // takes rational, sets integer
    mpq_set_z(*rsc, scratchVars.mpinteger);   // takes integer, sets rational
    scratchVars.mpinteger_in_use = 0;

    // Divide by 2^regionbitparam:
    mpq_div_2exp(*temp_q, *rsc, regionbitparam);  // shifts rational right by regionbitparam bits
    mpq_set(*rsc,*temp_q);
    scratchVars.mprational_in_use = 0;
}

#ifdef USE_CLIFFORD

static void round_sig_fig_clifford(clifford_t *rsc, const clifford_t sc)
{
   // use rsc to hold approximate value
   // all but real_coeffs[0] (and possibly imag_coeffs[0]) are 0/1
   convert_clifford_to_mprational_approx(*rsc,sc);

   // put rounded value into temp_q, copy over rsc afterward
   // note that round_sig_fig_mprational uses scratchVars.mprational
   if (scratchVars.mprational2_in_use)
       fprintf(stderr,"%s reusing scratchVars.mprational2!\n",__func__);
   scratchVars.mprational2_in_use = 1;

   mpq_t *temp_q = &scratchVars.mprational2;
   round_sig_fig_mprational(temp_q,(*rsc)->real_coeffs[0]);
   mpq_set((*rsc)->real_coeffs[0],*temp_q);
#ifdef IS_COMPLEX
   round_sig_fig_mprational(temp_q,(*rsc)->imag_coeffs[0]);
   mpq_set((*rsc)->imag_coeffs[0],*temp_q);
#endif // IS_COMPLEX
   scratchVars.mprational2_in_use = 0;
   return;
}

/*!
 * \brief Zeros coefficients of Clifford values if they fall below a threshold
 *
 * This routine is only called after round_sig_fig_clifford, so the input
 * Clifford value has all algebraic coefficients (except \sqrt{-1} for complex
 * types) set to zero.
 *
 * For real Clifford algebras, the multiprecision rational input is compared to
 * a threshold value and set to zero if smaller. For complex Clifford algebras,
 * the real and imaginary parts of the input are treated independently.
 *
 * \param rsc A pointer to the Clifford-type result of the function
 * \param sc The Clifford-valued input
 *
 */
static void collapse_near_zero_clifford(clifford_t *rsc, const clifford_t sc)
{
   // note that collapse_near_zero_mprational uses scratchVars.mprational
   if (scratchVars.mprational2_in_use)
       fprintf(stderr,"%s reusing scratchVars.mprational2!\n",__func__);
   scratchVars.mprational2_in_use = 1;

   mpq_t *temp_q = &scratchVars.mprational2;

   collapse_near_zero_mprational(temp_q,sc->real_coeffs[0]);
   mpq_set((*rsc)->real_coeffs[0],*temp_q);
#ifdef IS_COMPLEX
   collapse_near_zero_mprational(temp_q,sc->imag_coeffs[0]);
   mpq_set((*rsc)->imag_coeffs[0],*temp_q);
#endif // IS_COMPLEX
   scratchVars.mprational2_in_use = 0;
   return;
}
#endif // USE_CLIFFORD

#endif //# IS_RATIONAL

#if defined(USE_MPREAL) || defined(USE_MPCOMPLEX)

/*!
 * \ingroup larc
 * \brief Returns a scalar value that defines the SPR region or MAR tile containing the multiprecision floating point input (modified later if special treatment near zero).
 * \param rsc The multiprecision real result of the function
 * \param sc The multiprecision real input
 */
static void round_sig_fig_mpreal(mpfr_t *rsc, const mpfr_t sc)
{
    int regionbitparam = get_regionbitparam();

    //printf("in %s\n",__func__);
    //mpfr_printf("input value sc is %.128Rf\n", sc);

    // Multiply by 2^regionbitparam:
    if (scratchVars.mpreal_in_use)
        fprintf(stderr,"%s reusing scratchVars.mpreal!\n",__func__);
    scratchVars.mpreal_in_use = 1;

    mpfr_t *temp_q = &scratchVars.mpreal;
    mpfr_mul_2exp(*temp_q, sc, regionbitparam, MPFR_RNDN);  // shifts real left by regionbitparam bits
    mpfr_set(*rsc, *temp_q, MPFR_RNDN);
    //mpfr_printf("after multiply: value rsc is %.128Rf\n", *rsc);

    if (scratchVars.mpinteger_in_use)
        fprintf(stderr,"%s reusing scratchVars.mpinteger!\n",__func__);
    scratchVars.mpinteger_in_use = 1;

#ifdef ANOTHER_VERSION
    // round multiprecision real to nearest integer
    mpfr_get_z(scratchVars.mpinteger, *rsc, MPFR_RNDN);
#else
    // Change from Mark's original code: add \pm0.5, to change subsequent
    // truncation into round towards 0
    mpfr_set_d(*temp_q, 0.5, MPFR_RNDN);
    if (mpfr_cmp_ui(*rsc, 0) >= 0) {
        mpfr_add(*rsc, *rsc, *temp_q, MPFR_RNDN);
    } else {
        mpfr_sub(*rsc, *rsc, *temp_q, MPFR_RNDN);
    }
    //mpfr_printf("after add/sub: value rsc is %.128Rf\n", *rsc);

    // Truncate *rsc using RNDZ (round towards zero) rather than
    // RNDN (round towards nearest)
    mpfr_get_z(scratchVars.mpinteger, *rsc, MPFR_RNDZ);   // takes real, sets integer
#endif // ANOTHER_VERSION

    mpfr_set_z(*rsc, scratchVars.mpinteger, MPFR_RNDZ);   // takes integer, sets real
    scratchVars.mpinteger_in_use = 0;
    //mpfr_printf("after truncation: value rsc is %.128Rf\n", *rsc);

    // Divide by 2^regionbitparam:
    mpfr_div_2exp(*temp_q, *rsc, regionbitparam, MPFR_RNDN);  // shifts real right by regionbitparam bits
    mpfr_set(*rsc, *temp_q, MPFR_RNDN);
    //mpfr_printf("after rescale: value rsc is %.128Rf\n", *rsc);
    scratchVars.mpreal_in_use = 0;
}

/*!
 * \ingroup larc
 * \brief Compares the multiprecision real input to a threshold value; returns 0 if it is less than the threshold, returns the input if it is not
 * \param rsc The multiprecision real result of the function
 * \param sc The multiprecision real input
 */
static void collapse_near_zero_mpreal(mpfr_t *rsc, const mpfr_t sc)
{
  if (VERBOSE>DEBUG) printf("in %s\n",__func__);
  const mpfr_t* const zeromprealthresh = get_zerorealthresh();

  if (scratchVars.mpreal_in_use)
      fprintf(stderr,"%s reusing scratchVars.mpreal!\n",__func__);
  scratchVars.mpreal_in_use = 1;
  mpfr_t *abs_sc = &scratchVars.mpreal;
  mpfr_abs(*abs_sc, sc, MPFR_RNDN);
  scratchVars.mpreal_in_use = 0;

#ifdef DEBUG_SPR_C
  char *larc_sca_get_str_mpreal(const mpfr_t);
  char *temp1 = larc_sca_get_str_mpreal(*zeromprealthresh);
  char *temp2 = larc_sca_get_str_mpreal(*abs_sc);
  printf("threshold is %s\n",temp1);
  printf("abs(scalar) is %s\n",temp2);
  free(temp1);
  free(temp2);
#endif // DEBUG_SPR_C

  if (mpfr_cmp(*abs_sc, *zeromprealthresh) < 0) {
    if (VERBOSE>DEBUG) printf("\tsetting scalar to zero\n\n");
    mpfr_set_ui(*rsc, 0, MPFR_RNDN);
  } else {
    if (VERBOSE>DEBUG) printf("\tnot setting scalar to zero\n\n");
    mpfr_set(*rsc, sc, MPFR_RNDN);
  }
  if (VERBOSE>DEBUG) printf("exiting %s\n",__func__);
#undef DEBUG_SPR_C
}

#endif // USE_MPREAL or USE_MPCOMPLEX  

#if defined(USE_REAL) || defined(USE_COMPLEX)

/*!
 * \ingroup larc
 * \brief Returns a scalar value that defines the SPR region or MAR tile containing the double precision floating point input (modified later if special treatment near zero).
 * \param output The double precision floating point result of the function
 * \param input The double precision floating point input
 */
void round_sig_fig_real(long double *output, const long double input)
{
// This uses standard C functions (frexp, modf, and ldexp) to perform rounding
// in way that is independent of the specific floating point representation
// that is used by the underlying machine architecture.
  int regionbitparam = get_regionbitparam();

  if (isnan(input) ) {
    fprintf(stderr,"Inside %s a NaN showed up\n", __func__);
    raise(SIGSEGV);
    exit (1);
  }

  // Shift the input left by the number of significant bits to keep
  //   using the ldexp() function.
  long double shifted_real = ldexpl(input, regionbitparam);
  // Post-condition: shifted_real == input * (2 ^ regionbitparam)

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
  long double rounded_answer = ldexpl(shifted_real_int, - regionbitparam);
  // Post-condition: rounded_answer == shifted_real_int * (2 ^ (- regionbitparam))

  // printf("Hash debug:  input = %25.16la   output = %25.16la   regionbitparam = %d\n", input, rounded_answer, regionbitparam);

  *output = rounded_answer;
}

/*!
 * \ingroup larc
 * \brief Compares the double precision input to a threshold value; returns 0 if it is less than the threshold, returns the input if it is not
 * \param output The double precision floating point result of the function
 * \param input The double precision floating point input
 */
static void collapse_near_zero_real(long double *output, const long double input)
{
  if (VERBOSE>DEBUG) printf("in %s\n",__func__);
  long double threshold = get_zerorealthresh();
  long double abs_input = fabsl(input);

#ifdef DEBUG_SPR_C
  char *larc_sca_get_str_real(const long double);
  char *temp1 = larc_sca_get_str_real(threshold);
  char *temp2 = larc_sca_get_str_real(abs_input);
  printf("threshold is %s\n",temp1);
  printf("abs(scalar) is %s\n",temp2);
  free(temp1);
  free(temp2);
#endif // DEBUG_SPR_C

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

#endif // USE_REAL or USE_COMPLEX

/*!
 * \ingroup larc
 * \brief Returns a scalar value that defines the SPR region or MAR tile containing the scalarType input (modified later if special treatment near zero).
 * \param output The scalarType output of the function
 * \param input The scalarType input to the function
 */
static void round_sig_fig(scalarType *output, const scalarType input)
{
#if defined(USE_INTEGER) || defined(USE_BOOLEAN)
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
#elif defined(USE_CLIFFORD)
   round_sig_fig_clifford(output,input);
#else
   fprintf(stderr,"no valid scalarType for round_sig_fig!\n");
   raise(SIGSEGV);
   exit(1);
#endif
   return;
}

/*!
 * \ingroup larc
 * \brief Compares the input to a threshold value; returns 0 if the input is less than the threshold, or the input value if it is not. Integer types do not collapse.
 * \param output The scalarType result of the function
 * \param input The scalarType input to the function
 */
static void collapse_near_zero(scalarType *output, const scalarType input)
{
#if defined(USE_INTEGER) || defined(USE_BOOLEAN)
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
#elif defined(USE_CLIFFORD)
   collapse_near_zero_clifford(output,input);
#else
   fprintf(stderr,"no valid scalarType for collapse_near_zero!\n");
   raise(SIGSEGV);
   exit(1);
#endif
   return;
}

void return_SPR_region_center(scalarType *output, const scalarType input)
{
   round_sig_fig(output,input);
   // if the zero regions requested are larger than those
   // created by round_sig_fig then calculate these regions

   if (is_special_zeroregion()) {
     // ERROR CHECKING
     // fprintf(stderr,"We are case is_special_zero_region for scalar\n");
     // fprintf(stderr,"zerorealthresh is %lg\n",get_zerorealthresh());
     collapse_near_zero(output,*output);
   }
   return;
}

mats_ptr_t get_PTR_scalar_record(scalarType scalar)
{
  int verbose = 0;
  uint64_t hash = region_hash_from_scalarType(scalar, get_scalar_store_exp());
  if (verbose) {
    printf("The hash returned from region_hash_from_scalarType is %zd\n",hash);
  }

  mats_ptr_t s_ptr = find_scalar_SPRmode(hash, scalar);
  if (verbose) {
    printf("  called find_scalar_SPRmode, result (s_ptr) is %p\n", s_ptr);
  }
  
  // if unable to find, attempt to insert
  if (s_ptr == SCALAR_PTR_INVALID) {
    if (verbose) { 
      printf("Inside %s, SCALAR ptr not found by matrix_find\n", __func__);
      printf("   CALLING  insert_scalar_SPRmode\n");
    }
    s_ptr = insert_scalar_SPRmode(hash, scalar);
  }
  
  if (verbose) {   
    printf("    RETURNED from insert_scalar_SPRmode with scalar ptr %p\n", s_ptr);
  }
  return s_ptr;
}

int invalidate_recordPTR_in_indexTable(int64_t packedID);

/************************************************************************
 * Clean matrix store	                                                *
 *                                                                      *
 *  Removes all eligible matrices from the matrix store.                *
 *      The called function, remove_matrix_from_store, will         *
 *      recursively remove eligible children of any deleted matrix      *      	
 ***********************************************************************/
/*!
 * \ingroup larc
 * \brief removes all eligible matrices from the scalar store.
 *
 * This is the routine used to clean when LARC is run in SPRmode. A different
 * routine will be used when LARC is run in MARmode.
 *
 * \result Returns 1.
 */
int clean_scalar_matrices_SPR(struct matrix_store_t * matrix_store_ptr)
{
    int verbose = 0;
    uint64_t hash_value;

    // some scalars may not be in a larger matrix (e.g. those stored from
    // functions on matrices like norms or traces); we clean those too
    uint64_t max = (uint64_t)1 << get_scalar_store_exp();

    // loop though hash chains for each hashID without calling
    // clean_op_hash_chain so we don't have to repeat checks for invalid hashID
    hash_table_t *table_ptr = matrix_store_ptr->scalar_hash_table;
    for (hash_value = 0; hash_value < max; ++hash_value) {
        hash_node_t *node_ptr = table_ptr->heads[hash_value];

        while (node_ptr) {
            int64_t m_pID = node_ptr->packedID;
            if (m_pID == MATRIX_ID_INVALID)
            {
                fprintf(stderr,"ERROR in %s:\n",__func__);
                fprintf(stderr,"found scalarStore node with invalid packedID ");
                fprintf(stderr,"while traversing hash chain %lu\n",hash_value);
//                fprintf(stderr,"Removing the node, but ");
//                fprintf(stderr,"it's likely that cleaning is bugged.\n");
//                node_ptr = hash_node_remove_node(matrix_store.scalar_hash_table,
//                        node_ptr, hash_value);
//                continue;
		fprintf(stderr,"The scalar store is messed up, exiting.\n");
		exit(1);
            }
            record_ptr_t r_ptr = get_recordPTR_from_pID(m_pID,"",__func__,0);
            if (r_ptr == RECORD_PTR_INVALID)
            {
                fprintf(stderr,"ERROR in %s:\n",__func__);
                fprintf(stderr,"packedID %lu with null record ",m_pID);
                fprintf(stderr,"found while traversing scalarStore ");
                fprintf(stderr,"hash chain %lu\n", hash_value);
//                fprintf(stderr,"Removing the node, but ");
//                fprintf(stderr,"it's likely that cleaning is bugged.\n");
//                node_ptr = hash_node_remove_node(matrix_store_ptr->scalar_hash_table,
//                        node_ptr, hash_value);
//                continue;
		fprintf(stderr,"The scalar store is messed up, exiting.\n");
		exit(1);
            }
            else
            {
                if (verbose)
                {
                    fprintf(stderr,"remove_matrix_from_scalar_store_by_");
                    fprintf(stderr,"recordPTR_andor_nodePTR\n");
                    fprintf(stderr,"called from %s with pID %ld\n",
                            __func__, m_pID);
                    fprintf(stderr,"since this routine is NOT recursive, the");
                    fprintf(stderr,"next node pointer is n->next;n");
                }
                node_ptr = 
                remove_matrix_from_scalar_store_by_recordPTR_andor_nodePTR(
                    r_ptr, node_ptr, hash_value, m_pID);
            }
        }
    }

    return 1;
}

#endif

#endif
