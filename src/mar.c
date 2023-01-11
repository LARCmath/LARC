 //                         mar.c
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

#include "mar.h"

#ifdef MAR

MAR_tile_index_t *GLOBAL_primary_tile_index_PTR;
MAR_tile_index_t *GLOBAL_target_tile_index_PTR;
MAR_tile_index_t *GLOBAL_nhbr_tile_index_PTR;
MAR_tile_index_t *GLOBAL_tiny_tile_index_PTR;
char Gpti_in_use;
char Ggti_in_use;
char Gnti_in_use;
char Gtti_in_use;

void scratchVars_mar_init() {
  GLOBAL_primary_tile_index_PTR = (MAR_tile_index_t *)malloc(sizeof(MAR_tile_index_t));
  if (GLOBAL_primary_tile_index_PTR == NULL) { ALLOCFAIL(); }
  GLOBAL_target_tile_index_PTR = (MAR_tile_index_t *)malloc(sizeof(MAR_tile_index_t));
  if (GLOBAL_target_tile_index_PTR == NULL) { ALLOCFAIL(); }
  GLOBAL_nhbr_tile_index_PTR = (MAR_tile_index_t *)malloc(sizeof(MAR_tile_index_t));
  if (GLOBAL_nhbr_tile_index_PTR == NULL) { ALLOCFAIL(); }
  GLOBAL_tiny_tile_index_PTR = (MAR_tile_index_t *)malloc(sizeof(MAR_tile_index_t));
  if (GLOBAL_tiny_tile_index_PTR == NULL) { ALLOCFAIL(); }
#ifdef IS_COMPLEX
  mpz_init(GLOBAL_primary_tile_index_PTR->real_index);
  mpz_init(GLOBAL_primary_tile_index_PTR->imag_index);
  mpz_init(GLOBAL_target_tile_index_PTR->real_index);
  mpz_init(GLOBAL_target_tile_index_PTR->imag_index);
  mpz_init(GLOBAL_nhbr_tile_index_PTR->real_index);
  mpz_init(GLOBAL_nhbr_tile_index_PTR->imag_index);
  mpz_init(GLOBAL_tiny_tile_index_PTR->real_index);
  mpz_init(GLOBAL_tiny_tile_index_PTR->imag_index);
#else  
  mpz_init(GLOBAL_primary_tile_index_PTR->index);
  mpz_init(GLOBAL_target_tile_index_PTR->index);
  mpz_init(GLOBAL_nhbr_tile_index_PTR->index);
  mpz_init(GLOBAL_tiny_tile_index_PTR->index);
#endif // #ifdef IS_COMPLEX
  Gpti_in_use = Ggti_in_use = Gnti_in_use = Gtti_in_use = 0;
}

void scratchVars_mar_clear() {
  // This routine should only be called as part of a LARC shutdown
  // procedure.
  fprintf(stderr,"freeing memory for tile indices\n");
#ifdef IS_COMPLEX
  mpz_clear(GLOBAL_primary_tile_index_PTR->real_index);
  mpz_clear(GLOBAL_primary_tile_index_PTR->imag_index);
  mpz_clear(GLOBAL_target_tile_index_PTR->real_index);
  mpz_clear(GLOBAL_target_tile_index_PTR->imag_index);
  mpz_clear(GLOBAL_nhbr_tile_index_PTR->real_index);
  mpz_clear(GLOBAL_nhbr_tile_index_PTR->imag_index);
  mpz_clear(GLOBAL_tiny_tile_index_PTR->real_index);
  mpz_clear(GLOBAL_tiny_tile_index_PTR->imag_index);
#else
  mpz_clear(GLOBAL_primary_tile_index_PTR->index);
  mpz_clear(GLOBAL_target_tile_index_PTR->index);
  mpz_clear(GLOBAL_nhbr_tile_index_PTR->index);
  mpz_clear(GLOBAL_tiny_tile_index_PTR->index);
#endif // #ifdef IS_COMPLEX
  free(GLOBAL_primary_tile_index_PTR);
  free(GLOBAL_target_tile_index_PTR);
  free(GLOBAL_nhbr_tile_index_PTR);
  free(GLOBAL_tiny_tile_index_PTR);
}

// returns 1 if the tiles contain the same index values.
int tile_indices_equal(MAR_tile_index_t *tile1_index_PTR, MAR_tile_index_t *tile2_index_PTR) {
#ifdef IS_COMPLEX
  if ((mpz_cmp(tile1_index_PTR->real_index,tile2_index_PTR->real_index)==0) &&
      (mpz_cmp(tile1_index_PTR->imag_index,tile2_index_PTR->imag_index)==0)){
    return 1;
  }
#else
  if (mpz_cmp(tile1_index_PTR->index,tile2_index_PTR->index)==0) {
    return 1;
  }
#endif // #ifdef IS_COMPLEX
  return 0;
}


// returns 1 if primary tile index + tile_offset is the goal tile index.
int tile_indices_equal_after_offset(MAR_tile_index_t *primary_index_PTR,
                                    MAR_tile_index_t *target_index_PTR,
                                    unsigned int tile_offset_flag) {

    static mpz_t temp_mpz;
    static int temp_is_initialized = 0;
    if (!temp_is_initialized) {
        mpz_init(temp_mpz);
        temp_is_initialized = 1;
    }

    // Check if no offset, so just looking for primary.
    if (tile_offset_flag == 0) {
#ifdef IS_COMPLEX
        if (   (mpz_cmp(primary_index_PTR->real_index, target_index_PTR->real_index) == 0)
            && (mpz_cmp(primary_index_PTR->imag_index, target_index_PTR->imag_index) == 0) ) {
            return 1;
        }
#else
        if (mpz_cmp(primary_index_PTR->index, target_index_PTR->index)==0) {
            return 1;
        }
#endif // #ifdef IS_COMPLEX
        // Return 0 because tile indices do not match.
        return 0;

    } else {
#ifdef IS_COMPLEX
        // Translate offset flag to actual offset (for real).
        int real_correct_offset = 0;
        if ((tile_offset_flag & 3) == 1) {
            real_correct_offset = 1;
        } else if ((tile_offset_flag & 3) == 2) {
            real_correct_offset = -1;
        }

        // Translate offset flag to actual offset (for imaginary).
        int imag_correct_offset = 0;
        if ((tile_offset_flag >> 2) == 1) {
            imag_correct_offset = 1;
        } else if ((tile_offset_flag >> 2) == 2) {
            imag_correct_offset = -1;
        }

        // If tile indices match (both real and imaginary parts)
        //   after adjusting for offset, return 1.
        mpz_sub(temp_mpz, target_index_PTR->real_index, primary_index_PTR->real_index);
        if (mpz_cmp_si(temp_mpz, real_correct_offset) == 0) {
            mpz_sub(temp_mpz, target_index_PTR->imag_index, primary_index_PTR->imag_index);
            if (mpz_cmp_si(temp_mpz, imag_correct_offset) == 0) {
                return 1;
            }
        }
#else
        // Translate offset flag to actual offset.
        int correct_offset = 0;
        if (tile_offset_flag == 1) {
            correct_offset = 1;
        } else if (tile_offset_flag == 2) {
            correct_offset = -1;
        }

        // If tile indices match after adjusting for offset, return 1.
        mpz_sub(temp_mpz, target_index_PTR->index, primary_index_PTR->index);
        if (mpz_cmp_si(temp_mpz, correct_offset)==0) {
            return 1;
        }
#endif // #ifdef IS_COMPLEX

        // Return 0 because tile indices do not match.
        return 0;
    }
}

void print_tile_index(MAR_tile_index_t *tile_index_PTR)
{
    printf("the tile index is\n");
#ifdef IS_COMPLEX
    gmp_printf("\t (%Zd,",tile_index_PTR->real_index);
    gmp_printf("%Zd)\n",tile_index_PTR->imag_index);
#else
    gmp_printf("\t %Zd\n",tile_index_PTR->index);
#endif // #ifdef IS_COMPLEX
}

void set_tile_index(MAR_tile_index_t *tile_index_PTR, int i, int j)
{
#ifdef IS_COMPLEX
   mpz_set_si(tile_index_PTR->real_index,i);
   mpz_set_si(tile_index_PTR->imag_index,j);
#else
   mpz_set_si(tile_index_PTR->index,i);
#endif // #ifdef IS_COMPLEX
}


#ifdef IS_RATIONAL
static void get_tile_index_mprational(mpz_t tile_index, int tile_exp,
                               const mpq_t scalar)
{
    // Multiply by 2^regionbitparam (by shifting region left by tile_exp bits)
    if (scratchVars.mprational_in_use)
       fprintf(stderr,"%s reusing scratchVars.mprational!\n",__func__);
    scratchVars.mprational_in_use = 1;
    if (scratchVars.mprational2_in_use)
       fprintf(stderr,"%s reusing scratchVars.mprational2!\n",__func__);
    scratchVars.mprational2_in_use = 1;


    mpq_t *resulting_scalar = &scratchVars.mprational;
    mpq_t *temp_scalar = &scratchVars.mprational2;
    mpq_mul_2exp(*resulting_scalar, scalar, tile_exp);
    mpq_set(*temp_scalar,*resulting_scalar);

    // want to truncate toward negative infinity 
    //  but the following calls will truncate toward zero
    mpz_set_q(tile_index, *resulting_scalar); // takes rational, sets integer
    mpq_set_z(*resulting_scalar, tile_index); // takes integer, sets rational

    // check to see if truncated value is not equal to previous value
    // and if it started out negative, then subtract 1 to move toward
    // negative infinity
    if ((!mpq_equal(*temp_scalar,*resulting_scalar)) && (mpq_cmp_ui(*temp_scalar, 0, 1) < 0)) {
      mpz_sub_ui(tile_index, tile_index, 1);  // subtract 1
    }

    scratchVars.mprational_in_use = 0;
    scratchVars.mprational2_in_use = 0;
}
#endif   // end IS_RATIONAL

#if defined(USE_MPREAL) || defined(USE_MPCOMPLEX)
static void get_tile_index_mpreal(mpz_t tile_index,
                           int tile_exp, const mpfr_t scalar)
{
    //printf("in %s\n",__func__);
    //mpfr_printf("input value sc is %.128Rf\n", sc);
    if (scratchVars.mpreal_in_use)
       fprintf(stderr,"%s reusing scratchVars.mpreal!\n",__func__);
    scratchVars.mpreal_in_use = 1;

    // Multiply by 2^tile_exp (by shifting the real left by tile_exp)
    // then round toward negative infinity
    mpfr_t *temp_scalar = &scratchVars.mpreal;
    mpfr_mul_2exp(*temp_scalar, scalar, tile_exp, MPFR_RNDN); // shift left by tile_exp
    mpfr_get_z(tile_index, *temp_scalar, MPFR_RNDD);  // rounds toward neg
    scratchVars.mpreal_in_use = 0;
}

 
#endif // #if defined(USE_MPREAL) || defined(USE_MPCOMPLEX)
 
#if defined(USE_REAL) || defined(USE_COMPLEX)
static void get_tile_index_real(mpz_t tile_index, int tile_exp, 
                         const long double scalar)
{
#ifdef DEBUG_MAR_C
  if (isnanl(input) ) {
    fprintf(stderr,"Inside %s a NaN showed up\n", __func__);
    exit (1);
  }
#endif // #ifdef DEBUG_MAR_C
  
  // This uses standard C functions (frexp, modf, and ldexp) to perform rounding
  // in way that is independent of the specific floating point representation
  // that is used by the underlying machine architecture.

   if (scratchVars.mpreal_in_use)
      fprintf(stderr,"%s reusing scratchVars.mpreal!\n",__func__);
   scratchVars.mpreal_in_use = 1;

   // multiply by 2^tile_exp
   // Shift the input left by the number of significant bits to keep
   //   using the ldexp() function.
   // truncate toward negative infinity
   // printf("Hash debug:  input = %25.16la   output = %25.16la  tile_exp = %d\n", input, rounded_answer, tile_exp);
   // set the mpz_t tile_index to be the  ?????
   mpfr_t *temp_scalar = &scratchVars.mpreal;
   mpfr_set_ld(*temp_scalar,scalar,MPFR_RNDN);
   mpfr_mul_2exp(*temp_scalar, *temp_scalar, tile_exp, MPFR_RNDN); // shift left by tile_exp
   mpfr_get_z(tile_index, *temp_scalar, MPFR_RNDD);  // rounds toward neg

   scratchVars.mpreal_in_use = 0;
}
#endif   // REAL or COMPLEX


#ifdef IS_BOUNDING
static void get_tile_index_bounding(mpz_t tile_index, int tile_exp, 
                             const larc_exponent_scalar_t scalar)
{
    int start = 16 + scalar->is_zero + 2*scalar->is_one + 4*scalar->is_nan + 8*scalar->is_exact;
    mpz_set_ui(tile_index, start);
    for (int i = 0; i < NUM_EXPONENTS; ++i) {
        mpz_mul_2exp(tile_index, tile_index, 8*sizeof(bexp_type));
        mpz_add_ui(tile_index, tile_index, scalar->explist[i]);
    }
    // NOTE: It turns out that having unique tile indices is not enough. The
    // MAR algorithm will snap adjacent tile indices.  Thus, we need to scale
    // the tile indices far away from each other.
    // NEW NOTE: We fixed this in the MAR code so we are commenting out this
    // workaround:
    // mpz_mul_2exp(tile_index, tile_index, 8);
}
#endif   // IS_BOUNDING

void get_tile_index(MAR_tile_index_t *tile_index_PTR, int tile_exp, const scalarType scalar)
{
#if defined(USE_INTEGER) || defined(USE_BOOLEAN)
  // return the scalar
   mpz_set_si(tile_index_PTR->index, scalar);
#elif defined(USE_REAL)
   get_tile_index_real(tile_index_PTR->index,tile_exp,scalar);
#elif defined(USE_COMPLEX)
   get_tile_index_real(tile_index_PTR->real_index,tile_exp,creall(scalar));
   get_tile_index_real(tile_index_PTR->imag_index,tile_exp,cimagl(scalar));
#elif defined(USE_MPINTEGER)
   // return the scalar
   mpz_set(tile_index_PTR->index,scalar);
#elif defined(USE_MPRATIONAL)
   get_tile_index_mprational(tile_index_PTR->index,tile_exp,scalar);
#elif defined(USE_MPRATCOMPLEX)
   get_tile_index_mprational(tile_index_PTR->real_index,tile_exp,scalar->real);
   get_tile_index_mprational(tile_index_PTR->imag_index,tile_exp,scalar->imag);
#elif defined(USE_MPREAL)
   get_tile_index_mpreal(tile_index_PTR->index,tile_exp,scalar);
#elif defined(USE_MPCOMPLEX)
   get_tile_index_mpreal(tile_index_PTR->real_index,tile_exp,mpc_realref(scalar));
   get_tile_index_mpreal(tile_index_PTR->imag_index,tile_exp,mpc_imagref(scalar));
#elif defined(USE_CLIFFORD)
   if (scratchVars.clifford_in_use)
      fprintf(stderr,"%s reusing scratchVars.clifford!\n",__func__);
   scratchVars.clifford_in_use = 1;

   // first find a good rational approximation to the clifford variable
   clifford_t *temp = &scratchVars.clifford;
   convert_clifford_to_mprational_approx(*temp,scalar);
   #ifdef IS_COMPLEX
   get_tile_index_mprational(tile_index_PTR->real_index, tile_exp,
                             (*temp)->real_coeffs[0]);
   get_tile_index_mprational(tile_index_PTR->imag_index, tile_exp,
                             (*temp)->imag_coeffs[0]);
   #else
   get_tile_index_mprational(tile_index_PTR->index, tile_exp,
                             (*temp)->real_coeffs[0]);
   #endif // end IS_COMPLEX
   scratchVars.clifford_in_use = 0;
#elif defined(IS_BOUNDING)
   get_tile_index_bounding(tile_index_PTR->index,tile_exp,scalar);
#else
   fprintf(stderr,"no valid scalarType for get_tile_index!\n");
   exit(-1);
#endif // #if defined(USE_INTEGER) || defined(USE_BOOLEAN)
   return;
}


////////////////////////////////////////////////////////////////////////////////////////
// ORIGINAL
////////////////////////////////////////////////////////////////////////////////////////

mats_ptr_t retrieve_PTR_scalar_record(scalarType target_scalar)
{
  int verbose = 0;
  if (verbose) printf("entered %s\n",__func__);

  // go find out the value of the regionbitparam determining tile size = 1/2^regionbitparam
  int regionbitparam = get_regionbitparam();
  if (verbose) printf("regionbitparam = %d\n",regionbitparam);

  // Calculate the MAR tile index for the tile containing the target scalar
  // A tile will have a multiprecision index (i,j) for complex, (i) for real.
  // GLOBAL_target_tile_index_PTR is a global variable
  if (verbose) {
      printf("Before calling retrieve_PTR_scalar_record the value for target index was\n");
      printf("GLOBAL_target_tile_index_PTR = %ld\n",(long int)GLOBAL_target_tile_index_PTR);
  }

  if (Ggti_in_use)
      fprintf(stderr,"in %s reuse of GLOBAL_target_tile_index_PTR!\n",__func__);
  Ggti_in_use = 1;

  get_tile_index(GLOBAL_target_tile_index_PTR, regionbitparam, target_scalar);
//  printf("overriding returned value to debug\n");
//  set_tile_index(GLOBAL_target_tile_index_PTR, -237, 42);
  if (verbose) 
  {
    printf("returned from get_tile_index\n");
    char *tempStr = sca_get_readable_approx_str(target_scalar);
    printf("target_scalar is %s\n", tempStr);
    free(tempStr);
    print_tile_index(GLOBAL_target_tile_index_PTR);
    printf("is the new value for the target index.\n");
  }

  // Determine hash chain associated with the MAR tile index for target scalar.
  int64_t hash_chain = hash_tile_index(GLOBAL_target_tile_index_PTR);
  if (verbose)   printf("returned from hash_tile_index\n");

  // We traverse the hash chain looking for either a primary or nhbr record
  // associated with this tile (ie, the tile_index of the record matches the
  // tile_index of target) for record related to a previously stored scalar
  // (with its MAR region). If we find a previously stored tile, it could be
  // a primary tile, in which case we return its matrix PTR, or it could be
  // a nhbr tile, in which case we return the matrix PTR of the primary tile
  // for that nhbr. Either way we have returned the primary PTR for the unique
  // representative of a MAR region which contains the target_scalar.
  mats_ptr_t primary_ptr = find_record_from_tile_index(hash_chain,
      GLOBAL_target_tile_index_PTR);
  if (primary_ptr != NULL) {
#if ALERT_ON_SNAP
      if (!sca_eq(target_scalar, primary_ptr->scalar_value)) {
          handle_snap(target_scalar, primary_ptr->scalar_value, "MAR");
      }
#endif // #if ALERT_ON_SNAP
      Ggti_in_use = 0;
      return primary_ptr;
  }

  // If there is no previously stored tile, then we create one and save the
  // target scalar as the representative scalar.  Then we create a MAR tile
  // for this new primary tile by claiming any as yet unclaimed "coveted"
  // neighbors. Then we return the matrix PTR for this new primary tile.

  // When we store a scalar in LARC we will create an associated MARregion
  // from several tiles: the primary tile which
  // contains the unique scalar representing this region, and between
  // zero and three neighboring tiles which are closest to the scalar and had
  // not been previously claimed.

  // Our next step is to create a new primary tile matrix record
  primary_ptr = insert_primary_scalar_record(GLOBAL_target_tile_index_PTR, hash_chain, target_scalar);

  // This new primary record has everything in it except the list of nhbr tile
  // nodePTRs that will be claimed below.

  // For integer types and for the bounding types, the MARregion is just the
  // primary tile, we don't try to claim nhbrs
#if defined(IS_INTEGER) || defined(IS_BOUNDING)
  // return the ptr to this matrix record,
  Ggti_in_use = 0;
  return(primary_ptr);
#endif   //  INTEGER and BOUNDING types

   // For nonInteger types, MAR  will attempt to claim the closest nhbr tiles
   // to the scalar to form a MARregion with the primary tile. In order to
   // determine which nhbr tiles are closest to the scalar in the primary, we
   // check to see where the scalar lies in a grid of half the normal width.
   // Then choose nhbr tiles that are closest to this portion of the primary.
    if (Gtti_in_use)
        fprintf(stderr,"in %s reuse of GLOBAL_tiny_tile_index_PTR!\n",__func__);
    Gtti_in_use = 1;

    get_tile_index(GLOBAL_tiny_tile_index_PTR, (regionbitparam+1), target_scalar);
    if (verbose)
    {
      printf("returned from get_tile_index on smaller grid\n");
      print_tile_index(GLOBAL_tiny_tile_index_PTR);
    }

    if (scratchVars.mpinteger_in_use)
       fprintf(stderr,"%s reusing scratchVars.mpinteger!\n",__func__);
    scratchVars.mpinteger_in_use = 1;

    if (Gnti_in_use)
        fprintf(stderr,"in %s reuse of GLOBAL_nhbr_tile_index_PTR!\n",__func__);
    Gnti_in_use = 1;


#ifdef IS_COMPLEX
    //  When we are doing real case, below to figure out which half
    // of the primary tile the representative scalar (target_scalar)
    //  lies in side, we calculate  which half we are in:
    //          half = tiny_tile_index - 2* tile_ index
    // with half=0 meaning negative side, and half= 1 positive side.
    // Now for complex we do real and imaginary parts separately.
    // Calculate which half of the tile the scalar lies in for real and imag
    // axes with half=0 meaning negative side, and half= 1 positive side
    int real_half, imag_half;
    unsigned int real_tile_offset, imag_tile_offset;

    mpz_t *temp_mpz = &scratchVars.mpinteger;
    mpz_mul_ui(*temp_mpz,GLOBAL_target_tile_index_PTR->real_index,2);
    mpz_sub(*temp_mpz,GLOBAL_tiny_tile_index_PTR->real_index,*temp_mpz);
    real_half = mpz_get_si(*temp_mpz);
    mpz_mul_ui(*temp_mpz,GLOBAL_target_tile_index_PTR->imag_index,2);
    mpz_sub(*temp_mpz,GLOBAL_tiny_tile_index_PTR->imag_index,*temp_mpz);
    imag_half = mpz_get_si(*temp_mpz);

    // Say the primary tile had indices (i,j)
    // Depending on the value of (real_half,imag_half) we will covet different sets of
    // neighbors (the ones closest to the scalar in the primary tile).
    // if (real_half,imag_half) is:
    //        (0,0) then we covet the tiles with indices: (i-1,j-1), (i-1,j), and (i,j-1)
    //        (0,1) then we covet the tiles with indices: (i-1,j+1), (i,j+1), and (i-1,j)
    //        (1,0) then we covet the tiles with indices: (i+1,j-1), (i+1,j), and (i,j-1)
    //        (1,1) then we covet the tiles with indices: (i+1,j+1), (i,j+1), and (i+1,j)

    // Finding the index of a coveted nhbr_tile

    // FIRST COVETED NHBR CAN BE WRITTEN
    //        (i,j-1+(2*imag_half))
    // real part is still i

    mpz_set(GLOBAL_nhbr_tile_index_PTR->real_index,GLOBAL_target_tile_index_PTR->real_index);
    real_tile_offset = 0;
    // imag part is j-1 or j+1
    if (imag_half==0) {
       mpz_sub_ui(GLOBAL_nhbr_tile_index_PTR->imag_index,GLOBAL_target_tile_index_PTR->imag_index,1);
       imag_tile_offset = 2;
    }
    else {
       mpz_add_ui(GLOBAL_nhbr_tile_index_PTR->imag_index,GLOBAL_target_tile_index_PTR->imag_index,1);
       imag_tile_offset = 1;
    }
    if (verbose)
    {
       printf("have calculated a coveted neighbor tile:\n");
       print_tile_index(GLOBAL_nhbr_tile_index_PTR);
    }

    // now go try claim it
    uint64_t nhbr_hash = hash_tile_index(GLOBAL_nhbr_tile_index_PTR);
    mats_ptr_t nhbr_ptr = find_record_from_tile_index(nhbr_hash,GLOBAL_nhbr_tile_index_PTR);
    if (verbose) printf("checking hash chain for existing tile with coveted neighbors index\n");
    if (nhbr_ptr==NULL) {
      if (verbose) {
        printf("There is no existing record with the index for the coveted nbhr\n");
        printf("so will try to claim it.\n");
      }
      // create a neighbor node and add it to the list
      hash_node_t *nhbr_node_ptr =
      insert_nhbr_record(nhbr_hash, GLOBAL_nhbr_tile_index_PTR, primary_ptr, 4*imag_tile_offset+real_tile_offset);
      // add this nhbr_node_ptr into the list in primary record ptr
      // when only the imaginary is changed the position is 10 == index 2 
      primary_ptr->tile_node[2] = nhbr_node_ptr;
    }
    else {
      if (verbose) {
        printf("There is a tile already with the index for the coveted nbhr\n");
        printf("so we can not claim it.\n");
      }
    }

    // SECOND COVETED NHBR CAN BE WRITTEN
    //       (i-1+(2*real_half),j)
    // imag part is still j
    mpz_set(GLOBAL_nhbr_tile_index_PTR->imag_index,GLOBAL_target_tile_index_PTR->imag_index);
    imag_tile_offset = 0;
    // real part is i-1 or i+1
    if (real_half==0) {
       mpz_sub_ui(GLOBAL_nhbr_tile_index_PTR->real_index,GLOBAL_target_tile_index_PTR->real_index,1);
       real_tile_offset = 2;
    }
    else {
       mpz_add_ui(GLOBAL_nhbr_tile_index_PTR->real_index,GLOBAL_target_tile_index_PTR->real_index,1);
       real_tile_offset = 1;
    }
    if (verbose)
    {
       printf("have calculated a coveted neighbor tile:\n");
       print_tile_index(GLOBAL_nhbr_tile_index_PTR);
    }

    // now go try claim it
    nhbr_hash = hash_tile_index(GLOBAL_nhbr_tile_index_PTR);
    nhbr_ptr = find_record_from_tile_index(nhbr_hash,GLOBAL_nhbr_tile_index_PTR);
    if (verbose) printf("checking hash chain for existing tile with coveted neighbors index\n");
    if (nhbr_ptr==NULL) {
      if (verbose) {
        printf("There is no existing record with the index for the coveted nbhr\n");
        printf("so will try to claim it.\n");
      }
      // create a neighbor node and add it to the list
      hash_node_t *nhbr_node_ptr =
	insert_nhbr_record(nhbr_hash, GLOBAL_nhbr_tile_index_PTR, primary_ptr, 4*imag_tile_offset+real_tile_offset);
      // add this nhbr_node_ptr into the list in primary record ptr
      // when only the real is changed the position is 01 == index 1 
      primary_ptr->tile_node[1] = nhbr_node_ptr;
    }
    else {
      if (verbose) {
        printf("There is a tile already with the index for the coveted nbhr\n");
        printf("so we can not claim it.\n");
      }
    }

    // THE THIRD COVETED NEIGHBOR CAN BE WRITTEN
    //             (i-1+(2*real_half),j-1+(2*imag_half))
    // real part is i-1 or i+1
    if (real_half==0) {
       mpz_sub_ui(GLOBAL_nhbr_tile_index_PTR->real_index,GLOBAL_target_tile_index_PTR->real_index,1);
       real_tile_offset = 2;
    }
    else {
       mpz_add_ui(GLOBAL_nhbr_tile_index_PTR->real_index,GLOBAL_target_tile_index_PTR->real_index,1);
       real_tile_offset = 1;
    }	
    // imag part is j-1 or j+1
    if (imag_half==0) {
       mpz_sub_ui(GLOBAL_nhbr_tile_index_PTR->imag_index,GLOBAL_target_tile_index_PTR->imag_index,1);
       imag_tile_offset = 2;
    }
    else {
       mpz_add_ui(GLOBAL_nhbr_tile_index_PTR->imag_index,GLOBAL_target_tile_index_PTR->imag_index,1);
       imag_tile_offset = 1;
    }	

    if (verbose)
    {
       printf("have calculated a coveted neighbor tile:\n");
       print_tile_index(GLOBAL_nhbr_tile_index_PTR);
    }

    // now go try claim it
    nhbr_hash = hash_tile_index(GLOBAL_nhbr_tile_index_PTR);
    nhbr_ptr = find_record_from_tile_index(nhbr_hash,GLOBAL_nhbr_tile_index_PTR);
    if (verbose) printf("checking hash chain for existing tile with coveted neighbors index\n");
    if (nhbr_ptr==NULL) {
      if (verbose) {
        printf("There is no existing record with the index for the coveted nbhr\n");
        printf("so will try to claim it.\n");
      }
      // create a neighbor node and add it to the list
      hash_node_t *nhbr_node_ptr =
	insert_nhbr_record(nhbr_hash, GLOBAL_nhbr_tile_index_PTR, primary_ptr, 4*imag_tile_offset+real_tile_offset);
      // add this nhbr_node_ptr into the list in primary record ptr
      // when both the real and imag are changed the position is 11 == index 3 
      primary_ptr->tile_node[3] = nhbr_node_ptr;
    }
    else {
      if (verbose) {
        printf("There is a tile already with the index for the coveted nbhr\n");
        printf("so we can not claim it.\n");
      }
    }
//
// REAL types    
//
#else  
    //  To figure out which half of the primary tile the representative scalar (target_scalar)
    //  lies in side, we calculate  which half we are in:
    //          half = tiny_tile_index - 2* tile_ index
    // with half=0 meaning negative side, and half= 1 positive side
    int half;
    unsigned int tile_offset;
    mpz_t *temp_mpz = &scratchVars.mpinteger;
    mpz_mul_ui(*temp_mpz,GLOBAL_target_tile_index_PTR->index,2);
    mpz_sub(*temp_mpz,GLOBAL_tiny_tile_index_PTR->index,*temp_mpz);
    half = mpz_get_si(*temp_mpz);
    if (verbose) printf("have calculated which half of tile scalar lies on half=%d\n",half);
    
    //  The index for the nhbr tile that we covet depends on the value of half
    //  if half is 0 then we covet the nhbr (i-1) which is in the negative direction,
    //  if half is 1 then we covet the nhbr (i+1) which is in the positive direction.
    // We will use the global variable GLOBAL_nhbr_tile_index_PTR to store the coveted value
    if (half==0) {
       mpz_sub_ui(GLOBAL_nhbr_tile_index_PTR->index,GLOBAL_target_tile_index_PTR->index,1);
       tile_offset = 2;
    }
    else {
       mpz_add_ui(GLOBAL_nhbr_tile_index_PTR->index,GLOBAL_target_tile_index_PTR->index,1);
       tile_offset = 1;
    }
    if (verbose)
    {
       printf("have calculated coveted neighbor tile:\n");
       print_tile_index(GLOBAL_nhbr_tile_index_PTR);
    }
    // Now try to claim the nhbr_tile
    uint64_t nhbr_hash = hash_tile_index(GLOBAL_nhbr_tile_index_PTR);
    if (verbose) printf("checking hash chain for existing tile with coveted neighbors index\n");
    mats_ptr_t nhbr_ptr = find_record_from_tile_index(nhbr_hash, GLOBAL_nhbr_tile_index_PTR);
    if (verbose) printf("returned value of nhbr_ptr: %ld\n",(long int)nhbr_ptr);
    if (nhbr_ptr==NULL) {
        if (verbose) printf("inserting neighbor record\n");
	// create a neighbor node and add it to the list
        hash_node_t *nhbr_node_ptr =
	      insert_nhbr_record(nhbr_hash, GLOBAL_nhbr_tile_index_PTR, primary_ptr, tile_offset);
        // add this nhbr_node_ptr into the list in primary record ptr
        // when only the real is changed the position is 01 == index 1 
        primary_ptr->tile_node[1] = nhbr_node_ptr;
    }
#endif // #ifdef IS_COMPLEX

    // inserted the nhbr tiles into the MARregion and added there node ptrs to primary
    // done with global variables
    Gnti_in_use = Gtti_in_use = Ggti_in_use = scratchVars.mpinteger_in_use = 0;

    // We have created the full MARregion with a primary tile and possibly some
    // coveted nhbrs and we return the pointer to the primary matrix record
    return primary_ptr;
    
}  // end retrieve_PTR_scalar_record



/////////////////////////////////////////////////////////////


static uint64_t retrieve_hash_chain_index_for_nbhr_node(int ndir, scalarType scalar_value)
{
     int verbose = 0;
     if (verbose) printf("entered %s\n",__func__);

     // find  the value of the regionbitparam determining tile size = 1/2^regionbitparam
     int regionbitparam = get_regionbitparam();
     if (verbose) printf("regionbitparam = %d\n",regionbitparam);

     // Calculate the MAR tile index for the tile containing the target scalar */
     // A tile will have a multiprecision index (i,j) for complex, (i) for real. */

     // GLOBAL_target_tile_index_PTR is a global variable */
     if (Ggti_in_use)
         fprintf(stderr,"in %s reuse of GLOBAL_target_tile_index_PTR!\n",
                 __func__);
     Ggti_in_use = 1;

     get_tile_index(GLOBAL_target_tile_index_PTR,
		             regionbitparam, scalar_value);
     
      if (verbose) 
     {
          printf("returned from get_tile_index\n");
          char *tempStr = sca_get_readable_approx_str(scalar_value);
          printf("scalar to be deleted from store is %s\n", tempStr);
          free(tempStr);
          print_tile_index(GLOBAL_target_tile_index_PTR);
          printf("is the tile index associated with the scalar.\n");
      }

  // When we store a scalar in LARC we create an associated MARregion
  // from several tiles: the primary tile which
  // contains the unique scalar representing this region, and between
  // zero and three neighboring tiles which are closest to the scalar and had
  // not been previously claimed.
  // The exceptions to this are the INTEGER and BOUNDING types
  // where no neighbors are claimed.
  //
  // We will eventual delete the primary record for the scalar,
  // but first we need to delete all the hash node records for each
  // neighbor tile that is mentioned in the primary record.
  // The primary record is the argument scalar_PTR passed into this routine

// MAKE SURE WE HANDLE BELOW
/*   // For integer types and for the bounding types, the MARregion is just the */
/*   // primary tile, we don't try to claim nhbrs */
/* #if defined(IS_INTEGER) || defined(IS_BOUNDING) */
/*   // return the ptr to this matrix record, */
/*   return(primary_ptr); */
/* #endif   //  INTEGER and BOUNDING types */

   // The way MAR chooses regions depends on which quarter of the primary
   // tile the scalar lies inside off, we use a grid of half normal width to make
   // this determination.
    if (Gtti_in_use)
        fprintf(stderr,"in %s reuse of GLOBAL_tiny_tile_index_PTR!\n",__func__);
    Gtti_in_use = 1;

    get_tile_index(GLOBAL_tiny_tile_index_PTR, (regionbitparam+1), scalar_value);
    if (verbose)
    {
      printf("returned from get_tile_index on smaller grid\n");
      print_tile_index(GLOBAL_tiny_tile_index_PTR);
    }

#ifdef IS_COMPLEX

    // Say the primary tile had indices (i,j)
    // Depending on the value of (real_half,imag_half) we will covet different sets of
    // neighbors (the ones closest to the scalar in the primary tile).
    // if (real_half,imag_half) is:
    //        (0,0) then we covet the tiles with indices: (i-1,j-1), (i-1,j), and (i,j-1)
    //        (0,1) then we covet the tiles with indices: (i-1,j+1), (i,j+1), and (i-1,j)
    //        (1,0) then we covet the tiles with indices: (i+1,j-1), (i+1,j), and (i,j-1)
    //        (1,1) then we covet the tiles with indices: (i+1,j+1), (i,j+1), and (i+1,j)

    // If the tile_node[2] contains an entry then we will calculate
    // the tile index and hash chain   
	// delete the hash node from the hash chain by
	// we look for the record (remember to check to see if this points to scalar_PTR)
	// sewing up the prev and next records OR if at the beginning or end
	// changing LAST or FIRST to the right value
	// and then freeing the hash_node and all its contents

    // For complex we do real and imaginary parts of the tile index separately.
    // Calculate which half of the tile the scalar lies in for real and imag
    // axes with half=0 meaning negative side, and half= 1 positive side

    int real_half, imag_half;
    if (scratchVars.mpinteger_in_use)
       fprintf(stderr,"%s reusing scratchVars.mpinteger!\n",__func__);
    scratchVars.mpinteger_in_use = 1;

    mpz_t *temp_mpz = &scratchVars.mpinteger;
    mpz_mul_ui(*temp_mpz,GLOBAL_target_tile_index_PTR->real_index,2);
    mpz_sub(*temp_mpz,GLOBAL_tiny_tile_index_PTR->real_index,*temp_mpz);
    real_half = mpz_get_si(*temp_mpz);
    mpz_mul_ui(*temp_mpz,GLOBAL_target_tile_index_PTR->imag_index,2);
    mpz_sub(*temp_mpz,GLOBAL_tiny_tile_index_PTR->imag_index,*temp_mpz);
    imag_half = mpz_get_si(*temp_mpz);
    scratchVars.mpinteger_in_use = 0;

    if (ndir ==2) {
    // FIRST COVETED NHBR CAN BE WRITTEN
    //        (i,j-1+(2*imag_half))
    // real part is still i
    // and if it is claimed it is stored in tile_node[2]
    if (Gnti_in_use)
        fprintf(stderr,"in %s reuse of GLOBAL_nhbr_tile_index_PTR!\n",__func__);
    Gnti_in_use = 1;

    //  We want  (i,j-1+(2*imag_half))
    mpz_set(GLOBAL_nhbr_tile_index_PTR->real_index,
	           GLOBAL_target_tile_index_PTR->real_index);
    // imag part is j-1 or j+1
    if (imag_half==0) {
       mpz_sub_ui(GLOBAL_nhbr_tile_index_PTR->imag_index,
		           GLOBAL_target_tile_index_PTR->imag_index,1);
    }
    else {
       mpz_add_ui(GLOBAL_nhbr_tile_index_PTR->imag_index,
		           GLOBAL_target_tile_index_PTR->imag_index,1);
    }
    if (verbose)
    {
       printf("have calculated tile index of this nhbr:\n");
       print_tile_index(GLOBAL_nhbr_tile_index_PTR);
    }

    // calculate the hash chain nhbr_hash and the hash node ptr nhbr_ptr
    uint64_t nhbr_hash_index = hash_tile_index(GLOBAL_nhbr_tile_index_PTR);
    Gnti_in_use = Ggti_in_use = Gtti_in_use = 0;
    return(nhbr_hash_index);
    } // end of ndir == 2

    if (ndir == 1) {
    // SECOND COVETED NHBR CAN BE WRITTEN
    //       (i-1+(2*real_half),j)
    // imag part is still j
    // and if it is claimed it is stored in tile_node[1]
    if (Gnti_in_use)
        fprintf(stderr,"in %s reuse of GLOBAL_nhbr_tile_index_PTR!\n",__func__);
    Gnti_in_use = 1;

    mpz_set(GLOBAL_nhbr_tile_index_PTR->imag_index,
	           GLOBAL_target_tile_index_PTR->imag_index);
    // real part is i-1 or i+1
    if (real_half==0) {
       mpz_sub_ui(GLOBAL_nhbr_tile_index_PTR->real_index,
		           GLOBAL_target_tile_index_PTR->real_index,1);
    }
    else {
       mpz_add_ui(GLOBAL_nhbr_tile_index_PTR->real_index,
		           GLOBAL_target_tile_index_PTR->real_index,1);
    }

    if (verbose)
    {
       printf("have calculated tile index of this nhbr:\n");
       print_tile_index(GLOBAL_nhbr_tile_index_PTR);
    }

    // calculate the hash chain nhbr_hash and the hash node ptr nhbr_ptr
    uint64_t nhbr_hash_index = hash_tile_index(GLOBAL_nhbr_tile_index_PTR);
    Gnti_in_use = Ggti_in_use = Gtti_in_use = 0;
    return(nhbr_hash_index);
  } // end of ndir ===1

    if (ndir == 3) {
    // THE THIRD COVETED NEIGHBOR CAN BE WRITTEN
    //             (i-1+(2*real_half),j-1+(2*imag_half))
    // real part is i-1 or i+1
    // and if it is claimed it is stored in tile_node[3]

    if (Gnti_in_use)
        fprintf(stderr,"in %s reuse of GLOBAL_nhbr_tile_index_PTR!\n",__func__);
    Gnti_in_use = 1;

    if (real_half==0) {
       mpz_sub_ui(GLOBAL_nhbr_tile_index_PTR->real_index,
		           GLOBAL_target_tile_index_PTR->real_index,1);
    }
    else {
       mpz_add_ui(GLOBAL_nhbr_tile_index_PTR->real_index,
		            GLOBAL_target_tile_index_PTR->real_index,1);
    }	
    // imag part is j-1 or j+1
    if (imag_half==0) {
       mpz_sub_ui(GLOBAL_nhbr_tile_index_PTR->imag_index,
		           GLOBAL_target_tile_index_PTR->imag_index,1);
    }
    else {
       mpz_add_ui(GLOBAL_nhbr_tile_index_PTR->imag_index,
		           GLOBAL_target_tile_index_PTR->imag_index,1);
    }	

    if (verbose)
    {
       printf("have calculated tile index of this nhbr:\n");
       print_tile_index(GLOBAL_nhbr_tile_index_PTR);
    }

    // calculate the hash chain nhbr_hash and the hash node ptr nhbr_ptr
    uint64_t nhbr_hash_index = hash_tile_index(GLOBAL_nhbr_tile_index_PTR);
    Gnti_in_use = Ggti_in_use = Gtti_in_use = 0;
    return(nhbr_hash_index);
  } // end of ndir = 3


// REAL types    
#else
#ifndef IS_INTEGER     // we think this is handling BOOLEAN as off Nov 2021
#ifndef IS_BOUNDING    
    //  To figure out which half of the primary tile the representative scalar (target_scalar)
    //  lies in side, we calculate  which half we are in:
    //          half = tiny_tile_index - 2* tile_ index
    // with half=0 meaning negative side, and half= 1 positive side

    if (ndir == 1) {
    int half;
    if (scratchVars.mpinteger_in_use)
       fprintf(stderr,"%s reusing scratchVars.mpinteger!\n",__func__);
    scratchVars.mpinteger_in_use = 1;

    mpz_t *temp_mpz = &scratchVars.mpinteger;
    mpz_mul_ui(*temp_mpz,GLOBAL_target_tile_index_PTR->index,2);
    mpz_sub(*temp_mpz,GLOBAL_tiny_tile_index_PTR->index,*temp_mpz);
    half = mpz_get_si(*temp_mpz);
    scratchVars.mpinteger_in_use = 0;

    if (verbose) printf("have calculated which half of tile scalar lies on half=%d\n",half);
    
    if (Gnti_in_use)
        fprintf(stderr,"in %s reuse of GLOBAL_nhbr_tile_index_PTR!\n",__func__);
    Gnti_in_use = 1;

    //  The index for the nhbr tile that we covet depends on the value of half
    //  if half is 0 then we covet the nhbr (i-1) which is in the negative direction,
    //  if half is 1 then we covet the nhbr (i+1) which is in the positive direction.
    // We will use the global variable GLOBAL_nhbr_tile_index_PTR to  calculate
    // the tile index and hash chain we need to delete the hash node for this neighbor
    if (half==0) {
       mpz_sub_ui(GLOBAL_nhbr_tile_index_PTR->index,GLOBAL_target_tile_index_PTR->index,1);
    }
    else {
       mpz_add_ui(GLOBAL_nhbr_tile_index_PTR->index,GLOBAL_target_tile_index_PTR->index,1);
    }
    if (verbose)
    {
       printf("have calculated neighbor tile index:\n");
       print_tile_index(GLOBAL_nhbr_tile_index_PTR);
    }
    // Now calculate the hash chain and hash node pointer for this record.
    uint64_t nhbr_hash_index = hash_tile_index(GLOBAL_nhbr_tile_index_PTR);
    Gnti_in_use = Ggti_in_use = Gtti_in_use = 0;
    return(nhbr_hash_index);
  } // ndir = 1
    
# endif  // #ifndef IS_BOUNDING
# endif  // #ifndef IS_INTEGER
#endif  // #ifdef IS_COMPLEX
   // inserted the nhbr tiles into the MARregion and added there node ptrs to primary

   // should only get here if ndir is 0 or not valid
   fprintf(stderr,"in %s, ndir of %d - nothing done\n",__func__,ndir);

   Ggti_in_use = Gtti_in_use = 0;

   return -1;

}  // retrieve_nhbr_hash_chain

static void free_smatrix(mats_ptr_t primary_record_ptr)
{
#ifdef STORE_TILE_INDEX
    free(primary_record_ptr->tile);
#endif // #ifdef STORE_TILE_INDEX

    sca_clear(&(primary_record_ptr->scalar_value));
    free(primary_record_ptr);

}  // end free_smatrix


int invalidate_recordPTR_in_indexTable(int64_t packedID);


/************************************************************************
 * Clean matrix store	                                                *
 *                                                                      *
 *  Removes all eligible matrices from the matrix store.                *
 *      The called function, remove_matrix_from_store, will         *
 *      recursively remove eligible children of any deleted matrix      *      	
     WARNING: LARC does not toggle counters to show that a scalar is
    is needed for norms or traces, so these scalars will be cleaned
     locked, held, or used by some 2 by 2 matrix 

 ***********************************************************************/
int clean_scalar_matrices_MAR(struct matrix_store_t * matrix_store_ptr)
{
    int verbose = 0;

    // the sized the scalarStore hash table
    uint64_t max = (uint64_t)1 << (matrix_store_ptr->scalar_hash_table->exponent);
    uint64_t hash_index;

    // hash table for scalar store
    hash_table_t *table_ptr = matrix_store_ptr->scalar_hash_table;

    // traversing the hash chain looking for primary records
    for (hash_index = 0; hash_index < max; ++hash_index) {
        hash_node_t *node_ptr = table_ptr->heads[hash_index];

        if (verbose) {
            printf("In %s: Checking hash_index = %lu\n", __func__,hash_index);
        }

        // traverse the hash chain, looking for primaries, and if they
        // pass the removability test, remove their nbhr node, primary node and record
        while (node_ptr) {

            if (verbose) {
                printf("In %s: Checking node pointer = %p\n", __func__,node_ptr);
            }

            //  if hash node is not a primary node continue
            if (node_ptr->tile_offset_flag != 0) {
                node_ptr = node_ptr->next;
                continue;
            }

            // save off the primary_node_ptr
            hash_node_t *primary_node_ptr = node_ptr;

            // get the packed ID from the node
            int64_t primary_pID = primary_node_ptr->packedID;

            // check to see if the packedID corresponds to an existing record
            if (primary_pID == MATRIX_ID_INVALID) {
                fprintf(stderr,"ERROR in %s:\n",__func__);
                fprintf(stderr,"found scalarStore node with invalid packedID ");
                fprintf(stderr,"while traversing hash chain %lu\n",hash_index);
                fprintf(stderr,"The scalar store is messed up, exiting.\n");
                exit(1);
            }
            // check to see if the packedID corresponds to a scalar matrix
            if (!IS_SCALAR(primary_pID)) {
                fprintf(stderr,"in %s: input packedID %ld is for a nonscalar!\n",
                            __func__,primary_pID);
                exit(0);
            }

            if (verbose) {
                printf("In %s: Checking primary ID = %lu\n", __func__,primary_pID);
            }

            // get the primary record PTR from the packed ID  in the node
            mats_ptr_t primary_record_ptr = (mats_ptr_t)
                            get_recordPTR_from_pID(primary_pID,"",__func__,0);
            if (primary_record_ptr == (mats_ptr_t) RECORD_PTR_INVALID) {
                fprintf(stderr,"ERROR in %s:\n",__func__);
                fprintf(stderr,"packedID %lu with null record ",primary_pID);
                fprintf(stderr,"found while traversing scalarStore ");
                fprintf(stderr,"hash chain %lu\n", hash_index);
                fprintf(stderr,"The scalar store is messed up, exiting.\n");
                exit(1);
            }

            // carry out removability test for the primary:  no dependencies, no hold, no lock
            if (   (primary_record_ptr->appears_as_sub_count == 0)
                && (primary_record_ptr->hold == 0)
                && (primary_record_ptr->lock == 0)
               )
            {
                // remove the nbr tiles and sew up their hash chain */
                for (int ndir = 1; ndir < 4; ++ndir)   // 01 is imag offset, 10 is real only, 11 is both
                {
                    hash_node_t *nbhr_hash_node_ptr = primary_record_ptr->tile_node[ndir];
                    if (nbhr_hash_node_ptr == NULL) {continue;}

                    // check to make sure we really have the right node we are about to remove */
                    if (nbhr_hash_node_ptr->tile_offset_flag == 0) {
                        fprintf(stderr,"ERROR in %s:\n",__func__);
                        fprintf(stderr,"tile_node inside primary record points to a hash node\n");
                        fprintf(stderr,"with primary offset flag instead of nbhr offset flag\n");
                        fprintf(stderr,"The scalar store is messed up, exiting.\n");
                        exit(1);
                    }
                    if (nbhr_hash_node_ptr->packedID != primary_pID) {
                        fprintf(stderr,"ERROR in %s:\n",__func__);
                        fprintf(stderr,"tile_node inside primary record points to a hash node\n");
                        fprintf(stderr,"whose packedID is not equal to the primaryID.\n");
                        fprintf(stderr,"The scalar store is messed up, exiting.\n");
                        exit(1);
                    }

                    // Sew up the hash chain of the nbhr node
                    // SIMPLE METHOD - when nbhr node not at head or tail of the hash chain
                    if (nbhr_hash_node_ptr->prev != NULL)
                            nbhr_hash_node_ptr->prev->next = nbhr_hash_node_ptr->next;
                    if (nbhr_hash_node_ptr->next != NULL)
                            nbhr_hash_node_ptr->next->prev = nbhr_hash_node_ptr->prev;

                    // FANCY SEWING - when nbhr node is at the head or tail of the hash chain
                    // if either of the prev or next were NULL we need to find the hash_chain */
                    // for the nbhr so we can sew it up using HEAD and TAIL */
                    if (   (nbhr_hash_node_ptr->prev == NULL)    // nbhr node at head of chain
                        || (nbhr_hash_node_ptr->next == NULL)    // nbhr node at tail of chain
                       )
                    {
                        // find the nbhr hash chain index so we can sew up the ends of the chain
                        uint64_t nbr_hash_chain_index =
                                retrieve_hash_chain_index_for_nbhr_node(ndir,
                                        primary_record_ptr->scalar_value);
                        // SEW up the ends of the hash chain
                        if (nbhr_hash_node_ptr->prev == NULL)
                                table_ptr->heads[nbr_hash_chain_index]
                                        = nbhr_hash_node_ptr->next;
                        if (nbhr_hash_node_ptr->next == NULL)
                                table_ptr->tails[nbr_hash_chain_index]
                                        = nbhr_hash_node_ptr->prev;
                    }  // end FANCY SEWING

                    //  The sewing up of the nbhr hash chain is complete
                    //  and now we free  this hash_node
                    free(nbhr_hash_node_ptr);

                }  // loop through nbhr_nodes, ndir = 1, 2, 3

                // Remove the primary hash node

                // First sew up the hash chain for the primary hash node
                // SIMPLE METHOD - when primary node not at head or tail of the hash chain
                if (primary_node_ptr->prev != NULL)
                        primary_node_ptr->prev->next = primary_node_ptr->next;
                if (primary_node_ptr->next != NULL)
                        primary_node_ptr->next->prev = primary_node_ptr->prev;

                // FANCY SEWING - when primary node is at the head or tail of the hash chain
                // if either of the prev or next were NULL we need to repair the values */
                // in the head and/or tail of the hash chain */
                if (   (primary_node_ptr->prev == NULL)    // nbhr node at head of chain
                    || (primary_node_ptr->next == NULL)    // nbhr node at tail of chain
                   )
                {
                    // we know the hash chain for the primary it is hash_index
                    // SEW up the ends of the hash chain
                    if (primary_node_ptr->prev == NULL)
                            table_ptr->heads[hash_index]
                                    = primary_node_ptr->next;
                    if (primary_node_ptr->next == NULL)
                            table_ptr->tails[hash_index]
                                    = primary_node_ptr->prev;
                }  // end FANCY SEWING

                //  The sewing up of the primary hash chain is complete
                //  and now we free  this hash_node
                free(primary_node_ptr);

                // Remove the matrix pointer from the table indexed by matrixIDs */
                if (verbose) {
                     fprintf(stderr,"calling invalidate_recordPTR_in_indexTable(%ld)\n",primary_pID);
                }
                invalidate_recordPTR_in_indexTable(primary_pID);

                // Clean up required for ScalarRecord
                // update the statistics, and recordPTR table.
                // Decrement the histogram and counts of matrices and scalars
                matrix_store_ptr->hist[0][0]--;
                matrix_store_ptr->num_scalars--;

                // Now we can remove the primary record
                free_smatrix(primary_record_ptr);

            } // end PRIMARY RECORD IS REMOVABLE

            node_ptr = node_ptr->next;
        } // end LOOP THROUGH THIS HASH CHAIN
    } // end LOOP OVER HASH INDICES

    return 1;
}


#endif // MAR defined

