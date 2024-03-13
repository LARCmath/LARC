//                    exampleLARC.c
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
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <pthread.h>
#include <math.h>

// Our header files structures and functions
#include "global.h"
#include "larc.h"
#include "hash.h"
#include "io.h"
#include "json.h"
#include "matmath.h"
#include "matrix_store.h"
#include "op_store.h"
#include "organize.h"
#include "scalars.h"
#include "show_scalars.h"
#include "show_layout.h"
#include "exampleLARC.h"

/*!
 * \file exampleLARC.c
 * \brief Test program that runs through several LARC capabilities
 */

/*!
 * \ingroup larc
 * \brief A function that returns the norm squared of a scalar argument
 * \param y A pointer to the return value
 * \param x The input scalarType value
 */
static void square_scalar(scalarType *y, const scalarType x)
{
	if (scratchVars.calc_conj_in_use)
	    fprintf(stderr,"%s reusing scratchVars.calc_conj!\n",__func__);
	scratchVars.calc_conj_in_use = 1;
	
        sca_conj(&scratchVars.calc_conj,x);
        sca_mult(y,scratchVars.calc_conj,x);
	scratchVars.calc_conj_in_use = 0;
}

int64_t square_matrix_elements(int64_t input_pID, op_type_t op_memz)
{
        return apply_function_to_matrix_values(input_pID,
                square_scalar, op_memz);
}


#ifdef IS_BOUNDING
void larc_sca_init_bounding(larc_exponent_scalar_t *sc);
void larc_sca_clear_bounding(larc_exponent_scalar_t *sc);
char *larc_sca_get_str_bounding(const larc_exponent_scalar_t sc);
void larc_sca_add_bounding(larc_exponent_scalar_t *rsc, const larc_exponent_scalar_t sc1, const larc_exponent_scalar_t sc2);
void larc_sca_mult_bounding(larc_exponent_scalar_t *rsc, const larc_exponent_scalar_t sc1, const larc_exponent_scalar_t sc2);
void larc_sca_set_2ldoubles_bounding(larc_exponent_scalar_t *rsc, long double rv, long double iv);
void larc_sca_set_str_bounding(larc_exponent_scalar_t *rsc, const char *input_str);

static void set_bounding_test_var(larc_exponent_scalar_t *x,
    char **disp, int n)
{
    larc_sca_clear_bounding(x);
    (*x)->is_zero = 0;
    if (n == 0) {
        (*x)->is_zero = 1;
        (*disp) = "0";
    } else if (n == 1) {
        (*x)->is_one = 1;
        (*disp) = "1";
    } else if (n == 2) {
#ifdef USE_UPPER
        (*x)->is_one = 1;
        (*x)->is_exact = 0;
        (*disp) = "^1";
#else
        (*x)->is_nan = 1;
        (*disp) = "NaN";
#endif // #ifdef USE_UPPER
    } else if (n == 3) {
        (*x)->explist[0] = 2;
        (*x)->explist[1] = 4;
        (*x)->explist[2] = 5;
        (*disp) = "[2,4,5]";
    } else if (n == 4) {
        (*x)->is_exact = 0;
        (*x)->explist[0] = 2;
        (*x)->explist[1] = 4;
        (*x)->explist[2] = 5;
#ifdef USE_UPPER
        (*disp) = "^[2,4,5]";
#else
        (*disp) = "v[2,4,5]";
#endif // #ifdef USE_UPPER
    } else if (n == 5) {
        (*x)->explist[0] = 2;
        (*x)->explist[1] = 3;
        (*x)->explist[2] = 4;
        (*disp) = "[2,3,4]";
    } else if (n == 6) {
        (*x)->is_exact = 0;
        (*x)->explist[0] = 2;
        (*x)->explist[1] = 3;
        (*x)->explist[2] = 4;
#ifdef USE_UPPER
        (*disp) = "^[2,3,4]";
#else
        (*disp) = "v[2,3,4]";
#endif // #ifdef USE_UPPER
    } else {
        (*x)->is_zero = 1;
        (*disp) = "0";
    }
}

static void test_bounding_set_2ld()
{
    printf("\n");
    larc_exponent_scalar_t x1;
    larc_sca_init_bounding(&x1);
    long double a[] = {-0.1L, 0.0L, 0.1L, 0.2L, 0.25L, 0.3L, 0.4L, 0.5L, 0.6L, 0.75L, 0.8L, 0.9L, 1.0L, 1.1L,
                       0.1875L, 0.3125L, 0.375L, 0.46875L, 0.625L, 0.65625L, 0.6875L, -1.0L};
    for (int i = 0; a[i] != -1.0L; ++i) {
        larc_sca_set_2ldoubles_bounding(&x1, a[i], 0.0L);
        printf("LD: %Lf -> %s\n", a[i], larc_sca_get_str_bounding(x1));
    }
}

static void test_bounding_set_str()
{
#ifdef USE_LOWER
    char approx = 'v';
#else
    char approx = '^';
#endif // #ifdef USE_LOWER
    printf("\n");
    larc_exponent_scalar_t x1;
    larc_sca_init_bounding(&x1);
    char * a[] = {"0.0", "0", "0.5", "0.25", "0.4", "5/16", "1/2", "1/3", "1/65535", "1/4294967295",
                  "explist[1,2,3]", "explist[1]", "explist[5,10]", "explist[5,10,15,20]", "explist[2,4,6,8,10]", ""};
    for (int i = 0; a[i][0] != '\0'; ++i) {
        larc_sca_set_str_bounding(&x1, a[i]);
        printf("SS: %s -> %s\n", a[i], larc_sca_get_str_bounding(x1));
    }
    char buf[1024];
    for (int i = 0; a[i][0] != '\0'; ++i) {
        buf[0] = approx;
        buf[1] = '\0';
        strcat(buf, a[i]);
        larc_sca_set_str_bounding(&x1, buf);
        printf("SS: %s -> %s\n", buf, larc_sca_get_str_bounding(x1));
    }
}

static void test_bounding_code_all()
{
    printf("\nTesting bounding code ALL.\n");
    larc_exponent_scalar_t x1, x2, x3;
    char *dx1, *dx2;
    larc_sca_init_bounding(&x1);
    larc_sca_init_bounding(&x2);
    larc_sca_init_bounding(&x3);
    for (int n1 = 0; n1 <= 6; ++n1) {
        for (int n2 = 0; n2 <= 6; ++n2) {
            set_bounding_test_var(&x1, &dx1, n1);
            set_bounding_test_var(&x2, &dx2, n2);
            printf("\n==============================\n");
            printf("Testing %s op %s.\n", dx1, dx2);
            printf("x1 = %s\n", larc_sca_get_str_bounding(x1));
            printf("x2 = %s\n", larc_sca_get_str_bounding(x2));
            larc_sca_add_bounding(&x3, x1, x2);
            printf("x1 + x2 = %s\n", larc_sca_get_str_bounding(x3));
            larc_sca_mult_bounding(&x3, x1, x2);
            printf("x1 * x2 = %s\n", larc_sca_get_str_bounding(x3));
        }
    }
    test_bounding_set_2ld();
    test_bounding_set_str();
}

static void test_bounding_code()
{
    test_bounding_code_all();
    exit(1);
}
#endif

int main (int argc, char *argv[])
{

#ifdef IS_BOUNDING
  test_bounding_code();
#endif

  initialize_larc(25, 26, 10, -1, -1, 1);

// For debugging purposes, show layout of bit structures.
//#define SHOW_STRUCT_LAYOUT
#ifdef SHOW_STRUCT_LAYOUT
  show_layout_hash_node();
  show_layout_larc_matrix();
  printf("\n");
#endif // #ifdef SHOW_STRUCT_LAYOUT

  // read in a matrix and get its packedID. The matrix is size 8x8
  char* pathToMatrix = "../tests/dat/in/matrixMarketExchange1";
  int64_t A_pID = read_matrixMarketExchange_file(pathToMatrix);

  printf("\nWe compare two different methods for performing a nonlinear\n");
  printf("matrix operation: taking each element of a matrix and squaring\n");
  printf("it. The result of the two operations should be identical and\n");
  printf("therefore have the same packedID once stored in the matrixStore.\n");

  // We call the LARC function matrix_entrySquared with scale factor
  // set to 1. This should have the same result as the local functions above.
  int64_t A2_pID = matrix_entrySquared(A_pID,"1");

  // We call the function square_scalar_matID, which takes two arguments, the
  // input matrix and a type of operation, the latter of which determines how
  // the calculation is memoized. We must use one of the operation types set
  // aside for user-defined functions, because using any operation type for
  // more than one operation can result in wrong answers.
  int64_t A2x_pID = square_matrix_elements(A_pID,FUNC_A); 

  // We compare the two packedIDs, which should be identical
  if (A2_pID == A2x_pID) { printf("SUCCESS!\n"); }
  else
  {
     printf("ERROR in %s: the two ways of calculating the same ",__func__);
     printf("matrix return different results!\n");
  }
}
