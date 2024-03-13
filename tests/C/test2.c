//                     test2.c
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
#include <stdlib.h>
#include <assert.h>

#include <global.h>
#include <larc.h>
#include <hash.h>
#include <io.h>
#include <json.h>
#include <matmath.h>
#include <matrix_store.h>
#include <op_store.h>
#include <organize.h>

void initialize ( ) 
{
  // Default hash exponent size
  size_t hash_exponent_matrix = DEFAULT_MATRIX_STORE_EXPONENT;
  size_t hash_exponent_op     = DEFAULT_OP_STORE_EXPONENT;
  int    sighash              = DEFAULT_REGIONBITPARAM;
  int    zerobitthresh        = -1;
  mat_level_t max_level       = DEFAULT_MAX_LEVEL;
  int    verbose	      = 1;

  initialize_larc ( hash_exponent_matrix , 
		    hash_exponent_op     ,
		    max_level            ,
		    sighash              ,
		    zerobitthresh        ,
		    verbose              );
}

int main ( int argc , char *argv[] ) 
{
  // Initialize the larc library.
  initialize ( );
  
  scalarType one, counter, v;
  sca_init(&one);
  sca_init(&counter);
  sca_init(&v);

  // Work with matrices of different sizes
  for ( mat_level_t row_level = 0 ; row_level < 4 ; row_level++ ) {
    for ( mat_level_t col_level = 0 ; col_level < 4 ; col_level++ ) {

      mat_ptr_t m = get_zero_matrix_ptr ( row_level, col_level ); 
//      mat_ptr_t zPTR= get_zero_matrix_ptr (0,0);
//      int64_t zID = get_matID_from_matPTR(zPTR);

      sca_set_str(&one,"1");
      sca_set(&counter,one);

      for ( int r = 0 ; r < 1 << row_level ; r++ ) {
	for ( int c = 0 ; c < 1 << col_level ; c++ ) {

          get_valPTR_from_matPTR_and_coords(&v, m, r, c);
	  assert ( sca_eq((scalarType)0, v ) );

	  m = get_matPTR_from_oldMatPTR_newVal_and_coords (m, r, c, counter);

          sca_add(&counter,counter,one);
	}
      }

      sca_set(&counter,one);

      for ( int r = 0 ; r < 1 << row_level ; r++ ) {
	for ( int c = 0 ; c < 1 << col_level ; c++ ) {

          get_valPTR_from_matPTR_and_coords(&v, m, r, c);
	  assert ( sca_eq(counter,v) );
	  
          sca_add(&counter,counter,one);
	}
      }
    }
  }
  sca_clear(&one);
  sca_clear(&counter);
  sca_clear(&v);

  return 0;
}
