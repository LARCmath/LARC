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
  int    zerobitthresh        = ZEROBITTHRESH_DEFAULT;
  int    sighash              = SIGHASH_DEFAULT;
  mat_level_t max_level       = DEFAULT_MAX_LEVEL;

  initialize_larc ( hash_exponent_matrix , 
		    hash_exponent_op     ,
		    max_level            ,
		    sighash              ,
		    zerobitthresh        );
}

int main ( int argc , char *argv[] ) 
{
  // Initialize the larc library.
  initialize ( );
  
  // Work with matrices of different sizes
  for ( mat_level_t row_level = 0 ; row_level < 4 ; row_level++ ) {
    for ( mat_level_t col_level = 0 ; col_level < 4 ; col_level++ ) {

      mat_add_t m = get_zero_matrix_ptr ( row_level, col_level ); 

      ScalarType one     = (ScalarType) 1;
      ScalarType counter = one;

      for ( mat_level_t r = 0 ; r < 1 << row_level ; r++ ) {
	for ( mat_level_t c = 0 ; c < 1 << col_level ; c++ ) {

	  ScalarType v = matrix_get_value ( m , r , c );
	  assert ( 0 == v );

	  m = matrix_set_value ( m , r , c , counter );

	  counter = counter + one;
	}
      }

      counter = one;

      for ( mat_level_t r = 0 ; r < 1 << row_level ; r++ ) {
	for ( mat_level_t c = 0 ; c < 1 << col_level ; c++ ) {

	  ScalarType v = matrix_get_value ( m , r , c );

	  assert ( counter == v );
	  
	  counter = counter + one;
	}
      }
    }
  }

  return 0;
}
