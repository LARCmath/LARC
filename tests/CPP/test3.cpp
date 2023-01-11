#include <larc.hh>

#include <organize.h>
#include <matrix_store.h>

#include "support.hh"

void initialize ( ) 
{
  // Default hash exponent size
  size_t hash_exponent_matrix = DEFAULT_MATRIX_STORE_EXPONENT;
  size_t hash_exponent_op     = DEFAULT_OP_STORE_EXPONENT;
  int    zerobitthresh        = ZEROBITTHRESH_DEFAULT;
  int    sighash              = SIGHASH_DEFAULT;
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
  initialize();

  uint64_t start_count = num_matrices_in_store ( );

  int test_count = 0;

  for ( mat_level_t row_level = 0 ; row_level < 5 ; row_level++ ) {
    for ( mat_level_t col_level = 0 ; col_level < 5 ; col_level++ ) {

      LARC_t A ( row_level , col_level );
      LARC_t B ( col_level , row_level );
      LARC_t C ( row_level , row_level );
      
      for ( mat_level_t r = 0 ; r < 1 << row_level ; r++ ) 
	for ( mat_level_t c = 0 ; c < 1 << col_level ; c++ ) {
	  A[r][c] = r * 1001 + c + 1;
	  B[c][r] = c * 937  + r + 1;
	}
      
      C = A * B;

      // Ensure that LARC_t sets holds on matrices
      if ( 0 == test_count % 5 )
	clean_matrix_storage ( );

      for ( mat_level_t r = 0 ; r < 1 << row_level ; r++ ) 
	for ( mat_level_t c = 0 ; c < 1 << row_level ; c++ ) {

	  ScalarType sum = 0;
	  for ( mat_level_t k = 0 ; k < 1 << col_level ; k++ ) 
	    sum += A[r][k] * B[k][c];

	  assert_equal ( sum , C[r][c] , "verifying matrix multiplication" );
	}

      test_count++;
    }
  }

  // Ensure that LARC_t releases holds on matrices
  clean_matrix_store ( );
  uint64_t end_count = num_matrices_in_store ( );  
  assert_equal ( end_count , start_count , "matrix store clean up." );
}
