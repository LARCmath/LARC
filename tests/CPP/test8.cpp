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

  for ( mat_level_t row_level = 0 ; row_level < 3 ; row_level++ ) {
    for ( mat_level_t col_level = 0 ; col_level < 3 ; col_level++ ) {

      LARC_t A ( row_level , col_level );
      LARC_t B ( row_level , col_level );

      for ( int r = 0 ; r < 1 << row_level; r++ ) 
	for ( int c = 0 ; c < 1 << col_level; c++ ) {
	  A[r][c] = 123 * r + c + 1;
	  B[r][c] = r + 321 * c + 1;
	}
      
      LARC_t J = A.join  ( B );
      LARC_t S = A.stack ( B );

      for ( int r = 0 ; r < 1 << row_level; r++ ) 
	for ( int c = 0 ; c < 1 << col_level; c++ ) {
	  assert_equal ( J[r][c] , A[r][c] , "verifying the join operation(1)" );
	  assert_equal ( J[r][c + (1<<col_level)] , B[r][c] , "verifying the join opeation(2)" );
	  assert_equal ( S[r][c] , A[r][c] , "verifying the stack operation(1)" );
	  assert_equal ( S[r + (1<<row_level)][c] , B[r][c] , "verifying the stack operation(2)" );
	}      
    }
  }
  
  // Ensure that LARC_t releases holds on matrices
  clean_matrix_storage ( );
  uint64_t end_count = num_matrices_in_store ( );  
  assert_equal ( end_count , start_count , "verifying matrix store clean up" );
}

