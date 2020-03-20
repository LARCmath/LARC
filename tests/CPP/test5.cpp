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

  for ( mat_level_t row_level_A = 0 ; row_level_A < 3 ; row_level_A++ ) {
    for ( mat_level_t col_level_A = 0 ; col_level_A < 3 ; col_level_A++ ) {

      LARC_t A ( row_level_A , col_level_A );

      for ( int r = 0 ; r < 1 << row_level_A; r++ ) 
	for ( int c = 0 ; c < 1 << col_level_A; c++ ) 
	  A[r][c] = 123 * r + c + 1;
      
      for ( mat_level_t row_level_B = 0 ; row_level_B < 3 ; row_level_B++ ) {
	for ( mat_level_t col_level_B = 0 ; col_level_B < 3 ; col_level_B++ ) {
	  
	  LARC_t B ( row_level_B , col_level_B );
	  for ( int r = 0 ; r < 1 << row_level_B; r++ ) 
	    for ( int c = 0 ; c < 1 << col_level_B; c++ ) 
	      B[r][c] = 321 * r + c + 7;

	  
	  LARC_t C = A.kron ( B );

	  for ( int r = 0 ; r < 1 << (row_level_A + row_level_B) ; r++ ) 
	    for ( int c = 0 ; c < 1 << (col_level_A + col_level_B) ; c++ ) 
	      assert_equal ( A[r/(1<<row_level_B)][c/(1<<col_level_B)] * 
			     B[r%(1<<row_level_B)][c%(1<<col_level_B)] ,
			     C[r][c] , "verifying the kronecker product" );

	}
      }
    }
  }
  
  // Ensure that LARC_t releases holds on matrices
  clean_matrix_store ( );
  uint64_t end_count = num_matrices_in_store ( );  
  assert_equal ( end_count , start_count , "verifying matrix store clean up" );
}

