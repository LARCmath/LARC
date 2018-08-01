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

  initialize_larc ( hash_exponent_matrix , 
		    hash_exponent_op     ,
		    max_level            ,
		    sighash              ,
		    zerobitthresh        );
}

int main ( int argc , char *argv[] )
{
  initialize();

  uint64_t start_count = num_matrices_in_store ( );

  for ( mat_level_t level = 0 ; level < 5 ; level++ ) {

    LARC_t A ( level , level );

    for ( int r = 0 ; r < 1 << level; r++ ) 
      for ( int c = 0 ; c < 1 << level; c++ ) 
	A[r][c] = 123 * r + c + 1;
    
    ScalarType trace = A.trace();
    
    ScalarType sum = (ScalarType)0;
    for ( int d = 0 ; d < 1 << level; d++ ) 
      sum = sum + A[d][d];
    
    assert_equal ( trace , sum , "verifying traces" );
  }
  
  // Ensure that LARC_t releases holds on matrices
  clean_matrix_store ( );
  uint64_t end_count = num_matrices_in_store ( );  
  assert_equal ( end_count , start_count , "verifying matrix store clean up" );
}

