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

  for ( mat_level_t level = 0 ; level < 5 ; level++ ) {
    
    LARC_t A ( level , level );
    LARC_t B ( level , level );
    LARC_t C ( level , level );
    
    for ( mat_level_t r = 0 ; r < 1 << level ; r++ ) 
      for ( mat_level_t c = 0 ; c < 1 << level ; c++ ) {
	A[r][c] = r * 1001 + c + 1;
	B[r][c] = c * 937  + r + 1;
      }
    
    C = (A - 2*B)*(B - A*2);
    
    // Ensure that LARC_t sets holds on matrices
    clean_matrix_storage ( );
    
    for ( mat_level_t r = 0 ; r < 1 << level ; r++ ) 
      for ( mat_level_t c = 0 ; c < 1 << level ; c++ ) {
	
	ScalarType sum = 0;
	for ( mat_level_t k = 0 ; k < 1 << level ; k++ ) {
	  ScalarType two = 2;
	  sum = sum + (A[r][k] - two*B[r][k]) * (B[k][c] - two*A[k][c]);
	}
	assert_equal ( sum , C[r][c] , "matrix algebra" );
      }
  }
  
  // Ensure that LARC_t releases holds on matrices
  clean_matrix_storage ( );
  uint64_t end_count = num_matrices_in_store ( );  
  assert_equal ( end_count , start_count , "matrix store clean up" );
}
