#include <larc.hh>

#include <organize.h>
#include <matrix_store.h>

#include "support.hh"

using namespace std;

void initialize ( ) 
{
  // Default hash exponent size
  size_t hash_exponent_matrix = DEFAULT_MATRIX_STORE_EXPONENT;
  size_t hash_exponent_op     = DEFAULT_OP_STORE_EXPONENT;
  int    zerobitthresh        = ZEROBITTHRESH_DEFAULT;
  int    sighash              = SIGHASH_DEFAULT;
  mat_level_t max_level       = DEFAULT_MAX_LEVEL;
  int	 verbose	      = 1;

  initialize_larc ( hash_exponent_matrix , 
		    hash_exponent_op     ,
		    max_level            ,
		    sighash              ,
		    zerobitthresh        ,
		    verbose	 	 );
}

void copy_slice ( Const_Slice_t src_slice , Slice_t dst_slice )
{
  for ( mat_level_t col = 0 ; col < 1 << src_slice.get_level ( ) ; col++ )
    dst_slice[col] = src_slice[col];
}

void test_copy_by_slices ( const LARC_t &src , LARC_t &dst )
{
  for ( mat_level_t row = 0 ; row < 1 << src.get_row_level ( ) ; row++ ) {
    copy_slice ( src[row] , dst[row] );
  }
}

void test_copy ( const LARC_t &src , LARC_t &dst )
{
  for ( mat_level_t row = 0 ; row < 1 << src.get_row_level ( ) ; row++ ) 
    for ( mat_level_t col = 0 ; col < 1 << src.get_col_level ( ) ; col++ ) 
      dst[row][col] = src[row][col];
}

int main ( int argc , char *argv[] )
{
  // Initialize the larc library.
  initialize ( );

  uint64_t start_count = num_matrices_in_store ( );

  for ( mat_level_t row_level = 0 ; row_level < 5 ; row_level++ ) {
    for ( mat_level_t col_level = 0 ; col_level < 5 ; col_level++ ) {
      
      LARC_t M ( row_level , col_level );
      
      // Basic assignment and read
      for ( mat_level_t row = 0 ; row < 1 << row_level ; row++ ) 
	for ( mat_level_t col = 0 ; col < 1 << col_level ; col++ ) 
	  M[row][col] = (ScalarType)(row * 1001 + col);

      for ( mat_level_t row = 0 ; row < 1 << row_level ; row++ ) 
	for ( mat_level_t col = 0 ; col < 1 << col_level ; col++ ) 
	  assert_equal ( (ScalarType)(row * 1001 + col) , M[row][col] , 
			 "verifying element wise read/write" );

      // Assignment and read by slice
      for ( mat_level_t row = 0 ; row < 1 << row_level ; row++ ) {
	auto S = M[row];
	for ( mat_level_t col = 0 ; col < 1 << col_level ; col++ ) 
	  S[col] = (ScalarType)(row * 1001 + col);
      }      

      for ( mat_level_t row = 0 ; row < 1 << row_level ; row++ ) {
	auto S = M[row];
	for ( mat_level_t col = 0 ; col < 1 << col_level ; col++ ) 
	  assert_equal ( (ScalarType)(row * 1001 + col) , S[col] ,
			 "verifying element wise read/write for slices" );
      }

      // chained element-to-element assignment and read
      LARC_t N ( row_level , col_level );
      LARC_t P ( row_level , col_level );

      for ( mat_level_t row = 0 ; row < 1 << row_level ; row++ ) 
	for ( mat_level_t col = 0 ; col < 1 << col_level ; col++ ) 
	  P[row][col] = N[row][col] = M[row][col];

      for ( mat_level_t row = 0 ; row < 1 << row_level ; row++ ) 
	for ( mat_level_t col = 0 ; col < 1 << col_level ; col++ ) {
	  assert_equal ( N[row][col] , M[row][col] , "verifying chained assignment(1)" );
	  assert_equal ( P[row][col] , M[row][col] , "verifying chained assignment(2)" );
	}

      // Passing LARC_t by reference and updating it.
      LARC_t Q ( row_level , col_level );
      test_copy ( M , Q );

      for ( mat_level_t row = 0 ; row < 1 << row_level ; row++ ) 
	for ( mat_level_t col = 0 ; col < 1 << col_level ; col++ ) 
	  assert_equal ( M[row][col] , Q[row][col] , "verifying copy from reference" );

      // Passing LARC_t by reference and updating it slice by slice.
      LARC_t R ( row_level , col_level );
      test_copy_by_slices ( M , R );

      for ( mat_level_t row = 0 ; row < 1 << row_level ; row++ ) 
      	for ( mat_level_t col = 0 ; col < 1 << col_level ; col++ ) 
      	  assert_equal ( M[row][col] , R[row][col] , "verifying copy from reference by slice" );
      
      // Updating in place 
      for ( mat_level_t row = 0 ; row < 1 << row_level ; row++ ) 
      	for ( mat_level_t col = 0 ; col < 1 << col_level ; col++ ) 
	  (((R[row][col] /= 2) -= 3) *= 4) += 5;

      for ( mat_level_t row = 0 ; row < 1 << row_level ; row++ ) 
      	for ( mat_level_t col = 0 ; col < 1 << col_level ; col++ ) 
      	  assert_equal ( (M[row][col]/(ScalarType)2 - (ScalarType)3) 
			 * (ScalarType)4 + (ScalarType)5 ,  R[row][col] ,
			 "verifying chained update operators" );
    }
  }

  clean_matrix_store ( );
  uint64_t end_count = num_matrices_in_store ( );  
  assert_equal ( end_count , start_count , "verifying matrix store clean up" );

  return 0;
}
