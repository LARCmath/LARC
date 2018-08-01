#include <stdlib.h>

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

int  main ( int argc , char *argv[] ) 
{
  // Initialize the larc library.
  initialize ( );

  // Create a matrix
  ScalarType one = 1.0;
  ScalarType two = 2.0;

  int64_t mID_scalarOne = matrix_get_matrixID_from_scalar ( one );
  int64_t mID_scalarTwo = matrix_get_matrixID_from_scalar ( two );

  set_hold_matrix_from_matrixID ( mID_scalarOne );
  set_hold_matrix_from_matrixID ( mID_scalarTwo );

  int64_t mID_testMatrix = matrix_get_matrixID_from_panel ( mID_scalarOne , mID_scalarTwo , mID_scalarTwo , mID_scalarOne , 1 , 1 );
  set_hold_matrix_from_matrixID ( mID_testMatrix );
  
  print_matrix_naive_by_matrixID ( mID_testMatrix );

  return 0;
}
