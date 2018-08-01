#include <stdio.h>
#include <stdlib.h>

#include "matrix_store.h"
#include "matmath.h"

// Return the (row,col) entry in m.
ScalarType matrix_get_value ( mat_add_t m , mat_level_t row , mat_level_t col )
{
  mat_level_t rows = 1 << matrix_row_level ( m );
  mat_level_t cols = 1 << matrix_col_level ( m );
  
  if ( row >= rows ) {
    printf ( "ERROR: in matrix_get_value the row %d is too large for matrix %p\n",
	     row , m );
    exit ( 1 );
  }

  if ( col >= cols ) {
    printf ( "ERROR: in matrix_get_value the col %d is too large for matrix %p\n",
	     col , m );
    exit ( 1 );
  }

  if ( 1 == rows && 1 == cols ) 
    return matrix_trace ( m );

  if ( 2*row < rows ) {
    if ( 2*col < cols ) return matrix_get_value ( matrix_sub ( m , 0 ) , row , col          );
    else                return matrix_get_value ( matrix_sub ( m , 1 ) , row , col - cols/2 );
  } 

  if ( 2*col < cols ) return matrix_get_value ( matrix_sub ( m , 2 ) , row - rows/2 , col );
  return matrix_get_value ( matrix_sub ( m , 3 ) , row - rows/2 , col - cols/2 );
}

// Create a new matrix that is the same as m but with the entry at
// (row,col) set to v.
mat_add_t  matrix_set_value ( mat_add_t m , mat_level_t row , mat_level_t col , ScalarType v )
{ 
  if ( matrix_get_value ( m , row , col ) == v ) return m;

  mat_level_t rows = 1 << matrix_row_level ( m );
  mat_level_t cols = 1 << matrix_col_level ( m );
  
  if ( row >= rows ) {
    printf ( "ERROR: in matrix_set_value the row %d is too large for matrix %p\n",
	     row , m );
    exit ( 1 );
  }

  if ( col >= cols ) {
    printf ( "ERROR: in matrix_set_value the col %d is too large for matrix %p\n",
	     col , m );
    exit ( 1 );
  }

  mat_add_t panels[4];

  // Handle the scalar case

  if ( 1 == rows && 1 == cols ) 
    return matrix_get_ptr_scalar ( v );

  // Handle the one-row case

  if ( 1 == rows ) {

    panels[0] = matrix_sub ( m , 0 );
    panels[1] = matrix_sub ( m , 1 );
    
    if ( 2*col < cols ) panels[0] = matrix_set_value ( panels[0] , row , col          , v );
    else                panels[1] = matrix_set_value ( panels[1] , row , col - cols/2 , v );

    return join ( panels[0] , panels[1] );
  }

  // Handle the one-column case

  if ( 1 == cols ) {

    panels[0] = matrix_sub ( m , 0 );
    panels[1] = matrix_sub ( m , 2 );

    if ( 2*row < rows ) panels[0] = matrix_set_value ( panels[0] , row          , col , v );
    else                panels[1] = matrix_set_value ( panels[1] , row - rows/2 , col , v );

    return stack ( panels[0] , panels[1] );
  }

  // At this point matrices have at least two columns and at least two
  // rows, and the data structure at the next level down looks like four
  // copies of this one.

  for ( int quad = 0 ; quad < 4 ; quad++ ) 
    panels[quad] = matrix_sub ( m , quad );

  if ( 2*row < rows ) {
    if ( 2*col < cols ) panels[0] = matrix_set_value ( matrix_sub ( m , 0 ) , row          , col          , v );
    else                panels[1] = matrix_set_value ( matrix_sub ( m , 1 ) , row          , col - cols/2 , v );
  } else {
    if ( 2*col < cols ) panels[2] = matrix_set_value ( matrix_sub ( m , 2 ) , row - rows/2 , col          , v );
    else                panels[3] = matrix_set_value ( matrix_sub ( m , 3 ) , row - rows/2 , col - cols/2 , v );
  }

  return matrix_get_ptr_panel ( panels , matrix_row_level ( m ) , matrix_col_level ( m ) );
}

