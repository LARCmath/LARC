#include "larc.hh"

#include <algorithm>

#include <matrix_store.h>
#include <matmath.h>

LARC_t::LARC_t ( mat_level_t rlevel , mat_level_t clevel ) : 
  data_ref  ( get_zero_matrix_ptr ( rlevel , clevel ) ) ,
  row_level ( rlevel                                  ) , 
  col_level ( clevel                                  )
{
  set_hold_matrix ( this->data_ref );
}

LARC_t::LARC_t ( const LARC_t &rhs ) : 
  data_ref  ( rhs.data_ref  ) ,
  row_level ( rhs.row_level ) ,
  col_level ( rhs.col_level )
{
  set_hold_matrix ( this->data_ref );  
}

LARC_t::~LARC_t ( )
{
  release_hold_matrix ( this->data_ref );
}

LARC_t &LARC_t::operator+= ( const LARC_t &rhs )
{
  *this = *this + rhs;
  return *this;
}

LARC_t &LARC_t::operator-= ( const LARC_t &rhs )
{
  *this = *this - rhs;
  return *this;
}

LARC_t &LARC_t::operator*= ( const LARC_t &rhs )
{
  *this = *this * rhs;
  return *this;
}

LARC_t &LARC_t::operator*= ( ScalarType v )
{
  *this = *this * v;
  return *this;
}

LARC_t &LARC_t::operator= ( LARC_t rhs )
{
  std::swap ( this->data_ref , rhs.data_ref );
  return *this;
}

LARC_t LARC_t::kron ( const LARC_t &rhs ) const
{
  return LARC_t ( kronecker_product ( this->data_ref , rhs.data_ref ) );
}

LARC_t LARC_t::adjoint ( ) const
{
  return LARC_t ( matrix_adjoint ( this->data_ref ) );
}

ScalarType LARC_t::trace ( ) const
{
  return matrix_trace(this->data_ref);
}

static mat_add_t xjoin  ( mat_add_t x , mat_add_t y ) { return join  ( x , y ); }
LARC_t LARC_t::join ( const LARC_t &right ) const
{
  return LARC_t ( xjoin ( this->data_ref , right.data_ref ) );
}

static mat_add_t xstack ( mat_add_t x , mat_add_t y ) { return stack ( x , y ); }
LARC_t LARC_t::stack ( const LARC_t &lower ) const
{
  return LARC_t ( xstack ( this->data_ref , lower.data_ref ) );
}

ScalarType LARC_t::get ( uint64_t r , uint64_t c ) const
{
  return matrix_get_value ( this->data_ref , r , c );
}

void LARC_t::set ( uint64_t r , uint64_t c , ScalarType v )
{
  LARC_t temp = LARC_t ( get_matPTR_from_oldMatPTR_newVal_and_coords ( data_ref , r , c , v ) );
  std::swap ( this->data_ref , temp.data_ref );
}

mat_level_t LARC_t::get_row_level ( ) const 
{
  return row_level;
}

mat_level_t LARC_t::get_col_level ( ) const 
{
  return col_level;
}

LARC_t::LARC_t ( mat_add_t ptr ) : 
  data_ref  ( ptr                      ) , 
  row_level ( matrix_row_level ( ptr ) ) , 
  col_level ( matrix_col_level ( ptr ) )
{
  set_hold_matrix ( this->data_ref );  
}

LARC_t operator+ ( const LARC_t &op1 , const LARC_t &op2 )
{
  return LARC_t ( matrix_add ( op1.data_ref , op2.data_ref ) );
}

LARC_t operator- ( const LARC_t &op1 , const LARC_t &op2 )
{
  return op1 + (-1 * op2);
}

LARC_t operator* ( const LARC_t &op1 , const LARC_t &op2 )
{
  return LARC_t ( matrix_mult ( op1.data_ref , op2.data_ref ) );
}

LARC_t operator*  ( ScalarType    op1 , const LARC_t &op2 )
{
  return LARC_t ( scalar_mult ( get_valMatPTR_from_val ( op1 ) , op2.data_ref ) );
}

LARC_t operator*  ( const LARC_t &op1 , ScalarType    op2 )
{
  return operator* ( op2 , op1 );
}

std::ostream &operator<< ( std::ostream &out , const LARC_t &op )
{
  for ( uint64_t row = 0 ; row < ( (uint64_t)1 ) << op.row_level ; row++ ) {
    for ( uint64_t col = 0 ; col < ( (uint64_t)1 ) << op.col_level ; col++ ) {
      out << op[row][col] << " ";
    }
    out << std::endl;
  }

  out << std::endl;

  return out;

}

