#include "indexable.hh"

/*!
 *\cond
 */
Element_t::Element_t ( Indexable_t &rf , uint64_t r , uint64_t c ) : 
  ref ( rf ) , row ( r ) , col ( c ) 
{
}

ScalarType Element_t::operator= ( ScalarType v ) 
{
  ref.set ( row , col , v );
  return v;
}

ScalarType Element_t::operator= ( const Element_t &e )
{
  ScalarType v = e;
  ref.set ( row , col , v );
  return v;
}

Element_t::operator ScalarType ( ) const
{
  return ref.get ( row , col );
}

Element_t &Element_t::operator+= ( ScalarType v )
{
  *this = *this + v;
  return *this;
}

Element_t &Element_t::operator-= ( ScalarType v )
{
  *this = *this - v;
  return *this;
}

Element_t &Element_t::operator*= ( ScalarType v )
{
  *this = *this * v;
  return *this;
}

Element_t &Element_t::operator/= ( ScalarType v )
{
  *this = *this / v;
  return *this;
}

Const_Slice_t::Const_Slice_t ( const Indexable_t &rf , uint64_t r ) :
  ref ( rf ) , row ( r ) 
{
}

ScalarType Const_Slice_t::operator[] ( uint64_t c ) const
{
  return ref.get ( row , c );
}

mat_level_t Const_Slice_t::get_level ( ) const 
{
  return ref.get_col_level ( );
}

Slice_t::Slice_t ( Indexable_t &rf , uint64_t r ) :
  ref ( rf ) , row ( r ) 
{
}

ScalarType Slice_t::operator[] ( uint64_t c ) const
{
  return ref.get ( row , c );
}

Element_t Slice_t::operator[] ( uint64_t c ) 
{
  return Element_t ( ref , row , c );
}

mat_level_t Slice_t::get_level ( ) const 
{
  return ref.get_col_level ( );
}

Const_Slice_t Indexable_t::operator[] ( uint64_t r ) const
{
  return Const_Slice_t ( *this , r );
}

Slice_t Indexable_t::operator[] ( uint64_t r ) 
{
  return Slice_t ( *this , r );
}

bool operator== ( const Element_t &e1 , const Element_t &e2 )
{
  return e1.ref.get( e1.row , e1.col ) == e2.ref.get ( e2.row , e2.col );
}
bool operator== ( const Element_t &e1 , ScalarType e2 )
{
  return e1.ref.get( e1.row , e1.col ) == e2;
}
bool operator== ( ScalarType e1 , const Element_t &e2 )
{
  return e2 == e1;
}

bool operator!= ( const Element_t &e1 , const Element_t &e2 )
{
  return !(e1 == e2);
}
bool operator!= ( const Element_t &e1 , ScalarType e2 )
{
  return !(e1 == e2);
}
bool operator!= ( ScalarType e1 , const Element_t &e2 )
{
  return !(e1 == e2);
}

ScalarType operator+ ( const Element_t &e1 , const Element_t &e2 )
{
  return (ScalarType)e1 + (ScalarType)e2;
}
ScalarType operator+ ( const Element_t &e1 ,       ScalarType e2 )
{
  return (ScalarType)e1 + (ScalarType)e2;
}
ScalarType operator+ (       ScalarType e1 , const Element_t &e2 )
{
  return (ScalarType)e1 + (ScalarType)e2;
}

ScalarType operator* ( const Element_t &e1 , const Element_t &e2 )
{
  return (ScalarType)e1 * (ScalarType)e2;
}
ScalarType operator* ( const Element_t &e1 ,       ScalarType e2 )
{
  return (ScalarType)e1 * (ScalarType)e2;
}
ScalarType operator* (       ScalarType e1 , const Element_t &e2 )
{
  return (ScalarType)e1 * (ScalarType)e2;
}

ScalarType operator- ( const Element_t &e1 , const Element_t &e2 )
{
  return (ScalarType)e1 - (ScalarType)e2;
}
ScalarType operator- ( const Element_t &e1 ,       ScalarType e2 )
{
  return (ScalarType)e1 - (ScalarType)e2;
}
ScalarType operator- (       ScalarType e1 , const Element_t &e2 )
{
  return (ScalarType)e1 - (ScalarType)e2;
}

ScalarType operator/ ( const Element_t &e1 , const Element_t &e2 )
{
  return (ScalarType)e1 / (ScalarType)e2;
}
ScalarType operator/ ( const Element_t &e1 ,       ScalarType e2 )
{
  return (ScalarType)e1 / (ScalarType)e2;
}
ScalarType operator/ (       ScalarType e1 , const Element_t &e2 )
{
  return (ScalarType)e1 / (ScalarType)e2;
}

std::ostream &operator<< ( std::ostream &out , const Element_t &e )
{
  out << (ScalarType)e;
  return out;
}
/*!
 * \endcond
 */
