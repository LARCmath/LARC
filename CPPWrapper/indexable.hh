#ifndef __INDEXABLE_HH
#define __INDEXABLE_HH

#include <stdint.h>
#include <ostream>
#include <larc.h>

/*!
 * \cond
 */
class Indexable_t;

class Element_t {
public: 
  Element_t ( Indexable_t &rf , uint64_t r , uint64_t c );
  ScalarType operator= ( const Element_t &e );
  ScalarType operator= ( ScalarType v );
  operator ScalarType ( ) const;

  Element_t &operator+= ( ScalarType v );
  Element_t &operator-= ( ScalarType v );
  Element_t &operator*= ( ScalarType v );
  Element_t &operator/= ( ScalarType v );

private:
  Element_t ( const Element_t &rhs );
  Indexable_t &ref;
  uint64_t row , col;

  friend class Slice_t;

  friend bool       operator== ( const Element_t &e1 , const Element_t &e2 );
  friend bool       operator== ( const Element_t &e1 ,       ScalarType e2 );
  friend bool       operator== (       ScalarType e1 , const Element_t &e2 );

  friend bool       operator!= ( const Element_t &e1 , const Element_t &e2 );
  friend bool       operator!= ( const Element_t &e1 ,       ScalarType e2 );
  friend bool       operator!= (       ScalarType e1 , const Element_t &e2 );

  friend ScalarType operator+  ( const Element_t &e1 , const Element_t &e2 );
  friend ScalarType operator+  ( const Element_t &e1 ,       ScalarType e2 );
  friend ScalarType operator+  (       ScalarType e1 , const Element_t &e2 );

  friend ScalarType operator*  ( const Element_t &e1 , const Element_t &e2 );
  friend ScalarType operator*  ( const Element_t &e1 ,       ScalarType e2 );
  friend ScalarType operator*  (       ScalarType e1 , const Element_t &e2 );

  friend ScalarType operator-  ( const Element_t &e1 , const Element_t &e2 );
  friend ScalarType operator-  ( const Element_t &e1 ,       ScalarType e2 );
  friend ScalarType operator-  (       ScalarType e1 , const Element_t &e2 );

  friend ScalarType operator/  ( const Element_t &e1 , const Element_t &e2 );
  friend ScalarType operator/  ( const Element_t &e1 ,       ScalarType e2 );
  friend ScalarType operator/  (       ScalarType e1 , const Element_t &e2 );
  
  friend std::ostream &operator<< ( std::ostream &out , const Element_t &e );
};

class Const_Slice_t
{
public:
  Const_Slice_t ( const Indexable_t &rf , uint64_t r );
  ScalarType operator[] ( uint64_t c ) const;
  mat_level_t get_level ( ) const;

private:
  const Indexable_t &ref;
  uint64_t row;
};

class Slice_t
{
public:
  Slice_t ( Indexable_t &rf , uint64_t r );
  ScalarType operator[] ( uint64_t c ) const;
  Element_t  operator[] ( uint64_t c );
  mat_level_t get_level ( ) const;

private:
  Indexable_t &ref;
  uint64_t row;
};
/*!
 * \endcond
 */

/*!
 * \ingroup larcpp
 * \brief This class provides an index structure for the deriving
 * class.
 * 
 * This class provides rank-two indexing for the derving
 * class. Deriving classes that define appropriate get and set methods
 * can be accessed using array notation. 
 * 
 * There are two significant ways this structure differs from an array
 * of arrays.
 * 
 * 1. The assignment A = B, creates a copy of B and assigns it to
 * A. Setting an element of A does not change B.
 *
 * 2. Indexing is implemented with proxy objects, so a partial index
 * has the type of the proxy object. A useful type is auto, as in 
 *
 * \code
 * auto B = A[2];
 * \endcode
 *
 * A fully indexed object can be treated as a ScalarType.
 */
class Indexable_t
{
public:

  /*!
   * \ingroup larcpp
   * \brief The deriving class must implement the following getter.
   * \param r The row for which one would like the value.
   * \param c The column for which one would like the value.
   * \return The value at r,c
   */
  virtual ScalarType  get ( uint64_t r , uint64_t c                ) const = 0;

  /*!
   * \ingroup larcpp
   * \brief The deriving class must implement the following getter.
   * \param r The row for which one would like to set a value.
   * \param c The column for which one would like to set a value.
   * \param v The new value to set at r,c
   */
  virtual void        set ( uint64_t r , uint64_t c , ScalarType v )       = 0;
  
  /*!
   * \ingroup larcpp
   * \brief The deriving class must implememnt the following
   * \return The log-base-two of the number of rows.
   */
  virtual mat_level_t get_row_level ( ) const = 0;

  /*!
   * \ingroup larcpp
   * \brief The deriving class must implememnt the following
   * \return The log-base-two of the number of columns.
   */
  virtual mat_level_t get_col_level ( ) const = 0;
  /*!
   * \cond
   */
  Const_Slice_t operator[] ( uint64_t r ) const;
  Slice_t       operator[] ( uint64_t r );
  /*!
   * \endcond
   */
};

#endif // #ifndef __INDEXABLE_HH
