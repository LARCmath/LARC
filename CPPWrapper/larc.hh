#ifndef __LARC_HH
#define __LARC_HH

#include <ostream>

#include <larc.h>
#include <stdint.h>

#include "indexable.hh"

/*!
 * \ingroup larcpp
 * \brief A C++ wrapper for underlying LARC matrices.
 *
 * This class provides some C++ conveniences for working with LARC
 * matrices. These are
 *
 * 1. Setting and releasing holds in the constructors and
 * destructor. Clients of this code need not do this.
 *
 * 2. Overloading basic operations of +, - and *
 *
 * 3. By inheriting fron Indexable_t, providing a means to access
 * elements within the matrix with array notation. Note that this is
 * an expensive operation, setting a matrix element in particular
 * creates as many new matrices as there levels in the matrix. This is
 * intended for operations on small matrices.
 *
 * With this one can write code that looks like
 * \code
 * LARC_t A(1,1) 
 * LARC_t B(1,1)
 * LARC_t C(1,1)
 *
 * A[0][0] = 1; A[0][1] = 5; A[1][0] = 8; A[1][1] = 17;
 * B[0][0] = 2; B[0][1] = 4; B[1][0] = 6; B[1][1] = 8;
 *
 * C = (A + 2*B) * (B - 2*A)
 *
 * cout << C[1][1] << endl;
 * \endcode
 */
class LARC_t : public Indexable_t
{
public:
  
  /*!
   * \ingroup larcpp
   * \brief A constructor that specifies the size of the matrix. The
   * initial values are set to zero.
   * \param rlevel The number of row levels.
   * \param clevel The number of column levels.
   */ 
  LARC_t  ( mat_level_t rlevel , mat_level_t clevel );
  
  /*!
   * \ingroup larcpp
   * \brief A standard copy constructor.
   * \param rhs A reference to the LARC_t instance from which we
   * should initialize.
   */
  LARC_t  ( const LARC_t &rhs                       );
  
  /*!
   * \ingroup larcpp
   * \brief The destructor.
   */ 
  ~LARC_t (                                         );

  /*!
   * \ingroup larcpp
   * \brief Replace *this with the sum of *this and rhs.
   * \param rhs The object to add to *this.
   * \return The updated value of *this.
   */
  LARC_t &operator+= ( const LARC_t &rhs );

  /*!
   * \ingroup larcpp
   * \brief Replace *this with the difference of *this and rhs.
   * \param rhs The object to subtract from *this.
   * \return The updated value of *this.
   */
  LARC_t &operator-= ( const LARC_t &rhs );
  
  /*!
   * \ingroup larcpp
   * \brief Replace *this with the product of *this and rhs.
   * \param rhs The object to multiply on the right of *this.
   * \return The updated value of *this.
   */
  LARC_t &operator*= ( const LARC_t &rhs );

  /*!
   * \ingroup larcpp
   * \brief Replace *this with the product of *this and v.
   * \param v The scalar value to multiply by *this.
   * \return The updated value *this.
   */
  LARC_t &operator*= ( ScalarType    v   );
  
  /*!
   * \ingroup larcpp
   * \brief Replace *this with rhs.
   * \param rhs The new value.
   * \return The update value of *this.
   */
  LARC_t &operator=  ( LARC_t        rhs );

  /*!
   * \ingroup larcpp
   * \brief Compute a Kronecker product.
   * \param rhs The second operand in the Kronecker product.
   * \return The Kronecker product of *this and rhs.
   */
  LARC_t kron ( const LARC_t &rhs ) const;

  /*!
   * \ingroup larcpp
   * \brief Return the adjoint.
   * \return The adjoint of this matrix.
   */
  LARC_t adjoint ( ) const;

  /*!
   * \ingroup larcpp
   * \brief Return the trace.
   * \return The trace of this matrix.
   */
  ScalarType trace ( ) const;
  
  /*!
   * \ingroup larcpp
   * \brief Form a new matrix by appending a matrix to the right of
   * this one.
   * \param right The matrix to append to the right of this
   * matrix. Note that it must be the same size, both the number of
   * rows and the number of columns must agree.
   * \return A new matrix formed by joining a matrix on the right.
   */ 
  LARC_t join ( const LARC_t &right ) const;

  /*!
   * \ingroup larcpp
   * \brief Form a new matrix by appending a matrix below
   * this one.
   * \param lower The matrix to append below of this matrix. Note that
   * it must be the same size, both the number of rows and the number
   * of columns must agree.
   * \return A new matrix formed by joining a matrix on the bottom.
   */ 
  LARC_t stack ( const LARC_t &lower ) const;
  
  /*!
   * \ingroup larcpp
   * \brief Get a value in the underlying LARC matrix.
   * \param r The row of the value to access.
   * \param c The column of the value to access.
   * \return The value at r, c.
   */ 
  ScalarType get ( uint64_t r , uint64_t c                ) const;

  /*!
   * \ingroup larcpp
   * \brief Set a value in the underlying LARC matrix.
   * \param r The row of the value to set.
   * \param c The column of the value to set.
   * \param v The value to set at r, c.
   */
  void       set ( uint64_t r , uint64_t c , ScalarType v );

  /*!
   * \ingroup larcpp
   * \brief Get the row_level of the underlying LARC matrix.
   */
  mat_level_t get_row_level ( ) const;

  /*!
   * \ingroup larcpp
   * \brief Get the col_level of the underlying LARC matrix.
   */
  mat_level_t get_col_level ( ) const;

private:

  LARC_t ( mat_add_t ptr );

  mat_add_t data_ref;
  mat_level_t row_level, col_level;

  friend LARC_t        operator+  ( const LARC_t &op1 , const LARC_t &op2 );
  friend LARC_t        operator-  ( const LARC_t &op1 , const LARC_t &op2 );
  friend LARC_t        operator*  ( const LARC_t &op1 , const LARC_t &op2 );
  friend LARC_t        operator*  ( ScalarType    op1 , const LARC_t &op2 );
  friend LARC_t        operator*  ( const LARC_t &op1 , ScalarType    op2 );
  friend std::ostream &operator<< ( std::ostream &out , const LARC_t &op  );
};

/*!
 * \ingroup larcpp
 * \brief Add two LARC matrices.
 * \param op1 The first operand.
 * \param op2 The second operand.
 * \return The sum of op1 and op2.
 */
LARC_t        operator+  ( const LARC_t &op1 , const LARC_t &op2 );

/*!
 * \ingroup larcpp
 * \brief Subtract two LARC matrices.
 * \param op1 The first operand.
 * \param op2 The second operand.
 * \return The difference op1 - op2.
 */
LARC_t        operator-  ( const LARC_t &op1 , const LARC_t &op2 );

/*!
 * \ingroup larcpp
 * \brief Multiply two LARC matrices.
 * \param op1 The first operand.
 * \param op2 The second operand.
 * \return The product op1 * op2.
 */
LARC_t        operator*  ( const LARC_t &op1 , const LARC_t &op2 );

/*!
 * \ingroup larcpp
 * \brief Multiply a LARC matrix by a scalar.
 * \param op1 The scalar.
 * \param op2 The operand.
 * \return The product op1 * op2.
 */
LARC_t        operator*  ( ScalarType    op1 , const LARC_t &op2 );

/*!
 * \ingroup larcpp
 * \brief Multiply a LARC matrix by a scalar.
 * \param op1 The operand.
 * \param op2 The scalar.
 * \return The product op1 * op2.
 */
LARC_t        operator*  ( const LARC_t &op1 , ScalarType    op2 );

/*!
 * \ingroup larcpp
 * \brief Write a matrix in a niave form to an output stream.
 * \param out The stream to which to write.
 * \param op The matrix to write. 
 * \return The stream after the matrix has been written.
 */
std::ostream &operator<< ( std::ostream &out , const LARC_t &op  );

/*!
 *\cond
 */
void          swap       (       LARC_t &op1 ,       LARC_t &op2 );
/*!
 *\endcond
 */

#endif // #ifndef __LARC_HH
