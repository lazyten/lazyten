//
// Copyright (C) 2016-17 by the linalgwrap authors
//
// This file is part of linalgwrap.
//
// linalgwrap is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// linalgwrap is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with linalgwrap. If not, see <http://www.gnu.org/licenses/>.
//

#pragma once
#include "Base/Interfaces/MutableMemoryVector_i.hh"
#include "Base/Interfaces/Transposed.hh"
#include "TypeUtils/mat_vec_apply_enabled_t.hh"
#include "linalgwrap/LazyMatrixProduct.hh"
#include "linalgwrap/LazyMatrixSum.hh"
#include "linalgwrap/Matrix_i.hh"
#include "linalgwrap/StoredMatrix_i.hh"
#include <krims/GenMap.hh>

namespace linalgwrap {

// Forward declarations
template <typename StoredMatrix>
class LazyMatrixSum;

template <typename StoredMatrix>
class LazyMatrixProduct;

/** \brief Generic LazyMatrixExpression class
 *
 * Abstract base class for all lazy matrix expressions.
 * Implements a slightly modified paradigm of Expression templates.
 *
 * \tparam StoredMatrix The type to use for stored matrices.
 * */
template <typename StoredMatrix>
class LazyMatrixExpression : public Matrix_i<typename StoredMatrix::scalar_type> {
  static_assert(std::is_base_of<StoredMatrix_i<typename StoredMatrix::scalar_type>,
                                StoredMatrix>::value,
                "StoredMatrix is not a child of StoredMatrix_i of the same scalar "
                "type");

 public:
  typedef StoredMatrix stored_matrix_type;
  typedef typename stored_matrix_type::vector_type vector_type;
  typedef Matrix_i<typename StoredMatrix::scalar_type> base_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::size_type size_type;

  /** The pointer type to use for pointers to LazyMatrixExpressions */
  typedef std::unique_ptr<LazyMatrixExpression<StoredMatrix>>
        lazy_matrix_expression_ptr_type;

  /** \name Matrix properties
   */
  ///@{
  // TODO rename to has_transpose_mode
  /** Are operation modes Transposed::Trans and Transposed::ConjTrans
   *  supported for this matrix type
   *
   * These operation modes are important for the functions apply,
   * apply_inverse, mmult and extract_block
   **/
  virtual bool has_transpose_operation_mode() const = 0;

  /** Does this Matrix have an implemented inverse apply method?
   *
   * The idea is to allow special matrices to offer a method to take
   * advantage of their internal structure when one wants to apply
   * the inverse of such a matrix to a multivector. Possible examples
   * are diagonal and tridiagonal matrices.
   *
   * \note The inverse_apply should not be iterative or otherwise
   *       implicit. (With the exception of ImplicitlyInvertibleMatrix,
   *       which serves exactly the purpose to attach an implicit
   *       inverse to another matrix)
   **/
  virtual bool has_apply_inverse() const { return false; }

  // TODO has_element_access (i.e. operator() and extract_block
  //      has_mmult (i.e. has matrix-matrix multiplication
  //
  //      flags supported_operations() -> return flags
  //      bool has_support_for(flags)
  ///@}

  /** \name Data access
   */
  ///@{
  /** Extract a block of a matrix and (optionally) add it to
   * a different matrix.
   *
   *  Extract a block of values from this matrix, where the block
   *  size is given by the size of the Stored Matrix ``in`` and the
   *  element at which we start the extraction is given by
   *  ``start_row`` and ``start_col``. In other words if ``in`` is a
   *  2x2 matrix and ``start_row == 2`` and ``start_col==1`` we extract
   *  the four elements (2,1),(2,2), (3,1) and (3,2), given that
   *  mode == Normal
   *
   *  Loosely speaking we perform
   *  \[ M = c_M \cdot M + (A^{mode})_{rowrange,colrange} \]
   *  where
   *    - rowrange = [start_row, start_row+in.n_rows() ) and
   *    - colrange = [start_col, start_col+in.n_cols() )
   *
   *  There are two important details:
   *    - Mode is applied *before* forming the block of values,
   *      e.g. start_row and start_col are with respect to the
   *      transformed matrix if mode == Transp.
   *    - If $c_\text{this} == 0$ the function **may** simply scale $M$.
   *    - If $c_M == 0$ the function **must** overwrite $M$, since the
   *      values of M may be uninitialised.
   *  This is similar to the apply function in Matrix_i.
   *
   *  \note The function assumes that the sparsity pattern in ``in``
   *  already has the correct shape or is dynamically extended such
   *  that the addition does not hit non-existing elements.
   *
   *  \param M    Matrix to add the values to or to assign them into.
   *  \param start_row  The row index of the first element to extract
   *  \param start_col  The column index of the first element to extract
   *  \param c_this     The coefficient to multiply this matrix with
   *                    before extracting.
   *  \param c_M        Coefficient to multiply M with before adding.
   *
   */
  virtual void extract_block(
        stored_matrix_type& M, const size_type start_row, const size_type start_col,
        const Transposed mode = Transposed::None,
        const scalar_type c_this = Constants<scalar_type>::one,
        const scalar_type c_M = Constants<scalar_type>::zero) const = 0;

  /** \brief Convert the expression to a stored matrix
   *
   * Achieved by calling extract_block on the whole matrix.
   */
  virtual explicit operator stored_matrix_type() const {
    stored_matrix_type ret(this->n_rows(), this->n_cols(), false);
    extract_block(ret, 0, 0);
    return ret;
  }
  ///@}

  /** \name Matrix application and matrix products
   */
  ///@{
  /** \brief Compute the Matrix-Multivector application
   *
   * Loosely performs the operation
   * \[ y = c_this \cdot A^\text{mode} \cdot x + c_y \cdot y. \]
   * where mode gives the mode how the current object is to be
   * applied (normal, transposed, conjugate transposed)
   *
   * The implementation
   *  - **may** avoid the operator and just scale $y$ if $c_this == 0$.
   *  - **must** overwrite $y$ if $c_y == 0$ since the memory of
   *    y may be uninitialised and we want to avoid NaN values in y
   *    to be of effect)
   *
   *  \note This function is only implemented for stored vectors or
   *  modifyable vectors, whose storage can be accessed via a pointer.
   *
   * \note Whenever the virtual apply method is overwritten, this method
   * should be implemented as well as it assures that conversion to
   * MultiVector<MutableMemoryVector_i<scalar_type>> can actually occur
   * automatically.
   */
  template <typename VectorIn, typename VectorOut,
            mat_vec_apply_enabled_t<LazyMatrixExpression, VectorIn, VectorOut>...>
  void apply(const MultiVector<VectorIn>& x, MultiVector<VectorOut>& y,
             const Transposed mode = Transposed::None,
             const scalar_type c_this = Constants<scalar_type>::one,
             const scalar_type c_y = Constants<scalar_type>::zero) const {
    MultiVector<const MutableMemoryVector_i<scalar_type>> x_wrapped(x);
    MultiVector<MutableMemoryVector_i<scalar_type>> y_wrapped(y);
    apply(x_wrapped, y_wrapped, mode, c_this, c_y);
  }

  /** \brief Compute the Matrix-Multivector application
   *
   * Specialisation of the above function for MutableMemoryVector_i
   * This is the version child classes need to implement.
   */
  virtual void apply(const MultiVector<const MutableMemoryVector_i<scalar_type>>& x,
                     MultiVector<MutableMemoryVector_i<scalar_type>>& y,
                     const Transposed mode = Transposed::None,
                     const scalar_type c_this = Constants<scalar_type>::one,
                     const scalar_type c_y = Constants<scalar_type>::zero) const = 0;

  // TODO Split functions like apply and apply_inverse, ... into an outer interface
  // function and an apply_kernel, apply_inverse_kernel function which does the actual
  // work. This way we can do the assertions, simple cases and checks once in the outer
  // function and only need to do the real work in the kernel function.

  /** \brief Compute the application of the inverse of the matrix
   *  (or the inverse of the transpose of the matrix) to a MultiVector
   *
   * Loosely performs the operation
   * \[ y = c_this \cdot (A^{-1})^\text{mode} \cdot x + c_y \cdot y. \]
   * where mode gives the mode how the current object is to be
   * applied (normal, transposed, conjugate transposed)
   *
   * The implementation
   *  - **may** avoid the operator and just scale $y$ if $c_this == 0$.
   *  - **must** overwrite $y$ if $c_y == 0$ since the memory of
   *    y may be uninitialised and we want to avoid NaN values in y
   *    to be of effect)
   *
   *  \note This function is only implemented for stored vectors or
   *  modifyable vectors, whose storage can be accessed via a pointer.
   *
   * \note Whenever the virtual apply_inverse method is overwritten,
   * this method should be implemented as well as it assures that conversion to
   * MultiVector<MutableMemoryVector_i<scalar_type>> can actually occur
   * automatically.
   *
   */
  template <typename VectorIn, typename VectorOut,
            mat_vec_apply_enabled_t<LazyMatrixExpression, VectorIn, VectorOut>...>
  void apply_inverse(const MultiVector<VectorIn>& /*x*/, MultiVector<VectorOut>& /*y*/,
                     const Transposed /*mode*/ = Transposed::None,
                     const scalar_type /*c_this*/ = 1,
                     const scalar_type /*c_y*/ = 0) const {
    // In general there is no easy way to do an inverse:
    assert_throw(false, krims::ExcDisabled("The apply_inverse function is in general "
                                           "very expensive and is only implemented in "
                                           "some cases. Use the function "
                                           "has_apply_inverse() to check when."));
  }

  /** \brief Compute the application of the inverse of the matrix.
   * This is a specialisation of the above method.
   */
  virtual void apply_inverse(
        const MultiVector<const MutableMemoryVector_i<scalar_type>>& /*x*/,
        MultiVector<MutableMemoryVector_i<scalar_type>>& /*y*/,
        const Transposed /*mode */ = Transposed::None, const scalar_type /*c_this */ = 1,
        const scalar_type /*c_y */ = 0) const {
    // In general there is no easy way to do an inverse:
    assert_throw(false, krims::ExcDisabled("The apply_inverse function is in general "
                                           "very expensive and is only implemented in "
                                           "some cases. Use the function "
                                           "has_apply_inverse() to check when."));
  }

  /** Perform a matrix-matrix product.
   *
   * Loosely performs the operation
   * \[ out = c_this \cdot A^\text{mode} \cdot in + c_out \cdot out. \]
   * where mode gives the mode how the current object is to be
   * applied (normal, transposed, conjugate transposed)
   *
   * The implementation
   *  - **may** avoid the operator and just scale $y$ if $c_this == 0$.
   *  - **must** overwrite $y$ if $c_out == 0$, since the memory of out
   *    may be uninitialised.
   */
  virtual void mmult(const stored_matrix_type& in, stored_matrix_type& out,
                     const Transposed mode = Transposed::None,
                     const scalar_type c_this = Constants<scalar_type>::one,
                     const scalar_type c_out = Constants<scalar_type>::zero) const = 0;
  ///@}

  /** \brief Update the internal data of all objects in this expression
   * given the GenMap
   * */
  virtual void update(const krims::GenMap& map) = 0;

  /** \brief Clone the expression
   *
   * Return a clone of the current object as a pointer of type
   * lazy_matrix_expression_ptr_type
   */
  virtual lazy_matrix_expression_ptr_type clone() const = 0;
};

//@{
/** \brief struct representing a type (std::true_type, std::false_type) which
 *         indicates whether Matrix is a lazy matrix
 *
 * The definition is done using SFINAE, such that even for types not having a
 * typedef stored_matrix_type this expression is valid.
 *  */
template <typename Matrix, typename = void>
struct IsLazyMatrix : public std::false_type {};

template <typename Matrix>
struct IsLazyMatrix<Matrix, krims::VoidType<typename Matrix::stored_matrix_type>>
      : public std::is_base_of<LazyMatrixExpression<typename Matrix::stored_matrix_type>,
                               Matrix> {};
//@}

//
// Multiplication
//
/** Multiply two lazy matrix expressions */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator*(const LazyMatrixExpression<StoredMatrix>& lhs,
                                          const LazyMatrixExpression<StoredMatrix>& rhs) {

  // Construct product with one factor:
  LazyMatrixProduct<StoredMatrix> prod{lhs};

  // add the other factor and return the result.
  prod.push_factor(rhs);
  return prod;
}

/** Multiply lazy times stored */
template <typename StoredMatrix>
StoredMatrix operator*(const LazyMatrixExpression<StoredMatrix>& lhs,
                       const StoredMatrix& rhs) {
  StoredMatrix out(lhs.n_rows(), rhs.n_cols(), false);
  lhs.mmult(rhs, out);
  return out;
}

/** Multiply stored times lazy */
template <typename StoredMatrix>
StoredMatrix operator*(const StoredMatrix& lhs,
                       const LazyMatrixExpression<StoredMatrix>& rhs) {
  assert_dbg(false, krims::ExcDisabled(
                          "The operation \"StoredMatrix * LazyMatrixExpression\" is "
                          "disabled because it is usually pretty expensive. "
                          "Use \"StoredMatrix * LazyMatrixExpression * StoredMatrix\" "
                          "instead if possible. Otherwise you can enforce"
                          "\"StoredMatrix * LazyMatrixExpression\" by explicitly "
                          "converting the LazyMatrixExpression into a StoredMatrix, "
                          "i.e. \"lhs * static_cast<StoredMatrix>(rhs)\""));
  return lhs * static_cast<StoredMatrix>(rhs);
}

/** Scale a lazy matrix expression */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator*(const LazyMatrixExpression<StoredMatrix>& m,
                                          typename StoredMatrix::scalar_type s) {

  // Construct product with one factor:
  LazyMatrixProduct<StoredMatrix> prod{m};

  // and scale it:
  prod *= s;
  return prod;
}

/** Scale a lazy matrix expression */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator*(typename StoredMatrix::scalar_type s,
                                          const LazyMatrixExpression<StoredMatrix>& m) {
  return m * s;
}

/** Perform a Matrix-Vector product */
template <typename StoredMatrix, typename Vector,
          typename = typename std::enable_if<
                IsStoredVector<Vector>::value &&
                std::is_same<typename StoredMatrix::scalar_type,
                             typename Vector::scalar_type>::value>::type>
Vector operator*(const LazyMatrixExpression<StoredMatrix>& m, const Vector& v);

/** Perform a Matrix-MultiVector product */
template <typename StoredMatrix, typename Vector,
          typename = typename std::enable_if<
                IsStoredVector<Vector>::value &&
                std::is_same<typename StoredMatrix::scalar_type,
                             typename Vector::scalar_type>::value>::type>
MultiVector<typename std::remove_const<Vector>::type> operator*(
      const LazyMatrixExpression<StoredMatrix>& m, const MultiVector<Vector>& mv);

//
// Division by scalar
//
/** \brief Devide a lazy matrix by a scalar  */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator/(const LazyMatrixExpression<StoredMatrix>& m,
                                          typename StoredMatrix::scalar_type s) {

  typedef typename StoredMatrix::scalar_type scalar_type;
  scalar_type inverse = Constants<scalar_type>::one / s;
  return m * inverse;
}

//
// Addition
//
/** Add two lazy matrix expressions */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator+(const LazyMatrixExpression<StoredMatrix>& lhs,
                                      const LazyMatrixExpression<StoredMatrix>& rhs) {

  // Construct sum with one term:
  LazyMatrixSum<StoredMatrix> sum{lhs};

  // add the other factor and return the result.
  sum.push_term(rhs);
  return sum;
}

/** Add a lazy matrix and a stored matrix */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator+(const LazyMatrixExpression<StoredMatrix>& lhs,
                                      const StoredMatrix& rhs) {

  // Construct sum with one term:
  LazyMatrixSum<StoredMatrix> sum{lhs};

  // add the other factor and return the result.
  sum.push_term(rhs);
  return sum;
}

/** Add a lazy matrix and a stored matrix */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator+(const StoredMatrix& lhs,
                                      const LazyMatrixExpression<StoredMatrix>& rhs) {
  return rhs + lhs;
}

//
// Subtraction
//
/** Subtract two lazy matrix expressions */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator-(const LazyMatrixExpression<StoredMatrix>& lhs,
                                      const LazyMatrixExpression<StoredMatrix>& rhs) {

  // Construct sum with one term:
  LazyMatrixSum<StoredMatrix> sum{lhs};

  // typedef the scalar type
  typedef typename StoredMatrix::scalar_type scalar_type;

  // subtract the other factor and return the result.
  sum.push_term(rhs, -Constants<scalar_type>::one);
  return sum;
}

/** Subtract a stored matrix from a lazy matrix */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator-(const LazyMatrixExpression<StoredMatrix>& lhs,
                                      const StoredMatrix& rhs) {

  // Construct sum with one term:
  LazyMatrixSum<StoredMatrix> sum{lhs};

  // typedef the scalar type
  typedef typename StoredMatrix::scalar_type scalar_type;

  // add the other factor and return the result.
  sum.push_term(rhs, -Constants<scalar_type>::one);
  return sum;
}

/** Subtract a lazy matrix from a stored matrix */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator-(const StoredMatrix& lhs,
                                      const LazyMatrixExpression<StoredMatrix>& rhs) {
  // Construct sum with one term:
  LazyMatrixSum<StoredMatrix> sum{lhs};

  // typedef the scalar type
  typedef typename StoredMatrix::scalar_type scalar_type;

  // add the other factor and return the result.
  sum.push_term(rhs, -Constants<scalar_type>::one);
  return sum;
}

//
// Unary operator-
//
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator-(const LazyMatrixExpression<StoredMatrix>& mat) {
  typedef typename StoredMatrix::scalar_type scalar_type;
  return -Constants<scalar_type>::one * mat;
}

//
// ------------------------------------------------
//

//
// Multiplication
//
template <typename StoredMatrix, typename Vector, typename>
MultiVector<typename std::remove_const<Vector>::type> operator*(
      const LazyMatrixExpression<StoredMatrix>& m, const MultiVector<Vector>& mv) {
  assert_size(mv.n_elem(), m.n_cols());
  MultiVector<typename std::remove_const<Vector>::type> out(m.n_rows(), mv.n_vectors(),
                                                            false);
  m.apply(mv, out);
  return out;
}

template <typename StoredMatrix, typename Vector, typename>
Vector operator*(const LazyMatrixExpression<StoredMatrix>& m, const Vector& v) {
  assert_size(v.size(), m.n_cols());
  Vector out(m.n_rows(), false);
  MultiVector<const MutableMemoryVector_i<typename Vector::scalar_type>> v_wrapped(v);
  MultiVector<MutableMemoryVector_i<typename Vector::scalar_type>> out_wrapped(out);
  m.apply(v_wrapped, out_wrapped);
  return out;
}

}  // namespace linalgwrap
