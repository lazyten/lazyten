//
// Copyright (C) 2016 by the linalgwrap authors
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
#include "linalgwrap/Constants.hh"
#include "linalgwrap/DefaultMatrixIterator.hh"
#include "linalgwrap/Matrix_i.hh"
#include <krims/TypeUtils.hh>

namespace linalgwrap {

// TODO define a set of optional functions which make performance better
//      e.g. - an in-memory transpose() function
//           - transpose-multiply
//           - whatever else seems sensible

// Forward-declare the interface class
template <typename Scalar>
class Matrix_i;

/** \brief Interface class for a matrix which is actually stored in memory
 * in some way
 *
 * We expect any implementing class to also provide the following constructors:
 * - Construct matrix of fixed size and optionally fill with zeros or leave
 *   memory unassigned:
 *   ```
 *   StoredMatrix_i(n_rows, n_cols, fill_zero);
 *   ```
 * - Construct a matrix of the same size as a SmallMatrix and copy all entries
 *   from the SmallMatrix over, optionally providing a tolerance below which
 *   the entries are considered to be zero (The latter is useful for CRS
 *   matrices)
 *   ```
 *   StoredMatrix_i(const SmallMatrix&)
 *   StoredMatrix_i(const SmallMatrix&, scalar_type tolerance)
 *   ```
 *
 * All implementing classes should further provide the following types
 *   - vector_type  The type of the corresponding vector;
 *   - type_family  The type of the type family of corresponding
 *                  vector types (complex and real vectors)
 *
 * The following methods should be implemented:
 *
 * The following methods should be implemented between matrices of the
 * implementing type:
 * - Addition(+), Subtraction(-), Scalar multiplication and scalar division
 * - In-place addition(+=) and subtraction(-=)
 * - In-place scalar multiplication (*=)
 *
 * The following methods should be implemented:
 * - ```void extract_block(stored_matrix_type& M,
 *                         size_type start_row, size_type start_col,
 *                         Transposed mode = Transposed::None,
 *                         scalar_type c_this = Constants<scalar_type>::one,
 *                         scalar_type c_M = Constants<scalar_type>::zero);
 *   ```
 *   This method should extract a block of size M.n_rows() times M.n_cols() of
 *   values into the matrix M. Loosely speaking it should perform
 *   \[ M = c_\text{this} * A^\text{mode} + c_M * M \]
 *   where A is the matrix represented by this object. If the function
 *   has_transpose_operation_mode() returns false only normal operation mode
 *   (no transpose) needs to be supported, else all methods (normal,
 *   transpose, conjugate transpose). For more information
 *   see the documentation of this method in LazyMatrixExpression.
 *
 * - ```void mmult(const stored_matix_type& in, stored_matrix_type& out,
 *             const Transposed mode = Transposed::None,
 *             const scalar_type c_this = Constants<scalar_type>::one,
 *             const scalar_type c_out = Constants<scalar_type>::zero) const;
 *   ```
 *   This method should multiply two matrices and add the result to a third.
 *   (generalised matrix-matrix product, gemm), i.e. loosely needs to perform
 *   \[ out = c_out * out + c_this * A^\text{mode} * in \]
 *    See LazyMatrixExpression for more details.
 *
 * - ```template <typename VectorIn, typename VectorOut,
 *            mat_vec_apply_enabled_t<stored_matrix_type, VectorIn,
 * VectorOut>...>
 *  void apply(const MultiVector<VectorIn>& x, MultiVector<VectorOut>& y,
 *             const Transposed mode = Transposed::None,
 *             const scalar_type c_this = Constants<scalar_type>::one,
 *             const scalar_type c_y = Constants<scalar_type>::zero) const;
 *   ```
 *   Generalised matrix-vector application, performs
 *   \[ y = c_y * y + c_this * A^\text{mode} * x
 *   See LazyMatrixExpression for more details.
 *
 * Note that the operator() functions in derived classes are expected to return
 * zero even if an element is known to be zero by some sparsity pattern or
 * similar. Modification of a non-existing element should fail, however.
 *
 */
template <typename Scalar>
class StoredMatrix_i : public Matrix_i<Scalar> {
 public:
  typedef Matrix_i<Scalar> base_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::size_type size_type;

  //! The iterator type
  typedef DefaultMatrixIterator<StoredMatrix_i<Scalar>> iterator;

  //! The const_iterator type
  typedef typename base_type::const_iterator const_iterator;

  /** \name Matrix properties
   */
  ///@{
  /** Are operation modes Transposed::Trans and Transposed::ConjTrans
   *  supported for this matrix type
   *
   * These operation modes are important for the functions apply,
   * apply_inverse, mmult and extract_block
   *
   * \note It is highly recommended that stored matrices support
   * all operation modes!
   **/
  virtual bool has_transpose_operation_mode() const { return false; }

  /** Does this Matrix have an implemented inverse apply method?
   *
   * The idea is to allow special matrices to offer a method to take
   * advantage of their internal structure when one wants to apply
   * the inverse of such a matrix to a multivector. Possible examples
   * are diagonal and tridiagonal matrices.
   *
   * \note The inverse_apply should never be iterative.
   **/
  virtual bool has_apply_inverse() const { return false; }
  ///@}

  /** \name Data access */
  ///@{
  /** Read-write access to elements */
  virtual scalar_type& operator()(size_type row, size_type col) = 0;

  /** \brief Read-write access to vectorised matrix object
   *
   * Access the element in row-major ordering (i.e. the matrix is
   * traversed row by row)
   */
  virtual scalar_type& operator[](size_type i);

  /** \brief Read-only access to vectorised matrix object
   *
   * Access the element in row-major ordering (i.e. the matrix is
   * traversed row by row)
   * */
  scalar_type operator[](size_type i) const override { return base_type::operator[](i); }
  ///@}

  /** \name Standard operations */
  ///@{
  /** Set all elements to zero
   *
   * Overload this to get a more performant implementation.
   *
   * TODO generalise with a lambda and an arbitrary scalar
   * */
  virtual void set_zero() { std::fill(begin(), end(), Constants<scalar_type>::zero); }

  /** \brief Symmetrise the matrix, that is form (A + A^T)/2.
   *
   * Only sensible for square matrices.
   * */
  virtual void symmetrise();
  ///@}

  //
  // Iterators
  //
  /** Return an iterator to the beginning */
  iterator begin() { return iterator(*this, {0, 0}); }

  /** Return a const iterator to the beginning */
  const_iterator begin() const { return base_type::cbegin(); }

  /** Return an iterator to the end */
  iterator end() { return iterator(*this); }

  /** Return a const iterator to the end */
  const_iterator end() const { return base_type::cend(); }

  // TODO
  //   function to get actual number of non-zero entries
  //   function to get estimated/implicitly known number of non-zero entries
};

//@{
/** \brief struct representing a type (std::true_type, std::false_type) which
 *indicates
 *  whether T is a stored matrix
 *
 * The definition is done using SFINAE, such that even for types not having a
 *typedef
 * scalar_type this expression is valid.
 *  */
template <typename Matrix, typename = void>
struct IsStoredMatrix : public std::false_type {};

template <typename Matrix>
struct IsStoredMatrix<Matrix, krims::VoidType<typename Matrix::scalar_type>>
      : public std::is_base_of<StoredMatrix_i<typename Matrix::scalar_type>, Matrix> {};
//@}

//@{
/** \brief Convert a matrix to a stored matrix (in case it is none)
 *  else return a reference to the original object
 *
 *  \note Can be used in order to make sure that one is dealing
 *  with a stored matrix in any case without making an unnecessary
 *  copy in case it already is a stored matrix
 */
template <typename Matrix, typename StoredMatrix = typename Matrix::stored_matrix_type>
StoredMatrix as_stored(const Matrix& m) {
  return static_cast<StoredMatrix>(m);
}

template <typename StoredMatrix,
          typename = typename std::enable_if<IsStoredMatrix<StoredMatrix>::value>::type>
StoredMatrix& as_stored(StoredMatrix& m) {
  return m;
}

template <typename StoredMatrix,
          typename = typename std::enable_if<IsStoredMatrix<StoredMatrix>::value>::type>
const StoredMatrix& as_stored(const StoredMatrix& m) {
  return m;
}
//@}

//
// -------------------------------------------------------------
//

template <typename Scalar>
typename StoredMatrix_i<Scalar>::scalar_type& StoredMatrix_i<Scalar>::operator[](
      size_type i) {
  // Check that we do not overshoot.
  assert_range(0u, i, this->n_cols() * this->n_rows());

  const size_type i_row = i / this->n_cols();
  const size_type i_col = i % this->n_cols();
  return (*this)(i_row, i_col);
}

template <typename Scalar>
void StoredMatrix_i<Scalar>::symmetrise() {
  assert_dbg(this->n_rows() == this->n_cols(), ExcMatrixNotSquare());

  for (size_type i = 0; i < this->n_rows(); ++i) {
    for (size_type j = i + 1; j < this->n_rows(); ++j) {
      const scalar_type res = ((*this)(i, j) + (*this)(j, i)) / scalar_type(2.);
      (*this)(i, j) = (*this)(j, i) = res;
    }
  }
}

}  // namespace liblinalg
