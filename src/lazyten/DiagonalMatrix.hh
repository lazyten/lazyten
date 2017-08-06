//
// Copyright (C) 2016-17 by the lazyten authors
//
// This file is part of lazyten.
//
// lazyten is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// lazyten is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with lazyten. If not, see <http://www.gnu.org/licenses/>.
//

#pragma once
#include "LazyMatrixExpression.hh"
#include <krims/SubscriptionPointer.hh>

// TODO generalisation for non-square diagonal matrices

namespace lazyten {
/** \brief Make a diagonal matrix out of the vector of the diagonal elements */
template <typename StoredMatrix>
class DiagonalMatrix : public LazyMatrixExpression<StoredMatrix> {
 public:
  typedef LazyMatrixExpression<StoredMatrix> base_type;
  typedef typename base_type::stored_matrix_type stored_matrix_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::lazy_matrix_expression_ptr_type
        lazy_matrix_expression_ptr_type;
  typedef typename StoredMatrix::vector_type stored_vector_type;

  /** Construct from reference to diagonal elements */
  DiagonalMatrix(const stored_vector_type& diagonal)
        : m_diagonal_ptr{krims::make_subscription(diagonal, "DiagonalMatrix")} {}

  /** Construct by moving diagonal elements inside */
  DiagonalMatrix(stored_vector_type&& diagonal)
        : m_diagonal_ptr{std::make_shared<const stored_vector_type>(diagonal)} {}

  /** Number of rows */
  size_type n_rows() const override { return m_diagonal_ptr->size(); }

  /** Number of columns */
  size_type n_cols() const override { return m_diagonal_ptr->size(); }

  /** Element access */
  scalar_type operator()(size_type row, size_type col) const override {
    if (row == col) return (*m_diagonal_ptr)[row];
    return Constants<scalar_type>::zero;
  }

  //
  // LazyMatrixExpression interface
  //
  /** Are operation modes Transposed::Trans and Transposed::ConjTrans
   *  supported for this matrix type.
   **/
  bool has_transpose_operation_mode() const override { return true; }

  /** Is inverse_apply available for this matrix type */
  bool has_apply_inverse() const override { return true; }

  /** Extract a block of a matrix and (optionally) add it to
   * a different matrix.
   *
   *  Loosely speaking we perform
   *  \[ M = c_M \cdot M + (A^{mode})_{rowrange,colrange} \]
   *  where
   *    - rowrange = [start_row, start_row+in.n_rows() ) and
   *    - colrange = [start_col, start_col+in.n_cols() )
   *
   * More details can be found in the same function in
   * LazyMatrixExpression
   */
  void extract_block(stored_matrix_type& M, const size_type start_row,
                     const size_type start_col, const Transposed mode = Transposed::None,
                     const scalar_type c_this = Constants<scalar_type>::one,
                     const scalar_type c_M = Constants<scalar_type>::zero) const override;

  /** \brief Compute the Matrix-Multivector application -- generic version
   *
   * Loosely speaking we perform
   * \[ y = c_this \cdot A^\text{mode} \cdot x + c_y \cdot y. \]
   *
   * See LazyMatrixExpression for more details
   *
   * \note Whenever the virtual apply method is overwritten, this method
   * should be implemented as well as it assures that conversion to
   * MultiVector<MutableMemoryVector_i<scalar_type>> can actually occur
   * automatically.
   */
  template <typename VectorIn, typename VectorOut,
            mat_vec_apply_enabled_t<DiagonalMatrix, VectorIn, VectorOut>...>
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
   * Loosely speaking we perform
   * \[ y = c_this \cdot A^\text{mode} \cdot x + c_y \cdot y. \]
   *
   * See LazyMatrixExpression for more details
   */
  void apply(const MultiVector<const MutableMemoryVector_i<scalar_type>>& x,
             MultiVector<MutableMemoryVector_i<scalar_type>>& y,
             const Transposed mode = Transposed::None,
             const scalar_type c_this = Constants<scalar_type>::one,
             const scalar_type c_y = Constants<scalar_type>::zero) const override;

  /** \brief Compute the Inverse-Multivector application
   *
   * Loosely speaking we perform
   * \[ y = c_this \cdot (A^{-1})^\text{mode} \cdot x + c_y \cdot y. \]
   *
   * See LazyMatrixExpression for more details
   */
  template <typename VectorIn, typename VectorOut,
            mat_vec_apply_enabled_t<DiagonalMatrix, VectorIn, VectorOut>...>
  void apply_inverse(const MultiVector<VectorIn>& x, MultiVector<VectorOut>& y,
                     const Transposed mode = Transposed::None,
                     const scalar_type c_this = 1, const scalar_type c_y = 0) const {
    MultiVector<const MutableMemoryVector_i<scalar_type>> x_wrapped(x);
    MultiVector<MutableMemoryVector_i<scalar_type>> y_wrapped(y);
    apply_inverse(x_wrapped, y_wrapped, mode, c_this, c_y);
  }

  /** \brief Compute the Inverse-Multivector application
   *
   * Loosely speaking we perform
   * \[ y = c_this \cdot (A^{-1})^\text{mode} \cdot x + c_y \cdot y. \]
   *
   * See LazyMatrixExpression for more details
   */
  virtual void apply_inverse(
        const MultiVector<const MutableMemoryVector_i<scalar_type>>& x,
        MultiVector<MutableMemoryVector_i<scalar_type>>& y,
        const Transposed mode = Transposed::None, const scalar_type c_this = 1,
        const scalar_type c_y = 0) const override;

  /** Perform a matrix-matrix product.
   *
   * Loosely performs the operation
   * \[ out = c_this \cdot A^\text{mode} \cdot in + c_out \cdot out. \]
   *
   * See LazyMatrixExpression for more details
   */
  void mmult(const stored_matrix_type& in, stored_matrix_type& out,
             const Transposed mode = Transposed::None,
             const scalar_type c_this = Constants<scalar_type>::one,
             const scalar_type c_out = Constants<scalar_type>::zero) const override;

  /** Update method
   *
   * \note does nothing
   */
  void update(const krims::GenMap&) override {}

  /** Clone function */
  lazy_matrix_expression_ptr_type clone() const override {
    // return a copy enwrapped in the pointer type
    return lazy_matrix_expression_ptr_type(new DiagonalMatrix(*this));
  }

 private:
  krims::RCPWrapper<const stored_vector_type> m_diagonal_ptr;
};

template <typename StoredVector>
DiagonalMatrix<typename StoredVector::type_family::template matrix<
      typename StoredVector::scalar_type>>
make_diagmat(const StoredVector& v) {
  return DiagonalMatrix<typename StoredVector::type_family::template matrix<
        typename StoredVector::scalar_type>>(v);
}

template <typename StoredVector>
DiagonalMatrix<typename StoredVector::type_family::template matrix<
      typename StoredVector::scalar_type>>
make_diagmat(StoredVector&& v) {
  return DiagonalMatrix<typename StoredVector::type_family::template matrix<
        typename StoredVector::scalar_type>>(v);
}

//
// ----------------------------------------------------------------------
//

template <typename StoredMatrix>
void DiagonalMatrix<StoredMatrix>::extract_block(
      stored_matrix_type& M, const size_type start_row, const size_type start_col,
      const Transposed mode, const scalar_type c_this, const scalar_type c_M) const {
  assert_finite(c_this);
  assert_finite(c_M);
  // check that we do not overshoot the indices
  if (mode == Transposed::Trans || mode == Transposed::ConjTrans) {
    assert_greater_equal(start_row + M.n_rows(), n_cols());
    assert_greater_equal(start_col + M.n_cols(), n_rows());
  } else {
    assert_greater_equal(start_row + M.n_rows(), n_rows());
    assert_greater_equal(start_col + M.n_cols(), n_cols());
  }
  assert_sufficiently_tested(mode != Transposed::ConjTrans);

  // For empty matrices there is nothing to do
  if (M.n_rows() == 0 || M.n_cols() == 0) return;

  // Set elements of M to zero (if c_M == 0)
  // or scale them according to c_M.
  // This deals entirely with the coefficient c_M
  detail::scale_or_set(M, c_M);

  if (c_this == Constants<scalar_type>::zero) return;

  const size_type start = std::max(start_row, start_col);
  const size_type end = std::min(start_row + M.n_rows(), start_col + M.n_cols());

  for (size_type i = start; i < end; ++i) {
    switch (mode) {
      case Transposed::None:
      case Transposed::Trans:
        M(i - start_row, i - start_col) += c_this * (*m_diagonal_ptr)[i];
        break;
      case Transposed::ConjTrans:
        // A variant of std::conj, which does not return a complex
        // data type if scalar is real only.
        krims::ConjFctr conj;
        M(i - start_row, i - start_col) += c_this * conj((*m_diagonal_ptr)[i]);
        break;
    }  // mode
  }    // i
}

template <typename StoredMatrix>
void DiagonalMatrix<StoredMatrix>::apply(
      const MultiVector<const MutableMemoryVector_i<scalar_type>>& x,
      MultiVector<MutableMemoryVector_i<scalar_type>>& y, const Transposed mode,
      const scalar_type c_this, const scalar_type c_y) const {
  assert_finite(c_this);
  assert_finite(c_y);
  assert_size(x.n_vectors(), y.n_vectors());
  if (mode == Transposed::Trans || mode == Transposed::ConjTrans) {
    assert_size(x.n_elem(), this->n_cols());
    assert_size(y.n_elem(), this->n_rows());
  } else {
    assert_size(x.n_elem(), this->n_rows());
    assert_size(y.n_elem(), this->n_cols());
  }
  assert_sufficiently_tested(mode != Transposed::ConjTrans);

  // Scale the current values of out or set them to zero
  // (if c_y == 0): We are now done with c_y and do not
  // need to worry about it any more in this function
  for (auto& vec : y) detail::scale_or_set(vec, c_y);

  // if c_this == 0 we are done
  if (c_this == Constants<scalar_type>::zero) return;

  for (size_type vi = 0; vi < y.n_vectors(); ++vi) {
    for (size_type ei = 0; ei < y.n_elem(); ++ei) {
      const auto& vecin = x[vi];
      auto& vecout = y[vi];
      switch (mode) {
        case Transposed::None:
        case Transposed::Trans:
          vecout[ei] += c_this * vecin[ei] * (*m_diagonal_ptr)[ei];
          break;
        case Transposed::ConjTrans:
          // A variant of std::conj, which does not return a complex
          // data type if scalar is real only.
          krims::ConjFctr conj;
          vecout[ei] += c_this * vecin[ei] * conj((*m_diagonal_ptr)[ei]);
          break;
      }  // mode
    }    // ei
  }      // vi
}

template <typename StoredMatrix>
void DiagonalMatrix<StoredMatrix>::apply_inverse(
      const MultiVector<const MutableMemoryVector_i<scalar_type>>& x,
      MultiVector<MutableMemoryVector_i<scalar_type>>& y, const Transposed mode,
      const scalar_type c_this, const scalar_type c_y) const {
  assert_finite(c_this);
  assert_finite(c_y);
  assert_size(x.n_vectors(), y.n_vectors());
  assert_size(n_rows(), n_cols());  // only quadratic mat have inverses
  assert_size(x.n_elem(), this->n_rows());
  assert_size(y.n_elem(), this->n_rows());
  assert_sufficiently_tested(mode != Transposed::ConjTrans);

  // Scale the current values of out or set them to zero
  // (if c_y == 0): We are now done with c_y and do not
  // need to worry about it any more in this function
  for (auto& vec : y) detail::scale_or_set(vec, c_y);

  // if c_this == 0 we are done
  if (c_this == Constants<scalar_type>::zero) return;

  for (size_type vi = 0; vi < y.n_vectors(); ++vi) {
    for (size_type ei = 0; ei < y.n_elem(); ++ei) {
      const auto& vecin = x[vi];
      auto& vecout = y[vi];
      switch (mode) {
        case Transposed::None:
        case Transposed::Trans:
          vecout[ei] += c_this * vecin[ei] / (*m_diagonal_ptr)[ei];
          break;
        case Transposed::ConjTrans:
          // A variant of std::conj, which does not return a complex
          // data type if scalar is real only.
          krims::ConjFctr conj;
          vecout[ei] += c_this * vecin[ei] / conj((*m_diagonal_ptr)[ei]);
          break;
      }  // mode
    }    // ei
  }      // vi
}

template <typename StoredMatrix>
void DiagonalMatrix<StoredMatrix>::mmult(const stored_matrix_type& in,
                                         stored_matrix_type& out, const Transposed mode,
                                         const scalar_type c_this,
                                         const scalar_type c_out) const {
  assert_finite(c_this);
  assert_finite(c_out);
  assert_size(in.n_cols(), out.n_cols());
  if (mode == Transposed::Trans || mode == Transposed::ConjTrans) {
    assert_size(n_rows(), in.n_rows());
    assert_size(n_cols(), out.n_rows());
  } else {
    assert_size(n_cols(), in.n_rows());
    assert_size(n_rows(), out.n_rows());
  }
  assert_sufficiently_tested(mode != Transposed::ConjTrans);

  // Scale the current values of out or set them to zero
  // (if c_out == 0): We are now done with c_out.
  detail::scale_or_set(out, c_out);

  // if c_this == 0 we are done
  if (c_this == Constants<scalar_type>::zero) return;

  for (auto it = std::begin(in); it != std::end(in); ++it) {
    switch (mode) {
      case Transposed::None:
      case Transposed::Trans:
        out(it.row(), it.col()) += c_this * *it * (*m_diagonal_ptr)[it.row()];
        break;
      case Transposed::ConjTrans:
        // A variant of std::conj, which does not return a complex
        // data type if scalar is real only.
        krims::ConjFctr conj;
        out(it.row(), it.col()) += c_this * *it * conj((*m_diagonal_ptr)[it.row()]);
        break;
    }  // mode
  }    // it
}

}  // namespace lazyten
