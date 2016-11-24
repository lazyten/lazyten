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
#include "linalgwrap/LazyMatrixExpression.hh"
#include "linalgwrap/detail/GenericFunctionals.hh"

namespace linalgwrap {
/** \brief Interface of lazy matrices
 *
 * All classes which implement this interface perform lazy evaluation, e.g. the
 * product is not directly performed, but much rather it is stored as a
 * lazy_matrix_product expression. Copying the classes should not copy the
 * actual data, but only references to the data (i.e. the classes should contain
 * only small amounts of data that lives on the stack).
 *
 * Assumptions the library makes when using the matrices implementing this
 * interface:
 *    - Copying these type of objects is cheap.
 *    - Multiplication and addition of these objects is associative
 *
 * From the abstact class ``Matrix_i`` we expect abstract lazy matrices to
 * implement:
 *   - ``size_type n_rows() const``
 *   - ``size_type n_cols() const``
 *   - ``scalar_type operator() (size_type row, size_type col) const``
 *
 * From ``LazyMatrixExpression`` we expect the implementation of
 *   - ``lazy_matrix_expression_ptr_type clone() const``
 * see the appropriate classes for details what these methods are required to
 * do.
 *
 * For performance reasons it is highly recommended to overload the
 * following functions as well:
 *    - ```bool has_transpose_operation_mode() const;```
 *      Should return true iff Transposed::Trans and Transposed::ConjTrans
 *      operation modes are supported.
 *    - ```extract_block(stored_matrix_type& M, const size_type start_row,
 *        const size_type start_col, const Transposed mode,
 *        const scalar_type c_this, const scalar_type c_M) const;```
 *      Extract a block of values from the matrix and add them into
 *      a stored matrix.
 *    - ```void apply(const MultiVector<const
 * MutableMemoryVector_i<scalar_type>>& x,
 *        MultiVector<MutableMemoryVector_i<scalar_type>>& y,
 *        const Transposed mode, const scalar_type c_this,
 *        const scalar_type c_y) const;```
 *      Apply a matrix to a MultiVector (gemv style).
 *    - ```template <typename VectorIn, typename VectorOut,
 *            mat_vec_apply_enabled_t<LazyMatrix,
 *                                    VectorIn, VectorOut>...>
 *        void apply(const MultiVector<VectorIn>& x, MultiVector<VectorOut>& y,
 *             const Transposed mode, const scalar_type c_this,
 *             const scalar_type c_y) const;
 *      ```
 *      This method needs to be implement whenever the virtual method apply
 *      is implemented as the conversion to MultiVector<VectorIn>& cannot
 *      be made implicit. Also it allows for easier future transition
 *      to a virtual-free LazyMatrixProduct/LazyMatrixSum class.
 *    - ```void mmult(const stored_matrix_type& in, stored_matrix_type& out,
 *             const Transposed mode,const scalar_type c_this,
 *             const scalar_type c_out) const;```
 *     Perform a gemm-like matrix-matrix product.
 *
 * \tparam StoredMatrix   The type of stored matrix to use
 */
template <typename StoredMatrix>
class LazyMatrix_i : public LazyMatrixExpression<StoredMatrix> {

 public:
  typedef LazyMatrixExpression<StoredMatrix> base_type;
  typedef typename base_type::stored_matrix_type stored_matrix_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::lazy_matrix_expression_ptr_type
        lazy_matrix_expression_ptr_type;

  //
  // Partial implementation of the interface of a LazyMatrixExpression
  //
  /** Are operation modes Transposed::Trans and Transposed::ConjTrans
   *  supported for this matrix type.
   **/
  bool has_transpose_operation_mode() const override;

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
   *
   * \note This default version only uses operator() and is
   * hence very slow.
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
            mat_vec_apply_enabled_t<LazyMatrix_i, VectorIn, VectorOut>...>
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
   *
   * \note This default version only uses operator() and is
   * hence very slow.
   */
  void apply(const MultiVector<const MutableMemoryVector_i<scalar_type>>& x,
             MultiVector<MutableMemoryVector_i<scalar_type>>& y,
             const Transposed mode = Transposed::None,
             const scalar_type c_this = Constants<scalar_type>::one,
             const scalar_type c_y = Constants<scalar_type>::zero) const override;

  /** Perform a matrix-matrix product.
   *
   * Loosely performs the operation
   * \[ out = c_this \cdot A^\text{mode} \cdot in + c_out \cdot out. \]
   *
   * See LazyMatrixExpression for more details
   *
   * \note This default version only uses operator() and is
   * hence very slow.
   */
  void mmult(const stored_matrix_type& in, stored_matrix_type& out,
             const Transposed mode = Transposed::None,
             const scalar_type c_this = Constants<scalar_type>::one,
             const scalar_type c_out = Constants<scalar_type>::zero) const override;

  /** \brief Update the internal data of all objects in this expression
   *         given the ParameterMap
   *
   *   This function raises an ExcNotImplemented exception.
   * */
  virtual void update(const krims::ParameterMap&) override {
    assert_dbg(false, krims::ExcNotImplemented());
  }
};

//
// --------------------------------------------------------------
//

template <typename StoredMatrix>
bool LazyMatrix_i<StoredMatrix>::has_transpose_operation_mode() const {
  // We default to false here, but nevertheless implement
  // them below.
  //
  // This way one can easily enable the branches by overloading
  // this method, but if one forgets to implement
  // Transposed::Trans and Transposed::ConjTrans in an
  // overload of apply, mmult or extract_block
  // it is still fine.
  return false;
}

template <typename StoredMatrix>
void LazyMatrix_i<StoredMatrix>::extract_block(
      stored_matrix_type& M, const size_type start_row, const size_type start_col,
      const Transposed mode, const scalar_type c_this, const scalar_type c_M) const {
  assert_dbg(mode == Transposed::None || has_transpose_operation_mode(),
             ExcUnsupportedOperationMode(mode));
  assert_finite(c_this);
  assert_finite(c_M);
  // check that we do not overshoot the indices
  if (mode == Transposed::Trans || mode == Transposed::ConjTrans) {
    assert_greater_equal(start_row + M.n_rows(), this->n_cols());
    assert_greater_equal(start_col + M.n_cols(), this->n_rows());
  } else {
    assert_greater_equal(start_row + M.n_rows(), this->n_rows());
    assert_greater_equal(start_col + M.n_cols(), this->n_cols());
  }
  assert_sufficiently_tested(mode != Transposed::ConjTrans);

  // For empty matrices there is nothing to do
  if (M.n_rows() == 0 || M.n_cols() == 0) return;

  if (c_this == Constants<scalar_type>::zero) {
    detail::scale_or_set(M, c_M);
    return;
  }  // c_this == 0

  if (c_M == Constants<scalar_type>::zero) {
    // If c_out is zero we are not allowed to read
    // from out's actual memory as the values could
    // be NaN, hence overwrite it with zeros here
    // such that we know for sure.
    // (see LazyMatrixExpression for Details)
    M.set_zero();
  }

  for (size_type row = 0; row < M.n_rows(); ++row) {
    for (size_type col = 0; col < M.n_cols(); ++col) {
      switch (mode) {
        case Transposed::None:
          M(row, col) =
                c_this * (*this)(start_row + row, start_col + col) + c_M * M(row, col);
          break;
        case Transposed::Trans:
          M(row, col) =
                c_this * (*this)(start_col + col, start_row + row) + c_M * M(row, col);
          break;
        case Transposed::ConjTrans:
          // A variant of std::conj, which does not return a complex
          // data type if scalar is real only.
          detail::ConjFctr mconj;
          M(row, col) = c_this * mconj((*this)(start_col + col, start_row + row)) +
                        c_M * M(row, col);
          break;
      }  // mode
    }    // col
  }      // row
}

template <typename StoredMatrix>
void LazyMatrix_i<StoredMatrix>::apply(
      const MultiVector<const MutableMemoryVector_i<scalar_type>>& x,
      MultiVector<MutableMemoryVector_i<scalar_type>>& y, const Transposed mode,
      const scalar_type c_this, const scalar_type c_y) const {
  assert_dbg(mode == Transposed::None || has_transpose_operation_mode(),
             ExcUnsupportedOperationMode(mode));
  assert_finite(c_this);
  assert_finite(c_y);
  assert_size(x.n_vectors(), y.n_vectors());
  if (mode == Transposed::Trans || mode == Transposed::ConjTrans) {
    assert_size(x.n_elem(), this->n_rows());
    assert_size(y.n_elem(), this->n_cols());
  } else {
    assert_size(x.n_elem(), this->n_cols());
    assert_size(y.n_elem(), this->n_rows());
  }
  assert_sufficiently_tested(mode != Transposed::ConjTrans);

  // Scale the current values of out or set them to zero
  // (if c_y == 0): We are now done with c_y and do not
  // need to worry about it any more in this function
  for (auto& vec : y) detail::scale_or_set(vec, c_y);

  // if c_this == 0 we are done
  if (c_this == Constants<scalar_type>::zero) return;

  // The size of the resulting vectors
  const size_type size = mode == Transposed::None ? this->n_rows() : this->n_cols();
  for (size_type veci = 0; veci < x.n_vectors(); ++veci) {
    const auto& vec = x[veci];
    auto& out = y[veci];
    for (size_type row = 0; row < size; ++row) {
      scalar_type sum = 0;
      for (size_type k = 0; k < x.n_elem(); ++k) {
        switch (mode) {
          case Transposed::None:
            sum += (*this)(row, k) * vec(k);
            break;
          case Transposed::Trans:
            sum += (*this)(k, row) * vec(k);
            break;
          case Transposed::ConjTrans:
            // A variant of std::conj, which does not return a
            // complex data type if scalar is real only.
            detail::ConjFctr mconj;
            sum += mconj((*this)(k, row)) * vec(k);
            break;
        }  // mode
      }    // k
      out(row) += c_this * sum;
    }  // row
  }    // veci
}

template <typename StoredMatrix>
void LazyMatrix_i<StoredMatrix>::mmult(const stored_matrix_type& in,
                                       stored_matrix_type& out, const Transposed mode,
                                       const scalar_type c_this,
                                       const scalar_type c_out) const {
  assert_dbg(mode == Transposed::None || has_transpose_operation_mode(),
             ExcUnsupportedOperationMode(mode));
  assert_finite(c_this);
  assert_finite(c_out);
  assert_size(in.n_cols(), out.n_cols());
  if (mode == Transposed::Trans || mode == Transposed::ConjTrans) {
    assert_size(this->n_rows(), in.n_rows());
    assert_size(this->n_cols(), out.n_rows());
  } else {
    assert_size(this->n_cols(), in.n_rows());
    assert_size(this->n_rows(), out.n_rows());
  }
  assert_sufficiently_tested(mode != Transposed::ConjTrans);

  // Scale the current values of out or set them to zero
  // (if c_out == 0): We are now done with c_out.
  detail::scale_or_set(out, c_out);

  // if c_this == 0 we are done
  if (c_this == Constants<scalar_type>::zero) return;

  // The number of rows in the result matrix
  const size_type rows = mode == Transposed::None ? this->n_rows() : this->n_cols();
  // Note: The number of columns is always in.n_cols()
  // and the contraction index always runs over in.n_rows() in all modes!

  // Note: We run over k slowest, since in and out are for sure stored
  // matrices
  // in row-major format, i.e. column indices should run fastest.
  for (size_type k = 0; k < in.n_rows(); ++k) {
    for (size_type row = 0; row < rows; ++row) {
      for (size_type col = 0; col < in.n_cols(); ++col) {
        switch (mode) {
          case Transposed::None:
            out(row, col) += c_this * (*this)(row, k) * in(k, col);
            break;
          case Transposed::Trans:
            out(row, col) += c_this * (*this)(k, row) * in(k, col);
            break;
          case Transposed::ConjTrans:
            // A variant of std::conj, which does not return a
            // complex data type if scalar is real only.
            detail::ConjFctr mconj;
            out(row, col) += c_this * mconj((*this)(k, row)) * in(k, col);
            break;
        }  // mode
      }    // col
    }      // row
  }        // k
}

}  // namespace linalgwrap
