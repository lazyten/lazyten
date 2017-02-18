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
#include "LazyMatrixExpression.hh"
#include "detail/ProxyBase.hh"
#include "solve.hh"

namespace linalgwrap {

DefExceptionMsg(ExcMatrixHasNoInverse,
                "The matrix you passed to InverseProxy does not support an "
                "apply_inverse (has_apply_inverse() returned false).");

/** Class which represents an inverse of a matrix.
 *
 * Only supported for matrices with the apply_inverse function implemented.
 */
template <typename Matrix>
class InverseProxy : public detail::ProxyBase<Matrix> {
 public:
  typedef detail::ProxyBase<Matrix> base_type;
  typedef typename base_type::stored_matrix_type stored_matrix_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::lazy_matrix_expression_ptr_type
        lazy_matrix_expression_ptr_type;

  /* \brief Construct an inverse from an inner matrix, not taking ownership
   * of the passed object. */
  explicit InverseProxy(matrix_type& inner) : base_type{inner} {
    assert_dbg(base_type::inner_matrix().has_apply_inverse(), ExcMatrixHasNoInverse());
  }

  /* \brief Construct an inverse from an inner matrix, taking ownership
   * of the passed object. */
  explicit InverseProxy(matrix_type&& inner) : base_type{std::move(inner)} {
    assert_dbg(base_type::inner_matrix().has_apply_inverse(), ExcMatrixHasNoInverse());
  }

  //
  // Matrix_i interface
  //
  /** \brief Number of rows of the matrix */
  size_type n_rows() const override { return base_type::inner_matrix().n_cols(); }

  /** \brief Number of columns of the matrix  */
  size_type n_cols() const override { return base_type::inner_matrix().n_rows(); }

  /** \brief return an element of the matrix */
  scalar_type operator()(size_type /*row*/, size_type /*col*/) const override {
    assert_throw(false, krims::ExcDisabled("Obtaining elements from InverseProxy is "
                                           "extremely constly and hence disabled."));
    return 0.;
  }

  //
  // LazyMatrixExpression interface
  //
  /** Are operation modes Transposed::Trans and Transposed::ConjTrans
   *  supported for this matrix type.
   **/
  bool has_transpose_operation_mode() const override {
    return base_type::inner_matrix().has_transpose_operation_mode();
  }

  bool has_apply_inverse() const override {
    return true;  // By construction
  }

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
  void extract_block(stored_matrix_type& /*M*/, const size_type /*start_row*/,
                     const size_type /*start_col*/,
                     const Transposed /*mode = Transposed::None*/,
                     const scalar_type /*c_this = 1*/,
                     const scalar_type /*c_M = 0*/) const override {
    assert_throw(false, krims::ExcDisabled("Obtaning elements from InverseProxy is "
                                           "extremely constly and hence disabled."));
  }

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
            mat_vec_apply_enabled_t<InverseProxy, VectorIn, VectorOut>...>
  void apply(const MultiVector<VectorIn>& x, MultiVector<VectorOut>& y,
             const Transposed mode,
             const scalar_type c_this = Constants<scalar_type>::one,
             const scalar_type c_y = Constants<scalar_type>::zero) const {
    base_type::inner_matrix().apply_inverse(x, y, mode, c_this, c_y);
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
             const scalar_type c_y = Constants<scalar_type>::zero) const override {
    base_type::inner_matrix().apply_inverse(x, y, mode, c_this, c_y);
  }

  /** Perform a matrix-matrix product.
   *
   * Loosely performs the operation
   * \[ out = c_this \cdot A^\text{mode} \cdot in + c_out \cdot out. \]
   *
   * See LazyMatrixExpression for more details
   */
  void mmult(const stored_matrix_type& /*in*/, stored_matrix_type& /*out*/,
             const Transposed /*mode = Transposed::None*/,
             const scalar_type /*c_this = 1*/,
             const scalar_type /*c_out = 0*/) const override {
    assert_throw(false, krims::ExcDisabled(
                              "Matrix-matrix multiplication is disabled for Inverses"));
  }

  /** \brief Compute the Inverse-Multivector application
   *
   * Loosely speaking we perform
   * \[ y = c_this \cdot (A^{-1})^\text{mode} \cdot x + c_y \cdot y. \]
   *
   * See LazyMatrixExpression for more details
   */
  template <typename VectorIn, typename VectorOut,
            mat_vec_apply_enabled_t<InverseProxy, VectorIn, VectorOut>...>
  void apply_inverse(const MultiVector<VectorIn>& x, MultiVector<VectorOut>& y,
                     const Transposed mode = Transposed::None,
                     const scalar_type c_this = 1, const scalar_type c_y = 0) const {
    base_type::inner_matrix().apply(x, y, mode, c_this, c_y);
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
        const scalar_type c_y = 0) const override {
    base_type::inner_matrix().apply(x, y, mode, c_this, c_y);
  }

  /** \brief Clone the expression */
  lazy_matrix_expression_ptr_type clone() const override {
    // return a copy enwrapped in the pointer type
    return lazy_matrix_expression_ptr_type(new InverseProxy(*this));
  }
};

template <typename Matrix>
class InverseProxy<InverseProxy<Matrix>> {
  static_assert(IsMatrix<Matrix>::value && false,
                "InverseProxy<InverseProxy<Matrix>> has been disabled "
                "since it makes no sense to have it.");
};

}  // namespace linalgwrap
