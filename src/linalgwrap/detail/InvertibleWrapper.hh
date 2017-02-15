//
// Copyright (C) 2017 by the linalgwrap authors
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
#include "linalgwrap/TypeUtils.hh"
#include "linalgwrap/solve.hh"
#include "linalgwrap/trans.hh"
#include <krims/GenMap.hh>

namespace linalgwrap {
namespace detail {

// TODO See same comment as ProxyBase ... probably we want to use dynamic polymorphism
// here some day

/** Wrapper class which takes a matrix by const reference and makes a matrix object out of
 * it which can be inverted either directly by using an analytical inverse (as implemented
 * in apply_inverse) or implicitly using an iterative linear solver.
 *
 * The parameters for the linear solver are selected upon construction of this
 * class.
 */
template <typename Matrix>
class InvertibleWrapper
      : public LazyMatrixExpression<typename StoredTypeOf<Matrix>::type> {
 public:
  static_assert(IsMatrix<Matrix>::value, "Matrix needs to be a valid matrix type.");

  typedef LazyMatrixExpression<typename StoredTypeOf<Matrix>::type> base_type;
  typedef typename base_type::stored_matrix_type stored_matrix_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::lazy_matrix_expression_ptr_type
        lazy_matrix_expression_ptr_type;

  /* \brief Construct a wrapper object, which takes a const reference to an
   *        already existing matrix and a set of parameters, which are passed
   *        to the linear solver.
   */
  explicit InvertibleWrapper(const Matrix& inner, krims::GenMap params = krims::GenMap{})
        : m_inner_ptr("InvertibleWrapper", inner), m_params{std::move(params)} {
    assert_dbg(n_rows() == n_cols(), ExcMatrixNotSquare());

    // TODO Setup sensible preconditioner for the matrix ... we will probably
    //      use it a lot.
  }

  //
  // Matrix_i interface
  //
  /** \brief Number of rows of the matrix */
  size_type n_rows() const override { return m_inner_ptr->n_rows(); }

  /** \brief Number of columns of the matrix  */
  size_type n_cols() const override { return m_inner_ptr->n_cols(); }

  /** \brief return an element of the matrix */
  scalar_type operator()(size_type row, size_type col) const override {
    return (*m_inner_ptr)(row, col);
  }

  //
  // LazyMatrixExpression interface
  //
  /** Are operation modes Transposed::Trans and Transposed::ConjTrans
   *  supported for this matrix type.
   **/
  bool has_transpose_operation_mode() const override {
    return m_inner_ptr->has_transpose_operation_mode();
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
  void extract_block(stored_matrix_type& M, const size_type start_row,
                     const size_type start_col, const Transposed mode = Transposed::None,
                     const scalar_type c_this = 1,
                     const scalar_type c_M = 0) const override {
    m_inner_ptr->extract_block(M, start_row, start_col, mode, c_this, c_M);
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
            mat_vec_apply_enabled_t<InvertibleWrapper, VectorIn, VectorOut>...>
  void apply(const MultiVector<VectorIn>& x, MultiVector<VectorOut>& y,
             const Transposed mode,
             const scalar_type c_this = Constants<scalar_type>::one,
             const scalar_type c_y = Constants<scalar_type>::zero) const {
    m_inner_ptr->apply(x, y, mode, c_this, c_y);
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
    m_inner_ptr->apply(x, y, mode, c_this, c_y);
  }

  /** Perform a matrix-matrix product.
   *
   * Loosely performs the operation
   * \[ out = c_this \cdot A^\text{mode} \cdot in + c_out \cdot out. \]
   *
   * See LazyMatrixExpression for more details
   */
  void mmult(const stored_matrix_type& in, stored_matrix_type& out,
             const Transposed mode = Transposed::None, const scalar_type c_this = 1,
             const scalar_type c_out = 0) const override {
    m_inner_ptr->mmult(in, out, mode, c_this, c_out);
  }

  /** \brief Compute the Inverse-Multivector application
   *
   * Loosely speaking we perform
   * \[ y = c_this \cdot (A^{-1})^\text{mode} \cdot x + c_y \cdot y. \]
   *
   * See LazyMatrixExpression for more details
   */
  template <typename VectorIn, typename VectorOut,
            mat_vec_apply_enabled_t<InvertibleWrapper, VectorIn, VectorOut>...>
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
        const scalar_type c_y = 0) const override {

    if (m_inner_ptr->has_apply_inverse()) {
      // Call the analytical inverse of the inner matrix.
      return m_inner_ptr->apply_inverse(x, y, mode, c_this, c_y);
    }

    assert_dbg(mode == Transposed::None || has_transpose_operation_mode(),
               ExcUnsupportedOperationMode(mode));
    assert_finite(c_this);
    assert_finite(c_y);
    assert_size(x.n_vectors(), y.n_vectors());
    if (mode == Transposed::Trans || mode == Transposed::ConjTrans) {
      assert_size(x.n_elem(), n_rows());
      assert_size(y.n_elem(), n_cols());
    } else {
      assert_size(x.n_elem(), n_cols());
      assert_size(y.n_elem(), n_rows());
    }
    assert_sufficiently_tested(mode != Transposed::ConjTrans);

    if (c_this == Constants<scalar_type>::zero) {
      for (auto& vec : y) detail::scale_or_set(vec, c_y);
      return;
    }  // c_this == 0

    // Avoid copying y for later by explicitly disabling the case where we need to add a
    // fraction of y to the solution.
    assert_throw(c_y == 0,
                 krims::ExcDisabled("apply_inverse for the InvertibleWrapper is only "
                                    "reasonably efficient if c_y == 0."));

    // Solve the relevant linear system into y:
    const Matrix& A(*m_inner_ptr);
    switch (mode) {
      case Transposed::None:
        solve(A, y, x);  // Solve problem A y  = x
        break;
      case Transposed::Trans:
        solve(trans(A), y, x);  // Solve A^T y = x
        break;
      case Transposed::ConjTrans:
        solve(conjtrans(A), y, x);  // Solve A^H y = x
        break;
    }

    // We do not need to scale at all:
    if (c_this == 1) return;

    // Scale the result:
    for (auto& vec : y) detail::scale_or_set(vec, c_this);
  }

  /** \brief Clone the expression */
  lazy_matrix_expression_ptr_type clone() const override {
    // return a copy enwrapped in the pointer type
    return lazy_matrix_expression_ptr_type(new InvertibleWrapper(*this));
  }

  /** \brief Update the InvertibleWrapper */
  void update(const krims::ParameterMap& /* map */) override {
    // Cannot be done, since InvertibleWrapper is only
    // holding a const reference.
    assert_throw(false, krims::ExcDisabled("Cannot update an InvertibleWrapper."));
  }

 private:
  krims::SubscriptionPointer<const Matrix> m_inner_ptr;
  krims::GenMap m_params;
};

}  // namespace detail
}  // namespace linalgwrap
