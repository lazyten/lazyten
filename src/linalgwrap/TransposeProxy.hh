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
#include "LazyMatrixExpression.hh"
#include "TypeUtils.hh"
#include <krims/RCPWrapper.hh>

namespace linalgwrap {

DefExceptionMsg(ExcMatrixHasNoTransposeOperations,
                "The matrix you passed to TransposeProxy does not support a "
                "transpose operation mode (has_transpose_operation_mode() "
                "returned false). Hence it cannot be used with "
                "TransposeProxy.");

namespace detail {

/** TransposeProxyBase class for stored matrices */
template <typename Matrix,
          typename StoredMatrix = typename StoredTypeOf<Matrix>::type,
          bool isStoredMatrix = IsStoredMatrix<Matrix>::value,
          bool isReadonly = std::is_const<Matrix>::value>
class TransposeProxyBase : public LazyMatrixExpression<StoredMatrix> {
    static_assert(IsMatrix<Matrix>::value, "Matrix needs to be a matrix type");

  public:
    //! Matrix type underlying this TransposeProxyBase class
    typedef Matrix matrix_type;
    typedef StoredMatrix stored_matrix_type;

    /** \brief Update the internal data
     *
     *  In this case does nothing.
     * */
    void update(const krims::ParameterMap&) override {}

    /** Does this class own the inner matrix
     *
     * In other words are we the one who is responsible for
     * managing the inner storage or not.V
     */
    bool owns_inner_matrix() const { return m_inner_ptr.is_shared_ptr(); }

    /** Access to the inner object */
    Matrix& inner_matrix() { return *m_inner_ptr; }

    /** Const access to the inner object */
    const Matrix& inner_matrix() const { return *m_inner_ptr; }

  protected:
    TransposeProxyBase(Matrix& matrix)
          : m_inner_ptr(krims::make_subscription(matrix, "TransposeProxy")) {}

    TransposeProxyBase(Matrix&& matrix)
          : m_inner_ptr(std::make_shared<Matrix>(matrix)) {}

    krims::RCPWrapper<Matrix> m_inner_ptr;
};

/** TransposeProxyBase class for const lazy matrix classes */
template <typename Matrix, typename StoredMatrix>
class TransposeProxyBase<Matrix, StoredMatrix, /*isStored*/ false,
                         /*isReadonly*/ true>
      : public LazyMatrixExpression<StoredMatrix> {
    static_assert(IsMatrix<Matrix>::value, "Matrix needs to be a matrix type");

  public:
    //! Matrix type underlying this TransposeProxyBase class
    typedef Matrix matrix_type;
    typedef StoredMatrix stored_matrix_type;

    /** \brief Update the internal data
     *
     *  In this case does nothing, since the internal object is
     *  either a stored matrix or const.
     * */
    void update(const krims::ParameterMap&) override {
        assert_dbg(false,
                   krims::ExcDisabled("Update is not possible for this Matrix "
                                      "expression, since the matrix inside the "
                                      "TransposeProxy object is const."));
    }

    /** Does this class own the inner matrix
     *
     * In other words are we the one who is responsible for
     * managing the inner storage or not.V
     *
     * In this case it always owns it since it is a lazy matrix object.
     */
    constexpr bool owns_inner_matrix() const { return true; }

    /** Access to the inner object */
    Matrix& inner_matrix() { return m_inner; }

    /** Const access to the inner object */
    const Matrix& inner_matrix() const { return m_inner; }

  protected:
    TransposeProxyBase(const Matrix& matrix) : m_inner(matrix) {}
    TransposeProxyBase(Matrix&& matrix) : m_inner(std::move(matrix)) {}

    /** Matrix is lazy, so we can store a full copy */
    Matrix m_inner;
};

/** TransposeProxyBase class for mutable lazy matrix classes */
template <typename Matrix, typename StoredMatrix>
class TransposeProxyBase<Matrix, StoredMatrix, /*isStored*/ false,
                         /*isReadonly*/ false>
      : public LazyMatrixExpression<StoredMatrix> {
    static_assert(IsMatrix<Matrix>::value, "Matrix needs to be a matrix type");

  public:
    //! Matrix type underlying this TransposeProxyBase class
    typedef Matrix matrix_type;
    typedef StoredMatrix stored_matrix_type;

    /** \brief Update the internal data
     *
     *  In this case does nothing, since the internal object is
     *  either a stored matrix or const.
     * */
    void update(const krims::ParameterMap& map) override {
        m_inner.update(map);
    }

    /** Does this class own the inner matrix
     *
     * In other words are we the one who is responsible for
     * managing the inner storage or not.V
     *
     * In this case it always owns it since it is a lazy matrix object.
     */
    constexpr bool owns_inner_matrix() const { return true; }

    /** Access to the inner object */
    Matrix& inner_matrix() { return m_inner; }

    /** Const access to the inner object */
    const Matrix& inner_matrix() const { return m_inner; }

  protected:
    TransposeProxyBase(const Matrix& matrix) : m_inner(matrix) {}
    TransposeProxyBase(Matrix&& matrix) : m_inner(std::move(matrix)) {}

    /** Matrix is lazy, so we can store a full copy */
    Matrix m_inner;
};
}  // namespace detail

template <typename Matrix>
class TransposeProxy : public detail::TransposeProxyBase<Matrix> {
  public:
    typedef detail::TransposeProxyBase<Matrix> base_type;
    typedef typename base_type::stored_matrix_type stored_matrix_type;
    typedef typename base_type::matrix_type matrix_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::lazy_matrix_expression_ptr_type
          lazy_matrix_expression_ptr_type;

    /* \brief Construct from an inner matrix, not taking ownership of the passed
     * object. */
    explicit TransposeProxy(matrix_type& inner) : base_type{inner} {
        assert_dbg(base_type::inner_matrix().has_transpose_operation_mode(),
                   ExcMatrixHasNoTransposeOperations());
    }

    /** \brief Constructor from an inner matrix, taking ownership of the
     * passed object. */
    explicit TransposeProxy(matrix_type&& inner) : base_type{std::move(inner)} {
        assert_dbg(base_type::inner_matrix().has_transpose_operation_mode(),
                   ExcMatrixHasNoTransposeOperations());
    }

    //
    // Matrix_i interface
    //
    /** \brief Number of rows of the matrix */
    size_type n_rows() const override {
        return base_type::inner_matrix().n_cols();
    }

    /** \brief Number of columns of the matrix  */
    size_type n_cols() const override {
        return base_type::inner_matrix().n_rows();
    }

    /** \brief return an element of the matrix */
    scalar_type operator()(size_type row, size_type col) const override {
        return base_type::inner_matrix()(col, row);
    }

    //
    // LazyMatrixExpression interface
    //
    /** Are operation modes Transposed::Trans and Transposed::ConjTrans
     *  supported for this matrix type.
     **/
    bool has_transpose_operation_mode() const override {
        return true;  // By construction of this class
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
    void extract_block(
          stored_matrix_type& M, const size_type start_row,
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
              mat_vec_apply_enabled_t<TransposeProxy, VectorIn, VectorOut>...>
    void apply(const MultiVector<VectorIn>& x, MultiVector<VectorOut>& y,
               const Transposed mode,
               const scalar_type c_this = Constants<scalar_type>::one,
               const scalar_type c_y = Constants<scalar_type>::zero) const;

    /** \brief Compute the Matrix-Multivector application
     *
     * Loosely speaking we perform
     * \[ y = c_this \cdot A^\text{mode} \cdot x + c_y \cdot y. \]
     *
     * See LazyMatrixExpression for more details
     */
    void apply(
          const MultiVector<const MutableMemoryVector_i<scalar_type>>& x,
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
     */
    void mmult(const stored_matrix_type& in, stored_matrix_type& out,
               const Transposed mode = Transposed::None,
               const scalar_type c_this = Constants<scalar_type>::one,
               const scalar_type c_out =
                     Constants<scalar_type>::zero) const override;

    /** \brief Clone the expression */
    lazy_matrix_expression_ptr_type clone() const override {
        // return a copy enwrapped in the pointer type
        return lazy_matrix_expression_ptr_type(new TransposeProxy(*this));
    }
};

template <typename Matrix>
class TransposeProxy<TransposeProxy<Matrix>> {
    static_assert(IsMatrix<Matrix>::value && false,
                  "TransposeProxy<TransposeProxy<Matrix>> has been disabled "
                  "since it makes no sense to have it.");
};

//
// ---------------------------------------------------------
//

template <typename Matrix>
void TransposeProxy<Matrix>::extract_block(stored_matrix_type& M,
                                           const size_type start_row,
                                           const size_type start_col,
                                           const Transposed mode,
                                           const scalar_type c_this,
                                           const scalar_type c_M) const {
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

    switch (mode) {
        case Transposed::None:
            base_type::inner_matrix().extract_block(
                  M, start_row, start_col, Transposed::Trans, c_this, c_M);
            break;
        case Transposed::Trans:
            base_type::inner_matrix().extract_block(
                  M, start_row, start_col, Transposed::None, c_this, c_M);
            break;
        case Transposed::ConjTrans:
            // TODO Implement
            assert_dbg(false, krims::ExcNotImplemented());
            break;
    }  // mode
}

template <typename Matrix>
template <
      typename VectorIn, typename VectorOut,
      mat_vec_apply_enabled_t<TransposeProxy<Matrix>, VectorIn, VectorOut>...>
void TransposeProxy<Matrix>::apply(const MultiVector<VectorIn>& x,
                                   MultiVector<VectorOut>& y,
                                   const Transposed mode,
                                   const scalar_type c_this,
                                   const scalar_type c_y) const {
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

    switch (mode) {
        case Transposed::None:
            base_type::inner_matrix().apply(x, y, Transposed::Trans, c_this,
                                            c_y);
            break;
        case Transposed::Trans:
            base_type::inner_matrix().apply(x, y, Transposed::None, c_this,
                                            c_y);
            break;
        case Transposed::ConjTrans:
            // TODO Implement
            assert_dbg(false, krims::ExcNotImplemented());
            break;
    }  // mode
}

template <typename Matrix>
void TransposeProxy<Matrix>::apply(
      const MultiVector<const MutableMemoryVector_i<scalar_type>>& x,
      MultiVector<MutableMemoryVector_i<scalar_type>>& y, const Transposed mode,
      const scalar_type c_this, const scalar_type c_y) const {
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

    switch (mode) {
        case Transposed::None:
            base_type::inner_matrix().apply(x, y, Transposed::Trans, c_this,
                                            c_y);
            break;
        case Transposed::Trans:
            base_type::inner_matrix().apply(x, y, Transposed::None, c_this,
                                            c_y);
            break;
        case Transposed::ConjTrans:
            // TODO Implement
            assert_dbg(false, krims::ExcNotImplemented());
            break;
    }  // mode
}

template <typename Matrix>
void TransposeProxy<Matrix>::mmult(const stored_matrix_type& in,
                                   stored_matrix_type& out,
                                   const Transposed mode,
                                   const scalar_type c_this,
                                   const scalar_type c_out) const {
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

    switch (mode) {
        case Transposed::None:
            base_type::inner_matrix().mmult(in, out, Transposed::Trans, c_this,
                                            c_out);
            break;
        case Transposed::Trans:
            base_type::inner_matrix().mmult(in, out, Transposed::None, c_this,
                                            c_out);
            break;
        case Transposed::ConjTrans:
            // TODO Implement
            assert_dbg(false, krims::ExcNotImplemented());
            break;
    }  // mode
}

}  // namespace linalgwrap
