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
#include <krims/RCPWrapper.hh>

namespace linalgwrap {

/** \brief Class for enabeling lazy evaluation of matrix operations
 *  onto another matrix class, which is not already lazy.
 *
 * \tparam StoredMatrix:  The type of the stored matrix to make lazy.
 */
template <typename StoredMatrix>
class LazyMatrixWrapper : public LazyMatrixExpression<StoredMatrix> {
  public:
    typedef LazyMatrixExpression<StoredMatrix> base_type;
    typedef typename base_type::stored_matrix_type stored_matrix_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::lazy_matrix_expression_ptr_type
          lazy_matrix_expression_ptr_type;

    //
    // Construction and destruction
    //
    /** \brief Constructor from shared pointer
     *
     * Make a copy of the RCP pointer and use that as the inner matrix
     * object. All modifications done via this interface of cause are also
     * noticed from the Pointer provided here.
     */
    explicit LazyMatrixWrapper(
          krims::RCPWrapper<const stored_matrix_type> inner)
          : m_inner{std::move(inner)} {}

    /** \brief Constructor from stored_matrix_type, taking ownership of
     * the passed object. */
    explicit LazyMatrixWrapper(stored_matrix_type&& inner)
          : m_inner{std::make_shared<const stored_matrix_type>(
                  std::move(inner))} {}

    /** \brief Constructor from stored_matrix_type not taking
     *         ownership
     *
     *  \note This will not store a copy, but just a reference
     *  to the passed object inside, hence the passed object
     *  has to live longer than this class.
     */
    explicit LazyMatrixWrapper(const stored_matrix_type& inner)
          : m_inner{krims::make_subscription(inner, "LazyMatrixWrapper")} {}

    //
    // Access to inner matrix
    //
    /** Const access to the inner matrix object */
    const stored_matrix_type& inner_matrix() const { return *m_inner; }

    //
    // Matrix_i interface
    //
    /** \brief Number of rows of the matrix */
    size_type n_rows() const override { return m_inner->n_rows(); }

    /** \brief Number of columns of the matrix  */
    size_type n_cols() const override { return m_inner->n_cols(); }

    /** \brief return an element of the matrix */
    scalar_type operator()(size_type row, size_type col) const override {
        return (*m_inner)(row, col);
    }

    //
    // LazyMatrixExpression interface
    //
    /** Are operation modes Transposed::Trans and Transposed::ConjTrans
     *  supported for this matrix type.
     **/
    bool has_transpose_operation_mode() const override {
        return m_inner->has_transpose_operation_mode();
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
          const scalar_type c_M = Constants<scalar_type>::zero) const override {
        m_inner->extract_block(M, start_row, start_col, mode, c_this, c_M);
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
    template <
          typename VectorIn, typename VectorOut,
          mat_vec_apply_enabled_t<LazyMatrixWrapper, VectorIn, VectorOut>...>
    void apply(const MultiVector<VectorIn>& x, MultiVector<VectorOut>& y,
               const Transposed mode = Transposed::None,
               const scalar_type c_this = Constants<scalar_type>::one,
               const scalar_type c_y = Constants<scalar_type>::zero) const {
        m_inner->apply(x, y, mode, c_this, c_y);
    }

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
          const scalar_type c_y = Constants<scalar_type>::zero) const override {
        m_inner->apply(x, y, mode, c_this, c_y);
    }

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
                     Constants<scalar_type>::zero) const override {
        m_inner->mmult(in, out, mode, c_this, c_out);
    }

    /** \brief Update the internal data
     *
     *  In this case does nothing.
     * */
    void update(const krims::ParameterMap&) override {
        // Do nothing.
    }

    /** \brief Clone the expression */
    lazy_matrix_expression_ptr_type clone() const override {
        // return a copy enwrapped in the pointer type
        return lazy_matrix_expression_ptr_type(new LazyMatrixWrapper(*this));
    }

  private:
    krims::RCPWrapper<const stored_matrix_type> m_inner;
};

}  // namespace linalgwrap
