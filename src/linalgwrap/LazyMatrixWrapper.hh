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

#ifndef LINALG_ABSTRACT_MATRIX_HPP_
#define LINALG_ABSTRACT_MATRIX_HPP_

#include "linalgwrap/LazyMatrix_i.hh"
#include <memory>

namespace linalgwrap {

/** \brief Class for enabeling lazy evaluation of matrix operations
 *  onto another matrix class, which is not already lazy.
 *
 * \tparam InnerMatrix:   The type of the inner matrix that is made lazy
 * \tparam StoredMatrix:  The type of the stored matrix to use
 */
template <typename StoredMatrix, typename InnerMatrix>
class LazyMatrixWrapper : public LazyMatrixExpression<StoredMatrix> {
    static_assert(
          std::is_same<typename InnerMatrix::scalar_type,
                       typename StoredMatrix::scalar_type>::value,
          "InnerMatrix and StoredMatrix need to have the same scalar type");

    static_assert(
          std::is_base_of<StoredMatrix_i<typename InnerMatrix::scalar_type>,
                          InnerMatrix>::value,
          "InnerMatrix must be a child class of StoredMatrix_i");

  public:
    typedef LazyMatrixExpression<StoredMatrix> base_type;
    typedef typename base_type::stored_matrix_type stored_matrix_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::lazy_matrix_expression_ptr_type
          lazy_matrix_expression_ptr_type;

    typedef InnerMatrix inner_matrix_type;

    //
    // Construction and destruction
    //
    /** \brief Constructor from shared pointer
     *
     * Make a copy of the shared pointer and use that as the inner matrix
     * object. All modifications done via this interface of cause are also
     * noticed from the Pointer provided here.
     */
    explicit LazyMatrixWrapper(std::shared_ptr<inner_matrix_type> inner)
          : m_inner{std::move(inner)} {}

    /** \brief Constructor from inner_matrix_type (for implicit conversion) */
    explicit LazyMatrixWrapper(inner_matrix_type&& inner)
          : m_inner{std::make_shared<inner_matrix_type>(inner)} {}

    /** \brief Constructor from inner_matrix_type (for explicit conversion)
     *
     * Warning: This copies the full content of inner.
     */
    explicit LazyMatrixWrapper(const inner_matrix_type& inner)
          : m_inner{std::make_shared<inner_matrix_type>(inner)} {}

    //
    // Access to inner matrix
    //
    /** Const access to the inner matrix object */
    const inner_matrix_type& inner_matrix() const { return *m_inner; }

    /** Non-const access to the inner data object */
    inner_matrix_type& inner_matrix() { return *m_inner; }

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
    /** \brief Extract a block of values out of the matrix and
     *         return it as a stored matrix of the appropriate size
     *
     * For more details of the interface see the function of the same
     * name in ``LazyMatrixExpression``.
     *
     * \param row_range   The Range object representing the range of rows
     *                    to extract. Note that it is a half-open interval
     *                    i.e. the LHS is inclusive, but the RHS not.
     *                    The Range may not be empty.
     * \param col_range   The Range object representing the range of
     *                    columns to extract.
     *                    The Range may not be empty.
     */
    stored_matrix_type extract_block(
          Range<size_type> row_range,
          Range<size_type> col_range) const override {
        return m_inner->extract_block(row_range, col_range);
    }

    /** \brief Add a block of values of the matrix to the stored matrix
     *         provided by reference.
     *
     * For more details of the interface see the function of the same name
     * in ``LazyMatrixExpression``.
     *
     *  \param in   Matrix to add the values to. It is assumed that it
     *              already has the correct sparsity structure to take
     *              all the values. Its size defines the size of the
     *              block. May not be an empty matrix.
     *  \param start_row  The row index of the first element to extract
     *  \param start_col  The column index of the first element to extract
     *  \param c_this     The coefficient to multiply this matrix with
     *                    before extracting.
     */
    void add_block_to(
          stored_matrix_type& in, size_type start_row, size_type start_col,
          scalar_type c_this = Constants<scalar_type>::one) const override {
        m_inner->add_block_to(in, start_row, start_col, c_this);
    }

    /** \brief Update the internal data
     *
     *  In this case does nothing.
     * */
    void update(const ParameterMap&) override {
        // Do nothing.
    }

    /** \brief Multiplication with a stored matrix */
    stored_matrix_type operator*(const stored_matrix_type& m) const override {
        // TODO assume that stored_matrix_type and inner_matrix_type
        //      can be multiplied together.
        return (*m_inner) * m;
    }

    /** \brief Print the expression tree to this outstream
     * */
    void print_tree(std::ostream& o) const override {
        // We are the leaf, just print the name of inner:
        o << m_inner->name();
    }

    /** \brief Clone the expression */
    lazy_matrix_expression_ptr_type clone() const override {
        // return a copy enwrapped in the pointer type
        return lazy_matrix_expression_ptr_type(new LazyMatrixWrapper(*this));
    }

  private:
    std::shared_ptr<inner_matrix_type> m_inner;
};

/** \brief Convenience function to construct a LazyMatrixWrapper */
template <typename StoredMatrix, typename InnerMatrix, typename... Args>
LazyMatrixWrapper<StoredMatrix, InnerMatrix> make_lazy_matrix(Args&&... args) {
    return LazyMatrixWrapper<StoredMatrix, InnerMatrix>(
          std::move(std::make_shared<InnerMatrix>(args...)));
}

}  // namespace linalg

#endif
