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
#include "VectorOf.hh"
#include <krims/SubscriptionPointer.hh>

namespace linalgwrap {
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
    typedef VectorOf<stored_matrix_type> small_vector_type;

    /** Construct from reference to diagonal elements */
    DiagonalMatrix(const small_vector_type& diagonal)
          : m_diagonal_ptr{"DiagonalMatrix", diagonal} {}

    /** Number of rows */
    size_type n_rows() const override { return m_diagonal_ptr->size(); }

    /** Number of columns */
    size_type n_cols() const override { return m_diagonal_ptr->size(); }

    /** Element access */
    scalar_type operator()(size_type row, size_type col) const override {
        if (row == col) return (*m_diagonal_ptr)[row];
        return Constants<scalar_type>::zero;
    }

    /** \brief Extract_block function
     *
     * For more details see LazyMatrixExpression.hh
     */
    stored_matrix_type extract_block(
          Range<size_type> row_range,
          Range<size_type> col_range) const override {
        // At least one range is empty -> no work to be done:
        if (row_range.empty() || col_range.empty()) {
            return stored_matrix_type{row_range.length(), col_range.length()};
        }
        // Assertive checks:
        assert_greater_equal(row_range.last(), n_rows());
        assert_greater_equal(col_range.last(), n_cols());

        stored_matrix_type res(row_range.length(), col_range.length());

        const size_type start = std::max(row_range.first(), col_range.first());
        const size_type end = std::min(row_range.last(), col_range.last());

        for (size_type i = start; i < end; ++i) {
            res(i - row_range.first(), i - col_range.first()) =
                  (*m_diagonal_ptr)[i];
        }
        return res;
    }

    /** \brief Add a block to in.
     *
     * See LazyMatrixExpression.hh for more details
     */
    void add_block_to(
          stored_matrix_type& in, size_type start_row, size_type start_col,
          scalar_type c_this = Constants<scalar_type>::one) const override {
        assert_greater(0, in.n_rows());
        assert_greater(0, in.n_cols());

        // check that we do not overshoot the indices
        assert_greater_equal(start_row + in.n_rows(), this->n_rows());
        assert_greater_equal(start_col + in.n_cols(), this->n_cols());

        const size_type start = std::max(start_row, start_col);
        const size_type end =
              std::min(start_row + in.n_rows(), start_col + in.n_cols());

        for (size_t i = start; i < end; ++i) {
            in(i - start_row, i - start_col) += c_this * (*m_diagonal_ptr)[i];
        }
    }

    /** \brief Update the internal data of all objects in this expression
     * given the ParameterMap
     * */
    void update(const krims::ParameterMap& map) override {
        const std::string diagonal_key = "diagonal";

        if (map.exists(diagonal_key)) {
            auto new_diagonal_ptr = map.at_ptr<small_vector_type>(diagonal_key);
            assert_size(m_diagonal_ptr->size(), new_diagonal_ptr->size());
            m_diagonal_ptr = new_diagonal_ptr;
        }
    }

    /** \brief Multiplication with a stored matrix */
    stored_matrix_type operator*(const stored_matrix_type& m) const override {
        assert_size(n_cols(), m.n_rows());

        stored_matrix_type res(n_rows(), m.n_cols(), true);
        for (auto it = std::begin(m); it != std::end(m); ++it) {
            res(it.row(), it.col()) = *it * (*m_diagonal_ptr)[it.row()];
        }

        return res;
    }

    /** Clone function */
    lazy_matrix_expression_ptr_type clone() const override {
        // return a copy enwrapped in the pointer type
        return lazy_matrix_expression_ptr_type(new DiagonalMatrix(*this));
    }

  private:
    krims::SubscriptionPointer<const small_vector_type> m_diagonal_ptr;
};

template <typename StoredVector>
DiagonalMatrix<typename StoredVector::matrix_type> make_diagonalmatrix(
      const StoredVector& v) {
    return DiagonalMatrix<typename StoredVector::matrix_type>(v);
}

}  // namespace linalgwrap
