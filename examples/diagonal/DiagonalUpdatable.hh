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
#include <initializer_list>
#include <linalgwrap/Constants.hh>
#include <linalgwrap/Exceptions.hh>
#include <linalgwrap/LazyMatrix_i.hh>
#include <linalgwrap/SmallVector.hh>

using namespace linalgwrap;

/** Matrix to represent a diagonal matrix by a vector.
 *
 * The vector can be changed via the update() function */
template <typename StoredMatrix>
class DiagonalUpdatable : public LazyMatrix_i<StoredMatrix> {
  public:
    typedef LazyMatrix_i<StoredMatrix> base_type;
    typedef typename base_type::stored_matrix_type stored_matrix_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::lazy_matrix_expression_ptr_type
          lazy_matrix_expression_ptr_type;

    /** Construct form diagonal vector */
    DiagonalUpdatable(SmallVector<scalar_type> diagonal)
          : m_diagonal(diagonal) {}

    /** Construct from initialiser list giving the diagonal */
    DiagonalUpdatable(std::initializer_list<scalar_type> diagonal)
          : m_diagonal(diagonal.size(), 1) {
        size_t i = 0;
        for (auto it = std::begin(diagonal); it != std::end(diagonal);
             ++it, ++i) {
            m_diagonal(i, 0) = *it;
        }
    }

    /** Number of rows */
    size_type n_rows() const override { return m_diagonal.size(); }

    /** Number of colums */
    size_type n_cols() const override { return m_diagonal.size(); }

    /** Element access */
    scalar_type operator()(size_type row, size_type col) const override {
        if (row == col) return m_diagonal[row];
        return Constants<scalar_type>::zero;
    }

    /** Clone function */
    lazy_matrix_expression_ptr_type clone() const override {
        // return a copy enwrapped in the pointer type
        return lazy_matrix_expression_ptr_type(new DiagonalUpdatable(*this));
    }

    /** \brief Overriding default extract_block function
     *
     * Not strictly speaking required to get a lazy matrix, but
     * improves preformance here.
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

        size_type start = std::max(row_range.first(), col_range.first());
        size_type end = std::min(row_range.last(), col_range.last());

        for (size_t i = start; i < end; ++i) {
            res(i - start, i - start) = m_diagonal[0];
        }
        return res;
    }

    /** \brief Overriding default multiplication function
     *
     * Not strictly speaking required to get a lazy matrix, but
     * improves preformance.
     */
    stored_matrix_type operator*(const stored_matrix_type& m) const override {
        assert_size(n_cols(), m.n_rows());

        stored_matrix_type res(n_rows(), m.n_cols(), true);
        for (auto it = std::begin(m); it != std::end(m); ++it) {
            res(it.row(), it.col()) = *it * m_diagonal[it.row()];
        }

        return res;
    }

    /** \brief Provide an update mechanism */
    void update(const ParameterMap& m) override {
        auto diagonal = m.at<SmallVector<scalar_type>>("diagonal");

        // check that we have no size changes
        assert_size(diagonal.size(), m_diagonal.size());

        // assign:
        m_diagonal = std::move(diagonal);
    }

  private:
    // A n_rows times 1 vector of diagonal shifts
    SmallVector<scalar_type> m_diagonal;
};
