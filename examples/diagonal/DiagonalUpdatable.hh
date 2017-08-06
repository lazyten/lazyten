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
#include <initializer_list>
#include <krims/GenMap.hh>
#include <lazyten/Constants.hh>
#include <lazyten/LazyMatrix_i.hh>
#include <lazyten/SmallVector.hh>

using namespace lazyten;

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
  DiagonalUpdatable(SmallVector<scalar_type> diagonal) : m_diagonal(diagonal) {}

  /** Construct from initialiser list giving the diagonal */
  DiagonalUpdatable(std::initializer_list<scalar_type> diagonal)
        : m_diagonal(diagonal.size(), 1) {
    size_t i = 0;
    for (auto it = std::begin(diagonal); it != std::end(diagonal); ++it, ++i) {
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

  /** \brief Provide an update mechanism */
  void update(const krims::GenMap& m) override {
    const auto& diagonal = m.at<SmallVector<scalar_type>>("diagonal");

    // check that we have no size changes
    assert_size(diagonal.size(), m_diagonal.size());

    // assign:
    m_diagonal = diagonal;
  }

  // Note: One could override apply, mmult and extract_block for better
  // performance in those matrix operations, but this is not required
  // to get the basic functionality going.

 private:
  // A n_rows times 1 vector of diagonal shifts
  SmallVector<scalar_type> m_diagonal;
};
