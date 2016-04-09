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
#include "linalgwrap/BlockDiagonalMatrixBase.hh"
#include <algorithm>
#include <array>
#include <utility>

namespace linalgwrap {

// Forward definition of base class:
template <typename... Matrices>
class BlockDiagonalMatrixBase;

/** Class for representing block diagonal matrices. */
template <typename... Matrices>
class BlockDiagonalMatrix : public BlockDiagonalMatrixBase<Matrices...> {
  public:
    typedef BlockDiagonalMatrixBase<Matrices...> base_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;

    static_assert(sizeof...(Matrices) > 0,
                  "Empty BlockDiagonalMatrix objects with zero matrices along "
                  "the diagonal are not permitted.");

    /** Construct a block diagonal matrix from a tuple of matrices
     *
     * The matrices will be arranged as blocks along the diagonal.
     *
     * Matrices have to be square matrices!
     * */
    explicit BlockDiagonalMatrix(std::tuple<Matrices...> blocks);

    /** Construct a block diagonal matrix from a bunch of matrices
     *
     * The matrices will be arranged as blocks along the diagonal.
     *
     * Matrices have to be square matrices!
     * */
    explicit BlockDiagonalMatrix(Matrices... blocks);

    /** \brief Const access to the matrix blocks */
    std::tuple<const Matrices&...> blocks() const override;

    //
    // Matrix_i interface
    //
    /** \brief Number of rows of the matrix */
    size_type n_rows() const override;

    /** \brief Number of columns of the matrix */
    size_type n_cols() const override;

    /** \brief return an element of the matrix    */
    scalar_type operator()(size_type row, size_type col) const override;

    // TODO implement interface of StoredMatrix_i

  private:
    /** The Matrix blocks along the diagonal */
    std::tuple<Matrices...> m_blocks;

    /** Cache for the size of the matrix: */
    size_type m_size;
};

/** \name make_block_diagonal helper functions */
///@{
/** \brief Construct a block diagonal matrix from a tuple of matrices
 *  which will be arranged along the diagonal.
 */
template <typename... Matrices>
BlockDiagonalMatrix<Matrices...> make_block_diagonal(
      std::tuple<Matrices...>&& blocks);

/** \brief Construct a block diagonal matrix from matrices
 *  which will be arranged along the diagonal.
 */
template <typename... Matrices>
BlockDiagonalMatrix<Matrices...> make_block_diagonal(Matrices&&... blocks);
///@}

//
// ------------------------------------------
//

template <typename... Matrices>
BlockDiagonalMatrix<Matrices...>::BlockDiagonalMatrix(
      std::tuple<Matrices...> blocks)
      : m_blocks{std::move(blocks)}, m_size{0} {
    // Assert that all matrix blocks are quadratic.
    tuple_for_each([](auto& mat) { assert_size(mat.n_rows(), mat.n_cols()); },
                   tuple_ref(m_blocks));

    // Determine size of this matrix:
    // Lambda to add size of current block in:
    auto add_size = [&](auto& mat) { m_size += mat.n_rows(); };

    // Run over all blocks:
    tuple_for_each(add_size, tuple_ref(m_blocks));
}

template <typename... Matrices>
BlockDiagonalMatrix<Matrices...>::BlockDiagonalMatrix(Matrices... blocks)
      : BlockDiagonalMatrix<Matrices...>{std::forward_as_tuple(blocks...)} {}

template <typename... Matrices>
typename BlockDiagonalMatrix<Matrices...>::size_type
BlockDiagonalMatrix<Matrices...>::n_rows() const {
    return m_size;
}

template <typename... Matrices>
typename BlockDiagonalMatrix<Matrices...>::size_type
BlockDiagonalMatrix<Matrices...>::n_cols() const {
    // By construction this matrix is quadratic:
    return n_rows();
}

template <typename... Matrices>
typename BlockDiagonalMatrix<Matrices...>::scalar_type
BlockDiagonalMatrix<Matrices...>::operator()(size_type row,
                                             size_type col) const {
    assert_greater(row, n_rows());
    assert_greater(col, n_cols());

    // In this function we iterate over the tuple of matrix blocks,
    // i.e. down the block diagonal.

    scalar_type res = Constants<scalar_type>::zero;
    size_type cur_block_start_index = 0;
    size_type cur_block_end_index = 0;

    // Predicate that asserts that the current matrix block contains
    // the element we are interested in.
    auto block_of_interest = [&](auto& mat) {
        // current start index is the old end index;
        cur_block_start_index = cur_block_end_index;

        // Update the end index: Matrices are by assumption quadratic.
        cur_block_end_index += mat.n_rows();

        // Row index is in the appropriate range:
        // => This block either contains the element,
        //    or the result is zero
        const bool row_in_range =
              cur_block_start_index <= row && row < cur_block_end_index;
        return row_in_range;
    };

    // End recursion (since we passed the block of interest)
    auto extract_value = [&](auto& mat) {
        const bool col_in_range =
              cur_block_start_index <= col && col < cur_block_end_index;

        if (col_in_range) {
            res = mat(row - cur_block_start_index, col - cur_block_start_index);
        }
    };

    // Extract the value res from the first block of interest:
    tuple_for_first(block_of_interest, extract_value, blocks());

    return res;
}

template <typename... Matrices>
inline std::tuple<const Matrices&...> BlockDiagonalMatrix<Matrices...>::blocks()
      const {
    return m_blocks;
}

//
// Out of place functions
//

template <typename... Matrices>
inline BlockDiagonalMatrix<Matrices...> make_block_diagonal(
      std::tuple<Matrices...>&& blocks) {
    return BlockDiagonalMatrix<Matrices...>(std::move(blocks));
}

template <typename... Matrices>
inline BlockDiagonalMatrix<Matrices...> make_block_diagonal(
      Matrices&&... blocks) {
    return BlockDiagonalMatrix<Matrices...>{std::move(blocks)...};
}
}  // namespace linalgwrap
