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
// TODO Review and improve this class

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

    /** \name Constructors */
    ///@{
    /** Construct a block diagonal matrix from a bunch of matrices
     *
     * The matrices will be arranged as blocks along the diagonal.
     *
     * Matrices have to be square matrices!
     * */
    explicit BlockDiagonalMatrix(Matrices&&... blocks)
          : BlockDiagonalMatrix{std::tuple<Matrices&&...>{
                  std::forward<Matrices>(blocks)...}} {}

    /** Construct a block diagonal matrix from a tuple of matrices
     *
     * The matrices will be arranged as blocks along the diagonal.
     *
     * Matrices have to be square matrices!
     * */
    explicit BlockDiagonalMatrix(std::tuple<Matrices&&...> blocks);

    /** Construct a block diagonal matrix from a bunch of matrix
     *  references. The matrices will be copied, so be careful.
     *
     * The matrices will be arranged as blocks along the diagonal.
     *
     * Matrices have to be square matrices!
     * */
    explicit BlockDiagonalMatrix(const Matrices&... blocks)
          : BlockDiagonalMatrix{std::move(Matrices{blocks})...} {}
    ///@}

    /** \brief Const access to the matrix blocks */
    std::tuple<const Matrices&...> blocks() const override { return m_blocks; }

    //
    // Matrix_i interface
    //
    /** \brief Number of rows of the matrix */
    size_type n_rows() const override { return m_size; }

    /** \brief Number of columns of the matrix */
    size_type n_cols() const override {
        // By construction this matrix is quadratic:
        return n_rows();
    }

    /** \brief return an element of the matrix    */
    scalar_type operator()(size_type row, size_type col) const override;

    // TODO implement interface of StoredMatrix_i

    // TODO implement extract_block, mmult and apply

  private:
    /** The Matrix blocks along the diagonal */
    std::tuple<Matrices...> m_blocks;

    /** The cached total size of the Matrix */
    size_type m_size;

    /** \brief Map which maps the index past-the-end of a block to the
     * appropriate block.
     *
     * Using this map and the map::upper_bound function we can quickly find
     * the block which contains a certain element of which we only know
     * its index.
     *
     * E.g. consider a BlockDiagonalMatrix with the block sizes 4,5 and 7.
     * The content of this map would be
     * 4  -> block1  (has indices 0 to 3)
     * 9  -> block2  (has indices 4 to 8)
     * 16 -> block3  (has indices 9 to 15)
     */
    std::map<size_type, Matrix_i<scalar_type>*> m_start_index_to_block;
};

//@{
template <typename... Matrices>
BlockDiagonalMatrix<typename std::decay<Matrices>::type...> make_block_diagonal(
      Matrices&&... blocks) {
    return BlockDiagonalMatrix<typename std::decay<Matrices>::type...>(
          std::forward<Matrices>(blocks)...);
}

template <typename... Matrices>
BlockDiagonalMatrix<typename std::decay<Matrices>::type...> make_block_diagonal(
      std::tuple<Matrices...>&& blocks) {
    return BlockDiagonalMatrix<typename std::decay<Matrices>::type...>{
          std::move(blocks)};
}

//@}
///@}

//
// ------------------------------------------
//

template <typename... Matrices>
BlockDiagonalMatrix<Matrices...>::BlockDiagonalMatrix(
      std::tuple<Matrices&&...> blocks)
      : m_blocks{std::move(blocks)}, m_size{0} {
    auto build_sizemap = [&](Matrix_i<scalar_type>& mat) {
        assert_dbg(mat.n_rows() == mat.n_cols(), ExcMatrixNotSquare());

        m_size += mat.n_rows();
        m_start_index_to_block[m_size] = &mat;
    };
    tuple_for_each(build_sizemap, m_blocks);
}

template <typename... Matrices>
typename BlockDiagonalMatrix<Matrices...>::scalar_type
BlockDiagonalMatrix<Matrices...>::operator()(size_type row,
                                             size_type col) const {
    assert_greater(row, n_rows());
    assert_greater(col, n_cols());

    // Find the block which contains the row index:
    auto offset_blockptr = m_start_index_to_block.upper_bound(row);
    assert_dbg(offset_blockptr != m_start_index_to_block.end(),
               krims::ExcInternalError());

    // Get past-the-end index of the block and reference to the block
    const size_type end = offset_blockptr->first;
    const Matrix_i<scalar_type>& block = *(offset_blockptr->second);

    // Compute start index:
    const size_type start = end - block.n_rows();
    assert_dbg(start <= row, krims::ExcInternalError());
    assert_dbg(row < end, krims::ExcInternalError());

    if (col < start || col >= end) {
        // Col index is not in the same block, hence we are in an
        // off-diagonal zero block
        return Constants<scalar_type>::zero;
    } else {
        return block(row - start, col - start);
    }
}
}  // namespace linalgwrap
