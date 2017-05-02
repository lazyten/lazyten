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
#include "LazyMatrixExpression.hh"
#include "TypeUtils.hh"
#include "detail/GenericFunctionals.hh"
#include <array>

namespace linalgwrap {

template <typename Matrix, size_t N,
          typename Stored = typename StoredTypeOf<Matrix>::type>
class BlockDiagonalMatrix : public LazyMatrixExpression<Stored> {
 public:
  static_assert(IsLazyMatrix<Matrix>::value, "Matrix needs to be a lazy matrix");
  static_assert(N > 0, "We assume at least one block to be present");

  typedef LazyMatrixExpression<Stored> base_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::stored_matrix_type stored_matrix_type;
  typedef typename base_type::lazy_matrix_expression_ptr_type
        lazy_matrix_expression_ptr_type;

  typedef Matrix block_type;
  static constexpr size_t n_blocks = N;

  typedef std::array<Matrix, N> block_container_type;
  typedef typename block_container_type::iterator block_iterator_type;
  typedef typename block_container_type::const_iterator const_block_iterator_type;

  /** Construct from a bunch of square matrices */
  explicit BlockDiagonalMatrix(std::array<Matrix, N> diag_blocks);

  /** Access the array of matrix diagonal blocks */
  //@{
  const block_container_type& diag_blocks() const { return m_blocks; }
  block_container_type& diag_blocks() { return m_blocks; }
  //@}

  /** Return an iterator to the block, which contains this
   *  matrix element.
   *
   *  Checks that row and col are in the allowed range of indices.
   *
   *  Returns std::end(diag_blocks()) if no block contains this element,
   *  i.e. the element is not in any of the diagonal block.
   *  This implies the element is zero.
   **/
  const_block_iterator_type containing_block(size_t row, size_t col) const;

  //
  // Matrix_i interface
  //
  size_t n_rows() const override { return m_accu_sizes[N - 1]; }

  size_t n_cols() const override {
    return n_rows();  // Matrix is by construction quadratic
  }

  scalar_type operator()(size_t row, size_t col) const override;

  OperatorProperties properties() const override {
    // Merge the properties of all blocks
    // using the logical & (and) operator
    return std::accumulate(
          std::begin(m_blocks), std::end(m_blocks), m_blocks[0].properties(),
          [](OperatorProperties p, const block_type& b) { return p & b.properties(); });
  }

  //
  // Lazy matrix expression interface
  //
  bool has_transpose_operation_mode() const override {
    return std::all_of(std::begin(m_blocks), std::end(m_blocks),
                       [](const Matrix& m) { return m.has_transpose_operation_mode(); });
  }

  bool has_apply_inverse() const override {
    return std::all_of(std::begin(m_blocks), std::end(m_blocks),
                       [](const Matrix& m) { return m.has_apply_inverse(); });
  }

  void extract_block(stored_matrix_type& M, const size_t start_row,
                     const size_t start_col, const Transposed mode = Transposed::None,
                     const scalar_type c_this = 1,
                     const scalar_type c_M = 0) const override;

  template <typename VectorIn, typename VectorOut,
            mat_vec_apply_enabled_t<BlockDiagonalMatrix, VectorIn, VectorOut>...>
  void apply(const MultiVector<VectorIn>& x, MultiVector<VectorOut>& y,
             const Transposed mode = Transposed::None,
             const scalar_type c_this = Constants<scalar_type>::one,
             const scalar_type c_y = Constants<scalar_type>::zero) const {
    MultiVector<const MutableMemoryVector_i<scalar_type>> x_wrapped(x);
    MultiVector<MutableMemoryVector_i<scalar_type>> y_wrapped(y);
    apply(x_wrapped, y_wrapped, mode, c_this, c_y);
  }

  void apply(const MultiVector<const MutableMemoryVector_i<scalar_type>>& x,
             MultiVector<MutableMemoryVector_i<scalar_type>>& y,
             const Transposed mode = Transposed::None, const scalar_type c_this = 1,
             const scalar_type c_y = 0) const override;

  template <typename VectorIn, typename VectorOut,
            mat_vec_apply_enabled_t<BlockDiagonalMatrix, VectorIn, VectorOut>...>
  void apply_inverse(const MultiVector<VectorIn>& x, MultiVector<VectorOut>& y,
                     const Transposed mode = Transposed::None,
                     const scalar_type c_this = 1, const scalar_type c_y = 0) const {
    MultiVector<const MutableMemoryVector_i<scalar_type>> x_wrapped(x);
    MultiVector<MutableMemoryVector_i<scalar_type>> y_wrapped(y);
    apply_inverse(x_wrapped, y_wrapped, mode, c_this, c_y);
  }

  void apply_inverse(const MultiVector<const MutableMemoryVector_i<scalar_type>>& x,
                     MultiVector<MutableMemoryVector_i<scalar_type>>& y,
                     const Transposed mode = Transposed::None,
                     const scalar_type c_this = 1,
                     const scalar_type c_y = 0) const override;

  void mmult(const stored_matrix_type& in, stored_matrix_type& out,
             const Transposed mode = Transposed::None, const scalar_type c_this = 1,
             const scalar_type c_out = 0) const override;

  void update(const krims::GenMap& map) override {
    std::for_each(std::begin(m_blocks), std::end(m_blocks),
                  [&map](Matrix& m) { m.update(map); });
  }

  lazy_matrix_expression_ptr_type clone() const override {
    return lazy_matrix_expression_ptr_type(new BlockDiagonalMatrix(*this));
  }

 protected:
  /** Lookup the block, which contains this index
   * (Since the blocks are quadratic this works
   * equivalently for row and column indices)
   *
   * Throws an error if the index is out of bounds
   */
  const_block_iterator_type block_containing_index(size_t idx) const;

  /** Compute the range covered by a block */
  krims::Range<size_t> index_range_of(const_block_iterator_type block) const;

  /** Return a view into a multivector, which only looks at those elements of the
   *  vector, which are in the range of the current block */
  template <typename Vector>
  MultiVector<IteratorVector<scalar_type*>> block_view(
        const_block_iterator_type block, const MultiVector<Vector>& v) const;

 private:
  /** Contains the accumulated sizes,
   *  i.e. {blocks[0].size(), blocks[0].size() + blocks[1].size(), ...}
   */
  std::array<size_t, N> m_accu_sizes;

  /** Contains the actual blocks */
  std::array<block_type, N> m_blocks;
};

//@{
/** \brief helper struct to detect whether the scalar type
 *         underlying the matrix or matrix reference is complex.
 **/
template <typename Object, typename = void>
struct IsBlockDiagonalMatrix : public std::false_type {};

template <typename T>
struct IsBlockDiagonalMatrix<T, krims::VoidType<typename T::block_type>>
      : public std::is_base_of<BlockDiagonalMatrix<typename T::block_type, T::n_blocks>,
                               T> {};

//@}

/** Helper function to construct a BlockDiagonal matrix */
template <typename Matrix, size_t N>
BlockDiagonalMatrix<Matrix, N> make_block_diagonal(std::array<Matrix, N> ms) {
  return BlockDiagonalMatrix<Matrix, N>(std::move(ms));
}

/** \name BlockDiagonalMatrix operations */
///@{

// The following three macros ease the definition of operations with block diagonal
// matrices. They exploit that due to the block-diagonal structure, the operations
// can be done block by block.
//
// They each define a result array (assuming that the result type of the blockwise
// operation can be default constructed) and use the std::transform function
// to apply a blockwise operation (defined as a lambda) to each block or each pair
// of blocks. The resulting array is enwrapped as a new BlockDiagonalMatrix and returned.
#define BLOCK_DIAGONAL_BINARY_OP(op)                                                  \
  template <typename Mlhs, typename Mrhs, size_t N>                                   \
  auto operator op(const BlockDiagonalMatrix<Mlhs, N>& lhs,                           \
                   const BlockDiagonalMatrix<Mrhs, N>& rhs)                           \
        ->BlockDiagonalMatrix<decltype(lhs.diag_blocks()[0] op rhs.diag_blocks()[0]), \
                              N> {                                                    \
    typedef decltype(lhs.diag_blocks()[0] op rhs.diag_blocks()[0]) res_type;          \
    std::array<res_type, N> res;                                                      \
                                                                                      \
    auto fct = [](const Mlhs& lhs, const Mrhs& rhs) { return lhs op rhs; };           \
    std::transform(std::begin(lhs.diag_blocks()), std::end(lhs.diag_blocks()),        \
                   std::begin(rhs.diag_blocks()), std::begin(res), std::move(fct));   \
    return BlockDiagonalMatrix<res_type, N>{std::move(res)};                          \
  }

#define BLOCK_DIAGONAL_UNARY_OP(op)                                        \
  template <typename M, size_t N>                                          \
  auto operator op(const BlockDiagonalMatrix<M, N>& m)                     \
        ->BlockDiagonalMatrix<decltype(op m.diag_blocks()[0]), N> {        \
    typedef decltype(op m.diag_blocks()[0]) res_type;                      \
    std::array<res_type, N> res;                                           \
                                                                           \
    auto fct = [](const M& m) { return op m; };                            \
    std::transform(std::begin(m.diag_blocks()), std::end(m.diag_blocks()), \
                   std::begin(res), std::move(fct));                       \
    return BlockDiagonalMatrix<res_type, N>{std::move(res)};               \
  }

#define BLOCK_DIAGONAL_SCALAR_OP(op)                                              \
  template <typename M, size_t N>                                                 \
  auto operator op(const BlockDiagonalMatrix<M, N>& m, typename M::scalar_type s) \
        ->BlockDiagonalMatrix<decltype(m.diag_blocks()[0] op s), N> {             \
    typedef decltype(m.diag_blocks()[0] op s) res_type;                           \
    std::array<res_type, N> res;                                                  \
                                                                                  \
    auto fct = [s](const M& m) { return m op s; };                                \
    std::transform(std::begin(m.diag_blocks()), std::end(m.diag_blocks()),        \
                   std::begin(res), std::move(fct));                              \
    return BlockDiagonalMatrix<res_type, N>{std::move(res)};                      \
  }

//! Add two block diagonal matrices
BLOCK_DIAGONAL_BINARY_OP(+);

//! Subtract two block diagonal matrices
BLOCK_DIAGONAL_BINARY_OP(-);

//! Multiply two block diagonal matrices
BLOCK_DIAGONAL_BINARY_OP(*);

//! Unary minus on a block diagonal matrix
BLOCK_DIAGONAL_UNARY_OP(-);

//! Scale a block diagonal matrix
BLOCK_DIAGONAL_SCALAR_OP(*);

//! Divide block diagonal matrix elementwise by a scalar:
BLOCK_DIAGONAL_SCALAR_OP(/);

//! Scale a block diagonal matrix
template <typename M, size_t N>
auto operator*(typename M::scalar_type s, const BlockDiagonalMatrix<M, N>& m)
      -> decltype(m * s) {
  return m * s;
}

// Undo macro definition
#undef BLOCK_DIAGONAL_BINARY_OP
#undef BLOCK_DIAGONAL_UNARY_OP
#undef BLOCK_DIAGONAL_SCALAR_OP
///@}

//
// ---------------------------------------------------
//

template <typename Matrix, size_t N, typename Stored>
BlockDiagonalMatrix<Matrix, N, Stored>::BlockDiagonalMatrix(
      std::array<Matrix, N> diag_blocks)
      : m_accu_sizes{}, m_blocks(std::move(diag_blocks)) {
  size_t accu_size = 0;

  auto itsize = std::begin(m_accu_sizes);
  for (const Matrix& b : m_blocks) {
    assert_dbg(b.n_rows() == b.n_cols(), ExcMatrixNotSquare());
    accu_size += b.n_rows();
    *itsize = accu_size;
    ++itsize;
  }
}

template <typename Matrix, size_t N, typename Stored>
typename BlockDiagonalMatrix<Matrix, N, Stored>::const_block_iterator_type
BlockDiagonalMatrix<Matrix, N, Stored>::block_containing_index(size_t idx) const {
  assert_greater(idx, n_rows());

  // accu_sizes contains the accumulated sizes of the blocks.
  // So the upper bound of this gives the index of the block,
  // which contains this index
  const auto it_sizes =
        std::upper_bound(std::begin(m_accu_sizes), std::end(m_accu_sizes), idx);
  assert_dbg(it_sizes != std::end(m_accu_sizes), krims::ExcInternalError());

  return m_blocks.cbegin() + (it_sizes - std::begin(m_accu_sizes));
}

template <typename Matrix, size_t N, typename Stored>
krims::Range<size_t> BlockDiagonalMatrix<Matrix, N, Stored>::index_range_of(
      const_block_iterator_type block) const {
  const ptrdiff_t idx = block - m_blocks.cbegin();
  auto it = std::begin(m_accu_sizes);
  return krims::Range<size_t>{idx == 0 ? 0 : *(it + (idx - 1)), *(it + idx)};
}

template <typename Matrix, size_t N, typename Stored>
template <typename Vector>
MultiVector<IteratorVector<typename BlockDiagonalMatrix<Matrix, N, Stored>::scalar_type*>>
BlockDiagonalMatrix<Matrix, N, Stored>::block_view(const_block_iterator_type block,
                                                   const MultiVector<Vector>& v) const {
  // Make sure that there is no index overshooting in
  // the loop below:
  assert_equal(v.n_elem(), n_cols());

  krims::Range<size_t> range = index_range_of(block);
  MultiVector<IteratorVector<scalar_type*>> view;
  for (size_t vec = 0; vec < v.n_vectors(); ++vec) {
    auto begin = std::begin(v[vec]) + range.lower_bound();
    auto end = std::begin(v[vec]) + range.upper_bound();

    // TODO This is a little ugly
    scalar_type* pbegin = const_cast<scalar_type*>(begin);
    scalar_type* pend = const_cast<scalar_type*>(end);
    view.push_back(IteratorVector<scalar_type*>(pbegin, pend));
  }

  return view;
}

template <typename Matrix, size_t N, typename Stored>
typename BlockDiagonalMatrix<Matrix, N, Stored>::const_block_iterator_type
BlockDiagonalMatrix<Matrix, N, Stored>::containing_block(size_t row, size_t col) const {
  assert_greater(row, n_rows());
  assert_greater(col, n_cols());

  const auto row_block = block_containing_index(row);
  if (index_range_of(row_block).contains(col)) {
    return row_block;
  } else {
    return m_blocks.cend();
  }
}

template <typename Matrix, size_t N, typename Stored>
typename BlockDiagonalMatrix<Matrix, N, Stored>::scalar_type
BlockDiagonalMatrix<Matrix, N, Stored>::operator()(size_t row, size_t col) const {
  assert_greater(row, n_rows());
  assert_greater(col, n_cols());

  // Block which contains the row index and its index range
  const auto row_block = block_containing_index(row);
  const krims::Range<size_t> index_range = index_range_of(row_block);

  if (index_range.contains(col)) {
    size_t r_shifted = row - index_range.front();
    size_t c_shifted = col - index_range.front();
    return (*row_block)(r_shifted, c_shifted);
  } else {
    // Col index is not in the same block, hence we are in an
    // off-diagonal zero block
    return 0;
  }
}

template <typename Matrix, size_t N, typename Stored>
void BlockDiagonalMatrix<Matrix, N, Stored>::extract_block(
      stored_matrix_type& M, const size_t start_row, const size_t start_col,
      const Transposed mode, const scalar_type c_this, const scalar_type c_M) const {
  assert_dbg(mode == Transposed::None || has_transpose_operation_mode(),
             ExcUnsupportedOperationMode(mode));
  assert_finite(c_this);
  assert_finite(c_M);
  // check that we do not overshoot the indices
  if (mode == Transposed::Trans || mode == Transposed::ConjTrans) {
    assert_greater_equal(start_row + M.n_rows(), n_cols());
    assert_greater_equal(start_col + M.n_cols(), n_rows());
  } else {
    assert_greater_equal(start_row + M.n_rows(), n_rows());
    assert_greater_equal(start_col + M.n_cols(), n_cols());
  }

  // For empty matrices there is nothing to do
  if (M.n_rows() == 0 || M.n_cols() == 0) return;

  // Scale or set the full matrix => done with c_M
  detail::scale_or_set(M, c_M);

  // The last row and column index we care about
  // (i.e. *not* past the last, but the actual last)
  const size_t last_row = start_row + M.n_rows() - 1;
  const size_t last_col = start_col + M.n_cols() - 1;

  // Determine the blocks we need, i.e. the first containing both row and col range
  // up the one after the last containing both
  // (i.e. the begin and end we need to iterate over)
  const auto itbegin =
        std::min(block_containing_index(start_row), block_containing_index(start_col));
  const auto itend =
        1 + std::max(block_containing_index(last_row), block_containing_index(last_col));
  assert_dbg(itbegin < itend, krims::ExcInternalError());
  assert_dbg(itend <= std::end(m_blocks), krims::ExcInternalError());
  assert_dbg(itbegin < itend, krims::ExcInternalError());

  for (auto it = itbegin; it != itend; ++it) {
    // Ranges of actual indices for row and column to extract
    // using this block. We use min and max to assure that we never
    // overshoot the matrix.
    const krims::Range<size_t> block_range = index_range_of(it);
    const krims::Range<size_t> row_range{
          std::max(start_row, block_range.front()),
          std::min(last_row + 1, block_range.upper_bound())};
    const krims::Range<size_t> col_range{
          std::max(start_col, block_range.front()),
          std::min(last_col + 1, block_range.upper_bound())};
    if (row_range.empty() || col_range.empty()) continue;

    assert_dbg(block_range.front() <= row_range.front(), krims::ExcInternalError());
    assert_dbg(block_range.front() <= col_range.front(), krims::ExcInternalError());

    // Shift indices such that we are relative to the start of the current block
    const size_t block_row_start = row_range.front() - block_range.front();
    const size_t block_col_start = col_range.front() - block_range.front();

    stored_matrix_type extract(row_range.length(), col_range.length(), false);
    it->extract_block(extract, block_row_start, block_col_start, mode, c_this, 0);

    // TODO Would be nice to avoid the copy here!
    for (auto r : row_range) {
      for (auto c : col_range) {
        M(r - start_row, c - start_col) +=
              extract(r - row_range.front(), c - col_range.front());
      }
    }
  }
}

template <typename Matrix, size_t N, typename Stored>
void BlockDiagonalMatrix<Matrix, N, Stored>::apply(
      const MultiVector<const MutableMemoryVector_i<scalar_type>>& x,
      MultiVector<MutableMemoryVector_i<scalar_type>>& y, const Transposed mode,
      const scalar_type c_this, const scalar_type c_y) const {
  assert_finite(c_this);
  assert_finite(c_y);
  assert_size(x.n_vectors(), y.n_vectors());
  assert_size(x.n_elem(), n_rows());
  assert_size(y.n_elem(), n_cols());

  // Since we are block diagonal we can apply block by block:
  for (auto it = std::begin(m_blocks); it != std::end(m_blocks); ++it) {
    auto view_y = block_view(it, y);
    it->apply(block_view(it, x), view_y, mode, c_this, c_y);
  }
}

template <typename Matrix, size_t N, typename Stored>
void BlockDiagonalMatrix<Matrix, N, Stored>::apply_inverse(
      const MultiVector<const MutableMemoryVector_i<scalar_type>>& x,
      MultiVector<MutableMemoryVector_i<scalar_type>>& y, const Transposed mode,
      const scalar_type c_this, const scalar_type c_y) const {
  assert_throw(has_apply_inverse(),
               krims::ExcDisabled("The apply_inverse function is disabled, since "
                                  "has_apply_inverse() returns false."));
  assert_finite(c_this);
  assert_finite(c_y);
  assert_size(x.n_vectors(), y.n_vectors());
  assert_size(x.n_elem(), n_rows());
  assert_size(y.n_elem(), n_cols());

  // Since we are block diagonal we can apply block by block:
  for (auto it = std::begin(m_blocks); it != std::end(m_blocks); ++it) {
    auto view_y = block_view(it, y);
    it->apply_inverse(block_view(it, x), view_y, mode, c_this, c_y);
  }
}

template <typename Matrix, size_t N, typename Stored>
void BlockDiagonalMatrix<Matrix, N, Stored>::mmult(const stored_matrix_type& in,
                                                   stored_matrix_type& out,
                                                   const Transposed mode,
                                                   const scalar_type c_this,
                                                   const scalar_type c_out) const {
  assert_finite(c_this);
  assert_finite(c_out);
  assert_size(in.n_cols(), out.n_cols());
  assert_size(n_rows(), in.n_rows());
  assert_size(n_rows(), out.n_rows());  // since this is quadratic

  // This algorithm uses the fact that even though only this matrix
  // is block diagonal, the tiling of it can be applied to the input
  // and output matrices as well.
  //
  // In other words if we had 3 Blocks in A and multiply with a matrix B,
  // which we split up into 3 sheets with the same number of rows
  // as the individual blocks and with the number of columns equal to the
  // number of columns of B:
  //    A11  0   0              B1         A11.B1
  //     0  A12  0     \times   B2    =    A22.B2
  //     0   0  A13             B3         A33.B3

  if (c_this == Constants<scalar_type>::zero) {
    detail::scale_or_set(out, c_out);
    return;
  }

  for (auto ita = std::begin(m_blocks); ita != std::end(m_blocks); ++ita) {
    // Rows of the current sheets
    const krims::Range<size_t> rows = index_range_of(ita);

    // Number of columns in the sheets
    const size_t n_cols = in.n_cols();

    // TODO: This is bad copying. Use views instead.
    // Extract the parts, we care about in this iteration of the loop
    stored_matrix_type insheet(rows.length(), n_cols, false);
    in.extract_block(insheet, rows.front(), 0);

    stored_matrix_type outsheet(rows.length(), n_cols, false);
    out.extract_block(outsheet, rows.front(), 0);

    ita->mmult(insheet, outsheet, mode, c_this, c_out);

    // TODO This is bad copying
    // Copy back the results:
    for (size_t r : rows) {
      for (size_t c = 0; c < n_cols; ++c) {
        out(r, c) = outsheet(r - rows.front(), c);
      }  // c
    }    // r
  }      // ita
}

}  // namespace linalgwrap
