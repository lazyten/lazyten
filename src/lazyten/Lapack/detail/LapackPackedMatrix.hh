//
// Copyright (C) 2017 by the lazyten authors
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
#include "lazyten/config.hh"
#ifdef LAZYTEN_HAVE_LAPACK

#include "LapackSymmetricMatrix.hh"
#include "lazyten/LazyMatrixExpression.hh"
#include "lazyten/StoredMatrix_i.hh"
#include <vector>

namespace lazyten {
namespace detail {

/** Data structure to represent a packed matrix
 *
 * From the point of view of Fortran this *always* represents
 * a lower triangle of the symmetric matrix provided
 * upon construction.
 */
template <typename Scalar>
struct LapackPackedMatrix {
  typedef Scalar scalar_type;

  //! The number of rows and the number of columns
  size_t n;

  std::vector<Scalar> elements;

  /** Default construct */
  LapackPackedMatrix() = default;

  /** Construct from a usual symmetric lazyten stored matrix by
   * copying the values in */
  explicit LapackPackedMatrix(const StoredMatrix_i<Scalar>& m);

  /** Construct from a usual symmetric lazyten matrix expression by
   * copying the values in */
  template <typename Stored, typename = krims::enable_if_t<IsStoredMatrix<Stored>::value>>
  explicit LapackPackedMatrix(const LazyMatrixExpression<Stored>& m);

  /** Export to a lazyten matrix making it a full symmetric matrix
   *  again (i.e. both triangles will be set by this function) */
  void copy_symmetric_to(StoredMatrix_i<Scalar>& m);
};

//
// --------------------------------------------------
//

template <typename Scalar>
template <typename Stored, typename>
LapackPackedMatrix<Scalar>::LapackPackedMatrix(const LazyMatrixExpression<Stored>& m)
      : n(m.n_rows()), elements(n * (n + 1) / 2) {
  assert_dbg(m.is_symmetric(100 * lazyten::Constants<scalar_type>::default_tolerance),
             lazyten::ExcMatrixNotSymmetric());
  // For lazy matrices retrieving individual elements is extremely slow.
  // Therefore we tile the whole matrix into 4x4 block which we extract
  // using extract block.
  // TODO Check that 4 is sensible - we chose this value arbitrarily
  constexpr size_t bs = 4;

  // Convert indices i and j in row-major ordering (as the lazyten matrices use it)
  // into the packed column-major lower triagonal ordering
  // Lapack expects from us in this data structure.
  //
  // Effectively this means that we traverse the lower triangle of indices (i,j)
  // column-by-column. I.e. for a 5x5 matrix the indices of the packed array
  // Lapack expects map to the index of the lower triangle as follows:
  //    0 .  .  .  .
  //    1 5  .  .  .
  //    2 6  9  .  .
  //    3 7 10 12  .
  //    4 8 11 13 14
  // We notice that each column the difference between the first index values
  // is one less this results in the offset j*n - j*(j+1)/2 to get to the
  // index of the first value in column j. Overall we get:
  const size_t n = m.n_rows();
  auto toind = [&n](size_t i, size_t j) {
    assert_greater_equal(j, i);
    return j * n - j * (j + 1) / 2 + i;
  };

  // Use tiles of blocksize bs to extract the matrix data
  // Be careful not to overshoot by requiring roff+bs <= n
  // and coff <= roff, such that we go maximally to the
  // diagonal block and not overcover the matrix at the
  // bottom. We deal with the remaining entries in a second
  // step below

  size_t roff = 0;  //< The row index offset
  for (; roff + bs <= n; roff += bs) {
    for (size_t coff = 0; coff < 1 + roff; coff += bs) {
      Stored block(bs, bs, false);
      m.extract_block(block, roff, coff);

      for (size_t i = 0; i < bs; ++i) {
        // Block on the diagonal => only lower half is needed
        const size_t jend = roff == coff ? i + 1 : bs;

        for (size_t j = 0; j < jend; ++j) {
          elements[toind(roff + i, coff + j)] = block(i, j);
        }  // j
      }    // i
    }      // coff
  }        // roff

  // If something is missing => deal with it.
  if (roff < n) {
    assert_internal(n - roff > 0 && n - roff < bs);
    Stored slice(n - roff, n);
    m.extract_block(slice, roff, 0);
    for (size_t i = 0; i < slice.n_rows(); ++i) {
      const size_t iabs = roff + i;
      for (size_t j = 0; j < 1 + iabs; ++j) {
        elements[toind(iabs, j)] = slice(i, j);
      }  // i
    }    // j
  }      // roff
}

template <typename Scalar>
LapackPackedMatrix<Scalar>::LapackPackedMatrix(const StoredMatrix_i<Scalar>& m)
      : n(m.n_rows()), elements(n * (n + 1) / 2) {
  assert_dbg(m.is_symmetric(100 * lazyten::Constants<scalar_type>::default_tolerance),
             lazyten::ExcMatrixNotSymmetric());

  // Note: Since Fortran is column-major, but all our matrices are row-major,
  // we need to use the lower triangle-order when the upper triangle is requested
  // and vice versa.

  // TODO Use matrix iterator
  for (size_t i = 0, c = 0; i < n; ++i) {
    for (size_t j = i; j < n; ++j, ++c) {
      elements[c] = m(i, j);
    }
  }
}

template <typename Scalar>
void LapackPackedMatrix<Scalar>::copy_symmetric_to(StoredMatrix_i<Scalar>& m) {
  assert_size(m.n_rows(), n);
  assert_size(m.n_cols(), n);

  // TODO Improve cache misses
  for (size_t i = 0, c = 0; i < n; ++i) {
    for (size_t j = i; j < n; ++j, ++c) {
      m(i, j) = m(j, i) = elements[c];
    }
  }
}

}  // namespace detail
}  // namespace lazyten
#endif  // LAZYTEN_HAVE_LAPACK
