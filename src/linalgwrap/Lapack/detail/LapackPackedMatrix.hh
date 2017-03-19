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
#ifdef LINALGWRAP_HAVE_LAPACK
#include "linalgwrap/StoredMatrix_i.hh"
#include <vector>

namespace linalgwrap {
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

  //! The linearised packed matrix elements.
  std::vector<Scalar> elements;

  /** Default construct */
  LapackPackedMatrix() = default;

  /** Construct from a usual symmetric linalgwrap matrix by copying the
   *  values in
   */
  LapackPackedMatrix(const Matrix_i<Scalar>& m);

  /** Export to a linalgwrap matrix making it a full symmetric matrix
   *  again (i.e. both triangles will be set by this function) */
  void copy_symmetric_to(StoredMatrix_i<Scalar>& m);
};

//
// --------------------------------------------------
//

template <typename Scalar>
LapackPackedMatrix<Scalar>::LapackPackedMatrix(const Matrix_i<Scalar>& m)
      : n(m.n_rows()), elements(n * (n + 1) / 2) {
  assert_dbg(m.is_symmetric(100 * linalgwrap::Constants<scalar_type>::default_tolerance),
             linalgwrap::ExcMatrixNotSymmetric());

  // Note: Since Fortran is column-major, but all our matrices are row-major,
  // we need to use the lower triangle-order when the upper triangle is requested
  // and vice versa.

  // TODO Use extract_block
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
}  // namespace linalgwrap
#endif  // LINALGWRAP_HAVE_LAPACK
