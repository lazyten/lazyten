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
#include "linalgwrap/LazyMatrixExpression.hh"
#include "linalgwrap/StoredMatrix_i.hh"
#include <vector>

#ifdef LINALGWRAP_HAVE_ARMADILLO
#include "linalgwrap/Armadillo/ArmadilloMatrix.hh"
#endif  // LINALGWRAP_HAVE_ARMADILLO

namespace linalgwrap {
namespace detail {

/** Data structure to represent a matrix for Lapack
 *
 * From the point of view of Fortran this matrix contains
 * data in column-major form
 */
template <typename Scalar>
struct LapackSymmetricMatrix {
  typedef Scalar scalar_type;

  //! The number of rows and columns
  size_t n;

  //! The matrix elements in column-major form
  std::vector<Scalar> elements;

  /** Default construct */
  LapackSymmetricMatrix() = default;

  /** Construct from a usual *symmetric* linalgwrap stored matrix by
   * copying the values in */
  explicit LapackSymmetricMatrix(const StoredMatrix_i<Scalar>& m);

  /** Construct from an armadillo matrix */
  explicit LapackSymmetricMatrix(const LazyMatrixExpression<ArmadilloMatrix<Scalar>>& m);
  // TODO Ideally we want something like the above function for every stored matrix type.

  /** Construct from a usual *symmetric* linalgwrap matrix expression by
   * copying the values in */
  template <typename Stored, typename = krims::enable_if_t<IsStoredMatrix<Stored>::value>>
  explicit LapackSymmetricMatrix(const LazyMatrixExpression<Stored>& m)
        : LapackSymmetricMatrix(static_cast<Stored>(m)) {
    assert_dbg(false, krims::ExcDisabled("This version does two copies of the full "
                                         "matrix and is hence disabled. Implement "
                                         "similar things to the ArmadilloMatrix version "
                                         "above for other stored matrix types."));
  }
};

//
// --------------------------------------------------
//

template <typename Scalar>
LapackSymmetricMatrix<Scalar>::LapackSymmetricMatrix(const StoredMatrix_i<Scalar>& m)
      : n(m.n_cols()), elements(n * n) {
  // The matrix is symmetric so no need to take care about column-major / row-major
  // between linalgwrap matrices (column-major) and Fortran (row-major)
  assert_dbg(m.is_symmetric(100 * linalgwrap::Constants<scalar_type>::default_tolerance),
             linalgwrap::ExcMatrixNotSymmetric());
  std::copy(std::begin(m), std::end(m), std::begin(elements));
}

#ifdef LINALGWRAP_HAVE_ARMADILLO
template <typename Scalar>
LapackSymmetricMatrix<Scalar>::LapackSymmetricMatrix(
      const LazyMatrixExpression<ArmadilloMatrix<Scalar>>& m)
      : n(m.n_cols()), elements(n * n) {

  constexpr bool copy_into_arma = false;    // Do not copy the memory
  constexpr bool fixed_vector_size = true;  // No memory reallocation
  arma::Mat<Scalar> m_arma(elements.data(), n, n, copy_into_arma, fixed_vector_size);
  assert_dbg(m_arma.memptr() == elements.data(), krims::ExcInternalError());

  ArmadilloMatrix<Scalar> block(std::move(m_arma));
  assert_dbg(block.data().memptr() == elements.data(), krims::ExcInternalError());

  // Extract the block and write into this->elements
  m.extract_block(block, 0, 0);
}
#endif  // LINALGWRAP_HAVE_ARMADILLO

}  // namespace detail
}  // namespace linalgwrap
#endif  // LINALGWRAP_HAVE_LAPACK
