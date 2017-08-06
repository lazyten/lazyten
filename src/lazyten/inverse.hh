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
#include "InverseProxy.hh"
#include "detail/InvertibleWrapper.hh"

namespace lazyten {

/** Make a matrix invertible which is not invertible already by appending a linear solver
 * which applies the inverse by solving an appropriate linear system.
 *
 * \param params  Parameters for the linear solver.
 * */
template <typename Matrix>
detail::InvertibleWrapper<Matrix> make_invertible(
      const Matrix& m, krims::GenMap params = krims::GenMap()) {
  return detail::InvertibleWrapper<Matrix>(m, params);
}

//@{
/** Invert a matrix: Return an InverseProxy object, properly initialised
 *
 * \note The matrix needs to have an ``apply_inverse`` function for this
 *       to succeed.
 * */
template <typename Matrix,
          typename = krims::enable_if_t<
                IsMatrix<typename std::remove_reference<Matrix>::type>::value>>
InverseProxy<typename std::remove_reference<Matrix>::type> inverse(Matrix&& m) {
  return InverseProxy<typename std::remove_reference<Matrix>::type>(
        std::forward<Matrix>(m));
}

template <typename Matrix, krims::enable_if_t<IsMatrix<Matrix>::value>...>
Matrix inverse(InverseProxy<Matrix>&& mt) {
  if (mt.owns_inner_matrix()) {
    return Matrix(std::move(mt.inner_matrix()));
  }

  // TODO This is not yet implemented ... we need matrix views for that
  assert_implemented(false);
  return Matrix(mt.inner_matrix());
}

template <typename Matrix, krims::enable_if_t<IsMatrix<Matrix>::value>...>
Matrix& inverse(InverseProxy<Matrix>& mt) {
  // TODO This is not yet implemented ... we need matrix views for that
  assert_implemented(false);
  return mt.inner_matrix();
}

template <typename Matrix, krims::enable_if_t<IsMatrix<Matrix>::value>...>
const Matrix& inverse(const InverseProxy<Matrix>& mt) {
  // TODO This is not yet implemented ... we need matrix views for that
  assert_implemented(false);
  return mt.inner_matrix();
}
//@}
}  // namespace lazyten
