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
#include "TransposeProxy.hh"

namespace linalgwrap {
//
// trans
//
//@{
/** Transpose a matrix: Return a TransposeProxy object, properly initialised */
template <typename Matrix, typename = krims::enable_if_t<IsMatrix<
                                 typename std::remove_reference<Matrix>::type>::value>>
TransposeProxy<typename std::remove_reference<Matrix>::type> trans(Matrix&& m) {
  return TransposeProxy<typename std::remove_reference<Matrix>::type>(
        std::forward<Matrix>(m));
}

template <typename Matrix, krims::enable_if_t<IsMatrix<Matrix>::value>...>
Matrix trans(TransposeProxy<Matrix>&& mt) {
  if (mt.owns_inner_matrix()) {
    return Matrix(std::move(mt.inner_matrix()));
  }

  // TODO This is not yet implemented ... we need matrix views for that
  assert_dbg(false, krims::ExcNotImplemented());
  return Matrix(mt.inner_matrix());
}

template <typename Matrix, krims::enable_if_t<IsMatrix<Matrix>::value>...>
Matrix& trans(TransposeProxy<Matrix>& mt) {
  // TODO This is not yet implemented ... we need matrix views for that
  assert_dbg(false, krims::ExcNotImplemented());
  return mt.inner_matrix();
}

template <typename Matrix, krims::enable_if_t<IsMatrix<Matrix>::value>...>
const Matrix& trans(const TransposeProxy<Matrix>& mt) {
  // TODO This is not yet implemented ... we need matrix views for that
  assert_dbg(false, krims::ExcNotImplemented());
  return mt.inner_matrix();
}
//@}

//
// conjtrans
//
/** Take the conjugate transpose (Hermitian conjugate): Specialisation
 * for real scalar types */
template <typename Matrix,
          typename = krims::enable_if_t<
                IsMatrix<typename std::remove_reference<Matrix>::type>::value &&
                !HasComplexScalar<Matrix>::value>>
auto conjtrans(Matrix&& m) -> decltype(trans(std::forward<Matrix>(m))) {
  return trans(std::forward<Matrix>(m));
}

/** Take the conjugate transpose (Hermitian conjugate): Specialisation
 * for scalar data types */
template <
      typename Matrix,
      krims::enable_if_t<IsMatrix<typename std::remove_reference<Matrix>::type>::value &&
                               HasComplexScalar<Matrix>::value,
                         int> = 0>
auto conjtrans(Matrix&& m) -> decltype(trans(std::forward<Matrix>(m))) {
  static_assert(IsMatrix<Matrix>::value, "Matrix should be a matrix.");
  // TODO change once we have tested Transposed::ConjTrans on most objects
  // Also have a ConjugateTransposeProxy class
  static_assert(false && !HasComplexScalar<Matrix>::value,
                "conjtrans is currently only implemented for real scalar "
                "types.");
  return trans(std::forward<Matrix>(m));
}

}  // namespace linalgwrap
