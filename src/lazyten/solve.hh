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
#include "LinearSolver.hh"

namespace lazyten {

namespace detail {
template <bool cond>
using GenMap_if = typename std::enable_if<cond, krims::GenMap>::type;
}

//@{
/** Solve a linear system A x = b
 *
 * \param A     System matrix (assumed to be hermitian)
 * \param x     Lhs vector or vectors (solution)
 * \param b     Rhs vector or vectors (problem)
 * \param map   Specify some solver parameters
 *
 * \throws      Subclass of SolverException in case there is an error.
 **/
template <typename Matrix, typename Vector>
void solve(
      const Matrix& A, MultiVector<Vector>& x, const MultiVector<const Vector>& b,
      const detail::GenMap_if<IsMatrix<Matrix>::value && IsMutableVector<Vector>::value>&
            map = krims::GenMap()) {
  assert_size(x.n_vectors(), b.n_vectors());
  assert_size(A.n_rows(), b.n_elem());
  assert_size(A.n_cols(), x.n_elem());

  if (A.has_apply_inverse()) {
    // A has an analytic inverse implemented,
    // so we better use that instead of any implicit stuff
    A.apply_inverse(b, x);
    return;
  }

  typedef LinearProblem<Matrix, Vector> problem_type;
  problem_type problem{A, b};

  // Solve A x = b
  LinearSolver<problem_type>{map}.solve(problem, x);
}

template <typename Matrix, typename Vector>
void solve(
      const Matrix& A, Vector& x, const Vector& b,
      const detail::GenMap_if<IsMatrix<Matrix>::value && IsMutableVector<Vector>::value>&
            map = krims::GenMap()) {
  MultiVector<Vector> x_mv(x);
  MultiVector<const Vector> b_mv(b);
  solve(A, x_mv, b_mv, map);
}

}  // namespace lazyten
