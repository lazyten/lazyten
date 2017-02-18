//
// Copyright (C) 2016-17 by the linalgwrap authors
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
#include "EigensystemSolver.hh"

namespace linalgwrap {

/** \brief Solve a normal hermitian eigensystem
 *
 * The Parameter map \t map may be used to provide configurable parameters
 * to the underlying Eigensolver. The number of understood options and
 * their default values and supported values depend on the underlying
 * eigensolver which is used. The options which are supported for all
 * eigensolvers can be found in the documentation of EigensystemSolver and
 * in the struct EigensystemSolverKeys.
 **/
template <typename Matrix>
EigensolutionTypeFor<true, Matrix> eigensystem_hermitian(
      const Matrix& A,
      typename std::enable_if<IsMatrix<Matrix>::value, typename Matrix::size_type>::type
            n_ep = Constants<typename Matrix::size_type>::all,
      const krims::GenMap& map = krims::GenMap());

/** Solve a generalised hermitian eigensystem
 *
 * The Parameter map \t map may be used to provide configurable parameters
 * to the underlying Eigensolver. The number of understood options and
 * their default values and supported values depend on the underlying
 * eigensolver which is used. The options which are supported for all
 * eigensolvers can be found in the documentation of EigensystemSolver and
 * in the struct EigensystemSolverKeys.
 *
 * \param A     Matrix to diagonalise
 * \param B     Metric to use
 * \param n_ep  Number of eigenpairs to obtain
 *              (Constants<typename Matrix::size_type>::all for all)
 * \param map   Specify some solver parameters
 *
 * \throws      Subclass of SolverException in case there is an error.
 **/
template <typename MatrixA, typename MatrixB>
EigensolutionTypeFor<true, MatrixA> eigensystem_hermitian(
      const MatrixA& A, const MatrixB& B,
      typename std::enable_if<IsMatrix<MatrixA>::value && IsMatrix<MatrixB>::value,
                              typename MatrixA::size_type>::type n_ep =
            Constants<typename MatrixA::size_type>::all,
      const krims::GenMap& map = krims::GenMap());

/** Solve a normal eigensystem
 *
 * The Parameter map \t map may be used to provide configurable parameters
 * to the underlying Eigensolver. The number of understood options and
 * their default values and supported values depend on the underlying
 * eigensolver which is used. The options which are supported for all
 * eigensolvers can be found in the documentation of EigensystemSolver and
 * in the struct EigensystemSolverKeys.
 *
 * \param A     Matrix to diagonalise
 * \param n_ep  Number of eigenpairs to obtain
 *              (Constants<typename Matrix::size_type>::all for all)
 * \param map   Specify some solver parameters
 *
 * \throws      Subclass of SolverException in case there is an error.
 **/
template <typename Matrix>
EigensolutionTypeFor<false, Matrix> eigensystem(
      const Matrix& A,
      typename std::enable_if<IsMatrix<Matrix>::value, typename Matrix::size_type>::type
            n_ep = Constants<typename Matrix::size_type>::all,
      const krims::GenMap& map = krims::GenMap());

/** Solve a generalised eigensystem
 *
 * The Parameter map \t map may be used to provide configurable parameters
 * to the underlying Eigensolver. See the non-generalised version of
 * eigensolver_hermitian for which parameters are usually supported
 * and there default values are.
 *
 * \param A     Matrix to diagonalise
 * \param B     Metric to use
 * \param n_ep  Number of eigenpairs to obtain
 *              (Constants<typename Matrix::size_type>::all for all)
 * \param map   Specify some solver parameters
 *
 * \throws      Subclass of SolverException in case there is an error.
 **/
template <typename MatrixA, typename MatrixB>
EigensolutionTypeFor<false, MatrixA> eigensystem(
      const MatrixA& A, const MatrixB& B,
      typename std::enable_if<IsMatrix<MatrixA>::value && IsMatrix<MatrixB>::value,
                              typename MatrixA::size_type>::type n_ep =
            Constants<typename MatrixA::size_type>::all,
      const krims::GenMap& map = krims::GenMap());

//
// ------------------------------------------------------
//

template <typename Matrix>
EigensolutionTypeFor<true, Matrix> eigensystem_hermitian(
      const Matrix& A,
      typename std::enable_if<IsMatrix<Matrix>::value, typename Matrix::size_type>::type
            n_ep,
      const krims::GenMap& map) {
  typedef Eigenproblem<true, Matrix> problem_type;
  problem_type problem{A, n_ep};

  return EigensystemSolver<problem_type>{map}.solve(problem).eigensolution();
}

template <typename MatrixA, typename MatrixB>
EigensolutionTypeFor<true, MatrixA> eigensystem_hermitian(
      const MatrixA& A, const MatrixB& B,
      typename std::enable_if<IsMatrix<MatrixA>::value && IsMatrix<MatrixB>::value,
                              typename MatrixA::size_type>::type n_ep,
      const krims::GenMap& map) {
  typedef Eigenproblem<true, MatrixA, MatrixB> problem_type;
  problem_type problem{A, B, n_ep};
  return EigensystemSolver<problem_type>{map}.solve(problem).eigensolution();
}

template <typename Matrix>
EigensolutionTypeFor<false, Matrix> eigensystem(
      const Matrix& A,
      typename std::enable_if<IsMatrix<Matrix>::value, typename Matrix::size_type>::type
            n_ep,
      const krims::GenMap& map) {
  // TODO This code is untested!
  assert_sufficiently_tested(false);

  // Setup problem
  typedef Eigenproblem<false, Matrix> problem_type;
  problem_type problem{A, n_ep};
  return EigensystemSolver<problem_type>{map}.solve(problem).eigensolution();
}

template <typename MatrixA, typename MatrixB>
EigensolutionTypeFor<false, MatrixA> eigensystem(
      const MatrixA& A, const MatrixB& B,
      typename std::enable_if<IsMatrix<MatrixA>::value && IsMatrix<MatrixB>::value,
                              typename MatrixA::size_type>::type n_ep,
      const krims::GenMap& map) {
  // TODO This code is untested!
  assert_sufficiently_tested(false);

  // Setup problem
  typedef Eigenproblem<false, MatrixA, MatrixB> problem_type;
  problem_type problem{A, B, n_ep};
  return EigensystemSolver<problem_type>{map}.solve(problem).eigensolution();
}

}  // namespace linalgwrap
