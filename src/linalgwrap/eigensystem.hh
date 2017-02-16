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
#include "Armadillo/ArmadilloEigensolver.hh"
#include "Arpack/ArpackEigensolver.hh"
#include "Base/Solvers.hh"
#include <krims/ParameterMap.hh>

namespace linalgwrap {

#ifndef LINALGWRAP_HAVE_ARMADILLO
static_assert(false,
              "We need armadillo in order to have at least a working fallback "
              "eigensolver.");
#endif

namespace detail {
template <typename Eigenproblem>
EigensolutionTypeFor<Eigenproblem::hermitian, typename Eigenproblem::matrix_diag_type>
eigensystem_with_method(const std::string& method, const Eigenproblem problem,
                        const krims::ParameterMap& map = krims::ParameterMap());
}  // namespace detail

struct EigensystemKeys : public EigensolverBaseKeys {
  /** The solver method to use. Type: string */
  static const std::string method;
};

// TODO document default and supported keys (method, tolerance and which)
// valid method names : "arpack" or "armadillo" or "auto"

/** \brief Solve a normal hermitian eigensystem
 *
 * The Parameter map \t map may be used to provide configurable parameters
 * to the underlying Eigensolver. The number of understood options and
 * their default values and supported values depend on the underlying
 * eigensolver which is used. The options which are supported for all
 * eigensolvers are:
 *   - method:   Enforce that a particular eigensolver method should be
 *               used. Allowed values:
 *       - "auto"   Auto-select solver due to hardcoded criteria
 *                  (tolerance, number of eigenpairs, dimensionality /
 *                   sparsity of matrix, ...)
 *       - "arpack"   Use ARPACK
 *       - "armadillo"   Use Armadillo
 *   - which:    Which eigenvalues to target. Default: "SR";
 *     allowed values (for all eigensolvers):
 *       - "SM"   Smallest magnitude
 *       - "LM"   Largest magnitude
 *       - "SR"   Smallest real
 *       - "LR"   Largest real
 *       - "SI"   Smallest imaginary
 *                (only for complex scalar types)
 *       - "LI"   Largest imaginary
 *                (only for complex scalar types)
 *   - tolerance: Tolerance for eigensolver. Default: Default numeric
 *                tolerance (as in Constants.hh)
 *
 * \param A     Matrix to diagonalise
 * \param n_ep  Number of eigenpairs to obtain
 *              (Constants<typename Matrix::size_type>::all for all)
 * \param map   Specify some solver parameters
 *
 * \throws      Subclass of SolverException in case there is an error.
 **/
template <typename Matrix>
EigensolutionTypeFor<true, Matrix> eigensystem_hermitian(
      const Matrix& A,
      typename std::enable_if<IsMatrix<Matrix>::value, typename Matrix::size_type>::type
            n_ep = Constants<typename Matrix::size_type>::all,
      const krims::ParameterMap& map = krims::ParameterMap());

/** Solve a generalised hermitian eigensystem
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
EigensolutionTypeFor<true, MatrixA> eigensystem_hermitian(
      const MatrixA& A, const MatrixB& B,
      typename std::enable_if<IsMatrix<MatrixA>::value && IsMatrix<MatrixB>::value,
                              typename MatrixA::size_type>::type n_ep =
            Constants<typename MatrixA::size_type>::all,
      const krims::ParameterMap& map = krims::ParameterMap());

/** Solve a normal eigensystem
 *
 * The Parameter map \t map may be used to provide configurable parameters
 * to the underlying Eigensolver. See the non-generalised version of
 * eigensolver_hermitian for which parameters are usually supported
 * and there default values are.
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
      const krims::ParameterMap& map = krims::ParameterMap());

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
      const krims::ParameterMap& map = krims::ParameterMap());

//
// ------------------------------------------------------
//

namespace detail {
template <typename Eigenproblem>
EigensolutionTypeFor<Eigenproblem::hermitian, typename Eigenproblem::matrix_diag_type>
eigensystem_with_method(const std::string& method, const Eigenproblem problem,
                        const krims::ParameterMap& map) {

#ifdef LINALGWRAP_HAVE_ARPACK
  // We need a hack here, since Arpack is in fact only available for real
  // Hermitian non-general Eigenproblems. So we define a constexpr, which
  // tells us whether this is the case and a conditional type which is only
  // evaluating to Arpack in case we actually have it availabe (to make
  // the compiler happy)

  // Is Arpack available for this Eigenproblem?
  constexpr bool isArpackAvailable =
        Eigenproblem::hermitian && !Eigenproblem::generalised && Eigenproblem::real;

  typedef typename std::conditional<isArpackAvailable, ArpackEigensolver<Eigenproblem>,
                                    ArmadilloEigensolver<Eigenproblem>>::type
        conditional_arpack_type;

  if (method == "arpack") {
    // No method is supported!
    assert_throw(isArpackAvailable, ExcInvalidSolverParametersEncountered(
                                          "Arpack is currently only available for real "
                                          "symmetric non-general eigenproblems."));

    return conditional_arpack_type{map}.solve(std::move(problem)).eigensolution();
  }
#endif
#ifdef LINALGWRAP_HAVE_ARMADILLO
  if (method == "armadillo") {
    return ArmadilloEigensolver<Eigenproblem>{map}
          .solve(std::move(problem))
          .eigensolution();
  }
#endif

  // No method is supported!
  assert_throw(false, ExcInvalidSolverParametersEncountered(
                            "The eigensolver method " + method + "(set via the key " +
                            EigensystemKeys::method +
                            ") is not available. Either you spelled it wrong "
                            "or this method has not ben compiled into this "
                            "version of linalgwrapj"));

  return EigensolutionTypeFor<Eigenproblem::hermitian,
                              typename Eigenproblem::matrix_diag_type>{};
}
}  // namespace detail

template <typename Matrix>
EigensolutionTypeFor<true, Matrix> eigensystem_hermitian(
      const Matrix& A,
      typename std::enable_if<IsMatrix<Matrix>::value, typename Matrix::size_type>::type
            n_ep,
      const krims::ParameterMap& map) {
  // Setup problem
  typedef Eigenproblem<true, Matrix> problem_type;
  problem_type problem{A, n_ep};

  // Select method (auto or user-defined)
  const std::string method = map.at("method", std::string("auto"));
  if (method != "auto") {
    return detail::eigensystem_with_method(method, problem, map);
  }

#ifdef LINALGWRAP_HAVE_ARPACK
  if (problem.n_ep() < problem.dim() / 2 &&
      map.at(EigensolverBaseKeys::which, std::string("SR")) != std::string("SM")) {
    // Use Arpack since we want "few" eigenpairs
    // TODO This is just a shot, no clue whether half the dimension
    // is a sensible value or not.
    return detail::eigensystem_with_method("arpack", problem, map);
  }
#endif

  // Fallback: Armadillo
  return detail::eigensystem_with_method("armadillo", problem, map);
}

template <typename MatrixA, typename MatrixB>
EigensolutionTypeFor<true, MatrixA> eigensystem_hermitian(
      const MatrixA& A, const MatrixB& B,
      typename std::enable_if<IsMatrix<MatrixA>::value && IsMatrix<MatrixB>::value,
                              typename MatrixA::size_type>::type n_ep,
      const krims::ParameterMap& map) {
  // Setup problem
  typedef Eigenproblem<true, MatrixA, MatrixB> problem_type;
  problem_type problem{A, B, n_ep};

  // Select method (auto or user-defined)
  const std::string method = map.at("method", std::string("auto"));
  if (method != "auto") {
    return detail::eigensystem_with_method(method, problem, map);
  }

  // TODO Do something better here some day (e.g. use Arpack)

  // Fallback: Armadillo
  return detail::eigensystem_with_method("armadillo", problem, map);
}

template <typename Matrix>
EigensolutionTypeFor<false, Matrix> eigensystem(
      const Matrix& A,
      typename std::enable_if<IsMatrix<Matrix>::value, typename Matrix::size_type>::type
            n_ep,
      const krims::ParameterMap& map) {
  // TODO This code is untested!
  assert_sufficiently_tested(false);

  // Setup problem
  typedef Eigenproblem<false, Matrix> problem_type;
  problem_type problem{A, n_ep};

  // Select method (auto or user-defined)
  const std::string method = map.at("method", std::string("auto"));
  if (method != "auto") {
    return detail::eigensystem_with_method(method, problem, map);
  }

  // TODO Do something better here some day (e.g. use Arpack)

  // Fallback: Armadillo
  return detail::eigensystem_with_method("armadillo", problem, map);
}

template <typename MatrixA, typename MatrixB>
EigensolutionTypeFor<false, MatrixA> eigensystem(
      const MatrixA& A, const MatrixB& B,
      typename std::enable_if<IsMatrix<MatrixA>::value && IsMatrix<MatrixB>::value,
                              typename MatrixA::size_type>::type n_ep,
      const krims::ParameterMap& map) {
  // TODO This code is untested!
  assert_sufficiently_tested(false);

  // Setup problem
  typedef Eigenproblem<false, MatrixA, MatrixB> problem_type;
  problem_type problem{A, B, n_ep};

  // Select method (auto or user-defined)
  const std::string method = map.at("method", std::string("auto"));
  if (method != "auto") {
    return detail::eigensystem_with_method(method, problem, map);
  }

  // TODO Do something better here some day (e.g. use Arpack)

  // Fallback: Armadillo
  return detail::eigensystem_with_method("armadillo", problem, map);
}

}  // namespace linalgwrap
