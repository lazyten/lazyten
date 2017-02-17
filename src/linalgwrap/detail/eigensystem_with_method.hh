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
#include "linalgwrap/Armadillo/ArmadilloEigensolver.hh"
#include "linalgwrap/Arpack/ArpackEigensolver.hh"
#include "linalgwrap/Base/Solvers.hh"
#include "linalgwrap/EigensystemKeys.hh"

namespace linalgwrap {
namespace detail {

/** Insert the default values into parameter map.
 * This step determines what the defaults are. */
template <typename Real>
void eigensystem_insert_defaults(const krims::ParameterMap& map);

/** Is it best to use Arpack solvers over others?
 *
 *
 * \note This metric is very much written out of the current case
 * where Arpack can only do mode 1 and 2 (normal and generalised
 * eigenproblems, no shift-and-invert) and only for real scalars.
 * */
template <typename MatrixA>
bool best_to_use_arpack(const MatrixA& A, typename MatrixA::size_type n_ep,
                        const krims::ParameterMap& map);

/** Solve the eigensystem with exactly the specified method */
template <typename Eigenproblem>
EigensolutionTypeFor<Eigenproblem::hermitian, typename Eigenproblem::matrix_diag_type>
eigensystem_with_method(const std::string& method, const Eigenproblem problem,
                        const krims::ParameterMap& map = krims::ParameterMap());

//
// ----------------------------------------------------------------------
//

template <typename Real>
void eigensystem_insert_defaults(const krims::ParameterMap& map) {
  map.insert_default(EigensystemKeys::method, "auto");
  map.insert_default(EigensystemKeys::which, "SR");
  map.insert_default(EigensystemKeys::tolerance, Constants<Real>::default_tolerance);
}

template <typename MatrixA>
bool best_to_use_arpack(const MatrixA& A, typename MatrixA::size_type n_ep,
                        const krims::ParameterMap& map) {

  // Arpack is best used if
  //    - we want a small number of eigenpairs
  //      in fact if we want more than this number of epairs it will fail
  //      TODO Note that this is just a shot and I have no clue whether
  //           half the dimension is a sensible value or not
  bool use_arpack = (n_ep < A.n_cols() / 2);

  //    - we want to get large magnitude eigenpairs
  //      (for now since mode 3 is not supported yet)
  use_arpack = use_arpack &&
               map.at<std::string>(EigensolverBaseKeys::which) != std::string("SM");

  //    - we want relatively high accuracy (not checked)
  //    - we have a non-dense problem (not checked)
  return use_arpack;
}

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
  constexpr bool isArpackAvailable = Eigenproblem::hermitian && Eigenproblem::real;
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
}  // namespace linalgwrap
