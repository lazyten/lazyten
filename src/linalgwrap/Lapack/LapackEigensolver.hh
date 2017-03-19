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
#include "detail/lapack.hh"
#include "linalgwrap/Base/Solvers.hh"
#include "linalgwrap/Exceptions.hh"
#include <krims/Algorithm.hh>

namespace linalgwrap {

DefSolverException2(ExcLapackInfo, std::string, function, int, info,
                    << "Lapack function " << function << " returned with info value "
                    << info
                    << ". Check Arpack documentation for the meaning of this value.");

template <typename Eigenproblem>
struct LapackEigensolverState : public EigensolverStateBase<Eigenproblem> {
  typedef EigensolverStateBase<Eigenproblem> base_type;
  typedef typename base_type::eproblem_type eproblem_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::size_type size_type;

  /** Setup the initial state from an eigenproblem to solve */
  LapackEigensolverState(const eproblem_type problem) : base_type(std::move(problem)) {}
};

/** Class which contains all GenMap keys which are understood
 *  by the LapackEigensolver update_control_params as static string
 *  members.
 *  See their doc strings for the types required. */
struct LapackEigensolverKeys : public EigensolverBaseKeys {};

/** \brief Lapack eigensolver class
 *
 * ## Control parameters and their default values
 *   - which:    Which eigenvalues to target. Default: "SR";
 *     allowed values:
 *       - "SR"   Smallest real value
 *       - "LR"   Largest real magnitude
 *       - "SM"   Smallest magnitude
 *       - "LM"   Largest magnitude
 *       - "SI"   Smallest imaginary
 *                (only for complex scalar types)
 *       - "LI"   Largest imaginary
 *                (only for complex scalar types)
 *
 * \note Currently only double precision scalar types can be used
 * \note Currently the tolerance parameter has no effect on the
 *       solver.
 * \note This solver implicitly solves for all eigenpairs
 *       and then truncates the result down to the actually requested
 *       number.
 *
 * \tparam Eigenproblem  The eigenproblem to solve.
 * \tparam State         The state type of the solver.
 */
template <typename Eigenproblem, typename State = LapackEigensolverState<Eigenproblem>>
class LapackEigensolver : public EigensolverBase<State> {
  static_assert(std::is_same<typename Eigenproblem::scalar_type, double>::value,
                "LapackEigensolver is only implemented for double precision scalar type");

  static_assert(Eigenproblem::hermitian,
                "Can only solve Hermitian eigenproblems at the moment.");

 public:
  //@{
  /** Forwarded types */
  typedef EigensolverBase<State> base_type;
  typedef typename base_type::state_type state_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::real_type real_type;
  typedef typename base_type::evalue_type evalue_type;
  typedef typename base_type::evector_type evector_type;
  typedef typename base_type::esoln_type esoln_type;
  typedef typename base_type::stored_matrix_type stored_matrix_type;
  //@}

  /** \name Constructor */
  //@{
  /** Construct an eigensolver with the default parameters */
  LapackEigensolver() {}

  /** Construct an eigensolver setting the parameters from the map */
  LapackEigensolver(const krims::GenMap& map) : LapackEigensolver() {
    base_type::update_control_params(map);
  }
  //@}

  /** Implementation of the IterativeSolver method */
  virtual void solve_state(state_type& state) const override;

 private:
  /** Check that the control parameters are in a valid state */
  void assert_valid_control_params(state_type& state) const;

  /** Order the eigenvalues and copy the appropriate ones (selected by which)
   * into the solution data structure */
  void copy_to_solution(size_type n_ep, const std::vector<evalue_type>& eval,
                        const std::vector<typename evector_type::scalar_type>& evec,
                        esoln_type& soln) const;
};

template <typename Eigenproblem, typename State>
void LapackEigensolver<Eigenproblem, State>::assert_valid_control_params(
      state_type& state) const {
  // TODO I do not like this, since some checks could be performed
  // before the call of solve ... separate this?
  //
  // Some stuff is more general and affects multiple solver classes
  // this is not reflected here
  //
  // Here is code duplication with ArmadilloEigensolver

  const Eigenproblem& problem = state.eigenproblem();

  //
  // which
  //
  const std::string& which = base_type::which;
  /*
  const bool complexok =
        krims::IsComplexNumber<evalue_type>::value && (which == "SI" || which == "LI");
        */
  solver_assert(
        which == "SM" || which == "LM" || which == "LR" || which == "SR" /*|| complexok*/,
        state, ExcInvalidSolverParametersEncountered(
                     "The value " + which + " for which is not allowed in an "
                                            "Lapack solver call (only SR, LR, SM and LM "
                                            "are accepted for problems with real "
                                            "eigenvalues and LI and SI are additionally "
                                            "accepted for complex eigenvalues."));

  //
  // A and Diag
  //
  // note: We compare memory addresses
  solver_assert(&problem.A() == &problem.Diag(), state,
                ExcInvalidSolverParametersEncountered("The matrices A and Diag need "
                                                      "to be the same objects."));
}

template <typename Eigenproblem, typename State>
void LapackEigensolver<Eigenproblem, State>::solve_state(state_type& state) const {
  assert_dbg(!state.is_failed(), krims::ExcInvalidState("Cannot solve a failed state"));

  esoln_type& soln = state.eigensolution();
  const Eigenproblem& problem = state.eigenproblem();

  assert_valid_control_params(state);
  static_assert(std::is_same<typename Eigenproblem::scalar_type, double>::value,
                "LapackEigensolver is only implemented for double precision scalar type");
  static_assert(Eigenproblem::hermitian,
                "Can only solve Hermitian eigenproblems at the moment.");

  // Run hermitian generalised or normal hermitian
  // eigenproblem for double precision:
  int info;
  std::vector<double> evals;
  std::vector<double> evecs;
  detail::LapackPackedMatrix<double> Ap{problem.A()};

  if (Eigenproblem::generalised) {
    detail::LapackPackedMatrix<double> Bp{problem.B()};
    auto Bp_chol = detail::run_dspgv(std::move(Ap), std::move(Bp), evals, evecs, info);
    // TODO Bp_chol contains the choleski decomposition of B!
    //      Expose this result!
    (void)Bp_chol;
  } else {
    detail::run_dspev(std::move(Ap), evals, evecs, info);
  }
  solver_assert(info == 0, state, ExcLapackInfo("dspgv", info));

  // Check that we get the solution in the expected format.
  assert_dbg(evals.size() * evals.size() == evecs.size(), krims::ExcInternalError());
  assert_dbg(evals.size() >= problem.n_ep(), krims::ExcInternalError());

  // Copy the results over to solution data structure:
  copy_to_solution(problem.n_ep(), evals, evecs, soln);
}

template <typename Eigenproblem, typename State>
void LapackEigensolver<Eigenproblem, State>::copy_to_solution(
      size_type n_ep, const std::vector<evalue_type>& eval,
      const std::vector<typename evector_type::scalar_type>& evec,
      esoln_type& soln) const {

  // Lapack always sorts its eigenvalues canonically
  //    TODO check for complex eigenvalues!
  const bool sorted_canonically = true;

  // (Sorted) indices of eigenvalues to keep
  const std::vector<size_t> idcs = select_eigenvalues(
        std::begin(eval), std::end(eval), base_type::which, n_ep, sorted_canonically);

  // Reserve space to copy the values in
  soln.evalues().clear();
  soln.evectors().clear();
  soln.evalues().reserve(n_ep);
  soln.evectors().reserve(n_ep);

  // The size of the eigenproblem:
  const size_t N = eval.size();
  for (const auto& idx : idcs) {
    // The iterator range describing the current vector
    // This works, since Fortran arrays are column-major,
    // i.e. column-by-column, vector-by-vector
    auto vbegin = std::begin(evec) + idx * N;
    auto vend = std::begin(evec) + (idx + 1) * N;

    soln.evectors().emplace_back(vbegin, vend);
    soln.evalues().push_back(eval[idx]);
  }
}

}  // namespace linalgwrap
#endif  // LINALGWRAP_HAVE_LAPACK
