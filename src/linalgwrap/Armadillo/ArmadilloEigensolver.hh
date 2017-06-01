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
#include "linalgwrap/config.hh"
#ifdef LINALGWRAP_HAVE_ARMADILLO

#include "ArmadilloMatrix.hh"
#include "detail.hh"
#include "linalgwrap/Base/Solvers.hh"
#include <armadillo>
#include <krims/Algorithm.hh>

namespace linalgwrap {

DefSolverException0(ExcArmadilloEigensolverFailed, << "Armadillo eigensolver failed");

template <typename Eigenproblem>
struct ArmadilloEigensolverState : public EigensolverStateBase<Eigenproblem> {
  typedef EigensolverStateBase<Eigenproblem> base_type;
  typedef typename base_type::eproblem_type eproblem_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::size_type size_type;

  /** Setup the initial state from an eigenproblem to solve */
  ArmadilloEigensolverState(const eproblem_type problem)
        : base_type(std::move(problem)) {}
};

/** Class which contains all GenMap keys which are understood
 *  by the ArmadilloSolver update_control_params as static string
 *  members.
 *  See their doc strings for the types required. */
struct ArmadilloEigensolverKeys : public EigensolverBaseKeys {};

/** \brief Armadillo eigensolver class
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
template <typename Eigenproblem, typename State = ArmadilloEigensolverState<Eigenproblem>>
class ArmadilloEigensolver : public EigensolverBase<State> {
  static_assert(std::is_same<typename Eigenproblem::scalar_type, double>::value,
                "ArmadilloEigensolver has only been tested for real problems at "
                "double precision at the moment.");

  static_assert(std::is_same<typename Eigenproblem::stored_matrix_type,
                             ArmadilloMatrix<typename Eigenproblem::scalar_type>>::value,
                "This eigensolver only works for if AmadilloMatrix is the underlying "
                "stored matrix type");

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
  ArmadilloEigensolver() {}

  /** Construct an eigensolver setting the parameters from the map */
  ArmadilloEigensolver(const krims::GenMap& map) : ArmadilloEigensolver() {
    base_type::update_control_params(map);
  }
  //@}

  /** Implementation of the IterativeSolver method */
  virtual void solve_state(state_type& state) const override;

 private:
  struct eigenpairsorter {
    template <typename T>
    bool operator()(typename std::enable_if<krims::IsComplexNumber<T>::value, T>::type i,
                    T j) const {
      return std::abs(i) < std::abs(j);
    }

    template <typename T>
    bool operator()(typename std::enable_if<!krims::IsComplexNumber<T>::value, T>::type i,
                    T j) const {
      return i < j;
    }
  };

  /** Check that the control parameters are in a valid state */
  void assert_valid_control_params(state_type& state) const;

  /** Order the eigenvalues and copy the appropriate ones (selected by which)
   * into the solution data structure */
  void copy_to_solution(size_type n_ep, const arma::Col<evalue_type>& eval_arma,
                        const arma::Mat<typename evector_type::scalar_type>& evec_arma,
                        esoln_type& soln) const;
};

//
// ----------------------------------------------------------
//

template <typename Eigenproblem, typename State>
void ArmadilloEigensolver<Eigenproblem, State>::copy_to_solution(
      size_type n_ep, const arma::Col<evalue_type>& eval_arma,
      const arma::Mat<typename evector_type::scalar_type>& evec_arma,
      esoln_type& soln) const {

  // Armadillo eigenpairs are only sorted canonically if they
  // result from a real Hermitian non-general problem.
  //    TODO check for complex eigenvalues!
  const bool sorted_canonically =
        !Eigenproblem::generalised && Eigenproblem::hermitian && Eigenproblem::real;

  // (Sorted) indices of eigenvalues to keep
  const std::vector<size_t> idcs =
        select_eigenvalues(std::begin(eval_arma), std::end(eval_arma), base_type::which,
                           n_ep, sorted_canonically);

  // Reserve space to copy the values in
  soln.evalues().clear();
  soln.evectors().clear();
  soln.evalues().reserve(n_ep);
  soln.evectors().reserve(n_ep);

  for (const auto& idx : idcs) {
    evector_type v(evec_arma.begin_col(idx), evec_arma.end_col(idx));
    soln.evectors().push_back(std::move(v));
    soln.evalues().push_back(eval_arma[idx]);
  }
}

template <typename Eigenproblem, typename State>
void ArmadilloEigensolver<Eigenproblem, State>::assert_valid_control_params(
      state_type& state) const {
  // TODO I do not like this, since some checks could be performed
  // before the call of solve ... separate this?
  //
  // Some stuff is more general and affects multiple solver classes
  // this is not reflected here
  //
  // Here is code duplication with LapackEigensolver

  const Eigenproblem& problem = state.eigenproblem();

  //
  // which
  //
  const std::string& which = base_type::which;
  const bool complexok =
        krims::IsComplexNumber<evalue_type>::value && (which == "SI" || which == "LI");
  solver_assert(
        which == "SM" || which == "LM" || which == "LR" || which == "SR" || complexok,
        state,
        ExcInvalidSolverParametersEncountered(
              "The value " + which + " for which is not allowed in an "
                                     "Armadillo solver call (only SR, LR, SM and LM "
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
void ArmadilloEigensolver<Eigenproblem, State>::solve_state(state_type& state) const {
  assert_dbg(!state.is_failed(), krims::ExcInvalidState("Cannot solve a failed state"));

  esoln_type& soln = state.eigensolution();
  const Eigenproblem& problem = state.eigenproblem();

  // TODO The eigensolver has not been tested for non-hermitian problems yet
  assert_sufficiently_tested(Eigenproblem::hermitian);

  assert_valid_control_params(state);

  // Solve eigenproblem
  arma::Col<evalue_type> eval_arma;
  arma::Mat<typename evector_type::scalar_type> evec_arma;
  const bool res =
        detail::ArmadilloEigWrapper<Eigenproblem>::eig(problem, eval_arma, evec_arma);
  solver_assert(res, state, ExcArmadilloEigensolverFailed());

  // Check that we get the solution in the expected format.
  assert_internal(eval_arma.n_cols == 1);
  assert_internal(eval_arma.n_rows == evec_arma.n_rows);
  assert_internal(eval_arma.n_rows == evec_arma.n_cols);
  assert_internal(eval_arma.n_rows >= problem.n_ep());

  // Copy the results over to solution data structure
  copy_to_solution(problem.n_ep(), eval_arma, evec_arma, soln);
}

}  // namespace linalgwrap
#endif  // LINALGWRAP_HAVE_ARMADILLO
