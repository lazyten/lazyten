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
#include "linalgwrap/Armadillo/ArmadilloEigensolver.hh"
#include "linalgwrap/Arpack/ArpackEigensolver.hh"
#include "linalgwrap/Base/Solvers.hh"
#include "linalgwrap/Lapack/LapackEigensolver.hh"

namespace linalgwrap {

namespace detail {
/** Run an eigensolver and update the passed state accordingly */
template <typename Solver>
struct RunSolver {
  template <typename State>
  void run(State& state, const krims::GenMap& params) const;
};

/** Specialisation of RunSolver for a disabled solver (does nothing) */
template <>
struct RunSolver<void> {
  template <typename State>
  void run(State&, const krims::GenMap&) const {}
};
}  // namespace detail

template <typename Eigenproblem>
class EigensystemSolverState final : public EigensolverStateBase<Eigenproblem> {
  // Use final to prevent overwriting from this class, see EigensystemSolver
  // below for the reasons why
 public:
  typedef EigensolverStateBase<Eigenproblem> base_type;
  typedef typename base_type::eproblem_type eproblem_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::size_type size_type;

  /** Setup the initial state from an eigenproblem to solve */
  EigensystemSolverState(const eproblem_type problem)
        : base_type(std::move(problem)), m_n_iter{0}, m_n_mtx_applies(0) {}

  /** Transfer data from an arbitrary other state by copying
   * it into this one. This is needed, since we need to get the
   * data back here once the inner solvers are done with
   * (partially) solving the eigenproblem.
   *
   * \note Do not use this to setup a guess state.
   * For this case the function obtain_guess_from exists.
   **/
  void push_intermediate_results(base_type&& other_state) {
    m_n_iter += other_state.n_iter();
    m_n_mtx_applies += other_state.n_mtx_applies();
    static_cast<base_type&>(*this) = other_state;
  }

  size_t n_iter() const override final { return m_n_iter; }
  size_t n_mtx_applies() const override final { return m_n_mtx_applies; }

 private:
  size_t m_n_iter;
  size_t m_n_mtx_applies;
};

struct EigensystemSolverKeys final : public EigensolverBaseKeys {
  /** The solver method to use. Type: string */
  static const std::string method;
};

/** \brief Envelope eigensolver that calls some
 *         inner Eigensolver depending on certain criteria
 *
 * ## Control parameters and their default values
 *
 * The GenMap map \t map may be used to provide configurable parameters
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
 *       - "lapack"      Use Lapack
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
 * \tparam Eigenproblem  The eigenproblem to solve.
 */
template <typename Eigenproblem>
class EigensystemSolver final
      : public EigensolverBase<EigensystemSolverState<Eigenproblem>> {
  // We have the final keyword because when someone overrides from this to hook into
  // the algorithm this will fail (the handlers are actually never called)
 public:
  //@{
  /** Forwarded types */
  typedef EigensystemSolverState<Eigenproblem> state_type;
  typedef EigensolverBase<state_type> base_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::real_type real_type;
  typedef typename base_type::evalue_type evalue_type;
  typedef typename base_type::evector_type evector_type;
  typedef typename base_type::esoln_type esoln_type;
  typedef typename base_type::stored_matrix_type stored_matrix_type;
  typedef typename state_type::eproblem_type eproblem_type;
  //@}

  /** \name Constructor */
  //@{
  /** Construct an eigensolver with the default parameters */
  EigensystemSolver() {}

  /** Construct an eigensolver setting the parameters from the map */
  EigensystemSolver(const krims::GenMap& map) : EigensystemSolver() {
    update_control_params(map);
  }
  //@}

  /** \name Iteration control */
  ///@{
  // The method to use to actually solve the underlying eigenproblem.
  std::string method = "auto";

  /** Update control parameters from GenMap map */
  void update_control_params(const krims::GenMap& map) {
    base_type::update_control_params(map);
    method = map.at(EigensystemSolverKeys::method, method);

    // Copy the map to the internal storage such that we
    // can pass it on to the actual eigensolvers.
    m_solver_params = map;
  }

  /** Get the current settings of all internal control parameters and
   *  update the GenMap accordingly. */
  void get_control_params(krims::GenMap& map) const {
    base_type::get_control_params(map);
    map.update(EigensystemSolverKeys::method, method);
  }
  ///@}

  virtual void solve_state(state_type& state) const override final;

 private:
  bool should_use_arpack(const Eigenproblem& problem) const;
  bool should_use_armadillo(const Eigenproblem& problem) const;
  bool should_use_lapack(const Eigenproblem& problem) const;

  /** Setup and solve using the method provided */
  void solve_with_method(const std::string& method, state_type& state) const;

  /** Cache of the parameters which will be passed to the eigensolver.
   *
   * A mixture of the original parameters, which the user passed to us
   * and the current state which is reflected in the member variables
   * in this class and the subclasses. Should be updated with
   * get_control_params *before* the actual inner eigensolver invocation.
   */
  mutable krims::GenMap m_solver_params;
};

//
// ----------------------------------------------------------
//

namespace detail {
template <typename Solver>
template <typename State>
void RunSolver<Solver>::run(State& state, const krims::GenMap& params) const {
  typedef typename Solver::state_type solver_state_type;

  // Setup the inner solver state:
  solver_state_type inner_state{state.eigenproblem()};

  // TODO move the inner eigensolution here -> gets rid of a copy!
  inner_state.obtain_guess_from(state);

  try {
    Solver{params}.solve_state(inner_state);
    state.push_intermediate_results(std::move(inner_state));
  } catch (SolverException& e) {
    // On exception still update the state reference
    state.push_intermediate_results(std::move(inner_state));
    throw;
  }
}
}  // namespace detail

template <typename Eigenproblem>
bool EigensystemSolver<Eigenproblem>::should_use_arpack(
      const Eigenproblem& problem) const {
#ifdef LINALGWRAP_HAVE_ARPACK
  // Arpack should be used if:

  //   - not a complex or non-hermitian problem:
  //     (TODO Not yet implemented)
  if (!Eigenproblem::hermitian || !Eigenproblem::real) return false;

  //   - not a generalised problem without apply_inverse in B
  if (Eigenproblem::generalised && !problem.B().has_apply_inverse()) return false;

  //   - not a large number of eigenvalues is desired
  //     (TODO Note that this is just a shot and I have no clue whether
  //      half the dimension is a sensible value or not, for sure ARPACK
  //      cannot do larger problems)
  if (problem.n_ep() >= problem.dim() / 2) return false;

  //   - not if we desire small magnitude eigenvalues
  //     (TODO mode 3 is not yet implemented)
  if (base_type::which == std::string("SM")) return false;
  if (base_type::which == std::string("SR")) return false;

  return true;
#else
  return false;
#endif  // LINALGWRAP_HAVE_ARPACK
}

template <typename Eigenproblem>
bool EigensystemSolver<Eigenproblem>::should_use_lapack(
      const Eigenproblem& /*problem*/) const {
#ifdef LINALGWRAP_HAVE_LAPACK
  // Currently we only have Lapack for real, hermitian eigenproblems
  // implemented.
  if (Eigenproblem::real && Eigenproblem::hermitian) return true;

  return false;
#else
  return false;
#endif  // LINALGWRAP_HAVE_LAPACK
}

template <typename Eigenproblem>
bool EigensystemSolver<Eigenproblem>::should_use_armadillo(
      const Eigenproblem& /*problem*/) const {
#ifdef LINALGWRAP_HAVE_ARMADILLO
  // Armadillo sucks with real hermitian generalised problems
  // since it insists to do it in complex arithmetic:
  if (Eigenproblem::real && Eigenproblem::hermitian && Eigenproblem::generalised) {
    return false;
  }

  // It is perfect for problems with armadillo matrices:
  if (std::is_same<typename Eigenproblem::matrix_diag_type,
                   ArmadilloMatrix<typename Eigenproblem::scalar_type>>::value) {
    return true;
  }

  return false;
#else
  return false;
#endif  // LINALGWRAP_HAVE_ARMADILLO
}

template <typename Eigenproblem>
void EigensystemSolver<Eigenproblem>::solve_with_method(const std::string& method,
                                                        state_type& state) const {
  const std::string errorstring = "The eigensolver method " + method +
                                  "(set via the key '" + EigensystemSolverKeys::method +
                                  "') is not compiled into this version of linalgwrap.";

  // Make sure the control parameters are up to date:
  get_control_params(m_solver_params);

  //
  // Arpack
  //
  if (method == std::string("arpack")) {
#ifdef LINALGWRAP_HAVE_ARPACK
    // Only instantiate the Arpack Eigensolver type in case
    // the problem is hermitian.
    typedef typename std::conditional<Eigenproblem::hermitian,
                                      ArpackEigensolver<Eigenproblem>, void>::type
          cond_arpack_type;
    detail::RunSolver<cond_arpack_type>{}.run(state, m_solver_params);
    return;
#else
    assert_throw(false, ExcInvalidSolverParametersEncountered(errorstring));
#endif  // LINALGWRAP_HAVE_ARPACK
  }

  //
  // Lapack
  //
  if (method == std::string("lapack")) {
#ifdef LINALGWRAP_HAVE_LAPACK
    // Only instantiate the Lapack Eigensolver type in case
    // the problem is hermitian.
    typedef typename std::conditional<Eigenproblem::hermitian,
                                      LapackEigensolver<Eigenproblem>, void>::type
          cond_lapack_type;
    detail::RunSolver<cond_lapack_type>{}.run(state, m_solver_params);
    return;
#else
    assert_throw(false, ExcInvalidSolverParametersEncountered(errorstring));
#endif  // LINALGWRAP_HAVE_LAPACK
  }

  //
  // Armadillo
  //
  if (method == std::string("armadillo")) {
#ifdef LINALGWRAP_HAVE_ARMADILLO
    detail::RunSolver<ArmadilloEigensolver<Eigenproblem>>{}.run(state, m_solver_params);
    return;
#else
    assert_throw(false, ExcInvalidSolverParametersEncountered(errorstring));
#endif  // LINALGWRAP_HAVE_ARMADILLO
  }

  //
  // No method is supported!
  //
  assert_throw(false, ExcInvalidSolverParametersEncountered(
                            "The eigensolver method " + method + "(set via the key " +
                            EigensystemSolverKeys::method +
                            ") is unknown. Did you spell it wrong?"));
}

template <typename Eigenproblem>
void EigensystemSolver<Eigenproblem>::solve_state(state_type& state) const {
  assert_dbg(!state.is_failed(), krims::ExcInvalidState("Cannot solve a failed state"));

  /** User-selected */
  if (method != std::string("auto")) {
    solve_with_method(method, state);
    return;
  }

  if (should_use_arpack(state.eigenproblem())) {
    solve_with_method("arpack", state);
    return;
  }
  if (should_use_armadillo(state.eigenproblem())) {
    solve_with_method("armadillo", state);
    return;
  }
  if (should_use_lapack(state.eigenproblem())) {
    solve_with_method("lapack", state);
    return;
  }

  assert_throw(false,
               ExcInvalidSolverParametersEncountered(
                     "Could not autodetermine an eigensolver for your eigenproblem. "
                     "This could mean that there are not enough solvers available in "
                     "your linalgwrap installation. Try forcing the use of an existing "
                     "solver via the \"method\" parameter in this case."));
}

}  // namespace linalgwrap
