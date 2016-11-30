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

#ifndef LINALGWRAP_HAVE_ARMADILLO
static_assert(false,
              "We need armadillo in order to have at least a working fallback "
              "eigensolver.");
#endif

namespace linalgwrap {

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
 * \tparam Eigenproblem  The eigenproblem to solve.
 * \tparam State         The state type of the solver.
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
  EigensystemSolver(const krims::ParameterMap& map) : EigensystemSolver() {
    update_control_params(map);
  }
  //@}

  /** \name Iteration control */
  ///@{
  // The method to use to actually solve the underlying eigenproblem.
  std::string method = "auto";

  /** Update control parameters from Parameter map */
  void update_control_params(const krims::ParameterMap& map) {
    base_type::update_control_params(map);
    method = map.at(EigensystemSolverKeys::method, method);

    // Copy the map to the internal storage such that we
    // can pass it on to the actual eigensolvers.
    m_solver_params = map;
  }

  /** Get the current settings of all internal control parameters and
   *  update the ParameterMap accordingly.
   */
  void get_control_params(krims::ParameterMap& map) const {
    base_type::get_control_params(map);
    map.update(EigensystemSolverKeys::method, method);
  }
  ///@}

  virtual void solve_state(state_type& state) const override final;

 private:
  /** Should we use Arpack to do the solving */
  bool should_use_arpack(const Eigenproblem& problem) const;

  /** Run the given solver and update the state accordingly */
  template <typename Solver>
  void run_solver(state_type& state) const;

  /** Setup and solve using the method provided */
  void solve_with_method(const std::string& method, state_type& state) const;

  /** Cache of the parameters which will be passed to the eigensolver.
   *
   * A mixture of the original parameters, which the user passed to us
   * and the current state which is reflected in the member variables
   * in this class and the subclasses. Should be updated with
   * get_control_params *before* the actual inner eigensolver invocation.
   */
  mutable krims::ParameterMap m_solver_params;
};

//
// ----------------------------------------------------------
//

template <typename Eigenproblem>
bool EigensystemSolver<Eigenproblem>::should_use_arpack(
      const Eigenproblem& problem) const {
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

  return true;
}

template <typename Eigenproblem>
template <typename Solver>
void EigensystemSolver<Eigenproblem>::run_solver(state_type& state) const {
  typedef typename Solver::state_type solver_state_type;

  // Setup the inner solver state:
  solver_state_type inner_state{state.eigenproblem()};
  inner_state.obtain_guess_from(state);

  // Make sure that control parameters are up to date:
  get_control_params(m_solver_params);

  try {
    Solver{m_solver_params}.solve_state(inner_state);
    state.push_intermediate_results(std::move(inner_state));
  } catch (SolverException& e) {
    // On exception still update the state reference
    state.push_intermediate_results(std::move(inner_state));
    throw;
  }
}

template <typename Eigenproblem>
void EigensystemSolver<Eigenproblem>::solve_with_method(const std::string& method,
                                                        state_type& state) const {
#define throw_not_compiled_in()                                                         \
  assert_throw(false, ExcInvalidSolverParametersEncountered(                            \
                            "The eigensolver method " + method + "(set via the key '" + \
                            EigensystemSolverKeys::method +                             \
                            "' is not compiled into this version of linalgwrap."))

  if (method == std::string("arpack")) {
#ifdef LINALGWRAP_HAVE_ARPACK
    run_solver<ArpackEigensolver<Eigenproblem>>(state);
#else
    throw_not_compiled_in();
#endif
  } else if (method == std::string("armadillo")) {
#ifdef LINALGWRAP_HAVE_ARMADILLO
    run_solver<ArmadilloEigensolver<Eigenproblem>>(state);
#else
    throw_not_compiled_in();
#endif
  } else {
    // No method is supported!
    assert_throw(false, ExcInvalidSolverParametersEncountered(
                              "The eigensolver method " + method + "(set via the key " +
                              EigensystemSolverKeys::method +
                              ") is unknown. Did you spell it wrong?"));
  }
}

template <typename Eigenproblem>
void EigensystemSolver<Eigenproblem>::solve_state(state_type& state) const {
  assert_dbg(!state.is_failed(), krims::ExcInvalidState("Cannot solve a failed state"));

  if (method != std::string("auto")) {
    solve_with_method(method, state);
  } else if (should_use_arpack(state.eigenproblem())) {
    solve_with_method("arpack", state);
  } else {
    solve_with_method("armadillo", state);
  }
}

}  // namespace linalgwrap
