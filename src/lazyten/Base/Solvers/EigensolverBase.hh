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
#include "EigensolverBaseKeys.hh"
#include "EigensolverStateBase.hh"
#include "SolverBase.hh"

namespace lazyten {

/** Each eigensolver should support the following:
 *
 * - Default construction
 * - Construction from parameter map
 * - A function update_control_params to update the control parameters.
 *
 */
template <typename State>
class EigensolverBase : public SolverBase<State> {
  static_assert(IsEigensolverState<State>::value,
                "State needs to be a type derived from EigensolverStateBase");

 public:
  typedef SolverBase<State> base_type;
  typedef typename base_type::state_type state_type;

  /** \name Types forwarded from the state */
  ///@{
  /** The type of the eigenproblem */
  typedef typename state_type::eproblem_type eproblem_type;

  /** The stored matrix type of the eigenproblem */
  typedef typename state_type::stored_matrix_type stored_matrix_type;

  /** The size type of the eigenproblem */
  typedef typename state_type::size_type size_type;

  /** The scalar type of the eigenproblem
   *  (but not necessarily of the solution) */
  typedef typename state_type::scalar_type scalar_type;

  /** The real type of the eigenproblem */
  typedef typename state_type::real_type real_type;

  /** The type used for the eigenvalues */
  typedef typename state_type::evalue_type evalue_type;

  /** The type used for the eigenvectors */
  typedef typename state_type::evector_type evector_type;

  /** The type used for the eigensolution */
  typedef typename state_type::esoln_type esoln_type;
  ///@}

  /** \name Iteration control */
  ///@{
  /** Which eigenvalues one is interested.
   *
   * Strings, which are accepted by all Eigensolvers:
   *   - "SM"   Smallest magnitude
   *   - "LM"   Largest magnitude
   *   - "SR"   Smallest real
   *   - "LR"   Largest real
   *   - "SI"   Smallest imaginary
   *            (only for complex scalar types)
   *   - "LI"   Largest imaginary
   *            (only for complex scalar types)
   */
  std::string which{"SR"};

  /** The convergence tolerance to use within the eigensolver */
  real_type tolerance = Constants<real_type>::default_tolerance;

  /** Bulk-update control parameters from a parameter map.
   *
   * For the list of available keys, see EigensolverBaseKeys.hh
   */
  void update_control_params(const krims::GenMap& map) {
    base_type::update_control_params(map);
    which = map.at(EigensolverBaseKeys::which, which);
    tolerance = map.at(EigensolverBaseKeys::tolerance, tolerance);
  }

  /** Get the current settings of all internal control parameters and
   *  update the GenMap accordingly.
   */
  void get_control_params(krims::GenMap& map) const {
    base_type::get_control_params(map);
    map.update(EigensolverBaseKeys::which, which);
    map.update(EigensolverBaseKeys::tolerance, tolerance);
  }
  ///@}

  /** \name Run an Eigensolver
   */
  ///@{

  /** \brief Run the eigensolver with the provided eigenproblem and return
   * the final state.
   *
   * If the solver does not manage to achieve convergence a
   * SolverException
   * is thrown and the state's fail_bit will be set accompanied with
   * an appropriate fail message.
   */
  virtual state_type solve(const eproblem_type problem) const {
    // Note: eproblem_type only contains pointers or references,
    // so copying is ok.
    state_type state{std::move(problem)};
    this->solve_state(state);
    return state;
  }

  /** \brief Run the solver on a problem, starting from a guess state
   *
   * Here we use the EigensolverStateBase in order to be able to use
   * states of potentially different state_type as well.
   */
  template <typename GuessState,
            typename = krims::enable_if_t<
                  std::is_base_of<EigensolverStateBase<eproblem_type>,
                                  krims::remove_reference_t<GuessState>>::value>>
  state_type solve_with_guess(const eproblem_type problem,
                              GuessState&& guess_state) const {
    // Create a new state and install the guess state:
    state_type state{std::move(problem)};
    state.obtain_guess_from(std::forward<GuessState>(guess_state));
    this->solve_state(state);
    return state;
  }

  /** \brief Run the solver on a problem, starting from a guess state
   *
   * Here we use the EigensolverStateBase in order to be able to use
   * states of potentially different state_type as well.
   */
  template <typename GuessState,
            typename = krims::enable_if_t<
                  std::is_base_of<EigensolverStateBase<eproblem_type>,
                                  krims::remove_reference_t<GuessState>>::value>>
  state_type solve_with_guess(GuessState&& guess_state) const {
    return solve_with_guess(guess_state.problem(), std::forward<GuessState>(guess_state));
  }

  ///@}
};

}  // namespace lazyten
