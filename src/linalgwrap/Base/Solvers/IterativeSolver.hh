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
#include "IterativeSolverKeys.hh"
#include "IterativeSolverState.hh"
#include "SolverExceptions.hh"
#include <krims/ParameterMap.hh>

namespace linalgwrap {

/** Extra structure the state of a truely iterative solver
 *  needs to be able to do on top of another Base class. */
template <typename Base>
class IterativeSolver : public Base {
 public:
  typedef Base base_type;
  typedef typename base_type::state_type state_type;
  typedef typename state_type::count_type count_type;

  static_assert(std::is_base_of<linalgwrap::IterativeSolverState, state_type>::value,
                "Base's state_type needs to be derived off "
                "linalgwrap::IterationState");

  /** \name Iteration control */
  ///@{
  /** \brief Function to check convergence before a new iteration begins.
   *  This default implementation does no checking and just returns
   *  false.
   */
  virtual bool is_converged(const state_type&) const { return false; }

  /** Maximum number of iterations
   *
   * The iteration should be considered as failed
   * once we go beyond this number of iterations.
   **/
  count_type max_iter = 100;

  /** Bulk-update control parameters from a parameter map.
   *
   * For the list of available keys, see IterativeSolverKeys.hh
   */
  void update_control_params(const krims::ParameterMap& map) {
    base_type::update_control_params(map);
    max_iter = map.at(IterativeSolverKeys::max_iter, max_iter);
  }
  ///@}

 protected:
  /** \name Handler functions
   * Various virtual handler functions, which are called when
   * certain events happen.
   */
  ///@{
  /** Handler which is called once the iteration has converged
   *
   * \note This may not mark the end of the iteration yet.
   */
  virtual void on_converged(state_type&) const {}
  ///@}

  /** \name Solver building blocks.
   * Various handler functions, which should be called
   * by child classes to perform these tasks.
   */
  ///@{
  /** \brief Start the next iteration step
   *
   * Increase the iteration count.
   **/
  void start_iteration_step(state_type& s) const;

  /** \brief Use is_converged to check whether convergence
   * has been achieved.
   *
   * Return true if yes, else false
   * Also call on_convergence if convergence.
   **/
  bool convergence_reached(state_type& s) const;
  ///@}
};

template <typename State>
void IterativeSolver<State>::start_iteration_step(state_type& s) const {
  // Assert that we are not beyond the iteration count already
  solver_assert(s.n_iter() < max_iter, s, ExcMaximumNumberOfIterationsReached(max_iter));

  s.increase_iteration_count();
  base_type::start_iteration_step(s);
}

template <typename State>
bool IterativeSolver<State>::convergence_reached(state_type& s) const {
  if (is_converged(s)) {
    on_converged(s);
    return true;
  }
  return false;
}

}  // namespace linalgwrap
