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
#include "IterativeStateWrapper.hh"
#include "IterativeWrapperKeys.hh"
#include "SolverBase.hh"
#include "SolverExceptions.hh"
#include <krims/GenMap.hh>

namespace lazyten {

/** Extra structure the state of a truely iterative solver
 *  needs to be able to do on top of another Base class. */
template <typename Base>
class IterativeWrapper : public Base {
 public:
  typedef Base base_type;
  typedef typename base_type::state_type state_type;

  static_assert(std::is_base_of<IterativeStateWrapper<typename state_type::it_base_type>,
                                state_type>::value,
                "Base's state_type needs to be derived off "
                "lazyten::IterationState");

  static_assert(std::is_base_of<SolverBase<state_type>, Base>::value,
                "Base needs to be derived off SoverBase");

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
  size_t max_iter = 100;

  /** Bulk-update control parameters from a parameter map.
   *
   * For the list of available keys, see IterativeSolverKeys.hh
   */
  void update_control_params(const krims::GenMap& map) {
    base_type::update_control_params(map);
    max_iter = map.at(IterativeWrapperKeys::max_iter, max_iter);
  }

  /** Get the current settings of all internal control parameters and
   *  update the GenMap accordingly.
   */
  void get_control_params(krims::GenMap& map) const {
    base_type::get_control_params(map);
    map.update(IterativeWrapperKeys::max_iter, max_iter);
  }
  ///@}

 protected:
  /** \name Handler functions
   * Various virtual handler functions, which are called when
   * certain events happen.
   */
  ///@{
  /* Handler which is called before an iteration step is performed
   *
   * The iteration count has already been incremented.
   * */
  virtual void before_iteration_step(state_type&) const {}

  /** Handler which is called once the iteration has converged
   *
   * \note This may not mark the end of the iteration yet.
   */
  virtual void on_converged(state_type&) const {}

  /** Handler which is called once an iteration step finishes
   *
   * This is the last thing called before the convergence and sanity
   * checks.
   * */
  virtual void after_iteration_step(state_type&) const {}
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

  /** End the current iteration step.
   *
   *  In principle it checks that we are not beyond
   *   max_iter
   *
   *   \param check_itercount   Enable or disable check for iteration count
   */
  void end_iteration_step(state_type& s) const { after_iteration_step(s); }
  ///@}
};

template <typename Base>
void IterativeWrapper<Base>::start_iteration_step(state_type& s) const {
  // Assert that we are not beyond the iteration count already
  solver_assert(s.n_iter() < max_iter, s, ExcMaximumNumberOfIterationsReached(max_iter));

  s.increase_iteration_count();
  before_iteration_step(s);
}

template <typename Base>
bool IterativeWrapper<Base>::convergence_reached(state_type& s) const {
  if (is_converged(s)) {
    on_converged(s);
    return true;
  }
  return false;
}

}  // namespace lazyten
