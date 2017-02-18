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
#include "SolverExceptions.hh"
#include "SolverStateBase.hh"
#include <krims/GenMap.hh>

namespace linalgwrap {
template <typename State>
class SolverBase {
  static_assert(std::is_base_of<linalgwrap::SolverStateBase, State>::value,
                "State needs to be derived off linalgwrap::SolverStateBase");

 public:
  typedef State state_type;

  /** \name Iteration control */
  ///@{
  /** Bulk-update control parameters from a parameter map.
   *
   * This version does nothing since there are no parameters*/
  void update_control_params(const krims::GenMap&) {}

  /** Get the current settings of all internal control parameters and
   *  update the GenMap accordingly.
   */
  void get_control_params(krims::GenMap&) const {}
  ///@}

  /** \brief Run the solver by advancing the provided
   * state until convergence has been reached.
   *
   * If the solver fails, i.e. does not manage to converge the
   * provided state, an SolverBaseException is thrown
   * explaining the problem.
   * The fail bit of the state will also be set and the
   * fail message will also provide useful information
   * why the solver failed.
   */
  virtual void solve_state(state_type& state) const = 0;

  /** Defaults for the big five */
  //@{
  virtual ~SolverBase() = default;
  SolverBase(const SolverBase&) = default;
  SolverBase(SolverBase&&) = default;
  SolverBase& operator=(const SolverBase&) = default;
  SolverBase& operator=(SolverBase&&) = default;
  SolverBase() = default;
  //@}

 protected:
  /** \name Handler functions
   * Various virtual handler functions, which are called when
   * certain events happen.
   */
  ///@{
  /** Handler which is called by iterative_solver_assert once the
   * iteration has failed.
   */
  virtual void on_failed(state_type&) const {}
  ///@}
};
}  // namespace linalgwrap
