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
#include "SolverExceptions.hh"
#include "SolverStateBase.hh"
#include <krims/ParameterMap.hh>

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
    void update_control_params(const krims::ParameterMap&) {}
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
    /* Handler which is called before an iteration step is performed
     *
     * The iteration count has already been incremented.
     * */
    virtual void before_iteration_step(state_type&) const {}

    /** Handler which is called once an iteration step finishes
     *
     * This is the last thing called before the convergence and sanity
     * checks.
     * */
    virtual void after_iteration_step(state_type&) const {}

    /** Handler which is called by iterative_solver_assert once the
     * iteration has failed.
     */
    virtual void on_failed(state_type&) const {}
    ///@}

    /** \name Solver building blocks.
     * Various virtual handler functions, which should be called
     * by child classes to perform these tasks.
     */
    ///@{
    /** \brief Start the next iteration step
     *
     * Increase the iteration count and call the before_iteration_step.
     **/
    void start_iteration_step(state_type& s) const;

    /** End the current iteration step.
     *
     *  In principle it checks that we are not beyond
     *   max_iter
     *
     *   \param check_itercount   Enable or disable check for iteration count
     */
    void end_iteration_step(state_type& s) const;
    ///@}
};

//
// --------------------------------------------------------
//

template <typename State>
void SolverBase<State>::start_iteration_step(state_type& s) const {
    // Call the pre step handler:
    before_iteration_step(s);
}

template <typename State>
void SolverBase<State>::end_iteration_step(state_type& s) const {
    // Call the handler:
    after_iteration_step(s);
}

}  // namespace linalgwrap
