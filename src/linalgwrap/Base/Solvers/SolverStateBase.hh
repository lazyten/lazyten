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
#include <string>

namespace linalgwrap {
class SolverStateBase {
 public:
  /** \brief Default constructor
   *
   * The failbit is unset and the iteration count is 0.
   * */
  SolverStateBase() : m_failed{false}, m_fail_reason{} {}

  virtual ~SolverStateBase() = default;
  SolverStateBase(const SolverStateBase&) = default;
  SolverStateBase(SolverStateBase&&) = default;
  SolverStateBase& operator=(const SolverStateBase&) = default;
  SolverStateBase& operator=(SolverStateBase&&) = default;

  /** Fail the iteraton and specify a reason why */
  void fail(std::string reason);

  /** Clear the failed flag and clear the fail reason. */
  void clear_failed();

  /** Is the failed bit set? */
  bool is_failed() const { return m_failed; }

  /** Get the fail reason
   *
   * \note Returns an empty string in case that the iteration has not failed.
   * An empty string does, however, not imply that the iteration has not
   * failed.
   */
  const std::string& fail_reason() const { return m_fail_reason; }

  /** Return the number of iterations needed up to this point
   *  (for non-iterative solvers this should return exactly 1 */
  virtual size_t n_iter() const { return 1; }

  /** Setup the guess of this state. from another state.
   * In this case this function does nothing */
  void obtain_guess_from(const SolverStateBase&) {}

 protected:
  //! Has the iteration failed?
  bool m_failed;

  //! Why has it failed?
  std::string m_fail_reason;
};

}  // namespace linalgwrap
