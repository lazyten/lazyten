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
#include <string>

namespace lazyten {
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
   *
   * \note Non-iterative solvers return exactly 1 here
   **/
  virtual size_t n_iter() const { return 1; }

  /** Return the number of problem matrix applies the solver
   *  required to solve the problem
   *  \note Dense solvers work on the memory and return exactly 0
   **/
  virtual size_t n_mtx_applies() const { return 0; }

 protected:
  //! Has the iteration failed?
  bool m_failed;

  //! Why has it failed?
  std::string m_fail_reason;
};

}  // namespace lazyten
