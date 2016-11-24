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

/** Extra structure the state of a truely iterative solver
 *  needs to satisfy on top of SolverStateBase. */
class IterativeSolverState {
 public:
  //! Type for counting iteration numbers.
  typedef size_t count_type;

  /** \brief Default constructor
   *
   * The failbit is unset and the iteration count is 0.
   * */
  IterativeSolverState() : m_count{0} {}

  /** Increase the iteration count and return the new count */
  count_type increase_iteration_count() { return ++m_count; }

  /** Return the current iteration count */
  count_type n_iter() const { return m_count; }

 protected:
  count_type m_count;
};

}  // namespace linalgwrap
