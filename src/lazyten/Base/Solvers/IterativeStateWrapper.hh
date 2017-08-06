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

/** Extra structure the state of a truely iterative solver
 *  needs to satisfy on top of another Base class */
template <typename Base>
class IterativeStateWrapper : public Base {
 public:
  typedef Base it_base_type;

  /** \brief Default constructor */
  IterativeStateWrapper(Base&& base) : it_base_type(std::move(base)), m_count{0} {}

  /** Increase the iteration count and return the new count */
  size_t increase_iteration_count() { return ++m_count; }

  /** Return the current iteration count
   * \note Shadows this function in Base
   **/
  size_t n_iter() const override { return m_count; }

 protected:
  size_t m_count;
};

}  // namespace lazyten
