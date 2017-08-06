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
#include "NumericConstants.hh"
#include <cstddef>
#include <rapidcheck.h>

namespace lazyten {
namespace gen {

/** \brief Generator for a numeric size.
 *
 * Dependant on the dimensionality the size becomes smaller.
 * The precise values are compiled in and are determined from
 * experience.
 *
 * \tparam dim  The dimensionality of the object.
 */
template <unsigned int dim>
rc::Gen<size_t> numeric_size() {
  static_assert(dim > 0, "dimension has to be larger than 0");

  constexpr size_t max =
        (dim == 1) ? 80                                           // dim==1 -> 80
                   : ((dim == 2) ? 10                             // dim==2 -> 10
                                 : ((dim == 3) ? 4                // dim==3 -> 4
                                               : ((dim == 4) ? 3  // dim==4 -> 3
                                                             : 2  // else -> 2
                                                  )));

  // rhs of inRange is exclusive
  return rc::gen::inRange<size_t>(1u, max + 1);
  // TODO this really should be 0,max+1
  // return rc::gen::inRange<size_t>(0u, max + 1);
}
}  // namespace gen
}  // namespace lazyten
