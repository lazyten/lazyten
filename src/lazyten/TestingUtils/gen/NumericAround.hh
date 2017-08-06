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
#include "Numeric.hh"

namespace lazyten {
namespace gen {

/** \brief Generator for a numeric value, which is distributed
 * around a provided value.
 *
 *
 * All entries are in the range [min_value,max_value]+value and
 * centered around zero, i.e. the shrink towards **small values**
 * which are not zero. If this is not desired, consider
 * numeric_around(value), which shrinks towards the given value.
 */
template <typename Value>
rc::Gen<Value> numeric_around(Value v) {
  return rc::gen::map(gen::numeric<Value>(), [v](Value x) { return x + v; });
}

}  // gen
}  // namespace lazyten
