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
#include "Numeric.hh"

namespace linalgwrap {
namespace gen {

/** \brief Generator for a numeric value, which is not zero.
 *
 * All entries are in the range [min_value,max_value] and
 * centered around zero, i.e. the shrink towards *small values**
 * which are not zero. If this is not desired, consider
 * numeric_around(value), which shrinks towards the given value.
 **/
template <typename Value>
rc::Gen<Value> numeric_nonZero() {
    return rc::gen::distinctFrom(gen::numeric<Value>(), 0.0);
}

}  // gen
}  // linalgwrap
