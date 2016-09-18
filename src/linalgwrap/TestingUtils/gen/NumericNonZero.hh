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
 * Honours the bounds max_n_elem, max_value and min_value. This means all
 * 1D-objects contain no more than 100 elements and all 2D objects have a
 * smaller dimensionality than 10x10. All entries are in the range
 * [min_value,max_value].
 */
template <typename Value>
rc::Gen<Value> numeric_nonZero() {
    return rc::gen::distinctFrom(gen::numeric<Value>(), 0.0);
}
}  // gen
}  // linalgwrap
