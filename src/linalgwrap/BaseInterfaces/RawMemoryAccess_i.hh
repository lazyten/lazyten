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
#include <type_traits>

namespace linalgwrap {
/** Marker interface to assert that the objects derived off this class are
 *  stored in a continuous stride of memory which can be accessed and modified.
 *
 * This implies that all objects implementing this interfaces have the operation
 * memptr() which either returns a const or non-const pointer to the memory
 * (in whichever layout it may be).
 */
struct RawMemoryAccess_i {};

/** Struct representing a true_type if the class' storage is available
 *  via memptr, else false_type */
template <typename T>
struct IsRawMemoryAccessible : public std::is_base_of<RawMemoryAccess_i, T> {};
}
