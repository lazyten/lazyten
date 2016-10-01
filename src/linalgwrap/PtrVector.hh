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
#include "IteratorVector.hh"

namespace linalgwrap {

/** Using statement for a vector that operates directly on raw memory
 *
 * This is is a small wrapper interface to give a vector-like interface
 * to any stride of memory represented by pointers.
 * */
template <typename T>
using PtrVector = IteratorVector<T*>;
}  // namespace linalgwrap
