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

namespace lazyten {
/** Marker interface to assert that the objects derived off this class are
 * stored, i.e. that they manage their own storage allocation and deallocation
 * internally and sit fully in memory.
 *
 * This has certain implications on the type of constructors which are expected
 * by the library. See for example IsStoredVector.hh for details.
 */
struct Stored_i {};
}
