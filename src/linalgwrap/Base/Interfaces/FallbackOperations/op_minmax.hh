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
#include "linalgwrap/Base/Interfaces/Indexable_i.hh"
#include "macro_defs.hh"

namespace linalgwrap {

/** Compute the minimum of all values of an indexable object */
template <typename Indexable>
ValidIndexableScalarT<Indexable> min(const Indexable& i) {
    linalgwrap_called_fallback();
    return *std::min_element(i.begin(), i.end());
}

/** Compute the maximum of all values of an indexable object */
template <typename Indexable>
ValidIndexableScalarT<Indexable> max(const Indexable& i) {
    linalgwrap_called_fallback();
    return *std::max_element(i.begin(), i.end());
}

}  // namespace linalgwrap
