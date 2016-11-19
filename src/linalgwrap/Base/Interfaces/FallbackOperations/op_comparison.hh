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

/** \name relational operators */
///@{
//! Check whether two Indexables are equal
template <typename Indexable1, typename Indexable2>
typename std::enable_if<IsIndexable<Indexable1>::value &&
                              IsIndexable<Indexable2>::value,
                        bool>::type
operator==(const Indexable1& A, const Indexable2& B) {
    linalgwrap_called_fallback();
    if (A.n_elem() != B.n_elem()) return false;
    return std::equal(A.begin(), A.end(), B.begin());
}

//! Check whether two Indexables are not equal
template <typename Indexable1, typename Indexable2>
typename std::enable_if<IsIndexable<Indexable1>::value &&
                              IsIndexable<Indexable2>::value,
                        bool>::type
operator!=(const Indexable1& A, const Indexable2& B) {
    linalgwrap_called_fallback();
    return !operator==(A, B);
}
///@}

}  // namespace linalgwrap
