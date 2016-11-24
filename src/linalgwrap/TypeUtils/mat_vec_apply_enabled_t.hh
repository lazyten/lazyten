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
#include "linalgwrap/Base/Interfaces.hh"

namespace linalgwrap {

/** Using statement which selectively enables/disables the generic form of
 *  Matrix-Vector applies.
 *
 *  The apply is only enabled if the Vector has the same scalar type and is
 *  a mutable memory vector
 */
template <typename Matrix, typename VectorIn, typename VectorOut>
using mat_vec_apply_enabled_t =
      typename std::enable_if<IsMutableMemoryVector<VectorIn>::value &&
                              IsMutableMemoryVector<VectorOut>::value &&
                              std::is_same<typename VectorIn::scalar_type,
                                           typename VectorOut::scalar_type>::value &&
                              std::is_same<typename VectorIn::scalar_type,
                                           typename Matrix::scalar_type>::value>::type;
}  // namespace linalgwrap
