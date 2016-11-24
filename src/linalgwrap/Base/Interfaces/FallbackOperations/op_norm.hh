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
#include "linalgwrap/Base/Interfaces/Vector_i.hh"
#include "macro_defs.hh"
#include "op_accumulate.hh"
#include "op_dot.hh"
#include "op_minmax.hh"

namespace linalgwrap {

/** Compute vector norms */
///@{
/** Calculate the l1 norm of the vector (sum of abs values of elements) */
template <typename Vector,
          typename std::enable_if<IsVector<Vector>::value, int>::type = 0>
typename Vector::real_type norm_l1(const Vector& v) {
  linalgwrap_called_fallback();
  return accumulate(abs(v));
}

/** Calculate the linf norm of the vector (abs. largest element) */
template <typename Vector,
          typename std::enable_if<IsVector<Vector>::value, int>::type = 0>
typename Vector::real_type norm_linf(const Vector& v) {
  linalgwrap_called_fallback();
  return max(abs(v));
}

/** Calculate the l2 norm squared of the vector. */
template <typename Vector,
          typename std::enable_if<IsVector<Vector>::value, int>::type = 0>
typename Vector::real_type norm_l2_squared(const Vector& v) {
  linalgwrap_called_fallback();
  return std::real(cdot(v, v));
}

/** Calculate the l2 norm of the vector */
template <typename Vector,
          typename std::enable_if<IsVector<Vector>::value, int>::type = 0>
typename Vector::real_type norm_l2(const Vector& v) {
  linalgwrap_called_fallback();
  return std::sqrt(norm_l2_squared(v));
}
///@}

}  // namespace linalgwrap
