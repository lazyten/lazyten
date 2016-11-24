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
#include "linalgwrap/Base/Interfaces/MutableMemoryVector_i.hh"
#include "linalgwrap/Constants.hh"
#include <iterator>
namespace linalgwrap {
namespace detail {

/** Function which scales each element of an object, unless
 * the factor is zero, in this case the elements are explicitly
 * set to zero (to avoid nans)
 */
template <typename Object, typename Scalar = typename Object::scalar_type>
void scale_or_set(Object& o, Scalar s) {
  assert_finite(s);
  if (s == Constants<Scalar>::zero) {
    o.set_zero();
  } else {
    o *= s;
  }
}

/** Specialisation of above */
template <typename Scalar>
void scale_or_set(MutableMemoryVector_i<Scalar>& o, Scalar s) {
  assert_finite(s);
  if (s == Constants<Scalar>::zero) {
    o.set_zero();
  } else {
    std::transform(std::begin(o), std::end(o), std::begin(o),
                   [s](Scalar& elem) { return s * elem; });
  }
}

}  // namespace detail
}  // namespace linalgwrap
