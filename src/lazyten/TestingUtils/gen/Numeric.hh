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
#include "NumericConstants.hh"
#include <complex>
#include <krims/TypeUtils.hh>
#include <limits>
#include <rapidcheck.h>

namespace lazyten {
namespace gen {

namespace detail {
// TODO Have implementations for Matrices and Vectors here as well?

/** Default implementation: Not allowed */
template <typename Value, typename Enable = void>
struct Numeric {
  static_assert(!std::is_same<Enable, void>::value,
                "Numeric is only available for floating point types.");
};

/** Construct a numeric stored vector */
template <typename Scalar>
struct Numeric<Scalar,
               typename std::enable_if<std::is_floating_point<Scalar>::value>::type> {
  typedef Scalar scalar_type;
  static rc::Gen<Scalar> numeric();
};

/** Construct a numeric complex value */
template <typename Real>
struct Numeric<std::complex<Real>> {
  typedef std::complex<Real> scalar_type;
  static rc::Gen<scalar_type> numeric() {
    return rc::gen::construct<scalar_type>(Numeric<Real>::numeric(),
                                           Numeric<Real>::numeric());
  }
};
}  // namespace detail

/** \brief Generator for a numeric value.
 *
 * Honours the bounds max_n_elem, max_value and min_value. This means all
 * 1D-objects contain no more than 100 elements and all 2D objects have a
 * smaller dimensionality than 10x10. All entries are in the range
 * [min_value,max_value].
 */
template <typename Value>
rc::Gen<Value> numeric() {
  return detail::Numeric<Value>::numeric();
}

//
// -------------------------------------------------------------------
//

template <typename Scalar>
rc::Gen<Scalar> detail::Numeric<
      Scalar,
      typename std::enable_if<std::is_floating_point<Scalar>::value>::type>::numeric() {
  // TODO better do this in terms of 2-based numbers
  //      and convert them to decimal later.

  // Define a lambda that returns a new matrix element:
  auto gen_element = [=] {
    typedef long gen_type;

    // Generate an arbitrary value
    const gen_type gen_value = *rc::gen::arbitrary<gen_type>();

    // Bring it to scale between -1 and 1:
    const scalar_type ratio =
          static_cast<scalar_type>(gen_value) / std::numeric_limits<gen_type>::max();

    // Compute value and return it if large enough, else 0.
    const scalar_type value = ratio * max_value;
    if (std::fabs(value) < min_value) {
      return static_cast<scalar_type>(0.);
    } else {
      return value;
    }
  };

  return rc::gen::exec(gen_element);
}

}  // namespace gen
}  // namespace lazyten
