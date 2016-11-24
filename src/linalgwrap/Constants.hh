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
#include <limits>
#include <type_traits>

namespace linalgwrap {

/** Struct to hold constants for various data types */
template <typename T>
struct Constants {
  typedef T type;

  //! Constant representing 0
  static constexpr type zero = type(0.);

  //! Constant representing 1
  static constexpr type one = type(1.);

  /** \brief Constant which should be used for invalid values
   *
   * \note This value may well be a valid ``type`` value, so
   * do not check against it in order to see if the ``type`` value
   * is invalid.
   */
  static constexpr type invalid = std::numeric_limits<type>::has_signaling_NaN
                                        ? std::numeric_limits<type>::signaling_NaN()
                                        : std::numeric_limits<type>::max();

  /** \brief Constant which should be used to signal
   *  that all available eigenvalues, vectors, whatever should be used
   *  or computed.
   */
  static constexpr type all = invalid;

  /** \brief Constant which gives the default numerical tolerance to use
   *  for a number of checks */
  static constexpr type default_tolerance = std::numeric_limits<type>::epsilon();
};

template <typename T>
constexpr typename Constants<T>::type Constants<T>::zero;

template <typename T>
constexpr typename Constants<T>::type Constants<T>::one;

template <typename T>
constexpr typename Constants<T>::type Constants<T>::invalid;

template <typename T>
constexpr typename Constants<T>::type Constants<T>::all;

template <typename T>
constexpr typename Constants<T>::type Constants<T>::default_tolerance;
}
