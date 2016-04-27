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

#ifndef LIBLINALG_CONSTANTS_HH_
#define LIBLINALG_CONSTANTS_HH_
#include <limits>

namespace linalgwrap {

template <typename Scalar>
struct Constants {
    typedef Scalar scalar_type;

    //! Constant representing 0
    static constexpr scalar_type zero = scalar_type(0.);

    //! Constant representing 1
    static constexpr scalar_type one = scalar_type(1.);

    /** \brief Constant which should be used for invalid values
     *
     * \note This value may well be a valid ``scalar_type`` value, so
     * do not check against it in order to see if the ``scalar_type`` value
     * is invalid.
     */
    static constexpr scalar_type invalid =
          std::numeric_limits<scalar_type>::has_signaling_NaN
                ? std::numeric_limits<scalar_type>::signaling_NaN()
                : std::numeric_limits<scalar_type>::max();

    /** \brief Constant which gives the default numerical tolerance to use
     *  for a number of checks */
    static constexpr scalar_type default_tolerance =
          std::numeric_limits<scalar_type>::epsilon();
};
}

#endif
