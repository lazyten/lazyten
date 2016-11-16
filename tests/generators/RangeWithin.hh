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
#include <krims/Range.hh>
#include <rapidcheck.h>
#include <sstream>

namespace rc {

template <typename T>
struct RangeWithin {
    typedef T value_type;

    //! Construct a range which at most runs over the interval [min:max)
    static Gen<linalgwrap::Range<value_type>> range_within(
          value_type min, value_type max, value_type min_length = 0) {

        auto gen_range = [=] {
            if (min > max) {
                std::stringstream ss;
                ss << "Minimum " << min << " not smaller or equal to maximum "
                   << max;
                throw rc::GenerationFailure(ss.str());
            }

            if (min_length > max - min) {
                std::stringstream ss;
                ss << "Minimal length " << min_length
                   << " must be smaller or equal to the difference between max "
                   << "and min, which is " << max - min;
                throw rc::GenerationFailure(ss.str());
            }

            if (min_length == max - min) {
                return krims::Range<value_type>{min, max};
            }

            auto length = *gen::inRange<value_type>(min_length, max - min + 1);
            auto start = *gen::scale(
                  2.0, gen::inRange<value_type>(min, max - length + 1));
            if (start + length > max) {
                std::stringstream ss;
                ss << "Something went wrong: Too large maximum "
                   << start + length << " where desired max was " << max;
                throw rc::GenerationFailure(ss.str());
            }
            return krims::Range<value_type>{start, start + length};
        };
        return gen::exec(gen_range);
    }
};

/** Return a range within the half-open interval [min:max),
 *  i.e. a range representing a half-open subinterval
 *  within said [min:max)
 */
namespace gen {
template <typename T>
Gen<krims::Range<T>> range_within(T min, T max, T min_length = 0) {
    return RangeWithin<T>::range_within(min, max, min_length);
}
}  // namespace gen
}  // namespace rc
