#pragma once
#include <linalgwrap/Range.hh>
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
                return linalgwrap::Range<value_type>{min, max};
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
            return linalgwrap::Range<value_type>{start, start + length};
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
Gen<linalgwrap::Range<T>> range_within(T min, T max, T min_length = 0) {
    return RangeWithin<T>::range_within(min, max, min_length);
}
}  // namespace gen
}  // namespace rc
