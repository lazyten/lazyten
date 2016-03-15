#pragma once
#include <rapidcheck.h>
#include <Range.hh>

namespace rc {

template <typename T>
struct RangeWithin {
    typedef T value_type;

    //! Construct a range which at most runs over the interval [min:max)
    static Gen<linalgwrap::Range<value_type>> range_within(value_type min,
                                                           value_type max) {

        auto gen_range = [=] {
            if (min == max) return linalgwrap::Range<value_type>{min, max};

            RC_PRE(min < max);

            auto length = *gen::inRange<value_type>(0, max - min);
            auto start = *gen::inRange<value_type>(min, max - length);

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
Gen<linalgwrap::Range<T>> range_within(T min, T max) {
    return RangeWithin<T>::range_within(min, max);
}
}  // namespace gen
}  // namespace rc
