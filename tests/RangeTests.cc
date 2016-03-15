#include <rapidcheck.h>
#include <catch.hpp>
#include <Range.hh>
#include <type_traits>

namespace linalgwrap {
namespace tests {
using namespace rc;

namespace range_tests {
template <typename T>
struct RangeTests {
    typedef Range<T> range_type;
    typedef typename range_type::size_type size_type;

    static void construction(T t1, T t2) {
        bool is_valid = t1 <= t2;
        RC_CLASSIFY(is_valid, "Range is valid");

        if (is_valid) {
            range_type r{t1, t2};
            RC_ASSERT(r.length() == static_cast<size_type>(t2 - t1));
        } else {
#ifdef DEBUG
            RC_ASSERT_THROWS_AS((Range<T>{t1, t2}), ExcBelowLowerBound<T>);
#endif
        }
    }

    static void element_access() {
        const auto length = *gen::positive<T>();

        RC_CLASSIFY(length == 0, "Zero range");

        // Start value
        const auto start = *gen::arbitrary<T>().as("Start index");

        // Exclude overflow cases:
        RC_PRE(start + length > length);

        // Construct the range:
        range_type r = range(start, start + length);

        // Access arbitrary element:
        const auto acc = *gen::inRange<T>(0, r.length());

        if (r.is_empty()) {
#ifdef DEBUG
            typedef typename range_type::ExcEmptyRange exc_type;
            RC_ASSERT_THROWS_AS(r[acc], exc_type);
#endif
        } else {
            RC_ASSERT(r[acc] == start + acc);
        }
    }

    static void access_to_empty_range(T value) {
        // Generate an empty range
        range_type r{{value, value}};

        // it should be empty
        RC_ASSERT(r.is_empty());

        // length should be zero
        RC_ASSERT(r.length() == size_type{0});
        RC_ASSERT(r.size() == size_type{0});

        // begin and and should be identical.
        RC_ASSERT(r.begin() == r.end());

#ifdef DEBUG
        // check that element access throws:

        typedef typename range_type::ExcEmptyRange exc_type;
        RC_ASSERT_THROWS_AS(r[0], exc_type);
#endif
    }

    static void iteration() {
        const auto length = *gen::inRange<T>(0, 101);

        RC_CLASSIFY(length == 0, "Zero range");

        // Start value
        const auto start = *gen::arbitrary<T>().as("Start index");

        // Exclude overflow cases:
        RC_PRE(start + length > length);

        // Construct the range:
        range_type r = range(start, start + length);

        // Iterator variable:
        auto i = start;
        typename range_type::size_type counter{0};

        // Use range-based iterator:
        for (const auto val : r) {
            RC_ASSERT(i == val);
            ++i;
            ++counter;
        }

        RC_ASSERT(counter == r.length());
    }

    static void access_to_past_the_end_iterator() {
        range_type r{0, 1};
        auto it = std::end(r);

#ifdef DEBUG
        CHECK_THROWS_AS(*it, ExcIteratorPastEnd);
        CHECK_THROWS_AS(++it, ExcIteratorPastEnd);
        CHECK_THROWS_AS(it++, ExcIteratorPastEnd);
#endif
    }
};

}  // Range tests

TEST_CASE("Range tests", "[range]") {
    using namespace range_tests;

    // Make sure that the program does not get aborted
    exceptions::assert_dbg_effect = exceptions::ExceptionEffect::THROW;

    SECTION("Range construction") {
        REQUIRE(rc::check("Range construction (size_t)",
                          RangeTests<size_t>::construction));
        REQUIRE(rc::check("Range construction (int)",
                          RangeTests<int>::construction));
    }

    SECTION("Element access") {
        REQUIRE(rc::check("Element access (size_t)",
                          RangeTests<size_t>::element_access));
        REQUIRE(rc::check("Element access (int)",
                          RangeTests<int>::element_access));
    }

    SECTION("Empty range properties") {
        REQUIRE(rc::check("Empty range properties (size_t)",
                          RangeTests<size_t>::access_to_empty_range));
        REQUIRE(rc::check("Empty range properties (int)",
                          RangeTests<int>::access_to_empty_range));
    }

    SECTION("Iteraton") {
        REQUIRE(rc::check("Iteration (size_t)", RangeTests<size_t>::iteration));
        REQUIRE(rc::check("Iteration (int)", RangeTests<int>::iteration));
    }

    SECTION("Past-the-end-iterator properties (size_t)") {
        RangeTests<size_t>::access_to_past_the_end_iterator();
    }

    SECTION("Past-the-end-iterator properties (int)") {
        RangeTests<int>::access_to_past_the_end_iterator();
    }

}  // TEST_CASE
}  // namespace test
}  // namespace linalgwrap
