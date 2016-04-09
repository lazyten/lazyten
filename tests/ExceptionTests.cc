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

#include <catch.hpp>
#include <linalgwrap/Exceptions.hh>
#include <rapidcheck.h>

namespace linalgwrap {
namespace tests {
using namespace rc;

TEST_CASE("Exception system", "[exception]") {
    // Make sure that the program does not get aborted
    AssertDbgEffect::set(ExceptionEffect::THROW);

    SECTION("Test throw or raise mechanism") {
#ifdef DEBUG
        // Test if an ExcNotImplemented exception would be thrown:
        REQUIRE_THROWS_AS(assert_dbg(false, ExcNotImplemented()),
                          ExcNotImplemented);

        // Test if an ExcIO exception would be thrown:
        REQUIRE_THROWS_AS(assert_dbg(false, ExcIO()), ExcIO);
#else
        // Test that indeed nothing happens here:
        REQUIRE_NOTHROW(assert_dbg(false, ExcNotImplemented()));
        REQUIRE_NOTHROW(assert_dbg(false, ExcIO()));
#endif

        // Test if an ExcNotImplemented exception would be thrown:
        REQUIRE_THROWS_AS(assert_throw(false, ExcNotImplemented()),
                          ExcNotImplemented);

        // Test if an ExcIO exception would be thrown:
        REQUIRE_THROWS_AS(assert_throw(false, ExcIO()), ExcIO);
    }

    //
    // ----------------------------------------------------------------
    //

    SECTION("Test derived assertion macros") {
        // TODO template this in the type of value and bound
        //      use a testing namespace above for that
        auto test_assert_range = [](long lower_bound, long value) {
            // Arbitrary size:
            auto size = *gen::positive<decltype(value)>();

            // Calculate upper bound
            auto upper_bound = lower_bound + size + 1;

            // Should this assertion fail?
            const bool should_catch_something =
                  (value < lower_bound) || (upper_bound <= value);

            // Classify according to upper property
            RC_CLASSIFY(should_catch_something, "Assertion failed");

#ifdef DEBUG
            if (should_catch_something) {
                // sometimes we catch something in debug mode.
                RC_ASSERT_THROWS_AS(
                      assert_range(lower_bound, value, upper_bound),
                      ExcOutsideRange<decltype(value)>);
            } else {
                // Assert the lower bound. If error, throw
                assert_range(lower_bound, value, upper_bound);
            }
#else
            // we should never catch anything in release mode.
            assert_range(lower_bound, value, upper_bound);
#endif
        };
        CHECK(rc::check("Test assert_range", test_assert_range));

        //
        // ---------------------------------------------------------
        //

        // TODO template this in the type of value and bound
        //      use a testing namespace above for that
        auto test_assert_greater_equal = [](long value1, long value2) {
            // Should this assertion fail?
            const bool should_catch_something = !(value1 <= value2);

            // Classify according to upper property
            RC_CLASSIFY(should_catch_something, "Assertion failed");

#ifdef DEBUG
            if (should_catch_something) {
                // sometimes we catch something in debug mode.
                RC_ASSERT_THROWS_AS(assert_greater_equal(value1, value2),
                                    ExcTooLarge<decltype(value1)>);
            } else {
                // Assert that value2 greater or equal than value 2
                // If error throw.
                assert_greater_equal(value1, value2)
            }
#else
            // we should never catch anything in release mode.
            assert_greater_equal(value1, value2)
#endif
        };
        CHECK(rc::check("Test assert_greater_equal",
                        test_assert_greater_equal));

        //
        // ---------------------------------------------------------
        //

        // TODO template this in the type of value and bound
        //      use a testing namespace above for that
        auto test_assert_greater = [](long value1, long value2) {
            // Should this assertion fail?
            const bool should_catch_something = !(value1 < value2);

            // Classify according to upper property
            RC_CLASSIFY(should_catch_something, "Assertion failed");

#ifdef DEBUG
            if (should_catch_something) {
                // sometimes we catch something in debug mode.
                RC_ASSERT_THROWS_AS(assert_greater(value1, value2),
                                    ExcTooLargeOrEqual<decltype(value1)>);
            } else {
                // Assert the lower bound. If error, throw
                assert_greater(value1, value2);
            }
#else
            // we should never catch anything in release mode.
            assert_greater(value1, value2);
#endif
        };
        CHECK(rc::check("Test assert_greater", test_assert_greater));

        //
        // ---------------------------------------------------------
        //

        auto test_assert_size = [] {
            // The signed data type equivalent to size_t
            typedef typename std::make_signed<size_t>::type ssize;

            // Generate a first size and a difference
            const ssize ssize1 = *gen::positive<ssize>();
            const ssize diff = *gen::arbitrary<int>();

            // Assert that difference is not too large
            RC_PRE(-diff < ssize1);

            // Build second size:
            const ssize ssize2 = ssize1 + diff;

            // Assert that size2 is positive:
            RC_PRE(ssize2 > 0);

            // Convert both to unsigned size_t
            const size_t size1 = static_cast<size_t>(ssize1);
            const size_t size2 = static_cast<size_t>(ssize2);

            // Should this assertion fail?
            const bool should_catch_something = (size1 != size2);

            // Classify according to upper property
            RC_CLASSIFY(should_catch_something, "Assertion failed");

#ifdef DEBUG
            if (should_catch_something) {
                // sometimes we catch something in debug mode.
                RC_ASSERT_THROWS_AS(assert_size(size1, size2), ExcSizeMismatch);
            } else {
                // Assert the sizes, if they do not match, throw
                assert_size(size1, size2);
            }
#else
            // we should never catch anything in release mode.
            assert_size(size1, size2);
#endif
        };
        CHECK(rc::check("Test assert_size", test_assert_size));

        //
        // ---------------------------------------------------------
        //
        //
        // TODO test for assert_element_sizes
    }
}
}  // tests
}  // linalgwrap
