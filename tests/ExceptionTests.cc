#include <rapidcheck.h>
#include <catch.hpp>
#include <Exceptions.hh>

namespace linalgwrap {
namespace tests {
using namespace rc;

TEST_CASE("Exception system", "[exception]") {
    // Make sure that the program does not get aborted
    linalgwrap::exceptions::assert_dbg_effect =
          linalgwrap::exceptions::ExceptionEffect::THROW;

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
                          exceptions::ExcNotImplemented);

        // Test if an ExcIO exception would be thrown:
        REQUIRE_THROWS_AS(assert_throw(false, ExcIO()), exceptions::ExcIO);
    }

    //
    // ----------------------------------------------------------------
    //

    SECTION("Test derived assertion macros") {
        // TODO template this in the type of value and bound
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
        auto test_assert_lower_bound = [](long lower_bound, long value) {
            // Should this assertion fail?
            const bool should_catch_something = (value < lower_bound);

            // Classify according to upper property
            RC_CLASSIFY(should_catch_something, "Assertion failed");

#ifdef DEBUG
            if (should_catch_something) {
                // sometimes we catch something in debug mode.
                RC_ASSERT_THROWS_AS(assert_lower_bound(lower_bound, value),
                                    ExcBelowLowerBound<decltype(value)>);
            } else {
                // Assert the lower bound. If error, throw
                assert_lower_bound(lower_bound, value);
            }
#else
            // we should never catch anything in release mode.
            assert_lower_bound(value, lower_bound);
#endif
        };
        CHECK(rc::check("Test assert_lower_bound", test_assert_lower_bound));

        //
        // ---------------------------------------------------------
        //

        // TODO template this in the type of value and bound
        auto test_assert_upper_bound = [](long upper_bound, long value) {
            // Should this assertion fail?
            const bool should_catch_something = (value >= upper_bound);

            // Classify according to upper property
            RC_CLASSIFY(should_catch_something, "Assertion failed");

#ifdef DEBUG
            if (should_catch_something) {
                // sometimes we catch something in debug mode.
                RC_ASSERT_THROWS_AS(assert_upper_bound(value, upper_bound),
                                    ExcAboveUpperBound<decltype(value)>);
            } else {
                // Assert the lower bound. If error, throw
                assert_upper_bound(value, upper_bound);
            }
#else
            // we should never catch anything in release mode.
            assert_upper_bound(value, upper_bound);
#endif
        };
        CHECK(rc::check("Test assert_upper_bound", test_assert_upper_bound));

        //
        // ---------------------------------------------------------
        //

        auto test_assert_size = [](size_t size1, size_t size2) {
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
            assert_upper_bound(value, upper_bound);
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
