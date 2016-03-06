#include <rapidcheck.h>
#include <catch.hpp>
#include <ParameterMap.hh>
#include <SmallMatrix.hh>
#include <SubscriptionPointer.hh>

namespace linalgwrap {
namespace tests {
using namespace rc;

TEST_CASE("ParameterMap tests", "[parametermap]") {
    // Make sure that the program does not get aborted
    exceptions::assert_dbg_effect = exceptions::ExceptionEffect::THROW;

    // Some data:
    int i = 5;
    std::string s{"test"};
    SmallMatrix<double> mat(2, 2, false);
    mat(0, 0) = -1;
    mat(1, 1) = -5;
    mat(1, 0) = mat(0, 1) = 2;

    SECTION("Can add data to ParameterMap") {
        // Insert some data into a map.
        ParameterMap m{};
        m.update_copy("string", s);
        m.update_copy("integer", i);
        m.update("matrix", make_subscription(mat, "ParameterMap_FirstTest"));

        // See if we get it back:
        REQUIRE(m.at<int>("integer") == i);
        REQUIRE(m.at<std::string>("string") == s);
        REQUIRE(m.at<SmallMatrix<double>>("matrix") == mat);
    }

//
// ---------------------------------------------------------------
//

#ifdef DEBUG
    SECTION("Check that type safety is assured") {
        // Add data to map.
        ParameterMap m{};
        m.update_copy("s", s);
        m.update_copy("i", i);
        m.update("mat", make_subscription(mat, "ParameterMap_SecondTest"));

        // Extract using the wrong type
        REQUIRE_THROWS_AS(m.at<double>("i"),
                          ParameterMap::ExcWrongTypeRequested);
        REQUIRE_THROWS_AS(m.at<double>("s"),
                          ParameterMap::ExcWrongTypeRequested);
        REQUIRE_THROWS_AS(m.at<double>("mat"),
                          ParameterMap::ExcWrongTypeRequested);
    }

    //
    // ---------------------------------------------------------------
    //

    SECTION("Check for UnknownKey") {
        // Add data to map.
        ParameterMap m{};
        m.update_copy("i", i);

        REQUIRE_THROWS_AS(m.at<bool>("blubber"), ParameterMap::ExcUnknownKey);
        REQUIRE_THROWS_AS(m.at<int>("blubb"), ParameterMap::ExcUnknownKey);
        REQUIRE(m.at<int>("i") == i);
    }
#endif

    //
    // ---------------------------------------------------------------
    //

    SECTION("Check that data can be erased") {
        // Add data to map.
        ParameterMap m{};
        m.update_copy("s", s);
        m.update_copy("i", i);
        m.update("mat", make_subscription(mat, "ParameterMap_SecondTest"));

        // check it is there:
        REQUIRE(m.exists("i"));
        REQUIRE(m.exists("s"));
        REQUIRE(m.exists("mat"));

        // remove a few:
        m.erase("i");
        m.erase("mat");

        // check they are there (or not)
        REQUIRE(!m.exists("i"));
        REQUIRE(m.exists("s"));
        REQUIRE(!m.exists("mat"));
    }

}  // TEST_CASE
}  // namespace test
}  // namespace linalgwrap
