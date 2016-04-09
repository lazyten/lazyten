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
#include <linalgwrap/ParameterMap.hh>
#include <linalgwrap/SmallMatrix.hh>
#include <linalgwrap/SubscriptionPointer.hh>
#include <rapidcheck.h>

namespace linalgwrap {
namespace tests {
using namespace rc;

TEST_CASE("ParameterMap tests", "[parametermap]") {
    // Make sure that the program does not get aborted
    AssertDbgEffect::set(ExceptionEffect::THROW);

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
