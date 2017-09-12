//
// Copyright (C) 2017 by the lazyten authors
//
// This file is part of lazyten.
//
// lazyten is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// lazyten is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with lazyten. If not, see <http://www.gnu.org/licenses/>.
//

#include "stored_vector_tests.hh"
#include <catch.hpp>
#include <lazyten/Bohrium/BohriumVector.hh>
#include <rapidcheck.h>

#ifdef LAZYTEN_HAVE_BOHRIUM
namespace lazyten {
namespace tests {

template <typename S>
using genlib = vector_tests::GeneratorLibrary<BohriumVector<S>>;

TEST_CASE("BohriumVector", "[BohriumVector]") {
  typedef BohriumVector<double> vector_type;
  SECTION("Test initializer_list construction") {
    vector_type v{3., 4.};

    REQUIRE(v[0] == 3.);
    REQUIRE(v[1] == 4.);
    REQUIRE(v.size() == 2);
  }

  SECTION("Test vector norm") {
    vector_type v{-3., 4.};

    REQUIRE(norm_l2_squared(v) == 25.);
    REQUIRE(norm_l2(v) == 5.);
    REQUIRE(norm_l1(v) == 7.);
    REQUIRE(norm_linf(v) == 4.);
    REQUIRE(accumulate(v) == 1.);
  }

  SECTION("Check element access and matipulation") {
    vector_type zeros(10, true);
    zeros[4] = 6;
    CHECK(zeros[4] == 6.);
    for (size_t i = 0; i < 10; ++i) {
      if (i == 4) continue;
      CHECK(zeros[i] == 0);
    }
  }

  SECTION("Default stored vector tests") {
    SECTION("Test double vectors") {
      // Decrease tolerance to require a more accurate numerical agreement
      // for passing.
      auto lowertol = NumCompConstants::change_temporary(
            0.1 * krims::NumCompConstants::default_tolerance_factor);

      stored_vector_tests::run_with_generator(genlib<double>::testgenerator(),
                                              "BohriumVector<double>: ");
    }

    SECTION("Test complex vectors") {
      // Decrease tolerance to require a more accurate numerical agreement
      // for passing.
      // auto lowertol = NumCompConstants::change_temporary(
      //      0.1 * krims::NumCompConstants::default_tolerance_factor);

      stored_vector_tests::run_with_generator(
            genlib<std::complex<double>>::testgenerator(),
            "BohriumVector<complex double>: ");
    }
  }

}  // BohriumVector

}  // namespace tests
}  // namespace lazyten
#endif  // LAZYTEN_HAVE_BOHRIUM
