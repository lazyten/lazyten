//
// Copyright (C) 2016-17 by the lazyten authors
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
#include <lazyten/Armadillo/ArmadilloVector.hh>
#include <rapidcheck.h>

namespace lazyten {
namespace tests {
using namespace rc;

#ifdef LAZYTEN_HAVE_ARMADILLO
template <typename S>
using genlib = vector_tests::GeneratorLibrary<ArmadilloVector<S>>;

TEST_CASE("ArmadilloVector class", "[ArmadilloVector]") {
  SECTION("Test initializer_list construction") {
    typedef ArmadilloVector<double> vector_type;
    vector_type v{3., 4.};

    REQUIRE(v[0] == 3.);
    REQUIRE(v[1] == 4.);
    REQUIRE(v.size() == 2);
  }

  SECTION("Test vector norm") {
    typedef ArmadilloVector<double> vector_type;
    vector_type v{-3., 4.};

    REQUIRE(norm_l2_squared(v) == 25.);
    REQUIRE(norm_l2(v) == 5.);
    REQUIRE(norm_l1(v) == 7.);
    REQUIRE(norm_linf(v) == 4.);
    REQUIRE(accumulate(v) == 1.);
  }

  SECTION("Default stored vector tests") {
    {
      // Decrease tolerance to require a more accurate numerical agreement
      // for passing.
      auto lowertol = NumCompConstants::change_temporary(
            0.1 * krims::NumCompConstants::default_tolerance_factor);

      stored_vector_tests::run_with_generator(genlib<double>::testgenerator(),
                                              "ArmadilloVector<double>: ");
    }

    // TODO improve the accuracy for complex vectors:
    // run this at 0.1 * tolerance as well.
    auto highertol = NumCompConstants::change_temporary(
          10. * krims::NumCompConstants::default_tolerance_factor);

    stored_vector_tests::run_with_generator(genlib<std::complex<double>>::testgenerator(),
                                            "ArmadilloVector<complex double>: ");
  }
}
#endif

}  // namespace tests
}  // namespace lazyten
