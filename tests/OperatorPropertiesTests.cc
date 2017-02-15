//
// Copyright (C) 2017 by the linalgwrap authors
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

#include "linalgwrap/Base/Interfaces/OperatorProperties.hh"
#include <catch.hpp>

namespace linalgwrap {
namespace tests {

TEST_CASE("OperatorProperties enum", "[OperatorProperties]") {
  // Shorter aliases
  OperatorProperties none = OperatorProperties::None;
  OperatorProperties real = OperatorProperties::Real;
  OperatorProperties herm = OperatorProperties::Hermitian;
  OperatorProperties real_sym = OperatorProperties::RealSymmetric;
  OperatorProperties pos_semi_dev = OperatorProperties::PositiveSemiDefinite;
  OperatorProperties pos_dev = OperatorProperties::PositiveDefinite;
  OperatorProperties anti_herm = OperatorProperties::AntiHermitian;

  SECTION("Basic test") {
    // None is contained in everything:
    CHECK(props_contained_in(none, none));
    CHECK(props_contained_in(none, real));
    CHECK(props_contained_in(none, herm));
    CHECK(props_contained_in(none, real_sym));
    CHECK(props_contained_in(none, pos_semi_dev));
    CHECK(props_contained_in(none, pos_dev));
    CHECK(props_contained_in(none, anti_herm));

    // Hermitian / symmetric:
    CHECK_FALSE(props_contained_in(herm, none));
    CHECK_FALSE(props_contained_in(herm, real));
    CHECK(props_contained_in(herm, herm));
    CHECK(props_contained_in(herm, real_sym));
    CHECK(props_contained_in(herm, pos_semi_dev));
    CHECK(props_contained_in(herm, pos_dev));
    CHECK_FALSE(props_contained_in(herm, anti_herm));

    // Positive semi definite
    CHECK_FALSE(props_contained_in(pos_semi_dev, none));
    CHECK_FALSE(props_contained_in(pos_semi_dev, real));
    CHECK_FALSE(props_contained_in(pos_semi_dev, herm));
    CHECK_FALSE(props_contained_in(pos_semi_dev, real_sym));
    CHECK(props_contained_in(pos_semi_dev, pos_semi_dev));
    CHECK(props_contained_in(pos_semi_dev, pos_dev));
    CHECK_FALSE(props_contained_in(pos_semi_dev, anti_herm));

    // AntiHermitian
    CHECK_FALSE(props_contained_in(anti_herm, none));
    CHECK_FALSE(props_contained_in(anti_herm, real));
    CHECK_FALSE(props_contained_in(anti_herm, herm));
    CHECK_FALSE(props_contained_in(anti_herm, real_sym));
    CHECK_FALSE(props_contained_in(anti_herm, pos_semi_dev));
    CHECK_FALSE(props_contained_in(anti_herm, pos_dev));
    CHECK(props_contained_in(anti_herm, anti_herm));

  }  // Basic test
}  // OperatorProperties

}  // namespace tests
}  // namespace linalgwrap
