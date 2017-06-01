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

#include <catch.hpp>
#include <linalgwrap/version.hh>

namespace linalgwrap {
namespace tests {

TEST_CASE("Version informaton", "[version]") {
  SECTION("Version string") {
    const std::string expected = std::to_string(version::major) + "." +
                                 std::to_string(version::minor) + "." +
                                 std::to_string(version::patch);
    CHECK(version::version_string() == expected);
  }

  SECTION("Feature string") {
    const std::string feature_string = version::feature_string();
    for (auto kv : version::feature_availability) {
      const std::string prefix = kv.second ? "+" : "-";
      CHECK(feature_string.find(prefix + kv.first) != std::string::npos);

      CHECK(version::has_feature(kv.first) == kv.second);
    }

    CHECK_FALSE(version::has_feature("blabbel"));
  }
}  // version TEST_CASE
}  // namespace tests
}  // namespace linalgwrap
