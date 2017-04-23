//
// Copyright (C) 2016-17 by the linalgwrap authors
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

#include "stored_vector_tests.hh"
#include <linalgwrap/Builtin/BuiltinVector.hh>

namespace linalgwrap {
namespace tests {
namespace builtin_vector_tests {
using namespace rc;

template <typename S>
using genlib = vector_tests::GeneratorLibrary<BuiltinVector<S>>;

TEST_CASE("BuiltinVector class", "[BuiltinVector]") {
  SECTION("Default stored vector tests") {
    // Decrease tolerance to require a more accurate numerical agreement
    // for passing.
    auto lowertol = NumCompConstants::change_temporary(
          0.1 * krims::NumCompConstants::default_tolerance_factor);

    stored_vector_tests::run_with_generator(genlib<double>::testgenerator(),
                                            "BuiltinVector<double>: ");

    stored_vector_tests::run_with_generator(genlib<std::complex<double>>::testgenerator(),
                                            "BuiltinVector<complex double>: ");
  }
}

}  // namespace builtin_vector_tests
}  // namespace tests
}  // namespace linalgwrap
