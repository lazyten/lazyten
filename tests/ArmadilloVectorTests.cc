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

#include "stored_matrix_tests.hh"
#include <catch.hpp>
#include <linalgwrap/ArmadilloVector.hh>
#include <rapidcheck.h>

namespace linalgwrap {
namespace tests {
using namespace rc;

#ifdef LINALGWRAP_HAVE_ARMADILLO
TEST_CASE("ArmadilloVector class", "[ArmadilloVector]") {
    // Make sure that the program does not get aborted
    AssertDbgEffect::set(ExceptionEffect::THROW);

    // The type of vector we wish to test.
    typedef ArmadilloVector<double> vector_type;

    // TODO this is just a compile test, not a unit test!
    SECTION("Test initializer_list construction") {
        vector_type v{3., 4.};

        REQUIRE(v[0] == 3.);
        REQUIRE(v[1] == 4.);
        REQUIRE(v.size() == 2);
    }

    SECTION("Test vector norm") {
        vector_type v{-3., 4.};

        REQUIRE(v.norm_squared() == 25.);
        REQUIRE(v.l2_norm() == 5.);
        REQUIRE(v.l1_norm() == 7.);
        REQUIRE(v.linf_norm() == 4.);
    }

    // TODO extend in order to test:
    // l1, l2, linf norm
    // norm_squared
    // construction from initialiser list
    // +, -, +=, -=, *=, *, /, /=
    // size() function
    // vector operator() (with one size type)

    /*
    SECTION("Default stored vector tests") {
        typedef typename stored_matrix_tests::TestingLibrary<vector_type>
              testinglib;
        testinglib lib("ArmadilloVector: ",
                       0.01 * TestConstants::default_num_tol);
        lib.run_checks();
    }
    */
}
#endif

}  // tests
}  // linalgwrap
