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
#include <linalgwrap/ArmadilloMatrix.hh>
#include <rapidcheck.h>

// Generators for necessary matrices
#include "generators.hh"

namespace linalgwrap {
namespace tests {
using namespace rc;

#ifdef LINALGWRAP_HAVE_ARMADILLO
TEST_CASE("ArmadilloMatrix class", "[ArmadilloMatrix]") {
    // The type of matrix we wish to test.
    typedef ArmadilloMatrix<double> matrix_type;

    SECTION("Norm test") {
        matrix_type m = {{1.0},  //
                         {-2.0}};
        CHECK(m.n_cols() == 1);
        CHECK(m.n_rows() == 2);
        CHECK(m(0, 0) == 1.0);
        CHECK(m(1, 0) == -2.0);
        CHECK(norm_l1(m) == 3.);
        CHECK(norm_linf(m) == 2.);

        matrix_type n = {{1.0, -2.0}};
        CHECK(n.n_cols() == 2);
        CHECK(n.n_rows() == 1);
        CHECK(n(0, 0) == 1.0);
        CHECK(n(0, 1) == -2.0);
        CHECK(norm_l1(n) == 2.);
        CHECK(norm_linf(n) == 3.);

        matrix_type l{{1.0, 0.0},  //
                      {-2.0, 2.0}};
        CHECK(norm_l1(l) == 3.);
        CHECK(norm_linf(l) == 4.);

        matrix_type k{{1.0, -2.0},  //
                      {0.0, 0.0}};
        CHECK(norm_l1(k) == 2.);
        CHECK(norm_linf(k) == 3.);
    }

    SECTION("Default stored matrix tests") {
        typedef typename stored_matrix_tests::TestingLibrary<matrix_type>
              testinglib;

        // Run tests:
        testinglib("ArmadilloMatrix: ").run_checks();
    }
}
#endif

}  // tests
}  // linalgwrap
