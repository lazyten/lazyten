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

// Generators for neccessary matrices
#include "generators.hh"

// Numerical equality and comparison
#include "NumComp.hh"

namespace linalgwrap {
namespace tests {
using namespace rc;

#ifdef LINALGWRAP_HAVE_ARMADILLO
TEST_CASE("ArmadilloMatrix class", "[ArmadilloMatrix]") {
    // Make sure that the program does not get aborted
    AssertDbgEffect::set(ExceptionEffect::THROW);

    // The type of matrix we wish to test.
    typedef ArmadilloMatrix<double> matrix_type;

    SECTION("Default stored matrix tests") {
        typedef typename stored_matrix_tests::TestingLibrary<matrix_type>
              testinglib;
        testinglib lib("ArmadilloMatrix: ",
                       0.01 * TestConstants::default_num_tol);
        lib.run_checks();
    }
}
#endif

}  // tests
}  // linalgwrap
