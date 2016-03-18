#include <rapidcheck.h>
#include <catch.hpp>
#include "stored_matrix_tests.hh"

// Generators for neccessary matrices
#include "generators.hh"

// Numerical equality and comparison
#include "NumComp.hh"

namespace linalgwrap {
namespace tests {
using namespace rc;

TEST_CASE("SmallMatrix class", "[SmallMatrix]") {
    // Make sure that the program does not get aborted
    exceptions::assert_dbg_effect = exceptions::ExceptionEffect::THROW;

    // The type of matrix we wish to test.
    typedef SmallMatrix<double> small_matrix_type;

    SECTION("Default stored matrix tests") {
        typedef typename stored_matrix_tests::TestingLibrary<small_matrix_type>
              testinglib;
        testinglib lib("SmallMatrix: ", TestConstants::default_num_tol);
        REQUIRE(lib.run_checks());
    }
}

}  // tests
}  // linalgwrap
