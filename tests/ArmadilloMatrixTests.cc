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
    exceptions::assert_dbg_effect = exceptions::ExceptionEffect::THROW;

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
