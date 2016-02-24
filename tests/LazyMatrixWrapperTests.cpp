#include <rapidcheck.h>
#include <catch.hpp>
#include <LazyMatrixWrapper.hh>

// Generators for neccessary matrices
#include "generators.hh"

// Numerical equality and comparison
#include "NumComp.hh"

namespace linalgwrap {
namespace tests {
using namespace rc;

TEST_CASE("LazyMatrixWrapper class", "[LazyMatrixWrapper]") {
#ifdef DEBUG
    // Put this is some global context.
    // Make sure that the program does not get aborted
    exceptions::assert_dbg_effect = exceptions::ExceptionEffect::THROW;
#endif

    typedef double scalar_type;
    typedef SmallMatrix<scalar_type> small_matrix_type;

    SECTION("LazyMatrixWrapped SmallMatrix") {
        // Test if a wrapped small matrix is equivalent
        // to the unwrapped one.

        auto test_small_matrix_same = [](small_matrix_type m) {
            small_matrix_type inner{m};
            LazyMatrixWrapper<small_matrix_type, small_matrix_type> wrap{inner};

            // check that they are identical:
            RC_ASSERT(NumComp::is_equal_matrix(
                  m, wrap, std::numeric_limits<double>::epsilon()));
        };

        REQUIRE(rc::check("Equivalence of direct and lazy-wrapped SmallMatrix",
                          test_small_matrix_same));
    }

    // TODO - Test the other constructors
    //      - Test the other functionality of LazyMatrixWrapper.
}

}  // tests
}  // linalgwrap
