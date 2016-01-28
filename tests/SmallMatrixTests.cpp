#include <rapidcheck.h>
#include <catch.hpp>

// Generators for neccessary matrices
#include "generators.hh"

namespace linalgwrap {
namespace tests {
using namespace rc;

TEST_CASE("SmallMatrix class", "[SmallMatrix]") {
    // Test the basic interface first ...

    typedef double scalar_type;
    typedef typename Matrix_i<scalar_type>::size_type size_type;
    typedef SmallMatrix<scalar_type> small_matrix_type;

    SECTION("Copy") {
        // Test copy constructor

        auto of_small_matrix = [](small_matrix_type m1) {
            auto copy(m1);

            // check that they are identical:
            for (size_type i = 0; i < m1.n_rows() * m1.n_cols(); ++i) {
                RC_ASSERT(copy[i] == m1[i]);
            }
        };

        // can use CHECK macro to have execution continue in this test in case
        // stuff fails
        REQUIRE(rc::check("Copy of small matrices", of_small_matrix));
    }

    SECTION("Addition") {
        // tests both + and += operators

        auto with_small_matrix = [](small_matrix_type m1) {
            // generate another matrix of the same size:
            auto m2 = *FixedSize<small_matrix_type>::fixed_size(m1.n_rows(),
                                                                m1.n_cols());
            // add them
            auto res = m1 + m2;

            // check that it fits:
            for (size_type i = 0; i < res.n_rows() * res.n_cols(); ++i) {
                RC_ASSERT(res[i] == m1[i] + m2[i]);
            }
        };

        REQUIRE(rc::check("Addition of small matrices", with_small_matrix));
    }
}

}  // tests
}  // linalgwrap
