#include <rapidcheck.h>
#include <catch.hpp>

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

    typedef double scalar_type;
    typedef typename Matrix_i<scalar_type>::size_type size_type;
    typedef SmallMatrix<scalar_type> small_matrix_type;

    // TODO test setting elements
    //
    // TODO use matrix_test_utils and perform random test against a plain vector
    // (or similar) as model.

    SECTION("Copy") {
        // Test copy constructor

        auto of_small_matrix = [](small_matrix_type m1) {
            auto copy(m1);

            // check that they are identical:
            NumComp::is_equal_matrix(copy, m1,
                                     std::numeric_limits<double>::epsilon());
        };

        // can use CHECK macro to have execution continue in this test in case
        // stuff fails
        REQUIRE(rc::check("Copy of small matrices", of_small_matrix));
    }

    SECTION("Addition") {
        // tests both + and += operators

        auto with_small_matrix = [](small_matrix_type m1) {
            // generate another matrix of the same size:
            auto m2 = *FixedSizeMatrix<small_matrix_type>::fixed_size(
                            m1.n_rows(), m1.n_cols());
            // add them
            auto res = m1 + m2;

            // check that it fits:
            // TODO use matrix iterator
            for (size_type i = 0; i < res.n_rows() * res.n_cols(); ++i) {
                RC_ASSERT(NumComp::is_equal(res[i], m1[i] + m2[i]));
            }
        };

        REQUIRE(rc::check("Addition of small matrices", with_small_matrix));
    }

    SECTION("Subtraction") {
        // tests both - and -= operators

        auto with_small_matrix = [](small_matrix_type m1) {
            // generate another matrix of the same size:
            auto m2 = *FixedSizeMatrix<small_matrix_type>::fixed_size(
                            m1.n_rows(), m1.n_cols());
            // subtract them
            auto res = m1 - m2;

            // check that it fits:
            // TODO use matrix iterator
            for (size_type i = 0; i < res.n_rows() * res.n_cols(); ++i) {
                RC_ASSERT(NumComp::is_equal(res[i], m1[i] - m2[i]));
            }
        };

        REQUIRE(rc::check("Difference of small matrices", with_small_matrix));
    }

    SECTION("Multiplication") {
        // Tests * operator

        auto with_small_matrix = [](small_matrix_type m1) {

            // The size of the other matrix to multiply m1 with:
            // Note: the RHS of inRange is exclusive
            auto othersize = *gen::inRange<size_type>(
                                   1, TestConstants::max_matrix_size + 1);

            // generate another matrix of the appropriate size
            auto m2 = *FixedSizeMatrix<small_matrix_type>::fixed_size(
                            m1.n_cols(), othersize);

            // multipy them
            auto res = m1 * m2;

            // check that it is correct:
            for (size_type i = 0; i < m1.n_rows(); ++i) {
                for (size_type j = 0; j < m2.n_cols(); ++j) {
                    scalar_type sumk{0};
                    for (size_type k = 0; k < m1.n_cols(); ++k) {
                        sumk += m1(i, k) * m2(k, j);
                    }

                    RC_ASSERT(NumComp::is_equal(res(i, j), sumk));
                }
            }
        };

        REQUIRE(
              rc::check("Multiplication of small matrices", with_small_matrix));
    }

    // TODO test equality and inequality operators
}

}  // tests
}  // linalgwrap
