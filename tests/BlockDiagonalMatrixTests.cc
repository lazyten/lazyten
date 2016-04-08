#include <BlockDiagonalMatrix.hh>
#include <LazyMatrixWrapper.hh>
#include <SmallMatrix.hh>
#include <catch.hpp>
#include <rapidcheck.h>

namespace linalgwrap {
namespace tests {
using namespace rc;

TEST_CASE("BlockDiagonalMatrix class", "[BlockDiagonalMatrix]") {
    // Make sure that the program does not get aborted
    exceptions::assert_dbg_effect = exceptions::ExceptionEffect::THROW;

    typedef double scalar_type;
    typedef SmallMatrix<scalar_type> stored_matrix_type;
    typedef LazyMatrixWrapper<stored_matrix_type, stored_matrix_type>
          lazy_matrix_type;

    SECTION("Simple function test") {
        stored_matrix_type stored{{1.0, 2.0}, {3.0, 4.0}};
        lazy_matrix_type wrap{stored};

        auto diag1 = make_block_diagonal(stored_matrix_type{stored},
                                         stored_matrix_type{stored});

        auto diag2 = make_block_diagonal(std::move(wrap), std::move(stored));

        auto res = diag1 + diag2;

        auto res2 = diag1 * diag2;

        CHECK(diag1.n_rows() == 4u);
        CHECK(diag2.n_rows() == 4u);
        CHECK(res.n_rows() == 4u);

        CHECK(diag1(0, 0) == 1.);
        CHECK(diag1(0, 1) == 2.);
        CHECK(diag1(2, 2) == 1.);
        CHECK(diag1(2, 3) == 2.);
        CHECK(diag1(1, 3) == 0.);

        CHECK(diag2(0, 0) == 1.);
        CHECK(diag2(0, 1) == 2.);
        CHECK(diag2(2, 2) == 1.);
        CHECK(diag2(2, 3) == 2.);
        CHECK(diag2(1, 3) == 0.);

        CHECK(res(2, 2) == 2.);
        CHECK(res(0, 0) == 2.);
        CHECK(res(1, 3) == 0.);

        CHECK(res2(2, 2) == 7.);
        CHECK(res2(2, 3) == 10.);
    }
}

}  // tests
}  // linalgwrap
