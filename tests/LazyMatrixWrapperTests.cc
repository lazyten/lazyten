#include "lazy_matrix_tests.hh"
#include <catch.hpp>
#include <linalgwrap/LazyMatrixWrapper.hh>
#include <rapidcheck.h>

namespace linalgwrap {
namespace tests {
using namespace rc;

TEST_CASE("LazyMatrixWrapper class", "[LazyMatrixWrapper]") {
    // Test constructor
    // Test swapping

    // Make sure that the program does not get aborted
    exceptions::assert_dbg_effect = exceptions::ExceptionEffect::THROW;

    typedef double scalar_type;
    typedef SmallMatrix<scalar_type> stored_matrix_type;
    typedef LazyMatrixWrapper<stored_matrix_type, stored_matrix_type>
          lazy_matrix_type;
    typedef typename stored_matrix_type::size_type size_type;

    // Generator for the args
    auto args_generator = [] {
        stored_matrix_type m =
              *gen::arbitrary<stored_matrix_type>().as("Inner matrix");
        // TODO remove to test empty matrix case.
        RC_PRE(m.n_cols() > size_type{0} && m.n_rows() > size_type{0});
        return m;
    };

    // Generator for the model.
    auto model_generator = [](stored_matrix_type m) { return m; };

    // Generator for the sut
    auto lazy_generator = [](stored_matrix_type m) {
        lazy_matrix_type wrap{m};
        return wrap;
    };

    SECTION("Default lazy matrix tests") {
        typedef lazy_matrix_tests::TestingLibrary<lazy_matrix_type,
                                                  decltype(args_generator())>
              testlib;

        testlib lib{args_generator, lazy_generator, model_generator,
                    "LazyMatrixWrapper: ", TestConstants::default_num_tol};
        lib.run_checks();
    }
}

}  // tests
}  // linalgwrap
