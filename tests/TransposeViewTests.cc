#include <TransposeView.hh>
#include <ScaleView.hh>
#include <rapidcheck.h>
#include <catch.hpp>
#include "view_tests.hh"

namespace linalgwrap {
namespace tests {
using namespace rc;

TEST_CASE("TransposeView", "[TransposeView]") {
    // Make sure that the program does not get aborted
    exceptions::assert_dbg_effect = exceptions::ExceptionEffect::THROW;

    // Define the types we use for the test:
    struct TestTypes {
        typedef double scalar_type;
        typedef SmallMatrix<scalar_type> stored_matrix_type;
        typedef LazyMatrixWrapper<stored_matrix_type, stored_matrix_type>
              lazy_matrix_type;
        typedef view::ScaleView<stored_matrix_type> inner_scaleview_type;

        typedef view::TransposeView<const stored_matrix_type>
              view_of_stored_type;
        typedef view::TransposeView<inner_scaleview_type>
              view_of_scaleview_type;
        typedef view::TransposeView<lazy_matrix_type> view_of_lazy_type;
    };

    // Make types accessible
    typedef TestTypes::scalar_type scalar_type;
    typedef TestTypes::stored_matrix_type stored_matrix_type;
    typedef TestTypes::lazy_matrix_type lazy_matrix_type;
    typedef TestTypes::inner_scaleview_type inner_scaleview_type;

    // Generator for the args
    auto args_generator = []() {
        stored_matrix_type mat =
              *gen::arbitrary<stored_matrix_type>().as("Inner matrix");
        return std::make_pair(mat, std::make_tuple());
    };

    // Generator for the model:
    auto model_generator = [](std::pair<stored_matrix_type, std::tuple<>> p) {
        stored_matrix_type& m = p.first;

        stored_matrix_type transposed(m.n_cols(), m.n_rows());
        for (auto row : range(m.n_rows())) {
            for (auto col : range(m.n_cols())) {
                transposed(col, row) = m(row, col);
            }
        }
        return transposed;
    };

    // Generators for the test case views:
    typedef view_tests::StandardViewGenerators<TestTypes, std::tuple<>>
          standard_generators;

    SECTION("Default view tests on the stored view") {
        typedef view_tests::TestingLibrary<TestTypes::view_of_stored_type,
                                           decltype(args_generator())> testlib;

        auto make_view = [](const stored_matrix_type& sm,
                            std::tuple<>) { return view::transpose(sm); };

        // Generator for the stored view
        standard_generators::stored_view_generator svg(make_view);

        // Test library for the stored view
        testlib lib{args_generator, svg, model_generator,
                    "ScaleView(stored matrix): ",
                    0.1 * TestConstants::default_num_tol};

        // Run the tests:
        lib.run_checks();
    }

    SECTION("Default view tests on the lazy view") {
        typedef view_tests::TestingLibrary<TestTypes::view_of_lazy_type,
                                           decltype(args_generator())> testlib;

        auto make_view = [](lazy_matrix_type& sm,
                            std::tuple<>) { return view::transpose(sm); };

        // Generator for the lazy view
        standard_generators::lazy_view_generator lvg(make_view);

        // Test library for the lazy view
        testlib lib{args_generator, lvg, model_generator,
                    "ScaleView(lazy matrix): ",
                    0.1 * TestConstants::default_num_tol};

        // Run the tests:
        lib.disable_run_view_times_stored();
        lib.run_checks();
    }

    SECTION("Default view tests on the view view") {
        typedef view_tests::TestingLibrary<TestTypes::view_of_scaleview_type,
                                           decltype(args_generator())> testlib;

        auto make_view = [](inner_scaleview_type& sm,
                            std::tuple<>) { return view::transpose(sm); };

        // Generator for the scale-view view
        standard_generators::view_view_generator vvg(make_view);

        // Test library for the scale-view view
        testlib lib{args_generator, vvg, model_generator,
                    "ScaleView(inner stored view): ",
                    0.1 * TestConstants::default_num_tol};

        // Run the tests:
        lib.run_checks();
    }

}  // TEST_CASE
}  // tests
}  // linalgwrap
