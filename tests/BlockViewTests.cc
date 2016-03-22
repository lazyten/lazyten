#include <BlockView.hh>
#include <ScaleView.hh>
#include <rapidcheck.h>
#include <catch.hpp>
#include "view_tests.hh"

namespace linalgwrap {
namespace tests {
using namespace rc;

TEST_CASE("BlockView", "[BlockView]") {
    // Make sure that the program does not get aborted
    exceptions::assert_dbg_effect = exceptions::ExceptionEffect::THROW;

    // Define the types we use for the test:
    struct TestTypes {
        typedef double scalar_type;
        typedef SmallMatrix<scalar_type> stored_matrix_type;
        typedef LazyMatrixWrapper<stored_matrix_type, stored_matrix_type>
              lazy_matrix_type;
        typedef view::ScaleView<stored_matrix_type> inner_scaleview_type;

        typedef view::BlockView<const stored_matrix_type> view_of_stored_type;
        typedef view::BlockView<inner_scaleview_type> view_of_scaleview_type;
        typedef view::BlockView<lazy_matrix_type> view_of_lazy_type;
    };

    // Make types accessible
    typedef TestTypes::scalar_type scalar_type;
    typedef TestTypes::stored_matrix_type stored_matrix_type;
    typedef TestTypes::lazy_matrix_type lazy_matrix_type;
    typedef TestTypes::inner_scaleview_type inner_scaleview_type;
    typedef typename stored_matrix_type::size_type size_type;
    typedef Range<size_type> range_type;

    // Generator for the args
    auto args_generator = []() {
        stored_matrix_type mat =
              *gen::arbitrary<stored_matrix_type>().as("Inner matrix");

        // TODO remove later to test the case of an empty block explicitly
        RC_PRE(mat.n_rows() >= size_type{1} && mat.n_cols() >= size_type{1});

        // TODO Remove forcing a minimal length of 1 in the range
        //      in order to test the case of empty ranges / empty
        //      matrices
        range_type row_range =
              *gen::range_within<size_type>(0, mat.n_rows(), 1);
        range_type col_range =
              *gen::range_within<size_type>(0, mat.n_cols(), 1);

        return std::make_pair(mat, std::make_pair(row_range, col_range));
    };

    // Generator for the model:
    auto model_generator = [](
          std::pair<stored_matrix_type, std::pair<range_type, range_type>> p) {
        stored_matrix_type& m = p.first;

        range_type rows = p.second.first;
        range_type cols = p.second.second;

        stored_matrix_type ret(rows.length(), cols.length());

        if (rows.is_empty() || cols.is_empty()) {
            return ret;
        }

        for (auto row : rows) {
            for (auto col : cols) {
                ret(row - rows.first(), col - cols.first()) = m(row, col);
            }
        }

        return ret;
    };

    // Generators for the test case views:
    typedef view_tests::StandardViewGenerators<
          TestTypes, std::pair<range_type, range_type>> standard_generators;

    SECTION("Default view tests on the stored view") {
        typedef view_tests::TestingLibrary<TestTypes::view_of_stored_type,
                                           decltype(args_generator())> testlib;

        auto make_view = [](const stored_matrix_type& sm,
                            std::pair<range_type, range_type> p) {
            return view::block(sm, p.first, p.second);
            std::cout << p.first << "  " << p.second << "  " << sm << std::endl
                      << std::endl;
        };

        // Generator for the stored view
        standard_generators::stored_view_generator svg(make_view);

        // Test library for the stored view
        testlib lib{args_generator, svg, model_generator,
                    "BlockView(stored matrix): ",
                    0.1 * TestConstants::default_num_tol};

        // Run the tests:
        lib.run_checks();
    }

    SECTION("Default view tests on the lazy view") {
        typedef view_tests::TestingLibrary<TestTypes::view_of_lazy_type,
                                           decltype(args_generator())> testlib;

        auto make_view =
              [](lazy_matrix_type& sm, std::pair<range_type, range_type> p) {
            return view::block(sm, p.first, p.second);
        };

        // Generator for the lazy view
        standard_generators::lazy_view_generator lvg(make_view);

        // Test library for the lazy view
        testlib lib{args_generator, lvg, model_generator,
                    "BlockView(lazy matrix): ",
                    0.1 * TestConstants::default_num_tol};

        // Run the tests:
        lib.disable_run_view_times_stored();
        lib.run_checks();
    }

    SECTION("Default view tests on the view view") {
        typedef view_tests::TestingLibrary<TestTypes::view_of_scaleview_type,
                                           decltype(args_generator())> testlib;

        auto make_view = [](inner_scaleview_type& sm,
                            std::pair<range_type, range_type> p) {
            return view::block(sm, p.first, p.second);
        };

        // Generator for the scale-view view
        standard_generators::view_view_generator vvg(make_view);

        // Test library for the scale-view view
        testlib lib{args_generator, vvg, model_generator,
                    "BlockView(inner stored view): ",
                    0.1 * TestConstants::default_num_tol};

        // Run the tests:
        lib.disable_run_view_times_stored();
        lib.run_checks();
    }

}  // TEST_CASE
}  // tests
}  // linalgwrap
