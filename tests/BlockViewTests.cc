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

#include "view_tests.hh"
#include <catch.hpp>
#include <rapidcheck.h>

namespace linalgwrap {
namespace tests {
using namespace rc;

TEST_CASE("BlockView", "[BlockView]") {
    // Make sure that the program does not get aborted
    AssertDbgEffect::set(ExceptionEffect::THROW);

    // Define the types we use for the test:
    struct TestTypes {
        typedef double scalar_type;
        typedef SmallMatrix<scalar_type> stored_matrix_type;
        typedef LazyMatrixWrapper<stored_matrix_type, stored_matrix_type>
              lazy_matrix_type;
        typedef view::detail::ScaleView<stored_matrix_type>
              inner_scaleview_type;

        typedef view::detail::BlockView<const stored_matrix_type>
              view_of_stored_type;
        typedef view::detail::BlockView<inner_scaleview_type>
              view_of_scaleview_type;
        typedef view::detail::BlockView<lazy_matrix_type> view_of_lazy_type;
    };

    // Make types accessible
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

        if (rows.empty() || cols.empty()) {
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
          TestTypes, std::pair<range_type, range_type>>
          standard_generators;

    SECTION("Default view tests on the stored view") {
        typedef view_tests::TestingLibrary<TestTypes::view_of_stored_type,
                                           decltype(args_generator())>
              testlib;

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
                                           decltype(args_generator())>
              testlib;

        auto make_view = [](lazy_matrix_type& sm,
                            std::pair<range_type, range_type> p) {
            return view::block(sm, p.first, p.second);
        };

        // Generator for the lazy view
        standard_generators::lazy_view_generator lvg(make_view);

        // Test library for the lazy view
        testlib lib{args_generator, lvg, model_generator,
                    "BlockView(lazy matrix): ",
                    0.1 * TestConstants::default_num_tol};

        // Run the tests:
        lib.disable_run_matrix_times_stored();
        lib.run_checks();
    }

    SECTION("Default view tests on the view view") {
        typedef view_tests::TestingLibrary<TestTypes::view_of_scaleview_type,
                                           decltype(args_generator())>
              testlib;

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
        lib.disable_run_matrix_times_stored();
        lib.run_checks();
    }

}  // TEST_CASE
}  // tests
}  // linalgwrap
