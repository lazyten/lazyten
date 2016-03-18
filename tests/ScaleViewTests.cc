#include <ScaleView.hh>
#include <rapidcheck.h>
#include <catch.hpp>
#include "view_tests.hh"
#include <tuple>

namespace linalgwrap {
namespace tests {
using namespace rc;

TEST_CASE("ScaleView", "[ScaleView]") {
    // Make sure that the program does not get aborted
    exceptions::assert_dbg_effect = exceptions::ExceptionEffect::THROW;

    // Stored and lazy matrix types:
    typedef double scalar_type;
    typedef SmallMatrix<scalar_type> stored_matrix_type;
    typedef LazyMatrixWrapper<stored_matrix_type, stored_matrix_type>
          lazy_matrix_type;

    // View types to test
    typedef view::ScaleView<const stored_matrix_type> stored_view;
    typedef view::ScaleView<lazy_matrix_type> lazy_view;

    // Generator for the args
    auto args_generator = []() {
        scalar_type fac = *gen::arbitrary<scalar_type>().as("Scaling factor");
        stored_matrix_type mat =
              *gen::arbitrary<stored_matrix_type>().as("Inner matrix");
        return std::make_tuple(fac, mat);
    };

    // Generator for the model:
    auto model_generator = [](std::tuple<scalar_type, stored_matrix_type> t) {
        // Just multiply them:
        return std::get<0>(t) * std::get<1>(t);
    };

    // Generator for the stored view:
    struct stored_view_generator {
        stored_view operator()(std::tuple<scalar_type, stored_matrix_type> t) {
            // Store the stored matrix in permanent memory:
            m_stored_ptr.reset(
                  new stored_matrix_type{std::move(std::get<1>(t))});

            // Return a view to it:
            return view::scale(*m_stored_ptr, std::get<0>(t));
        }

      private:
        std::shared_ptr<const stored_matrix_type> m_stored_ptr;
    };

    // Generator for the lazy view:
    struct lazy_view_generator {
        lazy_view operator()(std::tuple<scalar_type, stored_matrix_type> t) {
            // Store the lazy matrix in permanent memory:
            m_lazy_ptr.reset(new lazy_matrix_type{std::move(std::get<1>(t))});

            // Return a view to it:
            return view::scale(*m_lazy_ptr, std::get<0>(t));
        }

      private:
        std::shared_ptr<lazy_matrix_type> m_lazy_ptr;
    };

    // Testing libraries
    typedef view_tests::TestingLibrary<stored_view, decltype(args_generator())>
          stored_testlib;
    typedef view_tests::TestingLibrary<lazy_view, decltype(args_generator())>
          lazy_testlib;

    SECTION("Default view tests on the stored view") {
        stored_view_generator svg{};
        stored_testlib lib{args_generator, svg, model_generator,
                           "ScaleView(stored matrix): ",
                           TestConstants::default_num_tol};
        REQUIRE(lib.run_checks());
    }

    SECTION("Default view tests on the lazy view") {
        lazy_view_generator lvg{};
        lazy_testlib lib{args_generator, lvg, model_generator,
                         "ScaleView(lazy matrix): ",
                         TestConstants::default_num_tol};
        REQUIRE(lib.run_checks());
    }

}  // TEST_CASE
}  // tests
}  // linalgwrap
