
// Copyright (C) 2016-17 by the lazyten authors
//
// This file is part of lazyten.
//
// lazyten is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// lazyten is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with lazyten. If not, see <http://www.gnu.org/licenses/>.
//

#include "generators.hh"
#include "lazy_matrix_tests_state.hh"
#include "rapidcheck_utils.hh"
#include <catch.hpp>
#include <lazyten/LazyMatrixProduct.hh>
#include <lazyten/SmallMatrix.hh>
#include <rapidcheck.h>

namespace lazyten {
namespace tests {
using namespace rc;

TEST_CASE("LazyMatrixProduct", "[LazyMatrixProduct]") {
  // TODO  Test swap function
  // TODO  Test constructors

  typedef double scalar_type;
  typedef SmallMatrix<scalar_type> stored_matrix_type;
  typedef stored_matrix_type model_matrix_type;

  SECTION("Test access into empty LazyMatrixProduct objects") {
    LazyMatrixProduct<stored_matrix_type> prod;
    CHECK(prod.n_rows() == 0);
    CHECK(prod.n_cols() == 0);

    auto random_access = [&]() {

      auto col = *rc::gen::arbitrary<size_t>().as("Column index");
      auto row = *rc::gen::arbitrary<size_t>().as("Row index");
      auto fac = *rc::gen::arbitrary<scalar_type>().as("Factor");

      LazyMatrixProduct<stored_matrix_type> copy(prod);
      copy *= fac;
      RC_ASSERT(fac == copy(col, col));
      RC_ASSERT(fac == copy(row, row));
      if (col != row) RC_ASSERT(0 == copy(row, col));
    };
    REQUIRE(rc::check("Random access into empty product", random_access));
  }

  SECTION("Test empty LazyMatrixProduct objects behave as diagonal matrices") {
    typedef lazy_matrix_tests::FunctionalityTests<stored_matrix_type,
                                                  LazyMatrixProduct<stored_matrix_type>>
          testing_lib;

    auto test_mult_with_empty = [](stored_matrix_type tomult) {
      LazyMatrixWrapper<stored_matrix_type> lazy{stored_matrix_type(tomult)};
      LazyMatrixProduct<stored_matrix_type> empty;

      LazyMatrixProduct<stored_matrix_type> res;
      res = res * lazy;
      testing_lib::run_all_tests(tomult, res);

      auto mult1 = empty * lazy;
      testing_lib::run_all_tests(tomult, mult1);

      auto mult2 = lazy * empty;
      testing_lib::run_all_tests(tomult, mult2);
    };

    // Enable hard tests in the testing library
    testing_lib::skip_easy_cases();

    // Increase numeric tolerance for this scope,
    // ie results need to be less exact for passing
    auto highertol = NumCompConstants::change_temporary(
          10. * krims::NumCompConstants::default_tolerance_factor);

    REQUIRE(rc::check("Test multiplying terms with empty product", test_mult_with_empty));
  }

  SECTION("Default lazy matrix tests") {
    typedef double scalar_type;
    typedef SmallMatrix<scalar_type> stored_matrix_type;
    typedef LazyMatrixWrapper<stored_matrix_type> lazy_matrix_type;
    typedef LazyMatrixProduct<stored_matrix_type> product_type;

    // Generator for the args
    auto args_generator = []() {
      stored_matrix_type mat1 =
            *gen::scale(0.9, gen::numeric_tensor<stored_matrix_type>())
                   .as("Matrix factor 1");

      // The size of the other matrix to multiply mat1 with:
      // Note: the RHS of inRange is exclusive
      auto othersize = *gen::numeric_size<2>().as("Number of columns of the 2nd matrix");

      stored_matrix_type mat2 = *gen::scale(0.9, gen::numeric_tensor<stored_matrix_type>(
                                                       mat1.n_cols(), othersize))
                                       .as("Matrix factor 2");

#ifdef LAZYTEN_TESTS_VERBOSE
      RC_CLASSIFY(norm_frobenius(mat1) == 0, "Factor 1 is zero");
      RC_CLASSIFY(norm_frobenius(mat2) == 0, "Factor 2 is zero");
#endif

      return std::make_pair(mat1, mat2);
    };

    // Generator for the model:
    auto model_generator = [](std::pair<stored_matrix_type, stored_matrix_type> p) {
      stored_matrix_type res(p.first.n_rows(), p.second.n_cols(), false);
      matrix_tests::matrix_product(p.first, p.second, res);
      return res;
    };

    // Generator for the sut:
    auto sut_generator = [](std::pair<stored_matrix_type, stored_matrix_type> p) {
      // Make two unit rectangles as extra factors
      RC_ASSERT(p.first.n_cols() == p.second.n_rows());
      stored_matrix_type unitrect1{p.first.n_cols(), p.first.n_cols() + 3};
      stored_matrix_type unitrect2{p.first.n_cols() + 3, p.first.n_cols()};
      for (typename stored_matrix_type::size_type i = 0; i < p.first.n_cols(); ++i) {
        unitrect1(i, i) = Constants<scalar_type>::one;
        unitrect2(i, i) = Constants<scalar_type>::one;
      }

      lazy_matrix_type lm1 = lazy_matrix_type{std::move(p.first)};
      lazy_matrix_type lm_unit1 = lazy_matrix_type{std::move(unitrect1)};
      lazy_matrix_type lm_unit2 = lazy_matrix_type{std::move(unitrect2)};
      lazy_matrix_type lm2 = lazy_matrix_type{std::move(p.second)};

      product_type prod{lm1};
      prod.push_factor(lm_unit1);
      prod.push_factor(lm_unit2);
      prod.push_factor(lm2);
      return prod;
    };

    typedef lazy_matrix_tests::TestingLibrary<product_type, decltype(args_generator())>
          testinglib;

    // Increase numeric tolerance for this scope,
    // ie results need to be less exact for passing
    // Smaller tolerance values do not work here,
    // since sometimes we have have rather in the matrices generated
    // by args_generator already.
    auto highertol = NumCompConstants::change_temporary(
          10. * krims::NumCompConstants::default_tolerance_factor);

    testinglib{args_generator, model_generator, sut_generator, "LazyMatrixProduct: "}
          .run_checks();
  }

  SECTION("Random function test") {
    // Increase numeric tolerance for this scope,
    // ie results need to be less exact for passing
    auto highertol = NumCompConstants::change_temporary(
          100. * krims::NumCompConstants::default_tolerance_factor);

    auto random_test = [] {
      model_matrix_type in(2, 3);
      in(0, 0) = 3;
      in(1, 1) = 2;
      in(0, 1) = -1;
      in(1, 0) = -4;
      in(0, 2) = 1;
      in(1, 2) = -1;

      //
      // The actual test
      //
      // The test library we use
      typedef lazy_matrix_tests::StatefulTestingLibrary<
            model_matrix_type, LazyMatrixProduct<stored_matrix_type>>
            test_library;

      // The commands we check
      auto gen_commands = state::gen::execOneOfWithArgs<
            typename test_library::op_MultiplyLazy, typename test_library::op_UnaryMinus,
            typename test_library::op_MultScalar, typename test_library::op_DivideScalar>;

      // Run the check:
      test_library().run_check(in, gen_commands(), 0.25);
    };

    REQUIRE(rc::check("Random function test of LazyMatrixProduct.", random_test));
  }  // Random function test
}

}  // namespace tests
}  // namespace lazyten
