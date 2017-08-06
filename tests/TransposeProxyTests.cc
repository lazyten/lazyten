//
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

#include "lazy_matrix_tests.hh"
#include <catch.hpp>
#include <lazyten/SmallMatrix.hh>
#include <lazyten/TransposeProxy.hh>
#include <rapidcheck.h>

namespace lazyten {
namespace tests {
using namespace rc;

TEST_CASE("TransposeProxy class", "[TransposeProxy]") {
  typedef double scalar_type;
  typedef SmallMatrix<scalar_type> stored_matrix_type;
  typedef TransposeProxy<stored_matrix_type> transpose_stored_type;
  typedef TransposeProxy<LazyMatrixWrapper<stored_matrix_type>> transpose_lazy_type;
  typedef typename stored_matrix_type::size_type size_type;

  // Generator for the args
  auto args_generator = [] {
    stored_matrix_type m = *gen::arbitrary<stored_matrix_type>().as("Transposed matrix");
    // RC_PRE(m.n_cols() > size_type{0} && m.n_rows() > size_type{0});
    return m;
  };

  // Generator for the model.
  auto model_generator = [](stored_matrix_type m) {
    stored_matrix_type mtrans(m.n_cols(), m.n_rows());
    for (size_type i = 0; i < m.n_rows(); ++i) {
      for (size_type j = 0; j < m.n_cols(); ++j) {
        mtrans(j, i) = m(i, j);
      }
    }
    return mtrans;
  };

  // Generator for the sut
  auto lazy_generator_stored = [](stored_matrix_type m) {
    transpose_stored_type proxy(std::move(m));
    RC_ASSERT(proxy.owns_inner_matrix());
    return proxy;
  };

  // Generator for the sut
  auto lazy_generator_lazy = [](stored_matrix_type m) {
    LazyMatrixWrapper<stored_matrix_type> wrap{std::move(m)};
    transpose_lazy_type proxy(std::move(wrap));
    RC_ASSERT(proxy.owns_inner_matrix());
    return proxy;
  };

  SECTION("Test class constructors") {
    stored_matrix_type test{{0, 1, 2}, {3, 4, 5}};
    TransposeProxy<stored_matrix_type> trans(test);
    REQUIRE(!trans.owns_inner_matrix());

    REQUIRE(trans.n_cols() == 2u);
    REQUIRE(trans.n_rows() == 3u);
    REQUIRE(trans(2, 0) == 2);
    REQUIRE(trans(1, 0) == 1);
    REQUIRE(trans(0, 1) == 3);

    //

    stored_matrix_type copy(test);
    TransposeProxy<stored_matrix_type> trans2(std::move(copy));
    REQUIRE(trans2.owns_inner_matrix());

    REQUIRE(trans2.n_cols() == 2u);
    REQUIRE(trans2.n_rows() == 3u);
    REQUIRE(trans2(2, 0) == 2);
    REQUIRE(trans2(1, 0) == 1);
    REQUIRE(trans2(0, 1) == 3);
  }

  SECTION("Default lazy matrix tests (with stored inner matrix)") {
    typedef lazy_matrix_tests::TestingLibrary<transpose_stored_type,
                                              decltype(args_generator())>
          testlib;

    testlib{args_generator, model_generator, lazy_generator_stored,
            "TransposeProxy (stored): "}
          .run_checks();
  }

  SECTION("Default lazy matrix tests (with lazy inner matrix)") {
    typedef lazy_matrix_tests::TestingLibrary<transpose_lazy_type,
                                              decltype(args_generator())>
          testlib;

    testlib{args_generator, model_generator, lazy_generator_lazy,
            "TransposeProxy (lazy): "}
          .run_checks();
  }
}  // TransposeProxy class

}  // namespace tests
}  // namespace lazyten
