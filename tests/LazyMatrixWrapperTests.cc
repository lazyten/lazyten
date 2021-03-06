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
#include <lazyten/LazyMatrixWrapper.hh>
#include <lazyten/SmallMatrix.hh>
#include <rapidcheck.h>

namespace lazyten {
namespace tests {
using namespace rc;

TEST_CASE("LazyMatrixWrapper class", "[LazyMatrixWrapper]") {
  // Test constructor
  // Test swapping

  typedef double scalar_type;
  typedef SmallMatrix<scalar_type> stored_matrix_type;
  typedef LazyMatrixWrapper<stored_matrix_type> lazy_matrix_type;
  typedef typename stored_matrix_type::size_type size_type;

  // Generator for the args
  auto args_generator = [] {
    stored_matrix_type m = *gen::arbitrary<stored_matrix_type>().as("Inner matrix");
    // TODO remove to test empty matrix case.
    RC_PRE(m.n_cols() > size_type{0} && m.n_rows() > size_type{0});
    return m;
  };

  // Generator for the model.
  auto model_generator = [](stored_matrix_type m) { return m; };

  // Generator for the sut
  auto lazy_generator = [](stored_matrix_type m) {
    lazy_matrix_type wrap{std::move(m)};
    return wrap;
  };

  SECTION("Default lazy matrix tests") {
    typedef lazy_matrix_tests::TestingLibrary<lazy_matrix_type,
                                              decltype(args_generator())>
          testlib;

    testlib{args_generator, model_generator, lazy_generator, "LazyMatrixWrapper: "}
          .run_checks();
  }
}

}  // namespace tests
}  // namespace lazyten
