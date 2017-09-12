//
// Copyright (C) 2017 by the lazyten authors
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

#include "stored_matrix_tests.hh"
#include <catch.hpp>
#include <lazyten/Bohrium/BohriumMatrix.hh>

#ifdef LAZYTEN_HAVE_BOHRIUM
namespace lazyten {
namespace tests {

TEST_CASE("BohriumMatrix", "[BohriumMatrix]") {
  typedef BohriumMatrix<double> matrix_type;

  SECTION("Check element access and matipulation") {
    matrix_type zeros(3, 4, true);
    zeros(1, 3) = 6;
    CHECK(zeros(1, 3) == 6.);
  }

  SECTION("Check equality") {
    matrix_type test = {{2, 3, 4}, {5, 6, 7}};
    matrix_type a = {{2, 2, 3}, {9, 8, 7}};
    CHECK(a != test);
    CHECK_FALSE(a == test);

    test.set_zero();
    a.set_zero();
    CHECK(a == test);
    CHECK_FALSE(a != test);
  }

  SECTION("Check matrix addition, multiplication and trace") {
    matrix_type A{{1, 2},  //
                  {3, 4}};
    matrix_type B{{6, 7},  //
                  {8, 9}};
    CHECK(trace(A) == 5);

    matrix_type res = A * B;
    matrix_type sum = A + B;

    CHECK(res(0, 0) == 22);
    CHECK(res(0, 1) == 25);
    CHECK(res(1, 0) == 50);
    CHECK(res(1, 1) == 57);

    CHECK(sum(0, 0) == 7);
    CHECK(sum(0, 1) == 9);
    CHECK(sum(1, 0) == 11);
    CHECK(sum(1, 1) == 13);

    CHECK(trace(A) == 5);
  }

  SECTION("Norm test") {
    matrix_type m = {{1.0},  //
                     {-2.0}};
    CHECK(m.n_cols() == 1);
    CHECK(m.n_rows() == 2);
    CHECK(m(0, 0) == 1.0);
    CHECK(m(1, 0) == -2.0);
    CHECK(norm_l1(m) == 3.);
    CHECK(norm_linf(m) == 2.);

    matrix_type n = {{1.0, -2.0}};
    CHECK(n.n_cols() == 2);
    CHECK(n.n_rows() == 1);
    CHECK(n(0, 0) == 1.0);
    CHECK(n(0, 1) == -2.0);
    CHECK(norm_l1(n) == 2.);
    CHECK(norm_linf(n) == 3.);

    matrix_type l{{1.0, 0.0},  //
                  {-2.0, 2.0}};
    CHECK(norm_l1(l) == 3.);
    CHECK(norm_linf(l) == 4.);

    matrix_type k{{1.0, -2.0},  //
                  {0.0, 0.0}};
    CHECK(norm_l1(k) == 2.);
    CHECK(norm_linf(k) == 3.);
  }

  SECTION("Default stored matrix tests") {
    typedef typename stored_matrix_tests::TestingLibrary<matrix_type> testinglib;

    // Run tests:
    const bool test_empty_matrices = false;  // Bohrium does not support empty matrices
    testinglib("BohriumMatrix: ", test_empty_matrices).run_checks();
  }

}  // BohriumMatrix
}  // namespace tests
}  // namespace lazyten
#endif  // LAZYTEN_HAVE_BOHRIUM
