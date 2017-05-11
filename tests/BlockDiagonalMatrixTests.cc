//
// Copyright (C) 2016-17 by the linalgwrap authors
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

#include "lazy_matrix_tests.hh"
#include <catch.hpp>
#include <linalgwrap/BlockDiagonalMatrix.hh>
#include <linalgwrap/LazyMatrixWrapper.hh>
#include <linalgwrap/SmallMatrix.hh>
#include <linalgwrap/SmallVector.hh>
#include <linalgwrap/TestingUtils.hh>
#include <rapidcheck.h>

namespace linalgwrap {
namespace tests {
using namespace rc;

/** TODO The next two functions are a little annoying and could be generalised using
 *  C++14 magic (use tuple_utils of krims for this! */
template <typename Stored>
std::array<LazyMatrixWrapper<Stored>, 0> make_lazy_array(std::array<Stored, 0> /* in */) {
  return {{}};
}
template <typename Stored>
std::array<LazyMatrixWrapper<Stored>, 1> make_lazy_array(std::array<Stored, 1> in) {
  return {{LazyMatrixWrapper<Stored>(std::move(in[0]))}};
}
template <typename Stored>
std::array<LazyMatrixWrapper<Stored>, 2> make_lazy_array(std::array<Stored, 2> in) {
  return {{LazyMatrixWrapper<Stored>(std::move(in[0])),
           LazyMatrixWrapper<Stored>(std::move(in[1]))}};
}
template <typename Stored>
std::array<LazyMatrixWrapper<Stored>, 3> make_lazy_array(std::array<Stored, 3> in) {
  return {{LazyMatrixWrapper<Stored>(std::move(in[0])),
           LazyMatrixWrapper<Stored>(std::move(in[1])),
           LazyMatrixWrapper<Stored>(std::move(in[2]))}};
}
template <typename Stored>
std::array<LazyMatrixWrapper<Stored>, 4> make_lazy_array(std::array<Stored, 4> in) {
  return {{LazyMatrixWrapper<Stored>(std::move(in[0])),
           LazyMatrixWrapper<Stored>(std::move(in[1])),
           LazyMatrixWrapper<Stored>(std::move(in[2])),
           LazyMatrixWrapper<Stored>(std::move(in[3]))}};
}

template <typename Array, size_t Size = std::tuple_size<Array>::value>
struct ArrayGen {
  static_assert(Size > 4, "The default case is not implemented.");
};
template <typename Array>
struct ArrayGen<Array, 0> {
  Array operator()(rc::Gen<typename Array::value_type> /*elementgen*/) { return Array{}; }
};
template <typename Array>
struct ArrayGen<Array, 1> {
  Array operator()(rc::Gen<typename Array::value_type> elementgen) {
    return {{*elementgen.as("Diagonal block 1")}};
  }
};
template <typename Array>
struct ArrayGen<Array, 2> {
  Array operator()(rc::Gen<typename Array::value_type> elementgen) {
    return {{*elementgen.as("Diagonal block 1"), *elementgen.as("Diagonal block 2")}};
  }
};
template <typename Array>
struct ArrayGen<Array, 3> {
  Array operator()(rc::Gen<typename Array::value_type> elementgen) {
    return {{*elementgen.as("Diagonal block 1"), *elementgen.as("Diagonal block 2"),
             *elementgen.as("Diagonal block 3")}};
  }
};
template <typename Array>
struct ArrayGen<Array, 4> {
  Array operator()(rc::Gen<typename Array::value_type> elementgen) {
    return {{*elementgen.as("Diagonal block 1"), *elementgen.as("Diagonal block 2"),
             *elementgen.as("Diagonal block 3"), *elementgen.as("Diagonal block 4")}};
  }
};

template <size_t N, typename Scalar = double>
struct Testfunctions {
  typedef Scalar scalar_type;
  typedef SmallMatrix<scalar_type> stored_matrix_type;
  typedef LazyMatrixWrapper<stored_matrix_type> lazy_matrix_type;
  typedef BlockDiagonalMatrix<lazy_matrix_type, N> bdmat_type;

  // Generator for the args
  static std::array<stored_matrix_type, N> args_generator() {
    auto gen_quadratic = gen::mapcat(gen::inRange(1, 21), [](size_t i) {
      return gen::numeric_tensor<stored_matrix_type>(i, i);
    });
    return ArrayGen<std::array<stored_matrix_type, N>>{}(std::move(gen_quadratic));
  }

  // Generator for the model.
  static stored_matrix_type model_generator(std::array<stored_matrix_type, N> blocks) {
    const size_t tot = std::accumulate(
          std::begin(blocks), std::end(blocks), 0ul,
          [](size_t v, const stored_matrix_type& m) { return v + m.n_rows(); });
    stored_matrix_type ret(tot, tot, true);

    size_t off = 0;
    for (auto& block : blocks) {
      for (size_t i = 0; i < block.n_rows(); ++i) {
        for (size_t j = 0; j < block.n_cols(); ++j) {
          ret(off + i, off + j) = block(i, j);
        }
      }
      off += block.n_rows();
    }

    return ret;
  }

  // Generator for the sut.
  static bdmat_type sut_generator(std::array<stored_matrix_type, N> blocks) {
    return bdmat_type(make_lazy_array(std::move(blocks)));
  }

  static void run_checks() {
    typedef lazy_matrix_tests::TestingLibrary<bdmat_type, decltype(args_generator())>
          testlib;
    testlib{args_generator, model_generator, sut_generator,
            "BlockDiagonal (" + std::to_string(N) + " blocks) "}
          .run_checks();
  }
};

//
// ----------------------------------------------------------------
//

TEST_CASE("BlockDiagonalMatrix class", "[BlockDiagonalMatrix]") {
  typedef double scalar_type;
  typedef SmallMatrix<scalar_type> stored_matrix_type;
  typedef SmallVector<scalar_type> vector_type;
  typedef LazyMatrixWrapper<stored_matrix_type> lazy_matrix_type;

  SECTION("Simple test case") {
    stored_matrix_type stored1{{1.0, 2.0},   // first row
                               {3.0, 4.0}};  // second row
    lazy_matrix_type lazy1{stored_matrix_type(stored1)};

    stored_matrix_type stored2{{5.0, 6.0},   // first row
                               {7.0, 8.0}};  // second row
    lazy_matrix_type lazy2{stored_matrix_type(stored2)};
    lazy2.add_properties(OperatorProperties::Real);

    BlockDiagonalMatrix<lazy_matrix_type, 2> diag1{{{lazy1, lazy2}}};

    lazy1.add_properties(OperatorProperties::Real);
    BlockDiagonalMatrix<lazy_matrix_type, 2> diag2{{{lazy2, lazy1}}};

    CHECK(diag1.properties() == OperatorProperties::None);
    CHECK(diag2.properties() == OperatorProperties::Real);

    CHECK(diag1.n_rows() == 4u);
    CHECK(diag2.n_rows() == 4u);
    CHECK(diag1.n_cols() == 4u);
    CHECK(diag2.n_cols() == 4u);

    auto sum = diag1 + diag2;
    auto prod = diag2 * diag1;

    CHECK(sum.n_rows() == 4u);
    CHECK(prod.n_rows() == 4u);
    CHECK(sum.n_cols() == 4u);
    CHECK(prod.n_cols() == 4u);

    CHECK(diag1(0, 0) == 1.);
    CHECK(diag1(0, 1) == 2.);
    CHECK(diag1(2, 2) == 5.);
    CHECK(diag1(2, 3) == 6.);
    CHECK(diag1(1, 3) == 0.);

    CHECK(diag2(0, 0) == 5.);
    CHECK(diag2(0, 1) == 6.);
    CHECK(diag2(2, 2) == 1.);
    CHECK(diag2(2, 3) == 2.);
    CHECK(diag2(1, 3) == 0.);

    CHECK(sum(2, 2) == 6.);
    CHECK(sum(0, 0) == 6.);
    CHECK(sum(3, 2) == 10.);
    CHECK(sum(1, 3) == 0.);

    CHECK(prod(0, 1) == 34.);
    CHECK(prod(2, 2) == 19.);
    CHECK(prod(2, 3) == 22.);
    CHECK(prod(1, 3) == 0.);

    //
    vector_type v{1, 2, 3, 4};
    vector_type vres = diag1 * v;
    CHECK(vres[0] == 5);
    CHECK(vres[1] == 11);
    CHECK(vres[2] == 39);
    CHECK(vres[3] == 53);
  }  // Simple test case

  // TODO test
  //   - properties()
  //   - has_apply_inverse()
  //   - has_transpose()
  //   - update()

  SECTION("Default lazy matrix tests for 1 block") {
    // Increase numeric tolerance for this scope,
    // ie results need to be less exact for passing
    auto highertol = NumCompConstants::change_temporary(
          10. * krims::NumCompConstants::default_tolerance_factor);
    Testfunctions<1>::run_checks();
  }

  SECTION("Default lazy matrix tests for 2 blocks") {
    // Increase numeric tolerance for this scope,
    // ie results need to be less exact for passing
    auto highertol = NumCompConstants::change_temporary(
          10. * krims::NumCompConstants::default_tolerance_factor);
    Testfunctions<2>::run_checks();
  }

  SECTION("Default lazy matrix tests for 4 blocks") {
    // Increase numeric tolerance for this scope,
    // ie results need to be less exact for passing
    auto highertol = NumCompConstants::change_temporary(
          10. * krims::NumCompConstants::default_tolerance_factor);
    Testfunctions<4>::run_checks();
  }
}  // BlockDiagonal tests

}  // namespace tests
}  // namespace linalgwrap
