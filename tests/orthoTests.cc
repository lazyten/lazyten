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

#include <catch.hpp>
#include <lazyten/SmallVector.hh>
#include <lazyten/TestingUtils.hh>
#include <lazyten/ortho.hh>
#include <lazyten/random.hh>

namespace lazyten {
namespace tests {
using namespace rc;

TEST_CASE("Test ortho function", "[ortho]") {
  typedef double scalar_type;
  typedef SmallVector<double> vector_type;
  typedef SmallMatrix<double> matrix_type;

  auto gen_id = []() {
    const size_t n_dim = *gen::map(gen::numeric_size<2>(), [](size_t i) {
                            return i + 3;
                          }).as("Dimensionality");

    // The identity matrix:
    matrix_type id(n_dim, n_dim);
    for (size_t i = 0; i < n_dim; ++i) id(i, i) = 1;

    return id;
  };

  auto gen_pos_dev_mat = [&gen_id]() {
    // TODO Use the WithProperties generator once it has
    //      support for arbitrary positive definite matrices.
    auto id = gen_id();
    for (size_t i = 0; i < id.n_cols(); ++i) {
      id(i, i) *= *gen::map(gen::numeric_around(1.5), [](scalar_type in) {
                     return std::abs(in);
                   }).as("Weighting factor w");
    }
    return id;
  };

  auto testable_pre = [](matrix_type m, bool is_id) {
    const size_t n_vecs = *gen::inRange<size_t>(2, m.n_cols()).as("Number of vectors");

    auto gen_non_zero_vec =
          gen::construct<vector_type>(gen::container<std::vector<scalar_type>>(
                m.n_cols(),
                gen::mapcat(gen::exec([] { return random<scalar_type>(); }),
                            [](scalar_type s) { return gen::numeric_around(s); })));
    const auto vecs =
          *gen::numeric_tensor<MultiVector<vector_type>>(
                 n_vecs, gen::with_l2_norm_in_range(0.2, 100, gen_non_zero_vec))
                 .as("Vectors to orthogonalise.");

    // Orthogonalise
    const auto res = is_id ? ortho(vecs) : ortho(vecs, m);

    // Build scalar products
    const matrix_type innprod = dot(res, m * res);

    // The identity matrix:
    matrix_type id(n_vecs, n_vecs);
    for (size_t i = 0; i < n_vecs; ++i) id(i, i) = 1;

    RC_ASSERT(numcomp(id) == innprod);
  };

  SECTION("Orthogonalise random real vectors") {
    CHECK(rc::check("Orthogonalise random real vectors",
                    [&testable_pre, &gen_id]() { testable_pre(gen_id(), true); }));
  }

  SECTION("M-orthogonalise random real vectors") {
    CHECK(rc::check("M-orthogonalise random real vectors",
                    [&testable_pre, &gen_pos_dev_mat]() {
                      testable_pre(gen_pos_dev_mat(), false);
                    }));
  }

}  // Test ortho function

}  // namespace tests
}  // namespace lazyten
