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

#include "rapidcheck_utils.hh"
#include <catch.hpp>
#include <linalgwrap/TestingUtils.hh>
#include <linalgwrap/trans.hh>

namespace linalgwrap {
namespace tests {
using namespace rc;

TEST_CASE("trans function", "[trans]") {
  typedef double scalar_type;
  typedef SmallMatrix<scalar_type> matrix_type;
  typedef typename matrix_type::size_type size_type;

  SECTION("trans by lvalue") {
    auto test = [] {
      auto mat = *gen::numeric_tensor<matrix_type>().as("Inner matrix");

      auto transp = trans(mat);
      RC_ASSERT(mat.n_rows() == transp.n_cols());
      RC_ASSERT(mat.n_cols() == transp.n_rows());
      for (size_type i = 0; i < mat.n_rows(); ++i) {
        for (size_type j = 0; j < mat.n_cols(); ++j) {
          RC_ASSERT(mat(i, j) == transp(j, i));
        }
      }

      RC_ASSERT(!transp.owns_inner_matrix());
    };
    CHECK(rc::check("trans called with lvalue", test));
  }

  SECTION("trans by rvalue") {
    auto test = [] {
      auto mat = *gen::numeric_tensor<matrix_type>().as("Inner matrix");

      std::unique_ptr<matrix_type> copyptr(new matrix_type{mat});
      auto transp = trans(std::move(*copyptr));

      RC_ASSERT(mat.n_rows() == transp.n_cols());
      RC_ASSERT(mat.n_cols() == transp.n_rows());
      for (size_type i = 0; i < mat.n_rows(); ++i) {
        for (size_type j = 0; j < mat.n_cols(); ++j) {
          RC_ASSERT(mat(i, j) == transp(j, i));
        }
      }

      RC_ASSERT(transp.owns_inner_matrix());
      // Check that move really happened, this should not raise an
      // exception:
      copyptr.reset();
    };
    CHECK(rc::check("trans called with rvalue", test));
  }

  SECTION("Transposing a trans") {
    auto test = [] {
      auto mat = *gen::numeric_tensor<matrix_type>().as("Inner matrix");

      matrix_type copy(mat);
      auto transp = trans(std::move(copy));
      auto backtransp = trans(std::move(transp));

      RC_ASSERT_NC(mat == numcomp(backtransp));
    };
    CHECK(rc::check("Transposing a trans", test));
  }
}  // Testing trans

TEST_CASE("conjtrans function", "[conjtrans]") {
  // TODO Extend to complex scalar types once these are implemented
  typedef double scalar_type;
  typedef SmallMatrix<scalar_type> matrix_type;
  typedef typename matrix_type::size_type size_type;

  SECTION("conjtrans by lvalue") {
    auto test = [] {
      auto mat = *gen::numeric_tensor<matrix_type>().as("Inner matrix");

      auto transp = conjtrans(mat);
      RC_ASSERT(mat.n_rows() == transp.n_cols());
      RC_ASSERT(mat.n_cols() == transp.n_rows());
      for (size_type i = 0; i < mat.n_rows(); ++i) {
        for (size_type j = 0; j < mat.n_cols(); ++j) {
          RC_ASSERT(mat(i, j) == transp(j, i));
        }
      }

      RC_ASSERT(!transp.owns_inner_matrix());
    };
    CHECK(rc::check("conjtrans called with lvalue", test));
  }

  SECTION("conjtrans by rvalue") {
    auto test = [] {
      auto mat = *gen::numeric_tensor<matrix_type>().as("Inner matrix");

      std::unique_ptr<matrix_type> copyptr(new matrix_type{mat});
      auto transp = conjtrans(std::move(*copyptr));

      RC_ASSERT(mat.n_rows() == transp.n_cols());
      RC_ASSERT(mat.n_cols() == transp.n_rows());
      for (size_type i = 0; i < mat.n_rows(); ++i) {
        for (size_type j = 0; j < mat.n_cols(); ++j) {
          RC_ASSERT(mat(i, j) == transp(j, i));
        }
      }

      RC_ASSERT(transp.owns_inner_matrix());
      // Check that move really happened, this should not raise an
      // exception:
      copyptr.reset();
    };
    CHECK(rc::check("conjtrans called with rvalue", test));
  }

  SECTION("Transposing a conjtrans") {
    auto test = [] {
      auto mat = *gen::numeric_tensor<matrix_type>().as("Inner matrix");

      matrix_type copy(mat);
      auto transp = conjtrans(std::move(copy));
      auto backtransp = conjtrans(std::move(transp));

      RC_ASSERT_NC(mat == numcomp(backtransp));
    };
    CHECK(rc::check("Transposing a conjtrans", test));
  }
}  // Testing trans

}  // tests
}  // linalgwrap
