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

#include "eigensolver_tests.hh"
#include <lazyten/Lapack/LapackEigensolver.hh>
#include <lazyten/LazyMatrixWrapper.hh>
#include <lazyten/TestingUtils.hh>
#include <rapidcheck.h>

#ifdef LAZYTEN_HAVE_LAPACK
namespace lazyten {
namespace tests {
using namespace rc;

/** Traits class needed for the tests */
struct LapackEigensolverTraits {
  template <typename Eigenproblem>
  using Solver = LapackEigensolver<Eigenproblem>;
};

TEST_CASE("LapackEigensolver", "[LapackEigensolver]") {
  using namespace eigensolver_tests;
  typedef double scalar_type;
  typedef ArmadilloMatrix<scalar_type> matrix_type;
  typedef ArmadilloVector<scalar_type> vector_type;

  SECTION("Test LapackPackedMatrix") {
    SECTION("Simple hard-coded test") {
      matrix_type mat{{1, 2, 3}, {2, 4, 5}, {3, 5, 6}};
      LazyMatrixWrapper<matrix_type> mat_wrap(mat);

      detail::LapackPackedMatrix<double> pack(mat);
      detail::LapackPackedMatrix<double> packwrap(mat_wrap);

      CHECK(pack.elements.size() == 6);
      CHECK(pack.elements[0] == 1);
      CHECK(pack.elements[1] == 2);
      CHECK(pack.elements[2] == 3);
      CHECK(pack.elements[3] == 4);
      CHECK(pack.elements[4] == 5);
      CHECK(pack.elements[5] == 6);

      CHECK(packwrap.elements.size() == 6);
      CHECK(packwrap.elements[0] == 1);
      CHECK(packwrap.elements[1] == 2);
      CHECK(packwrap.elements[2] == 3);
      CHECK(packwrap.elements[3] == 4);
      CHECK(packwrap.elements[4] == 5);
      CHECK(packwrap.elements[5] == 6);

      matrix_type Mback(3, 3);
      pack.copy_symmetric_to(Mback);
      CHECK(Mback == mat);
    }  // Simple hard-coded test

    auto test = []() {
      const size_t size = *gen::numeric_size<2>().as("Symmetric matrix size");
      matrix_type mat(size, size, false);
      vector_type comp(size * (size + 1) / 2, false);
      for (size_t j = 0, c = 0; j < size; ++j) {
        for (size_t i = j; i < size; ++i, ++c) {
          const auto val = *gen::numeric<scalar_type>().as(
                "Element " + std::to_string(i) + "," + std::to_string(j));
          comp[c] = mat(i, j) = mat(j, i) = val;
        }
      }
      LazyMatrixWrapper<matrix_type> mat_wrap(mat);
      detail::LapackPackedMatrix<scalar_type> pack(mat);
      detail::LapackPackedMatrix<scalar_type> packwrap(mat_wrap);

      // Check sizes and content:
      RC_ASSERT(pack.elements.size() == comp.size());
      RC_ASSERT(packwrap.elements.size() == comp.size());
      RC_ASSERT(comp == vector_type(pack.elements));
      RC_ASSERT(comp == vector_type(packwrap.elements));

      // Transform back and check identity
      matrix_type Mback(size, size, false);
      pack.copy_symmetric_to(Mback);
      RC_ASSERT(Mback == mat);
    };

    REQUIRE(rc::check("LapackPackedMatrix generation and unpacking", test));
  }  // LapackPackedMatrix

  krims::GenMap params1{{LapackEigensolverKeys::prefer_packed_matrices, false}};
  krims::GenMap params2{{LapackEigensolverKeys::prefer_packed_matrices, true}};

  SECTION("Real hermitian normal problems") {
    typedef EigensolverTestProblem<matrix_type, /* Hermitian= */ true> tprob_type;
    TestProblemRunner<tprob_type, DefaultSolveFunctor<LapackEigensolverTraits>> tr;

    SECTION("Run with first parameter set") {
      tr.solve_functor().extra_params = params1;
      tr.run_normal();
    }

    SECTION("Run with second parameter set") {
      tr.solve_functor().extra_params = params2;
      tr.run_normal();
    }
  }  // real hermitian normal problems

  SECTION("Real hermitian generalised problems") {
    typedef EigensolverTestProblem<matrix_type, /* Hermitian= */ true> tprob_type;
    TestProblemRunner<tprob_type, DefaultSolveFunctor<LapackEigensolverTraits>> tr;

    // Run all problems as generalised problems.
    tr.solve_functor().force_generalised = true;

    SECTION("Run with first parameter set") {
      tr.solve_functor().extra_params = params1;
      tr.run_all();
    }

    SECTION("Run with second parameter set") {
      tr.solve_functor().extra_params = params2;
      tr.run_all();
    }
  }  // real hermitian generalised problems

  //  TODO Not yet there
  //  SECTION("Real non-hermitian normal problems") {
  //    typedef EigensolverTestProblem<matrix_type, /* Hermitian= */ false> tprob_type;
  //    TestProblemRunner<tprob_type, DefaultSolveFunctor<LapackEigensolverTraits>> tr;
  //    tr.run_normal();
  //  }  // real hermitian normal problems
  //
  //  SECTION("Real non-hermitian generalised problems") {
  //    typedef EigensolverTestProblem<matrix_type, /* Hermitian= */ false> tprob_type;
  //    TestProblemRunner<tprob_type, DefaultSolveFunctor<LapackEigensolverTraits>> tr;
  //
  //    // Run all problems as generalised problems.
  //    tr.solve_functor().force_generalised = true;
  //    tr.run_all();
  //  }  // real hermitian generalised problems

}  // LapackEigensolver
}  // namespace tests
}  // namespace lazyten

#endif  // LAZYTEN_HAVE_LAPACK
