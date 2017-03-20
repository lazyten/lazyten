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

#ifdef LINALGWRAP_HAVE_LAPACK
#include "eigensolver_tests.hh"
#include <linalgwrap/Lapack/LapackEigensolver.hh>

namespace linalgwrap {
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
    matrix_type M{{1, 2, 3}, {2, 4, 5}, {3, 5, 6}};

    detail::LapackPackedMatrix<double> pack(M);

    CHECK(pack.elements.size() == 6);
    CHECK(pack.elements[0] == 1);
    CHECK(pack.elements[1] == 2);
    CHECK(pack.elements[2] == 3);
    CHECK(pack.elements[3] == 4);
    CHECK(pack.elements[4] == 5);
    CHECK(pack.elements[5] == 6);

    matrix_type Mback(3, 3);
    pack.copy_symmetric_to(Mback);
    CHECK(Mback == M);
  }

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
}  // tests
}  // linalgwrap

#endif  // LINALGWRAP_HAVE_LAPACK
