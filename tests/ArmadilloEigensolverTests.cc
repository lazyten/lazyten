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
#include <lazyten/Armadillo/ArmadilloEigensolver.hh>

#ifdef LAZYTEN_HAVE_ARMADILLO
namespace lazyten {
namespace tests {
using namespace rc;

/** Traits class needed for the tests */
struct ArmadilloEigensolverTraits {
  template <typename Eigenproblem>
  using Solver = ArmadilloEigensolver<Eigenproblem>;
};

TEST_CASE("ArmadilloEigensolver", "[ArmadilloEigensolver]") {
  using namespace eigensolver_tests;
  typedef ArmadilloMatrix<double> matrix_type;

  SECTION("Real hermitian normal problems") {
    typedef EigensolverTestProblem<matrix_type, /* Hermitian= */ true> tprob_type;
    TestProblemRunner<tprob_type, DefaultSolveFunctor<ArmadilloEigensolverTraits>> tr;
    tr.run_normal();
  }  // real hermitian normal problems

  SECTION("Real hermitian generalised problems") {
    typedef EigensolverTestProblem<matrix_type, /* Hermitian= */ true> tprob_type;
    TestProblemRunner<tprob_type, DefaultSolveFunctor<ArmadilloEigensolverTraits>> tr;

    // Run all problems as generalised problems.
    tr.solve_functor().force_generalised = true;
    tr.run_all();
  }  // real hermitian generalised problems

  SECTION("Real non-hermitian normal problems") {
    typedef EigensolverTestProblem<matrix_type, /* Hermitian= */ false> tprob_type;
    TestProblemRunner<tprob_type, DefaultSolveFunctor<ArmadilloEigensolverTraits>> tr;
    tr.run_normal();
  }  // real hermitian normal problems

  SECTION("Real non-hermitian generalised problems") {
    typedef EigensolverTestProblem<matrix_type, /* Hermitian= */ false> tprob_type;
    TestProblemRunner<tprob_type, DefaultSolveFunctor<ArmadilloEigensolverTraits>> tr;

    // Run all problems as generalised problems.
    tr.solve_functor().force_generalised = true;
    tr.run_all();
  }  // real hermitian generalised problems

}  // ArmadilloEigensolver
}  // namespace tests
}  // namespace lazyten

#endif
