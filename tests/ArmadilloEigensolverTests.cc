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

#ifdef LINALGWRAP_HAVE_ARMADILLO
#include "eigensolver_tests.hh"
#include <linalgwrap/Armadillo/ArmadilloEigensolver.hh>

namespace linalgwrap {
namespace tests {
using namespace rc;

TEST_CASE("ArmadilloEigensolver", "[ArmadilloEigensolver]") {
  using namespace eigensolver_tests;
  typedef ArmadilloMatrix<double> matrix_type;

  SECTION("Real hermitian normal problems") {
    typedef EigensolverTestProblem<matrix_type,
                                   /* Hermitian= */ true,
                                   /* general= */ false>
          tprob_type;
    typedef ArmadilloEigensolver<typename tprob_type::prob_type> solver_type;

    TestProblemRunner<tprob_type>(DefaultSolveFunctor<tprob_type, solver_type>())
          .run_all();
  }  // real hermitian normal problems

  SECTION("Real hermitian generalised problems") {
    typedef EigensolverTestProblem<matrix_type,
                                   /* Hermitian= */ true,
                                   /* general= */ true>
          tprob_type;
    typedef ArmadilloEigensolver<typename tprob_type::prob_type> solver_type;

    TestProblemRunner<tprob_type>(DefaultSolveFunctor<tprob_type, solver_type>())
          .run_all();
  }  // real hermitian generalised problems

  SECTION("Real non-hermitian normal problems") {
    typedef EigensolverTestProblem<matrix_type,
                                   /* Hermitian= */ false,
                                   /* general= */ false>
          tprob_type;
    typedef ArmadilloEigensolver<typename tprob_type::prob_type> solver_type;

    TestProblemRunner<tprob_type>(DefaultSolveFunctor<tprob_type, solver_type>())
          .run_all();
  }  // real non-hermitian normal problems

  SECTION("Real non-hermitian generalised problems") {
    typedef EigensolverTestProblem<matrix_type,
                                   /* Hermitian= */ false,
                                   /* general= */ true>
          tprob_type;
    typedef ArmadilloEigensolver<typename tprob_type::prob_type> solver_type;

    TestProblemRunner<tprob_type>(DefaultSolveFunctor<tprob_type, solver_type>())
          .run_all();
  }  // real non-hermitian generalised problems

}  // ArmadilloEigensolver
}  // tests
}  // linalgwrap

#endif
