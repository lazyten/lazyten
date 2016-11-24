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

#include "eigensolver_tests.hh"
#include <linalgwrap/SmallMatrix.hh>
#include <linalgwrap/SmallVector.hh>
#include <linalgwrap/eigensystem.hh>

namespace linalgwrap {
namespace tests {
using namespace rc;

TEST_CASE("eigensystem", "[eigensystem]") {
  using namespace eigensolver_tests;
  typedef SmallMatrix<double> matrix_type;

  /* The filter functor to filter out problems which make no sense for us
   * here*/
  auto filter = [](const EigensolverTestProblemBase<matrix_type>& problem) {
    // We cannot deal with cases where diag is different from A
    // from the eigensystem.hh framework
    if (problem.have_Diag()) return false;

    return true;
  };

  SECTION("Real hermitian normal problems") {
    typedef EigensolverTestProblem<matrix_type,
                                   /* Hermitian= */ true,
                                   /* general= */ false>
          tprob_type;

    /** The solver functor: just call eigensystem_hermitian */
    auto solver = [](const tprob_type& problem) {
      return eigensystem_hermitian(problem.A(), problem.evalues.size(), problem.params);
    };

    TestProblemRunner<tprob_type> runner(solver);
    runner.run_matching(filter);
  }  // real hermitian normal problems

  SECTION("Real hermitian generalised problems") {
    typedef EigensolverTestProblem<matrix_type,
                                   /* Hermitian= */ true,
                                   /* general= */ true>
          tprob_type;

    /** The solver functor: just call eigensystem_hermitian */
    auto solver = [](const tprob_type& problem) {
      return eigensystem_hermitian(problem.A(), problem.B(), problem.evalues.size(),
                                   problem.params);
    };

    TestProblemRunner<tprob_type> runner(solver);
    runner.run_matching(filter);
  }  // real hermitian generalised problems

  SECTION("Real non-hermitian normal problems") {
    typedef EigensolverTestProblem<matrix_type,
                                   /* Hermitian= */ true,
                                   /* general= */ false>
          tprob_type;

    /** The solver functor: just call eigensystem_hermitian */
    auto solver = [](const tprob_type& problem) {
      return eigensystem_hermitian(problem.A(), problem.evalues.size(), problem.params);
    };

    TestProblemRunner<tprob_type> runner(solver);
    runner.run_matching(filter);
  }  // real non-hermitian normal problems

  SECTION("Real non-hermitian generalised problems") {
    typedef EigensolverTestProblem<matrix_type,
                                   /* Hermitian= */ true,
                                   /* general= */ true>
          tprob_type;

    /** The solver functor: just call eigensystem_hermitian */
    auto solver = [](const tprob_type& problem) {
      return eigensystem_hermitian(problem.A(), problem.B(), problem.evalues.size(),
                                   problem.params);
    };

    TestProblemRunner<tprob_type> runner(solver);
    runner.run_matching(filter);
  }  // real non-hermitian generalised problems

}  // eigensystem

}  // tests
}  // linalgwrap
