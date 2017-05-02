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

#include "eigensolver_tests.hh"
#include <linalgwrap/SmallMatrix.hh>
#include <linalgwrap/SmallVector.hh>
#include <linalgwrap/eigensystem.hh>

namespace linalgwrap {
namespace tests {
using namespace rc;

struct HermitianSolveFunctor {
  template <typename Testproblem>
  typename Testproblem::esoln_type operator()(const Testproblem& problem) const {
    if (problem.generalised_only()) {
      return eigensystem_hermitian(problem.A(), problem.B(), problem.n_ep,
                                   problem.params);
    } else {
      return eigensystem_hermitian(problem.A(), problem.n_ep, problem.params);
    }
  }
};

struct SolveFunctor {
  template <typename Testproblem>
  typename Testproblem::esoln_type operator()(const Testproblem& problem) const {
    if (problem.generalised_only()) {
      return eigensystem(problem.A(), problem.B(), problem.n_ep, problem.params);
    } else {
      return eigensystem(problem.A(), problem.n_ep, problem.params);
    }
  }
};

TEST_CASE("eigensystem", "[eigensystem]") {
  using namespace eigensolver_tests;
  typedef SmallMatrix<double> matrix_type;

  /* The filter functor to filter out problems which make no sense for eigensystem() */
  auto filter = [](const EigensolverTestProblemBase<matrix_type>& problem) {
    // We cannot deal with cases where diag is different from A
    // from the eigensystem.hh framework
    return !problem.have_Diag();
  };

  SECTION("Real hermitian problems") {
    typedef EigensolverTestProblem<matrix_type, /* Hermitian= */ true> tprob_type;
    TestProblemRunner<tprob_type, HermitianSolveFunctor>{}.run_matching(filter);
  }  // real hermitian normal problems

  SECTION("Real non-hermitian problems") {
    typedef EigensolverTestProblem<matrix_type, /* Hermitian= */ false> tprob_type;
    TestProblemRunner<tprob_type, SolveFunctor>{}.run_matching(filter);
  }  // real hermitian normal problems
}  // eigensystem

}  // namespace tests
}  // namespace linalgwrap
