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

#ifdef LINALGWRAP_HAVE_ARPACK
#include "eigensolver_tests.hh"
#include <linalgwrap/Arpack.hh>
#include <linalgwrap/SmallMatrix.hh>

namespace linalgwrap {
namespace tests {
using namespace rc;

/** Traits class needed for the tests */
struct ArpackEigensolverTraits {
  template <typename Eigenproblem>
  using Solver = ArpackEigensolver<Eigenproblem>;
};

TEST_CASE("ArpackEigensolver", "[ArpackEigensolver]") {
  using namespace eigensolver_tests;
  typedef SmallMatrix<double> matrix_type;

  /* The filter functor to filter out problems which make no sense
   * for us here*/
  auto filter = [](const EigensolverTestProblemBase<matrix_type>& problem) {
    // Arpack cannot deal with cases where the dimensionality is equal
    // to the number of eigenpairs to be seeked
    if (problem.n_ep >= problem.dim) return false;
    //
    // Arpack is bad at finding the smallest magnitude eigenvalues
    if (problem.params.at<std::string>(EigensolverBaseKeys::which, "SR") ==
        std::string("SM"))
      return false;
    if (problem.params.at<std::string>(EigensolverBaseKeys::which, "SR") ==
        std::string("SR"))
      return false;

    return true;
  };

  SECTION("Real hermitian normal problems") {
    typedef EigensolverTestProblem<matrix_type, /* Hermitian= */ true> tprob_type;
    TestProblemRunner<tprob_type, DefaultSolveFunctor<ArpackEigensolverTraits>> tr;
    tr.run_normal_matching(filter);
  }  // real hermitian normal problems

  SECTION("Real hermitian generalised problems") {
    typedef EigensolverTestProblem<matrix_type, /* Hermitian= */ true> tprob_type;
    TestProblemRunner<tprob_type, DefaultSolveFunctor<ArpackEigensolverTraits>> tr;

    // Run all problems as generalised problems.
    tr.solve_functor().force_generalised = true;
    tr.run_matching(filter);
  }  // real hermitian generalised problems

  // TODO test real non-hermitian problems
  // TODO test real non-hermitian generalised problems

  SECTION("Check that providing a guess reduces the number of steps needed") {
    // TODO Get this into the general testing library

    typedef EigensolverTestProblem<matrix_type, /* Hermitian= */ true> tprob_type;
    auto allprobs = EigensolverTestProblemLibrary<tprob_type>::get_all();

    for (const tprob_type& testproblem : allprobs) {
      if (!filter(testproblem)) continue;

      INFO("#");
      INFO("# " + testproblem.description);
      INFO("#");

      // Assume the last problem is a generalised eigenproblem
      auto prob = testproblem.generalised_eigenproblem();
      ArpackEigensolverState<decltype(prob)> guess_state{prob};

      // Set the results to expect:
      guess_state.eigensolution().evalues() = testproblem.evalues;

      auto& evecs = guess_state.eigensolution().evectors();
      evecs.clear();
      evecs.reserve(testproblem.evectors.size());
      for (auto& vec : testproblem.evectors) {
        evecs.push_back(typename tprob_type::evector_type{vec});
      }

      ArpackEigensolver<decltype(prob)> solver{testproblem.params};
      auto ret = solver.solve_with_guess(prob, guess_state);
      size_t n_iter = ret.n_iter();

      // TODO or compare that we have less than normal solve
      CHECK(n_iter < 10);

      // Check eigenvalues
      typedef typename tprob_type::evalue_type evalue_type;
      SmallVector<evalue_type> evals(ret.eigensolution().evalues());
      SmallVector<evalue_type> evals_ref(testproblem.evalues);
      CHECK(evals == numcomp(evals_ref).tolerance(testproblem.tolerance));

      // Check eigenvectors
      for (size_t i = 0; i < testproblem.n_ep; ++i) {
        const auto& evec = ret.eigensolution().evectors()[i];
        auto evec_ref = testproblem.evectors[i];

        // Sign-normalise the reference
        // Note that this is necessary due to the ambiguity in the sign
        // in the eigenvectors.
        adjust_phase(evec, evec_ref);

        CHECK(evec == numcomp(evec_ref).tolerance(testproblem.tolerance));
      }
    }
  }

}  // ArpackEigensolver

}  // tests
}  // linalgwrap
#endif
