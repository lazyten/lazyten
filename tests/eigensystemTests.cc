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
#include <linalgwrap/LazyMatrixWrapper.hh>
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

  SECTION("Dummy test for a block-diagonal problem") {
    using linalgwrap::EigensystemSolverKeys;

    auto testable = [] {
      typedef LazyMatrixWrapper<matrix_type> lazy_type;

      auto addone = [](std::vector<double> v) {
        for (auto& e : v) e += 1;
        return v;
      };
      auto evals1 = *gen::map(gen::numeric_container<std::vector<double>>(), addone)
                           .as("The first diagonal block");
      auto evals2 = *gen::map(gen::numeric_container<std::vector<double>>(), addone)
                           .as("The second diagonal block");

      auto ones1 = std::vector<double>(evals1.size(), 1);
      auto ones2 = std::vector<double>(evals2.size(), 1);

      matrix_type mm1(evals1.size(), evals1.size());
      matrix_type mm2(evals2.size(), evals2.size());
      matrix_type bb1(ones1.size(), ones1.size());
      matrix_type bb2(ones2.size(), ones2.size());

      for (size_t i = 0; i < evals1.size(); ++i) mm1(i, i) = evals1[i];
      for (size_t i = 0; i < evals2.size(); ++i) mm2(i, i) = evals2[i];
      for (size_t i = 0; i < ones1.size(); ++i) bb1(i, i) = ones1[i];
      for (size_t i = 0; i < ones2.size(); ++i) bb2(i, i) = ones2[i];

      BlockDiagonalMatrix<lazy_type, 2> diag{
            {{lazy_type(std::move(mm1)), lazy_type(std::move(mm2))}}};
      BlockDiagonalMatrix<lazy_type, 2> bdiag{
            {{lazy_type(std::move(bb1)), lazy_type(std::move(bb2))}}};

      const size_t n_ep = std::min(evals1.size(), evals2.size());
      const size_t n_ep2 = n_ep / 2;
      const size_t n_ep1 = n_ep - n_ep2;

      krims::GenMap params{{EigensystemSolverKeys::method, "lapack"},
                           {EigensystemSolverKeys::which, "SR"}};

      const bool explict_n_per_blocks =
            *rc::gen::arbitrary<bool>().as("Explict n_ep_per_blocks");
      RC_CLASSIFY(explict_n_per_blocks, "Use explicit n per blocks");
      if (explict_n_per_blocks) {
        std::array<size_t, 2> npb{{n_ep1, n_ep2}};
        params.update("n_ep_per_block", std::move(npb));
      }

      typedef decltype(eigensystem_hermitian(diag, n_ep, params)) sol_type;
      sol_type sol;
      if (*gen::arbitrary<bool>().as("Solve pseudo-generalised problem")) {
        sol = eigensystem_hermitian(diag, bdiag, n_ep, params);
        RC_CLASSIFY(true, "Solve pseudo-general eigenproblem");
      } else {
        sol = eigensystem_hermitian(diag, n_ep, params);
        RC_CLASSIFY(true, "Solve non-general eigenproblem");
      }
      RC_ASSERT(sol.evalues().size() == n_ep);

      // Sort by size:
      std::vector<size_t> idcs1 = krims::argsort(evals1.begin(), evals1.end());
      std::vector<size_t> idcs2 = krims::argsort(evals2.begin(), evals2.end());

      std::stringstream ss;
      ss << "Got eigenvalues " << std::endl;
      for (const auto& val : sol.evalues()) {
        ss << "  " << std::setw(15) << val << std::endl;
      }
      ss << std::endl;
      ss << "Got eigenvectors " << std::endl;
      for (const auto& vec : sol.evectors()) {
        ss << "  " << std::setw(15) << vec << std::endl;
      }
      ss << std::endl;
      ss << "Expected eigenvalues (block 1)" << std::endl;
      for (size_t i = 0; i < idcs1.size(); ++i) {
        if (i == n_ep1) ss << "              ----" << std::endl;
        ss << "  " << std::setw(15) << evals1[idcs1[i]] << std::endl;
      }
      ss << std::endl;
      ss << "Expected eigenvalues (block 2)" << std::endl;
      for (size_t i = 0; i < idcs2.size(); ++i) {
        if (i == n_ep2) ss << "              ----" << std::endl;
        ss << "  " << std::setw(15) << evals2[idcs2[i]] << std::endl;
      }
      ss << std::endl;
      RC_LOG(ss.str());

      for (size_t i = 0; i < n_ep1; ++i) {
        // Check eval:
        RC_ASSERT(sol.evalues()[i] == evals1[idcs1[i]]);
      }
      for (size_t i = n_ep1; i < n_ep; ++i) {
        // Check eval:
        RC_ASSERT(sol.evalues()[i] == evals2[idcs2[i - n_ep1]]);
      }

      for (size_t i = 0; i < n_ep1; ++i) {
        // Find the non-zero element
        const auto it = std::find_if(sol.evectors()[i].begin(), sol.evectors()[i].end(),
                                     [](double& v) { return std::abs(v) > 1e-12; });
        RC_ASSERT(it != std::end(sol.evectors()[i]));
        const size_t idx = static_cast<size_t>(it - std::begin(sol.evectors()[i]));
        RC_ASSERT(idx < evals1.size());

        for (auto itt = sol.evectors()[i].begin(); itt != sol.evectors()[i].end();
             ++itt) {
          if (it == itt) continue;
          RC_ASSERT(*itt == 0);
        }

        // We cannot check a stronger property, since the
        // degenerate eigenspaces might be interchanged
        RC_ASSERT(evals1[idcs1[i]] == evals1[idx]);
      }

      for (size_t i = n_ep1; i < n_ep; ++i) {
        // Find the non-zero element (kind of code duplication)
        const auto it = std::find_if(sol.evectors()[i].begin(), sol.evectors()[i].end(),
                                     [](double& v) { return std::abs(v) > 1e-12; });
        RC_ASSERT(it != std::end(sol.evectors()[i]));
        const size_t idx =
              static_cast<size_t>(it - std::begin(sol.evectors()[i])) - evals1.size();
        RC_ASSERT(idx < evals2.size());

        for (auto itt = sol.evectors()[i].begin(); itt != sol.evectors()[i].end();
             ++itt) {
          if (it == itt) continue;
          RC_ASSERT(*itt == 0);
        }

        // We cannot check a stronger property, since the
        // degenerate eigenspaces might be interchanged
        RC_ASSERT(evals2[idcs2[i - n_ep1]] == evals2[idx]);
      }

      // TODO Better test for generalised version of these problems, too
    };

    // TODO Test version where n_ep_per_block is specified

    // TODO Test general block-diagonal problems
    //      (e.g. by combining two library problems together)

    CHECK(rc::check("Test block-diagonal problems", testable));
  }  //
}  // eigensystem

}  // namespace tests
}  // namespace linalgwrap
