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

#pragma once
#include "EigensolverTestProblemLibrary.hh"
#include "rapidcheck_utils.hh"
#include <catch.hpp>
#include <lazyten/Base/Solvers.hh>
#include <lazyten/SmallVector.hh>
#include <lazyten/TestingUtils.hh>

namespace lazyten {
namespace tests {

namespace eigensolver_tests {
using namespace krims;
using namespace rc;

/** A default functor for the solvefctn in TestProblemRunner, which is usually
 *  well-suited for eigensolvers */
template <typename SolverTraits>
struct DefaultSolveFunctor {

  template <typename Testproblem>
  typename Testproblem::esoln_type operator()(const Testproblem& problem) const {
    if (problem.generalised_only() || force_generalised) {
      // Solve as a generalised eigenproblem:
      typedef typename Testproblem::gen_prob_type prob_type;
      typedef typename SolverTraits::template Solver<prob_type> solver_type;
      krims::GenMap fparams{problem.params};
      fparams.update(extra_params);
      return solver_type{fparams}
            .solve(problem.generalised_eigenproblem())
            .eigensolution();
    } else {
      typedef typename Testproblem::prob_type prob_type;
      typedef typename SolverTraits::template Solver<prob_type> solver_type;
      krims::GenMap fparams{problem.params};
      fparams.update(extra_params);
      return solver_type{fparams}.solve(problem.eigenproblem()).eigensolution();
    }
  }

  /** Force all eigenproblems to be generalised problems including the
   *  once where the metric matrix is the identity */
  bool force_generalised = false;

  /** Extra parameters to use for this solver */
  krims::GenMap extra_params{};
};

/** Class which runs all or a selected subset of the eigensolver test problems
 * the EigensolverTestProblemLibrary yields for the provided type fo test
 * problems
 *
 * \tparam TestProblem     testproblem to solve
 * \tparam SolverFunctor   Functor which takes a testproblem and returns
 *                         a solution for it.
 * */
template <typename Testproblem, typename SolveFunctor>
class TestProblemRunner {
 public:
  typedef Testproblem testproblem_type;
  typedef typename testproblem_type::size_type size_type;
  typedef typename testproblem_type::scalar_type scalar_type;
  typedef typename testproblem_type::matrix_type matrix_type;
  typedef typename testproblem_type::evector_type evector_type;
  typedef typename testproblem_type::evalue_type evalue_type;
  typedef EigensolverTestProblemLibrary<testproblem_type> testproblemlib_type;

  /** Run all test problems which match the given predicate (i.e. for which
   * the given predicate returns true */
  void run_matching(std::function<bool(const EigensolverTestProblemBase<matrix_type>&)>
                          predicate) const;

  /** Run all test problems which the library yields */
  void run_all() const {
    run_matching([](const EigensolverTestProblemBase<matrix_type>&) { return true; });
  }

  /** Run all normal (i.e. non-general) test problems the library yields */
  void run_normal() const {
    run_matching([](const EigensolverTestProblemBase<matrix_type>& p) {
      return !p.generalised_only();
    });
  }

  /** Run matching normal (i.e. non-general) test problems the library yields */
  void run_normal_matching(
        std::function<bool(const EigensolverTestProblemBase<matrix_type>&)> predicate)
        const {
    run_matching([predicate](const EigensolverTestProblemBase<matrix_type>& p) {
      return !p.generalised_only() && predicate(p);
    });
  }

  /** Construct a TestProblemRunner from a functor which takes an
   * EigensolverTestProblem and returns the solution to it */
  TestProblemRunner(SolveFunctor solvefctn = SolveFunctor{})
        : m_solve_testproblem(solvefctn) {}

  /** Access to the inner solve functor which is used */
  SolveFunctor& solve_functor() { return m_solve_testproblem; }
  const SolveFunctor& solve_functor() const { return m_solve_testproblem; }

 private:
  /** Function which takes the test problem and returns the solution*/
  SolveFunctor m_solve_testproblem;
};

//
// ---------------------------------------------------------------------
//

template <typename Testproblem, typename SolveFunctor>
void TestProblemRunner<Testproblem, SolveFunctor>::run_matching(
      std::function<bool(const EigensolverTestProblemBase<matrix_type>&)> predicate)
      const {
  const std::vector<testproblem_type> problems = testproblemlib_type::get_all();

  for (auto& problem : problems) {
    // Skip those where predicate does not match.
    if (!predicate(problem)) continue;

    INFO("#");
    INFO("# " + problem.description);
    INFO("#");

    try {
      const auto solution = m_solve_testproblem(problem);

      // TODO do this better throwing a nicer exception which uses the
      // description string
      // and specifies more properly what is actually going wrong and
      // where!

      // Check sizes
      CHECK(problem.n_ep == solution.evalues().size());
      CHECK(problem.n_ep == solution.evectors().n_vectors());

      // Log eigenvalues and eigenvectors, such that we know
      // them in case things fail.
      std::stringstream ss;
      ss << "Got eigenvalues " << std::endl;
      for (const auto& val : solution.evalues()) {
        ss << "  " << val << std::endl;
      }
      ss << std::endl;
      ss << "Got eigenvectors " << std::endl;
      for (const auto& evec : solution.evectors()) {
        ss << "  " << evec << std::endl;
      }
      ss << std::endl;
      ss << "Expected eigenvalues " << std::endl;
      for (const auto& val : problem.evalues) {
        ss << "  " << val << std::endl;
      }
      ss << std::endl;
      ss << "Expected eigenvectors " << std::endl;
      for (const auto& evec : problem.evectors) {
        ss << "  " << evec << std::endl;
      }
      ss << std::endl;
      INFO(ss.str());

      // Check eigenvalues
      SmallVector<evalue_type> evals(solution.evalues());
      SmallVector<evalue_type> evals_ref(problem.evalues);
      CHECK(evals == numcomp(evals_ref).tolerance(problem.tolerance));

      // Check eigenvectors
      for (size_type i = 0; i < problem.n_ep; ++i) {
        const auto& evec = solution.evectors()[i];
        auto evec_ref = problem.evectors[i];

        // Sign-normalise the reference
        // Note that this is necessary due to the ambiguity in the sign
        // in the eigenvectors.
        adjust_phase(evec, evec_ref);

        CHECK(evec == numcomp(evec_ref).tolerance(problem.tolerance));
      }
    } catch (const SolverException& e) {
      FAIL("Solver failed with message" + std::string(e.what()));
    }
  }
}

}  // namespace eigensolver_tests
}  // namespace tests
}  // namespace lazyten
