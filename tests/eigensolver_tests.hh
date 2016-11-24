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

#pragma once
#include "EigensolverTestProblemLibrary.hh"
#include "rapidcheck_utils.hh"
#include <catch.hpp>
#include <linalgwrap/Base/Solvers.hh>
#include <linalgwrap/SmallVector.hh>

namespace linalgwrap {
namespace tests {

namespace eigensolver_tests {
using namespace krims;
using namespace rc;

/** Class which runs all or a selected subset of the eigensolver test problems
 * the EigensolverTestProblemLibrary yields for the provided type fo test
 * problems
 *
 * \tparam TestProblem  testproblem to solve
 * \tparam Solver       eigensolver to use
 * */
template <typename Testproblem>
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

  /** Construct a TestProblemRunner from a functor which takes an
   * EigensolverTestProblem and returns the solution to it */
  TestProblemRunner(
        std::function<typename testproblem_type::soln_type(const testproblem_type&)>
              solvefctn)
        : m_solve_testproblem(solvefctn) {}

 private:
  /** Function which takes the test problem and returns the solution*/
  std::function<typename testproblem_type::soln_type(const testproblem_type&)>
        m_solve_testproblem;

  /** \brief Check that the signs of the entries of the vectors from and to
   * agree
   *  if not multiply to by -1.
   *
   *  The idea behind this function is that eigensolvers can only compute
   *  the eigenvector up to the sign and this ambiguity has to be fixed before
   *  comparing the solver results and a reference
   *
   *  \param tolerance   Entries smaller than the numeric tolerance
   *                     (which compare to 0 under this tolerance)
   *                     will not be considered.
   *  \param from        Vector to read the sign from
   *  \param to          Vector to apply the sign to.
   */
  void sign_normalise(const evector_type& from, evector_type& to) const;
};

/** A default functor for the solvefctn in TestProblemRunner, which is usually
 *  well-suited for eigensolvers */
template <typename Testproblem, typename Solver>
struct DefaultSolveFunctor {
  typedef typename Solver::esoln_type esoln_type;
  typedef Testproblem testproblem_type;

  esoln_type operator()(const testproblem_type& problem) {
    Solver solver(problem.params);
    return solver.solve(problem.eigenproblem()).eigensolution();
  }
};

//
// ---------------------------------------------------------------------
//

template <typename Testproblem>
void TestProblemRunner<Testproblem>::sign_normalise(const evector_type& from,
                                                    evector_type& to) const {
  auto itfrom = std::begin(from);
  auto itto = std::begin(to);
  for (; itfrom != std::end(from); ++itfrom, ++itto) {
    // Skip if the element is numerically zero
    if (numcomp(*itfrom).failure_action(NumCompActionType::Return) ==
        Constants<typename evector_type::scalar_type>::zero) {
      continue;
    }

    // TODO The case of complex eigenvalues has not been properly thought
    // through when it comes to sign normalisation: It seems that one could
    // perform any rotation of real and imaginary part on the unit circle,
    // but I am not sure whether this is correct ... mfh
    assert_dbg(!IsComplexNumber<typename evector_type::scalar_type>::value,
               krims::ExcNotImplemented());

    // The sign of the first important element is different:
    if ((std::real(*itfrom) < 0. && std::real(*itto) > 0.) ||
        (std::real(*itfrom) > 0. && std::real(*itto) < 0.)) {
      to *= -1;
      return;
    }
  }
}

template <typename Testproblem>
void TestProblemRunner<Testproblem>::run_matching(
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
        sign_normalise(evec, evec_ref);

        CHECK(evec == numcomp(evec_ref).tolerance(problem.tolerance));
      }
    } catch (const SolverException& e) {
      FAIL("Solver failed with message" + std::string(e.what()));
    }
  }
}

}  // namespace eigensolver_tests
}  // namespace tests
}  // namespace linalgwrap
