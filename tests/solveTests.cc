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

#include "linalgwrap/SmallMatrix.hh"
#include "linalgwrap/SmallVector.hh"
#include "linalgwrap/TestingUtils.hh"
#include "linalgwrap/solve.hh"
#include "rapidcheck_utils.hh"
#include <catch.hpp>

namespace linalgwrap {
namespace tests {

namespace linear_solver_tests {
using namespace krims;
using namespace rc;

template <typename Matrix>
struct SolveTestProblem {
  typedef Matrix matrix_type;
  typedef SmallVector<typename Matrix::scalar_type> vector_type;

  std::string description;         //< A short description of the test
  matrix_type M;                   //< The system matrix
  vector_type rhs;                 //< The rhs vector
  vector_type ref_soln;            //< The expected solution
  NumCompAccuracyLevel tolerance;  //< The accuracy level
  ParameterMap params;             //< Parameters to run the solver with

  SolveTestProblem(std::string description_, matrix_type M_, vector_type rhs_,
                   vector_type ref_soln_)
        : description(description_), M(M_), rhs(rhs_), ref_soln(ref_soln_) {}
};

/** Get the list of all stored real hermitian problems */
template <typename Matrix>
static std::vector<SolveTestProblem<Matrix>> real_hermitian_problems() {
  typedef Matrix matrix_type;
  typedef SolveTestProblem<Matrix> prob_type;
  typedef typename prob_type::vector_type vector_type;

  std::vector<prob_type> res;

  res.push_back({"4x4 Diagonal problem",
                 // system matrix
                 matrix_type{{-3., 0, 0, 0},     //
                             {0, 1., 0, 0},      //
                             {0, 0, -3., 0},     //
                             {0, 0, 0, -100.}},  //
                 // RHS
                 vector_type{1., 2., 3., 4.},
                 // reference solution
                 vector_type{-1. / 3., 2, -1., -0.04}});
  //
  res.push_back({"2x2 Problem",
                 // system matrix
                 matrix_type{{32.05145843625442, -20.956703837816036},    //
                             {-20.956703837816036, 18.964996172649933}},  //
                 // RHS
                 vector_type{-31.948228195051612, -5.145090279241401},
                 // reference solution
                 vector_type{-4.2314116993688815, -4.947089428580436}});
  //
  // The system matrix for the 10x10 problem
  matrix_type M10{
        {-29.77414501699191, 5.023041160944544, -12.299250170229122, -23.004950731622138,
         -3.6892477692559567, 19.036437445495082, 17.064221469147384, -30.938276653698182,
         -31.455418095862406, -16.559794292926625},
        {5.023041160944544, 16.57963908376911, 14.596928349481004, -33.50979554240813,
         32.696998538453414, 2.292918092266511, 30.93791127888244, -4.4034350350539455,
         -29.605717882921688, 4.896970126919982},
        {-12.299250170229122, 14.596928349481004, -22.967147006294496, 22.96144654074645,
         -2.446537082695528, -12.140205299564322, -40.803186849031775, 6.9104926794307175,
         7.870014128430384, -16.687754857817247},
        {-23.004950731622138, -33.50979554240813, 22.96144654074645, 30.802824237784137,
         -8.336047174559965, 7.627125903208707, 12.923036985575564, 7.683114917488595,
         -12.708738923471074, 4.70605814521852},
        {-3.6892477692559567, 32.696998538453414, -2.446537082695528, -8.336047174559965,
         15.794311678824357, -13.667193086106408, -5.681969319870632,
         0.056116494717912246, 6.064236968987672, 15.660441002796148},
        {19.036437445495082, 2.292918092266511, -12.140205299564322, 7.627125903208707,
         -13.667193086106408, 19.154978567158565, -30.560786590359584, 32.269311830082515,
         -31.95760681639232, -16.020671463641037},
        {17.064221469147384, 30.93791127888244, -40.803186849031775, 12.923036985575564,
         -5.681969319870632, -30.560786590359584, 9.84485299717116, 12.882370745924277,
         -5.456425519293546, 11.593271341349407},
        {-30.938276653698182, -4.4034350350539455, 6.9104926794307175, 7.683114917488595,
         0.056116494717912246, 32.269311830082515, 12.882370745924277, 42.263316333443385,
         14.238796644516086, 23.983856059200335},
        {-31.455418095862406, -29.605717882921688, 7.870014128430384, -12.708738923471074,
         6.064236968987672, -31.95760681639232, -5.456425519293546, 14.238796644516086,
         36.038967743911115, -19.81723783017545},
        {-16.559794292926625, 4.896970126919982, -16.687754857817247, 4.70605814521852,
         15.660441002796148, -16.020671463641037, 11.593271341349407, 23.983856059200335,
         -19.81723783017545, -8.409379036253625}};

  vector_type rhs10{-26.368937616170868, -25.400337363622086, -13.173325291668903,
                    -46.995513797484165, -30.864863546306637, 23.247194086851778,
                    -20.759935374594562, -43.96548838519901,  -6.611790128003037,
                    21.34528557985766};

  vector_type sol10{0.9013244219509855,  -17.674562334753404, -14.677410251628974,
                    -5.211675107624927,  11.718173724192823,  0.6134289026917301,
                    -10.322090730773851, -4.220837757412481,  -4.6076227749253125,
                    16.847802793373393};

  res.push_back({"10x10 Problem", M10, rhs10, sol10});

  for (auto& prob : res) {
    prob.M.check_and_set_properties(OperatorProperties::RealSymmetric);
  }

  return res;
}

}  // namespace linear_solver_tests

using namespace linear_solver_tests;

TEST_CASE("solve", "[solve]") {
  typedef SmallMatrix<double> matrix_type;

  // Real hermitian problems
  std::vector<SolveTestProblem<matrix_type>> real_hermitian =
        real_hermitian_problems<matrix_type>();

  SECTION("solve with hermitian reference problems") {
    typedef typename SolveTestProblem<matrix_type>::vector_type vector_type;
    for (const auto& prb : real_hermitian) {
      vector_type res(prb.ref_soln.size());
      solve(prb.M, res, prb.rhs);
      REQUIRE(res == numcomp(prb.ref_soln));
    }
  }  // solve with hermitian reference problems

  SECTION("solve with random hermitian problems") {
    typedef SmallVector<typename matrix_type::scalar_type> vector_type;

    auto test = [&] {
      // Generate the matrix
      auto symm = [](matrix_type m) {
        m.symmetrise();

        // Add one on the diagonal to make matrix less singular:
        for (size_t i = 0; i < m.n_rows(); ++i) {
          m(i, i) += 1.;
        }

        m.check_and_set_properties(OperatorProperties::RealSymmetric);
        return m;
      };

      auto n = *gen::scale(0.8, gen::numeric_size<2>()).as("Matrix size");
      auto M =
            *gen::map(gen::numeric_tensor<matrix_type>(n, n), symm).as("Problem matrix");

      // Generate the vectors and solve
      auto rhs = *gen::numeric_tensor<vector_type>(M.n_rows());
      vector_type sol(M.n_rows());
      solve(M, sol, rhs);

      // Check that we agree:
      RC_ASSERT_NC(M * sol == numcomp(rhs).tolerance(NumCompAccuracyLevel::Sloppy));
    };

    CHECK(rc::check("solve with random hermitian problems", test));
  }  // solve_hermitian with random problems

}  // solve

}  // tests
}  // linalgwrap
