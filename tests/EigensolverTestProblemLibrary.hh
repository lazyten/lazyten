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
#include <linalgwrap/Base/Solvers.hh>
#include <linalgwrap/TestingUtils.hh>

namespace linalgwrap {
namespace tests {

namespace eigensolver_tests {
using namespace krims;

/** Base class for describing test problems for testing eigensolvers */
template <typename Matrix>
class EigensolverTestProblemBase {
 public:
  typedef Matrix matrix_type;
  typedef typename matrix_type::size_type size_type;

  std::string description;  //< A short description of the test

  //! The number of eigenpairs to be obtained
  size_type n_ep;

  //! The dimensionality of the eigenproblem (number of columns of A)
  size_type dim;

  //! The accuracy level to be used when comparing results
  NumCompAccuracyLevel tolerance;

  /** The parameters which should be used when
      running the solver. */
  ParameterMap params;

  /** Matrix getters */
  ///@{
  matrix_type& A() { return *A_ptr; }
  const matrix_type& A() const { return *A_ptr; }

  matrix_type& Diag() { return *Diag_ptr; }
  const matrix_type& Diag() const { return *Diag_ptr; }
  ///@}

  /** Do we have a special Diag matrix, which the eigensolver should
   * diagonalise instead of A */
  bool have_Diag() const { return Diag_ptr != nullptr; }

  EigensolverTestProblemBase(std::string description_, size_type n_ep_,
                             std::shared_ptr<matrix_type> A_ptr_,
                             std::shared_ptr<matrix_type> Diag_ptr_)
        : description(description_),
          n_ep(n_ep_),
          dim(A_ptr_->n_cols()),
          tolerance(NumCompAccuracyLevel::Default),
          params(),
          A_ptr(A_ptr_),
          Diag_ptr(Diag_ptr_) {}

 protected:
  //! The matrix A of $Ax =\lambda B x$ or $Ax = \lambda x$
  std::shared_ptr<matrix_type> A_ptr;

  //! The matrix to diagonalise (or nullptr if equal to A)
  std::shared_ptr<matrix_type> Diag_ptr;
};

/** Class describing a test problem for testing eigensolvers
 *
 * This is the case for generalised eigenproblems */
template <typename Matrix, bool isHermitian, bool isGeneral>
class EigensolverTestProblem : public EigensolverTestProblemBase<Matrix> {
  static_assert(isGeneral, "Something went wrong with the template specialisation");

 public:
  typedef Matrix matrix_type;

  typedef Eigenproblem<isHermitian, matrix_type, matrix_type, matrix_type> prob_type;

  typedef EigensolutionTypeFor<isHermitian, Matrix> soln_type;
  typedef typename soln_type::evector_type evector_type;
  typedef typename soln_type::evalue_type evalue_type;
  typedef typename prob_type::size_type size_type;
  typedef typename prob_type::scalar_type scalar_type;

  /** The eigenvalues ordered by magnitude (complex case)
   *  or actual value (real case) */
  std::vector<evalue_type> evalues;

  //! The according eigenvectors
  std::vector<evector_type> evectors;

  /** Get the eigenproblem object representing this test problem.
   */
  prob_type eigenproblem() const {
    if (base_type::Diag_ptr == nullptr) {
      return prob_type(base_type::A(), B(), base_type::n_ep);
    } else {
      return prob_type(base_type::A(), B(), base_type::n_ep, base_type::Diag());
    }
  }

  EigensolverTestProblem(std::string description_, matrix_type A_, matrix_type B_,
                         std::vector<evalue_type> evalues_,
                         std::vector<evector_type> evectors_, matrix_type Diag_)
        : base_type{description_, evalues_.size(),
                    std::make_shared<matrix_type>(std::move(A_)),
                    std::make_shared<matrix_type>(std::move(Diag_))},
          evalues(evalues_),
          evectors(evectors_),
          B_ptr{new matrix_type(std::move(B_))} {}

  EigensolverTestProblem(std::string description_, matrix_type A_, matrix_type B_,
                         std::vector<evalue_type> evalues_,
                         std::vector<evector_type> evectors_)
        : base_type{description_, evalues_.size(),
                    std::make_shared<matrix_type>(std::move(A_)), nullptr},
          evalues(evalues_),
          evectors(evectors_),
          B_ptr{new matrix_type(std::move(B_))} {}

  /** Matrix getters */
  ///@{
  matrix_type& B() { return *B_ptr; }
  const matrix_type& B() const { return *B_ptr; }
  ///@}

 protected:
  typedef EigensolverTestProblemBase<Matrix> base_type;

  //! The matrix B of $Ax =\lambda B x$ ( or nullptr if non-general problem)
  std::shared_ptr<matrix_type> B_ptr;
};

/** Class describing a test problem for testing eigensolvers
 *
 * This is the case for normal eigenproblems */
template <typename Matrix, bool isHermitian>
class EigensolverTestProblem<Matrix, isHermitian, false>
      : public EigensolverTestProblemBase<Matrix> {
 public:
  typedef Matrix matrix_type;

  typedef Eigenproblem<isHermitian, matrix_type, void, matrix_type> prob_type;

  typedef EigensolutionTypeFor<isHermitian, Matrix> soln_type;
  typedef typename soln_type::evector_type evector_type;
  typedef typename soln_type::evalue_type evalue_type;
  typedef typename prob_type::size_type size_type;
  typedef typename prob_type::scalar_type scalar_type;

  //! The eigenvalues ordered by magnitude (complex case)
  //  or ordered algebraically (real case)
  //  This number also gives the number of eigenvalues to be seeked
  std::vector<evalue_type> evalues;

  //! The according eigenvectors
  std::vector<evector_type> evectors;

  /** Get the eigenproblem object representing this test problem.
   */
  prob_type eigenproblem() const {
    if (base_type::Diag_ptr == nullptr) {
      return prob_type(base_type::A(), base_type::n_ep);
    } else {
      return prob_type(base_type::A(), base_type::n_ep, base_type::Diag());
    }
  }

  EigensolverTestProblem(std::string description_, matrix_type A_,
                         std::vector<evalue_type> evalues_,
                         std::vector<evector_type> evectors_, matrix_type Diag_)
        : base_type{description_, evalues_.size(),
                    std::make_shared<matrix_type>(std::move(A_)),
                    std::make_shared<matrix_type>(std::move(Diag_))},
          evalues(evalues_),
          evectors(evectors_) {}

  EigensolverTestProblem(std::string description_, matrix_type A_,
                         std::vector<evalue_type> evalues_,
                         std::vector<evector_type> evectors_)
        : base_type{description_, evalues_.size(),
                    std::make_shared<matrix_type>(std::move(A_)), nullptr},
          evalues(evalues_),
          evectors(evectors_) {}

 protected:
  typedef EigensolverTestProblemBase<Matrix> base_type;
};

/** Struct to obtain a list of library test problems */
template <typename EigensolverTestProblem,
          bool isComplex =
                IsComplexNumber<typename EigensolverTestProblem::scalar_type>::value>
struct EigensolverTestProblemLibrary {};

template <typename Matrix>
struct EigensolverTestProblemLibrary<
      EigensolverTestProblem<Matrix, /*Hermitian=*/true, /*generalised=*/false>,
      /*complex=*/false> {
  typedef EigensolverTestProblem<Matrix, /*Hermitian=*/true,
                                 /*generalised=*/false>
        testproblem_type;

  /** Get all Hermitian non-generalised test problems */
  static std::vector<testproblem_type> get_all();
};

template <typename Matrix>
struct EigensolverTestProblemLibrary<
      EigensolverTestProblem<Matrix, /*Hermitian=*/true, /*generalised=*/true>,
      /*complex=*/false> {
  typedef EigensolverTestProblem<Matrix, /*Hermitian=*/true,
                                 /*generalised=*/true>
        testproblem_type;

  /** Get all Hermitian generalised test problems */
  static std::vector<testproblem_type> get_all();
};

template <typename Matrix>
struct EigensolverTestProblemLibrary<EigensolverTestProblem<Matrix, /*Hermitian=*/false,
                                                            /*generalised=*/false>,
                                     /*complex=*/false> {
  typedef EigensolverTestProblem<Matrix, /*Hermitian=*/false,
                                 /*generalised=*/false>
        testproblem_type;

  /** Get all non-Hermitian test problems */
  static std::vector<testproblem_type> get_all();
};

template <typename Matrix>
struct EigensolverTestProblemLibrary<
      EigensolverTestProblem<Matrix, /*Hermitian=*/false, /*generalised=*/true>,
      /*complex=*/false> {
  typedef EigensolverTestProblem<Matrix, /*Hermitian=*/false,
                                 /*generalised=*/true>
        testproblem_type;

  /** Get all non-Hermitian generalised test problems */
  static std::vector<testproblem_type> get_all();
};

//
// ----------------------------------------------------------------------------------
//

template <typename Matrix>
std::vector<EigensolverTestProblem<Matrix, true, false>> EigensolverTestProblemLibrary<
      EigensolverTestProblem<Matrix, /*Hermitian=*/true, /*generalised=*/false>,
      /*complex=*/false>::get_all() {
  typedef typename testproblem_type::matrix_type matrix_type;
  typedef typename testproblem_type::evalue_type evalue_type;
  typedef typename testproblem_type::evector_type evector_type;

  std::vector<testproblem_type> res;

  //
  // Simple near-diagonal problem (all eigenvalues)
  //
  {
    matrix_type A{{-3, 2, 0, 0},  //
                  {2, 1, 0, 0},   //
                  {0, 0, -2, 0},  //
                  {0, 0, 0, -100}};

    std::vector<evalue_type> evalues{-100.00000000000000000, -3.8284271247461900976,
                                     -2.0000000000000000000, 1.8284271247461900976};
    std::vector<evector_type> evectors{
          evector_type{0, 0, 0, 1.0000000000000000000},
          evector_type{-0.92387953251128675613, 0.38268343236508977173, 0, 0},
          evector_type{0, 0, 1.0000000000000000000, 0},
          evector_type{-0.38268343236508977173, -0.92387953251128675613, 0, 0}};
    res.push_back(testproblem_type{"4x4 near-diagonal problem (all eigvals)", A, evalues,
                                   evectors});
  }

  //
  // Simple near-diagonal problem (2 largest eigenvalues)
  //
  {
    matrix_type A{{-3, 2, 0, 0},  //
                  {2, 1, 0, 0},   //
                  {0, 0, -2, 0},  //
                  {0, 0, 0, -100}};

    std::vector<evalue_type> evalues{-100.00000000000000000, -3.8284271247461900976};
    std::vector<evector_type> evectors{
          evector_type{0, 0, 0, 1.0000000000000000000},
          evector_type{-0.92387953251128675613, 0.38268343236508977173, 0, 0}};
    testproblem_type prob{"4x4 near-diagonal problem (2 eigvals)", A, evalues, evectors};
    prob.params.update({{EigensolverBaseKeys::which, "LM"}});
    res.push_back(std::move(prob));
  }

  //
  // Dense 4x4 problem (all eigenvalues)
  //
  {
    matrix_type A{{-31.295016419381554, 29.88921976956469, -22.819339045315687,
                   -2.9748571479165946},  //
                  {29.88921976956469, 12.956259716168347, -5.124711961884337,
                   37.53664340401812},  //
                  {-22.819339045315687, -5.124711961884337, -1.2738042733796533,
                   -19.060589257868173},  //
                  {-2.9748571479165946, 37.53664340401812, -19.060589257868173,
                   9.006913232655023}};  //

    std::vector<evalue_type> evalues{-62.03057355556985, -9.352391446701581,
                                     -0.35870676156019954, 61.13602401989404};

    std::vector<evector_type> evectors{
          evector_type{-0.7351645834700887, 0.4509543660677544, -0.35213494949439583,
                       -0.36355766018646946},
          evector_type{-0.5285029987514623, -0.47396413059156434, -0.20062255452199496,
                       0.6751245617220307},
          evector_type{0.31688501798857666, -0.3444942470129464, -0.8487577663019811,
                       -0.24600376682366906},
          evector_type{0.2824915943817493, 0.6732918211384983, -0.33964955778525674,
                       0.5928868362410665}};

    testproblem_type prob{"4x4 dense problem (all eigvals)", A, evalues, evectors};
    prob.tolerance = NumCompAccuracyLevel::Lower;
    res.push_back(std::move(prob));
  }

  //
  // Dense 4x4 problem (one eigenvalue)
  //
  {
    matrix_type A{
          {-31.295016419381554, 29.88921976956469, -22.819339045315687,
           -2.9748571479165946},
          {29.88921976956469, 12.956259716168347, -5.124711961884337, 37.53664340401812},
          {-22.819339045315687, -5.124711961884337, -1.2738042733796533,
           -19.060589257868173},
          {-2.9748571479165946, 37.53664340401812, -19.060589257868173,
           9.006913232655023}};

    std::vector<evalue_type> evalues{61.13602401989404};

    std::vector<evector_type> evectors{
          evector_type{0.2824915943817493, 0.6732918211384983, -0.33964955778525674,
                       0.5928868362410665}};

    testproblem_type prob{"4x4 dense problem (one eigval)", A, evalues, evectors};
    prob.params.update({{EigensolverBaseKeys::which, "LR"}});
    res.push_back(std::move(prob));
  }

  //
  // Dense 4x4 problem (one eigenvalue)
  //
  {
    matrix_type A{
          {-31.295016419381554, 29.88921976956469, -22.819339045315687,
           -2.9748571479165946},
          {29.88921976956469, 12.956259716168347, -5.124711961884337, 37.53664340401812},
          {-22.819339045315687, -5.124711961884337, -1.2738042733796533,
           -19.060589257868173},
          {-2.9748571479165946, 37.53664340401812, -19.060589257868173,
           9.006913232655023}};

    std::vector<evalue_type> evalues{-0.35870676156019954};

    std::vector<evector_type> evectors{
          evector_type{0.31688501798857666, -0.3444942470129464, -0.8487577663019811,
                       -0.24600376682366906}};

    testproblem_type prob{"4x4 dense problem (one eigval)", A, evalues, evectors};
    prob.params.update({{EigensolverBaseKeys::which, "SM"}});
    prob.tolerance = NumCompAccuracyLevel::Lower;
    res.push_back(std::move(prob));
  }

  //
  // Dense 10x10 problem (4 eigenvalues)
  //
  {
    matrix_type A{
          {-14.57424130732258, 4.860830701921898, 11.735331819442138, -24.293371826576752,
           24.05822629417608, -1.2204720017617205, -15.85018059434573, 3.4557346760715433,
           42.9044463934365, -5.9197963799792745},
          {4.860830701921898, -13.045142289081411, -15.01609004969778, -48.945232444559,
           31.458600992124246, 6.3661425453486515, -5.955683573936213, 3.5484322265099735,
           -16.203738772858628, -16.692386392871313},
          {11.735331819442138, -15.01609004969778, -7.861417280209906,
           -26.931440705592664, -41.12529185654063, -0.8852914039228494,
           -23.558876235231374, -8.994913532732781, -13.861997648927307,
           -22.245581697162187},
          {-24.293371826576752, -48.945232444559, -26.931440705592664, -44.46223926627806,
           11.886486301489185, -24.749366718628586, 30.045652088708863, 3.740636150863601,
           -6.772251960789731, -2.6132506338280805},
          {24.05822629417608, 31.458600992124246, -41.12529185654063, 11.886486301489185,
           7.41257485052347, -25.60373428454045, -14.195790040972724, 17.465798498724624,
           -11.321141191871085, 7.556424052704543},
          {-1.2204720017617205, 6.3661425453486515, -0.8852914039228494,
           -24.749366718628586, -25.60373428454045, 32.107733199027166,
           -19.316181762535123, 12.915914066842973, 3.5976205623419446,
           -10.252025228591862},
          {-15.85018059434573, -5.955683573936213, -23.558876235231374,
           30.045652088708863, -14.195790040972724, -19.316181762535123,
           -46.612962833455484, 5.454204246123318, -23.591724235569608,
           3.1674507141161925},
          {3.4557346760715433, 3.5484322265099735, -8.994913532732781, 3.740636150863601,
           17.465798498724624, 12.915914066842973, 5.454204246123318, 39.654346833529644,
           -14.160199907180676, 3.757319654298165},
          {42.9044463934365, -16.203738772858628, -13.861997648927307, -6.772251960789731,
           -11.321141191871085, 3.5976205623419446, -23.591724235569608,
           -14.160199907180676, 22.352153940761923, 20.801864534114983},
          {-5.9197963799792745, -16.692386392871313, -22.245581697162187,
           -2.6132506338280805, 7.556424052704543, -10.252025228591862,
           3.1674507141161925, 3.757319654298165, 20.801864534114983,
           -2.228736245292424}};

    std::vector<evalue_type> evalues{-105.85950499281432, -87.66911187239057,
                                     68.36805117811491, 83.28666165536679};

    std::vector<evector_type> evectors{
          evector_type{0.22403180380813342, 0.463087112893307, 0.0697485193937882,
                       0.711765589477168, -0.2809517658930841, 0.02540349394174031,
                       -0.32873633731346297, 0.0005293705063164483, -0.08490172855981404,
                       0.1703806648492464},
          evector_type{-0.2151910176148905, 0.24612449261040636, 0.5633618884166837,
                       0.23325910167946295, 0.35017995192436724, 0.19239521704944448,
                       0.4995110625599768, -0.023469756097665342, 0.32367723858606845,
                       0.0828138280444607},
          evector_type{0.4705794400224978, 0.322782767196029, -0.23107500311145912,
                       -0.24558510439353567, 0.47674486539774374, 0.09558640881000315,
                       -0.2495820142033457, 0.22251976212644745, 0.4410315332268127,
                       0.133869740356561},
          evector_type{0.18747075925995604, -0.03728759548794888, 0.4128610868173742,
                       -0.3210021970321167, -0.4622338545647539, 0.4427411053904194,
                       -0.25950637523118303, -0.3026407483904913, 0.3177288444600617,
                       -0.14283290746786131}};

    testproblem_type prob{"10x10 dense problem (4 eigenvalues)", A, evalues, evectors};
    prob.params.update({{EigensolverBaseKeys::which, "LM"}});
    res.push_back(std::move(prob));
  }

  return res;
}

template <typename Matrix>
std::vector<EigensolverTestProblem<Matrix, true, true>> EigensolverTestProblemLibrary<
      EigensolverTestProblem<Matrix, /*Hermitian=*/true, /*generalised=*/true>,
      /*complex=*/false>::get_all() {
  typedef typename testproblem_type::matrix_type matrix_type;
  typedef typename testproblem_type::evalue_type evalue_type;
  typedef typename testproblem_type::evector_type evector_type;

  std::vector<testproblem_type> res;

  //
  // Simple generalised 4x4 problem (all eigenpairs)
  //
  {
    matrix_type A{{-3, 2, 0, 0},  //
                  {2, 1, 0, 0},   //
                  {0, 0, -2, 0},  //
                  {0, 0, 0, -100}};

    matrix_type B{{2, 0, 0, 0},  //
                  {0, 1, 0, 0},  //
                  {0, 0, 3, -1},
                  {0, 0, -1, 6}};

    std::vector<evalue_type> evalues{
          -17.687810469499261500914145, -2.137458608817687011338649,
          -0.665130706971330476484638, 1.637458608817687455427858};

    std::vector<evector_type> evectors{
          evector_type{0, 0, 0.145505916561683712462383, 0.420065072649349724454027},
          evector_type{-0.644644510279216098602717, 0.410934192704548451047231, 0, 0},
          evector_type{0, 0, -0.575994101285883575158664, 0.003990360308374735556769},
          evector_type{-0.290574354282805646931820, -0.911665009346229515685422, 0, 0}};

    testproblem_type prob{"4x4 near-diagonal problem (all eigenpairs)", A, B, evalues,
                          evectors};
    res.push_back(std::move(prob));
  }

  {  // 4x4
    matrix_type A{
          {-72.64085387384407, -48.30700858824612, -19.3702010361321, 51.516503440727945},
          {-48.30700858824612, -87.11535639272842, 2.3738314209291502, 29.42327462447335},
          {-19.3702010361321, 2.3738314209291502, 55.21571891343558, 7.83398454828334},
          {51.516503440727945, 29.42327462447335, 7.83398454828334, -38.61950221926611}};

    matrix_type B{
          {96.83468091117089, -4.913791379700475, -18.339421312281036,
           -20.346322316982484},
          {-4.913791379700475, 88.59720090255685, 28.357509831861982, -33.19377288633707},
          {-18.339421312281036, 28.357509831861982, 34.13774454920027, -3.14768974712774},
          {-20.346322316982484, -33.19377288633707, -3.14768974712774,
           59.28831461068096}};

    std::vector<evalue_type> evalues_all{-1.4897627147607764, -0.9161984422453742,
                                         -0.029585986915773334, 2.0250925458971345};

    std::vector<evector_type> evectors_all{
          evector_type{0.06226711022185643, 0.09668084160249292, -0.013652970551548528,
                       0.01325439248149767},
          evector_type{-0.04233543038268334, 0.09457751570693788, -0.05511965993565457,
                       0.11198168964857742},
          evector_type{-0.08111587346126517, 0.009768215698732107, -0.01448140623123017,
                       -0.10753324172783782},
          evector_type{0.021283358651158414, -0.044129449484843375, 0.20582238784498755,
                       0.003988733734177379}};

    {
      //
      // Dense generalised 4x4 problem (all eigenpairs)
      //

      testproblem_type prob{"4x4 dense problem (all eigenpairs)", A, B, evalues_all,
                            evectors_all};
      res.push_back(std::move(prob));
    }

    //
    // Dense generalised 4x4 problem (2 largest real)
    //
    {
      std::vector<evalue_type> evalues{evalues_all[2], evalues_all[3]};
      std::vector<evector_type> evectors{evectors_all[2], evectors_all[3]};

      testproblem_type prob{"4x4 dense problem (2 LargestReal)", A, B, evalues, evectors};
      prob.params.update({{EigensolverBaseKeys::which, "LR"}});
      res.push_back(std::move(prob));
    }
  }  // 4x4 problem

  {  // 10x10 problems
    matrix_type A{
          {-53.935136233822675, -28.677589284881904, -7.56982046435148,
           -2.0722645098813643, 27.156455301373313, 0.8524176648965636,
           -28.520432253006334, 53.95615623932778, -7.317845285379093,
           -33.28607335580642},
          {-28.677589284881904, 89.2618557832887, 39.077977669512734, 30.032224630025866,
           48.15781713041957, -7.002016512685032, -33.188129930458075,
           -63.236049287823505, 66.57963444157375, -19.41260115004485},
          {-7.56982046435148, 39.077977669512734, -35.07163005380278, 34.291098587603926,
           -5.508454237261006, 50.975126813586314, -32.06122534107001, 18.69964061697776,
           8.843268928583427, -24.19973486253687},
          {-2.0722645098813643, 30.032224630025866, 34.291098587603926,
           -68.52378937802081, 56.51456700055036, 68.94141993752322, 18.853338230966727,
           -41.3845590368809, 1.5204993126399131, -73.0816884019255},
          {27.156455301373313, 48.15781713041957, -5.508454237261006, 56.51456700055036,
           28.3175145162096, -18.21880172780247, -30.90022374445236, -42.34432706810645,
           17.323538958056474, 5.0086116408434975},
          {0.8524176648965636, -7.002016512685032, 50.975126813586314, 68.94141993752322,
           -18.21880172780247, -95.66986405032512, 35.699988990153685, 19.17046465038132,
           -42.98093897084816, -33.367594431732485},
          {-28.520432253006334, -33.188129930458075, -32.06122534107001,
           18.853338230966727, -30.90022374445236, 35.699988990153685, 52.50278461230829,
           -57.34733617944396, -6.330261551024364, -69.71421300355263},
          {53.95615623932778, -63.236049287823505, 18.69964061697776, -41.3845590368809,
           -42.34432706810645, 19.17046465038132, -57.34733617944396, 29.36169015959274,
           -57.816883275177375, -29.084474879859044},
          {-7.317845285379093, 66.57963444157375, 8.843268928583427, 1.5204993126399131,
           17.323538958056474, -42.98093897084816, -6.330261551024364,
           -57.816883275177375, 69.86107748405686, -6.375991509337922},
          {-33.28607335580642, -19.41260115004485, -24.19973486253687, -73.0816884019255,
           5.0086116408434975, -33.367594431732485, -69.71421300355263,
           -29.084474879859044, -6.375991509337922, 50.66098665912392}};

    matrix_type B{
          {69.73012444971515, -11.430661933986705, 15.730350842775003, -24.52988856533893,
           -14.172605343131824, -7.0501253589043404, 7.27170645584755, -9.258484935937416,
           17.7668745034304, 9.41475873631007},
          {-11.430661933986705, 126.22683677796444, -13.639569018127501,
           28.91169798703286, -18.992190953527633, -28.25788991482182, -6.762052074665895,
           5.434877254425396, -47.88103619784025, 27.24610286717661},
          {15.730350842775003, -13.639569018127501, 112.87045255337152,
           12.826767727878451, -3.3473568601920713, 19.201723091765622,
           -4.406826909805092, 20.92719172336139, 30.606272253597652, 23.34013760044534},
          {-24.52988856533893, 28.911697987032852, 12.82676772787846, 156.9893092045161,
           -40.45460649944304, -19.25808118067705, -4.2945415882091975,
           44.862097365098904, 35.84853985200436, -16.874840344838177},
          {-14.172605343131824, -18.992190953527633, -3.347356860192068,
           -40.45460649944304, 94.03814128117055, 18.0470675071459, -26.24594967630993,
           -20.158240416988537, -1.4758063788071172, 28.77500164248446},
          {-7.0501253589043404, -28.25788991482182, 19.201723091765622,
           -19.258081180677046, 18.0470675071459, 130.36177605175888, 3.0608453179691413,
           -13.363776832528758, 40.28636700785224, -5.54509522032295},
          {7.2717064558475855, -6.7620520746659025, -4.406826909805086,
           -4.294541588209197, -26.24594967630993, 3.0608453179691413, 135.21757248783427,
           9.729839129580744, -19.26774548815396, 16.880710356279746},
          {-9.258484935937416, 5.434877254425396, 20.927191723361403, 44.86209736509889,
           -20.158240416988537, -13.363776832528755, 9.729839129580741,
           140.93585709877664, -5.332475110946257, -1.074412898097807},
          {17.766874503430504, -47.88103619784027, 30.606272253597663, 35.84853985200435,
           -1.4758063788071079, 40.286367007852235, -19.267745488153952,
           -5.332475110946267, 183.05436092158007, 15.347950944544298},
          {9.41475873631007, 27.246102867176596, 23.34013760044534, -16.874840344838177,
           28.77500164248446, -5.545095220322949, 16.880710356279746, -1.0744128980978052,
           15.347950944544298, 129.1949048477422}};

    // TODO properly symmetrise the matrix such that this is not needed any
    // more.
    B.symmetrise();

    std::vector<evalue_type> evalues_all{-1.7136036278526083,  -1.355211095728974,
                                         -0.8527445634266658,  -0.698706672984702,
                                         -0.16786604430660576, 0.2677240969637297,
                                         0.46961764077586304,  1.050204884045193,
                                         1.5178488974871367,   3.563158995708878};

    std::vector<evector_type> evectors_all{
          evector_type{-0.0787821494940039, -0.008269089783958577, 0.05916106731163469,
                       -0.04250664671974947, 0.027764377889736177, -0.05967689999629079,
                       0.018288710723297158, 0.011502320651175561, 0.010634561400675589,
                       -0.035186784107419766},
          evector_type{0.04628835863486506, -0.022812705148012785, 0.03144094229897967,
                       0.047644299749638734, -0.011316652441932586, -0.052069711745034344,
                       0.0020123289148461596, -0.03191501503743128, -0.019373672338860306,
                       0.015308083953496004},
          evector_type{0.047641816814945954, 0.013341632173566946, 0.03631061258988824,
                       -0.049384344909928275, -0.006653041703626039, 0.009519277252012221,
                       0.000871591361769514, -0.02978539517616976, -0.00847141629956872,
                       -0.020149474705948658},
          evector_type{-0.03474053572451757, -0.02750381674370305, -0.005709892260547733,
                       0.01326148068635601, -0.03490396970501798, 0.0015708666271720823,
                       -0.05399525128097285, -0.0452603280856187, -0.011379572786059316,
                       -0.023776889462778093},
          evector_type{0.012482448220842204, 0.038541519250349175, -0.042191174130216294,
                       -0.023260640196367983, -0.005224916129290883,
                       -0.030270158049210367, 0.005686573319119298, 0.013581275191611461,
                       -0.016835664376319462, -0.02663297767715056},
          evector_type{0.041535634931841585, -0.02413642886899026, -0.03893864113625138,
                       -0.014759218323871028, 0.021325031935225624, -0.03213609458943114,
                       0.009056739862796478, 0.020772751088208603, 0.05097203325318504,
                       -0.032173702104134436},
          evector_type{0.01347366865101061, -0.033355400819688956, -0.03201231867917881,
                       0.02824837865343379, 0.07048413978414196, -0.006038923914580567,
                       0.005366263261109896, -0.03133738408808794, -0.028989356051589962,
                       0.03099889686974452},
          evector_type{0.042061981027501166, 0.02396136921700857, 0.009085400417279834,
                       -0.0018875204010733834, 0.038313539286835126,
                       -0.005481967868511801, -0.06288687592457373, 0.03668595768990976,
                       -0.001106582679160941, 0.004452228153220459},
          evector_type{-0.055662642711443534, 0.026205441433511927, -0.009323313613061534,
                       -0.0423549853365753, -0.04098342788382772, -0.03379740469278036,
                       -0.012681039829628035, -0.024715277107019505, 0.050437266940886634,
                       0.044287991495435615},
          evector_type{-0.013970575918636721, -0.07846617601398624, -0.023615010088959235,
                       -0.019587155987521972, -0.07809374384529433, 0.006709492767634691,
                       -0.03977942570794509, 0.025812247968306334, -0.04010839910887779,
                       0.06676247421329834}};

    //
    // 10x10 dense problem (all eigenpairs)
    //
    {
      testproblem_type prob{"10x10 dense problem (all epairs)", A, B, evalues_all,
                            evectors_all};
      res.push_back(std::move(prob));
    }  // all

    //
    // Dense generalised 10x10 problem (2 smallest real)
    //
    {
      std::vector<evalue_type> evalues{evalues_all[0], evalues_all[1]};

      std::vector<evector_type> evectors{evectors_all[0], evectors_all[1]};

      testproblem_type prob{"10x10 dense problem (2 SmallestReal)", A, B, evalues,
                            evectors};
      prob.params.update({{EigensolverBaseKeys::which, "SR"}});
      res.push_back(std::move(prob));
    }  // 2 SR

    //
    // Dense generalised 10x10 problem (4 largest magnitude)
    //
    {
      std::vector<evalue_type> evalues{evalues_all[0], evalues_all[1], evalues_all[8],
                                       evalues_all[9]};

      std::vector<evector_type> evectors{evectors_all[0], evectors_all[1],
                                         evectors_all[8], evectors_all[9]};

      testproblem_type prob{"10x10 dense problem (4 LargestMagnitude)", A, B, evalues,
                            evectors};
      prob.params.update({{EigensolverBaseKeys::which, "LM"}});
      res.push_back(std::move(prob));
    }  // 4 LM

    //
    // Dense generalised 10x10 problem (3 largest real)
    //
    {
      std::vector<evalue_type> evalues{evalues_all[7], evalues_all[8], evalues_all[9]};

      std::vector<evector_type> evectors{evectors_all[7], evectors_all[8],
                                         evectors_all[9]};

      testproblem_type prob{"10x10 dense problem (3 LargestReal)", A, B, evalues,
                            evectors};
      prob.params.update({{EigensolverBaseKeys::which, "LR"}});
      res.push_back(std::move(prob));

    }  // 3 LR
  }    // 10x10
  return res;
}

template <typename Matrix>
std::vector<EigensolverTestProblem<Matrix, false, false>>
EigensolverTestProblemLibrary<EigensolverTestProblem<Matrix, /*Hermitian=*/false,
                                                     /*generalised=*/false>,
                              /*complex=*/false>::get_all() {
  /*typedef typename testproblem_type::matrix_type matrix_type;
  typedef typename testproblem_type::evalue_type evalue_type;
  typedef typename testproblem_type::evector_type evector_type;
  */

  std::vector<testproblem_type> res;

  // TODO

  return res;
}

template <typename Matrix>
std::vector<EigensolverTestProblem<Matrix, false, true>> EigensolverTestProblemLibrary<
      EigensolverTestProblem<Matrix, /*Hermitian=*/false, /*generalised=*/true>,
      /*complex=*/false>::get_all() {
  /*typedef typename testproblem_type::matrix_type matrix_type;
  typedef typename testproblem_type::evalue_type evalue_type;
  typedef typename testproblem_type::evector_type evector_type;
  */

  std::vector<testproblem_type> res;

  // TODO

  return res;
}
}  // namespace eigensolver_tests
}  // namespace tests
}  // namespace linalgwrap
