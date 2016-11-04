#pragma once
#include "Armadillo/ArmadilloEigensolver.hh"
#include "Arpack/ArpackEigensolver.hh"
#include "Base/Solvers.hh"
#include <krims/ParameterMap.hh>

namespace linalgwrap {

#ifndef LINALGWRAP_HAVE_ARMADILLO
static_assert(false,
              "We need armadillo in order to have at least a working fallback "
              "eigensolver.");
#endif

/** Solve a normal hermitian eigensystem
 *
 * \param A     Matrix to diagonalise
 * \param n_ep  Number of eigenpairs to obtain
 *              (Constants<typename Matrix::size_type>::all for all)
 * \param map   Specify some solver parameters
 *
 * \throws      Subclass of SolverException in case there is an error.
 **/
template <typename Matrix>
EigensolutionTypeFor<true, Matrix> eigensystem_hermitian(
      const Matrix& A,
      typename std::enable_if<IsMatrix<Matrix>::value,
                              typename Matrix::size_type>::type n_ep =
            Constants<typename Matrix::size_type>::all,
      const krims::ParameterMap& map = krims::ParameterMap());

/** Solve a generalised hermitian eigensystem
 *
 * \param A     Matrix to diagonalise
 * \param B     Metric to use
 * \param n_ep  Number of eigenpairs to obtain
 *              (Constants<typename Matrix::size_type>::all for all)
 * \param map   Specify some solver parameters
 *
 * \throws      Subclass of SolverException in case there is an error.
 **/
template <typename MatrixA, typename MatrixB>
EigensolutionTypeFor<true, MatrixA> eigensystem_hermitian(
      const MatrixA& A, const MatrixB& B,
      typename std::enable_if<IsMatrix<MatrixA>::value &&
                                    IsMatrix<MatrixB>::value,
                              typename MatrixA::size_type>::type n_ep =
            Constants<typename MatrixA::size_type>::all,
      const krims::ParameterMap& map = krims::ParameterMap());

/** Solve a normal eigensystem
 *
 * \param A     Matrix to diagonalise
 * \param n_ep  Number of eigenpairs to obtain
 *              (Constants<typename Matrix::size_type>::all for all)
 * \param map   Specify some solver parameters
 *
 * \throws      Subclass of SolverException in case there is an error.
 **/
template <typename Matrix>
EigensolutionTypeFor<false, Matrix> eigensystem(
      const Matrix& A,
      typename std::enable_if<IsMatrix<Matrix>::value,
                              typename Matrix::size_type>::type n_ep =
            Constants<typename Matrix::size_type>::all,
      const krims::ParameterMap& map = krims::ParameterMap());

/** Solve a generalised eigensystem
 *
 * \param A     Matrix to diagonalise
 * \param B     Metric to use
 * \param n_ep  Number of eigenpairs to obtain
 *              (Constants<typename Matrix::size_type>::all for all)
 * \param map   Specify some solver parameters
 *
 * \throws      Subclass of SolverException in case there is an error.
 **/
template <typename MatrixA, typename MatrixB>
EigensolutionTypeFor<false, MatrixA> eigensystem(
      const MatrixA& A, const MatrixB& B,
      typename std::enable_if<IsMatrix<MatrixA>::value &&
                                    IsMatrix<MatrixB>::value,
                              typename MatrixA::size_type>::type n_ep =
            Constants<typename MatrixA::size_type>::all,
      const krims::ParameterMap& map = krims::ParameterMap());

//
// ------------------------------------------------------
//

template <typename Matrix>
EigensolutionTypeFor<true, Matrix> eigensystem_hermitian(
      const Matrix& A,
      typename std::enable_if<IsMatrix<Matrix>::value,
                              typename Matrix::size_type>::type n_ep,
      const krims::ParameterMap& map) {
    typedef Eigenproblem<true, Matrix> problem_type;
    problem_type problem{A, n_ep};

#ifdef LINALGWRAP_HAVE_ARPACK
    if (problem.n_ep() <= problem.dim() / 2 &&
        map.at<std::string>(EigensolverBaseKeys::which, "SR") !=
              std::string("SM")) {
        // Use Arpack since we want "few" eigenpairs
        // TODO This is just a shot, no clue whether half the dimension
        // is a sensible value or not.
        return ArpackEigensolver<problem_type>{map}
              .solve(problem)
              .eigensolution();
    }
#endif

    // Fallback: Armadillo
    ArmadilloEigensolver<problem_type> solver{map};
    return solver.solve(problem).eigensolution();
}

template <typename MatrixA, typename MatrixB>
EigensolutionTypeFor<true, MatrixA> eigensystem_hermitian(
      const MatrixA& A, const MatrixB& B,
      typename std::enable_if<IsMatrix<MatrixA>::value &&
                                    IsMatrix<MatrixB>::value,
                              typename MatrixA::size_type>::type n_ep,
      const krims::ParameterMap& map) {
    typedef Eigenproblem<true, MatrixA, MatrixB> problem_type;
    problem_type problem{A, B, n_ep};

    // TODO Do something better here some day (e.g. use Arpack)
    ArmadilloEigensolver<problem_type> solver{map};
    return solver.solve(problem).eigensolution();
}

template <typename Matrix>
EigensolutionTypeFor<false, Matrix> eigensystem(
      const Matrix& A,
      typename std::enable_if<IsMatrix<Matrix>::value,
                              typename Matrix::size_type>::type n_ep,
      const krims::ParameterMap& map) {
    // TODO This code is untested!
    assert_sufficiently_tested(false);

    typedef Eigenproblem<false, Matrix> problem_type;
    problem_type problem{A, n_ep};

    // TODO Do something better here some day (e.g. use Arpack)
    ArmadilloEigensolver<problem_type> solver{map};
    return solver.solve(problem).eigensolution();
}

template <typename MatrixA, typename MatrixB>
EigensolutionTypeFor<false, MatrixA> eigensystem(
      const MatrixA& A, const MatrixB& B,
      typename std::enable_if<IsMatrix<MatrixA>::value &&
                                    IsMatrix<MatrixB>::value,
                              typename MatrixA::size_type>::type n_ep,
      const krims::ParameterMap& map) {
    // TODO This code is untested!
    assert_sufficiently_tested(false);

    typedef Eigenproblem<false, MatrixA, MatrixB> problem_type;
    problem_type problem{A, B, n_ep};

    // TODO Do something better here some day (e.g. use Arpack)
    ArmadilloEigensolver<problem_type> solver{map};
    return solver.solve(problem).eigensolution();
}

}  // namespace linalgwrap
