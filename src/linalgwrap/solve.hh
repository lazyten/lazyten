#pragma once
#include "Armadillo.hh"
#include "MultiVector.hh"
#include <krims/ParameterMap.hh>

namespace linalgwrap {

// TODO This is just a quick and dirty implementation to get it going.

#ifndef LINALGWRAP_HAVE_ARMADILLO
static_assert(false, "We need armadillo for this at the moment.");
#endif

namespace detail {
template <bool cond>
using ParameterMap_if =
      typename std::enable_if<cond, krims::ParameterMap>::type;
}

// TODO transpose solve?

//@{
/** Solve a linear system A x = b, where A is hermitian
 *
 * \param A     System matrix (assumed to be hermitian)
 * \param x     Lhs vector or vectors (solution)
 * \param b     Rhs vector or vectors (problem)
 * \param map   Specify some solver parameters
 *
 * \throws      Subclass of SolverException in case there is an error.
 **/
template <typename Matrix, typename Vector>
void solve_hermitian(
      const Matrix& A, MultiVector<Vector>& x, const MultiVector<Vector>& b,
      const detail::ParameterMap_if<IsMatrix<Matrix>::value &&
                                    IsMutableVector<Vector>::value>& map =
            krims::ParameterMap()) {
    assert_size(x.n_vectors(), b.n_vectors());
    for (size_t i = 0; i < x.n_vectors(); ++i) {
        solve_hermitian(A, x[i], b[i], map);
    }
}

template <typename Matrix, typename Vector>
void solve_hermitian(
      const Matrix& A, Vector& x, const Vector& b,
      const detail::ParameterMap_if<IsMatrix<Matrix>::value &&
                                    IsMutableVector<Vector>::value>& map =
            krims::ParameterMap()) {

    static_assert(std::is_same<typename Matrix::scalar_type,
                               typename Vector::scalar_type>::value,
                  "The scalar types need to agree.");
    static_assert(std::is_same<typename Matrix::scalar_type, double>::value,
                  "This simple implementation only works for double");

    assert_dbg(A.is_symmetric(), ExcMatrixNotSymmetric());
    assert_size(A.n_rows(), b.size());
    assert_size(A.n_cols(), x.size());

    // Assert some types:
    static_assert(std::is_same<typename StoredTypeOf<Matrix>::type,
                               ArmadilloMatrix<double>>::value,
                  "This hack only works for armadillo matrices");

    static_assert(IsStoredVector<Vector>::value,
                  "Currently we assume that the Vector is a stored vector");

    // Make arma matrices
    const arma::Mat<double>& m_arma =
          as_stored(A).data();  //.t() is skipped since m is symmetric
    arma::Col<double> b_arma(b.memptr(), b.size());
    arma::Col<double> x_arma(x.memptr(), x.size(), false);

    // Solve the linear system:
    bool result = arma::solve(x_arma, m_arma, b_arma);
    assert_throw(result, SolverException());

    // Fake-use the map.
    (void)map;
}
//@}

//
// ---------------------------------------------
//

//@{
/** Solve a linear system A x = b
 *
 * \param A     System matrix
 * \param x     Lhs vector or vectors (solution)
 * \param b     Rhs vector or vectors (problem)
 * \param map   Specify some solver parameters
 *
 * \throws      Subclass of SolverException in case there is an error.
 **/
template <typename Matrix, typename Vector>
void solve(const Matrix& A, MultiVector<Vector>& x,
           const MultiVector<Vector>& b,
           const detail::ParameterMap_if<IsMatrix<Matrix>::value &&
                                         IsMutableVector<Vector>::value>& map =
                 krims::ParameterMap()) {
    assert_size(x.n_vectors(), b.n_vectors());
    for (size_t i = 0; i < x.n_vectors(); ++i) {
        solve(A, x[i], b[i], map);
    }
}

template <typename Matrix, typename Vector>
void solve(const Matrix& A, Vector& x, const Vector& b,
           const detail::ParameterMap_if<IsMatrix<Matrix>::value &&
                                         IsMutableVector<Vector>::value>& map =
                 krims::ParameterMap()) {

    assert_dbg(false, krims::ExcNotImplemented());
}
//@}

}  // namespace linalgwrap
