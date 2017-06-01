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
#include "linalgwrap/config.hh"
#ifdef LINALGWRAP_HAVE_ARMADILLO

#include <armadillo>
#include <krims/ExceptionSystem.hh>
#include <krims/TypeUtils.hh>
#include <linalgwrap/Constants.hh>

// TODO: Performance analysis: Maybe it is better to first use transpose()
// to get the transpose and then use as_stored to get the stored matrix
// representation

namespace linalgwrap {

namespace detail {

/** In armadillo t() is the normal transpose for real matrices and the conjugate
 * transpose for complex ones.
 * Here we make sure that we actually get exactly what we want. */
template <typename Scalar, bool real = !krims::IsComplexNumber<Scalar>::value>
struct ArmaTranspose {
  static auto trans(const arma::Mat<Scalar>& m) -> decltype(m.t()) { return m.t(); }
  static auto conjtrans(const arma::Mat<Scalar>& m) -> decltype(m.t()) { return m.t(); }
};

/** Specialisation of above for complex numbers */
template <typename Scalar>
struct ArmaTranspose<Scalar, false> {
  static auto trans(const arma::Mat<Scalar>& m) -> decltype(m.st()) { return m.st(); }
  static auto conjtrans(const arma::Mat<Scalar>& m) -> decltype(m.t()) { return m.t(); }
};

/** Wrapper around the armadillo eig_sym, eig_gen and eig_pair
 *  functions calling whatever is appropriate
 **/
template <typename Eigenproblem, bool isHermitian = Eigenproblem::hermitian,
          bool isGeneralised = Eigenproblem::generalised,
          bool isReal = Eigenproblem::real>
struct ArmadilloEigWrapper {};

/** Hermitian, normal eigenproblems */
template <typename Eigenproblem, bool real>
struct ArmadilloEigWrapper<Eigenproblem, /* hermitian= */ true,
                           /* generalised= */ false,
                           /* real= */ real> {

  // Note: The arma::Mat object stored inside ArmadilloMatrix is actually
  // the *transpose* of the matrix the ArmadilloMatrix object represents.
  //
  // This is due to the fact that we work with row-major matrices,
  // but armadillo works with column-major matrices.
  //
  // For Hermitian problems this means that if we just use data(),
  // we actually solve for the eigenpairs of the complex conjugate of
  // the matrix we want to get the eigenpairs for. This means
  // that we need to take the complex conjugate of the eigenvalues
  // and eigenvectors we obtain before returning them.
  // But since Hermetian matrices have real eigenvalues, we can even
  // spare taking the complex conjugate on them.
  // (This is since Av = b v  =>  A^\dagger v = b v
  //                          =>  A^T conj(v) = conj(b) conj(v)
  //                          =>  A^T conj(v) = b conj(v)
  // and A^T is what we put in)

  typedef typename Eigenproblem::scalar_type scalar_type;
  typedef typename Eigenproblem::stored_matrix_type stored_matrix_type;

  static bool eig(const Eigenproblem& prb, arma::Col<scalar_type>& eval_arma,
                  arma::Mat<scalar_type>& evec_arma) {
    const bool res = arma::eig_sym(eval_arma, evec_arma, as_stored(prb.A()).data());

    static_assert(std::is_same<scalar_type, double>::value,
                  "Here we assume that the scalar type is real only");
    // For complex scalar types we needed:
    // evec_arma = arma::conj(evec_arma);

    return res;
  }
};

/** non-Hermitian normal eigenproblems */
template <typename Eigenproblem, bool real>
struct ArmadilloEigWrapper<Eigenproblem, /* hermitian= */ false,
                           /* generalised= */ false,
                           /* real= */ real> {
  // If not real: leave type as is, else make it complex.
  typedef
        typename std::conditional<real, std::complex<typename Eigenproblem::scalar_type>,
                                  typename Eigenproblem::scalar_type>::type scalar_type;
  typedef typename Eigenproblem::stored_matrix_type stored_matrix_type;

  static bool eig(const Eigenproblem& prb, arma::Col<scalar_type>& eval_arma,
                  arma::Mat<scalar_type>& evec_arma) {
    // See notes in ArmadilloEigWrapper<Eigenproblem,true,false,real>
    // why we need the .t() in the end
    return arma::eig_gen(eval_arma, evec_arma, as_stored(prb.A()).data().t());
  }
};

DefException3(ExcImaginaryPartNonZero, size_t, size_t, arma::cx_mat&,
              << "The imaginary part of the element (" << arg1 << "," << arg2
              << ") of a generalised isHermitian "
                 "eigenproblem was not zero."
              << "The object is " << std::endl
              << arg3);

/** real symmetric generalised eigenproblems */
template <typename Eigenproblem>
struct ArmadilloEigWrapper<Eigenproblem, /* hermitian= */ true,
                           /* generalised= */ true,
                           /* real= */ true> {
  typedef typename Eigenproblem::scalar_type scalar_type;
  typedef typename Eigenproblem::stored_matrix_type stored_matrix_type;

  static bool eig(const Eigenproblem& prb, arma::Col<scalar_type>& eval_arma,
                  arma::Mat<scalar_type>& evec_arma) {

    // Arma eig_pair insists on complex algebra
    arma::cx_mat inner_evecs;
    arma::cx_vec inner_evals;  // column vector

    // No need to transpose or take complex conjugate, since real symmetric
    // (see ArmadilloEigWrapper<Eigenproblem,true,false,real> for details)
    const arma::Mat<scalar_type> b_arma = as_stored(prb.B()).data();
    const bool res =
          arma::eig_pair(inner_evals, inner_evecs, as_stored(prb.A()).data(), b_arma);
    if (!res) return false;

    // Since eigenvectors are only unique up to a phase,
    // we may have eigenvectors which have an imaginary part, but
    // which can (hopefully) be rotated such that they are
    // entirely real.
    for (size_t coli = 0; coli < inner_evecs.n_cols; ++coli) {
      rotate_imaginary_part(inner_evecs.unsafe_col(coli));
    }

#ifdef DEBUG
    // Check imaginary parts are now really zero:
    for (size_t i = 0; i < inner_evals.n_elem; ++i) {
      assert_dbg(
            std::abs(inner_evals(i).imag()) < Constants<scalar_type>::default_tolerance,
            ExcImaginaryPartNonZero(i, 0, inner_evals));
    }

    for (size_t i = 0; i < inner_evecs.n_rows; ++i) {
      for (size_t j = 0; j < inner_evecs.n_cols; ++j) {
        assert_dbg(std::abs(inner_evecs(i, j).imag()) <
                         Constants<scalar_type>::default_tolerance,
                   ExcImaginaryPartNonZero(i, j, inner_evecs));
      }
    }
#endif
    // Copy results:
    eval_arma = arma::real(inner_evals);
    evec_arma = arma::real(inner_evecs);

    // Armadillo does not B-normalise the eigenvectors properly
    // so we need to do this here:
    for (size_t i = 0; i < evec_arma.n_cols; ++i) {
      double norm = arma::as_scalar(evec_arma.col(i).t() * b_arma * evec_arma.col(i));
      evec_arma.col(i) /= std::sqrt(norm);
    }

    return true;
  }

 private:
  /** Rotate the imaginary part of each element of this vectors
   * such that of this vector such that hopefully the vector
   * is entirely real afterwards
   */
  static void rotate_imaginary_part(arma::cx_vec v) {
    // Tolerance for comparisons against zero we use here
    const scalar_type tolerance = Constants<scalar_type>::default_tolerance;

    for (size_t i = 0; i < v.n_elem; ++i) {
      if (std::abs(v[i].imag()) >= tolerance) {
        // We have a non-zero imaginary part
        // => determine phase factor to rotate away
        const std::complex<double> fac = std::conj(v[i]) / std::sqrt(std::abs(v[i]));

        // Scale vector with this factor:
        v[i] *= fac;

        // We can only do this scaling once
        // so if this didn't do it, we have no
        // chance to rotate the other imaginary
        // parts away as well.
        return;
      }
    }
  }
};

/** complex hermitian generalised Eigenproblems */
template <typename Eigenproblem>
struct ArmadilloEigWrapper<Eigenproblem, /* hermitian= */ true,
                           /* generalised= */ true,
                           /* real= */ false> {
  typedef typename Eigenproblem::scalar_type scalar_type;
  typedef typename Eigenproblem::stored_matrix_type stored_matrix_type;

  static bool eig(const Eigenproblem& prb, arma::Col<scalar_type>& eval_arma,
                  arma::Mat<scalar_type>& evec_arma) {
    // Since this is a hermitian problem we can avoid the .t()
    // and take the complex conjugate of the eigenvectors instead
    // (see ArmadilloEigWrapper<Eigenproblem,true,false,real> for details)

    const bool res = arma::eig_pair(eval_arma, evec_arma, as_scalar(prb.A()).data(),
                                    as_scalar(prb.B()).data());
    static_assert(std::is_same<scalar_type, double>::value, "Untested code");
    evec_arma = arma::conj(evec_arma);
    return res;
  }
};

/** non-hermitian generalised Eigenproblems */
template <typename Eigenproblem, bool real>
struct ArmadilloEigWrapper<Eigenproblem, /* hermitian= */ false,
                           /* generalised= */ true,
                           /* real= */ real> {
  // If not real: leave type as is, else make it complex.
  typedef
        typename std::conditional<real, std::complex<typename Eigenproblem::scalar_type>,
                                  typename Eigenproblem::scalar_type>::type scalar_type;
  typedef typename Eigenproblem::stored_matrix_type stored_matrix_type;

  static bool eig(const Eigenproblem& prb, arma::Col<scalar_type>& eval_arma,
                  arma::Mat<scalar_type>& evec_arma) {
    return arma::eig_pair(eval_arma, evec_arma, as_stored(prb.A()).data().t(),
                          as_stored(prb.B()).data().t());
  }
};

}  // namespace detail

}  // namespace linalgwrap
#endif
