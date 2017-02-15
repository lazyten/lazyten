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

#pragma once
#include "linalgwrap/Exceptions.hh"
#include "linalgwrap/Matrix_i.hh"
#include "linalgwrap/TypeUtils.hh"
#include <krims/SubscriptionPointer.hh>

namespace linalgwrap {

// TODO Would it be a big problem if the interface here was more flexible,
//      i.e. that more parameters could be altered post-construction time?

/** \brief Base class for eigenproblems.
 *  The common stuff both general and non-general eigenproblems
 *  share.
 *
 * \tparam isHermitian   Is the problem an Hermitian problem
 * \tparam MatrixA       The type of the operator A in the eigenproblem
 *                       $Ax = \lambda Bx$ or $Ax = \lambda x$.
 * \tparam MatrixDiag    The type of the matrix to diagonalise
 *                       (may be different for example if
 *                       shift-and-invert is used)
 */
template <bool isHermitian, typename MatrixA, typename MatrixDiag>
class EigenproblemBase {
 public:
  /** \name Type definitions */
  ///@{
  /** The stored matrix type used within the eigenproblem */
  typedef typename StoredTypeOf<MatrixA>::type stored_matrix_type;

  /** The size type used in the matrices*/
  typedef typename stored_matrix_type::size_type size_type;

  /** The scalar type used in the matrices*/
  typedef typename stored_matrix_type::scalar_type scalar_type;

  /** The real type used in the matrices */
  typedef typename stored_matrix_type::real_type real_type;

  /** The type of the matrix A in $A x = \lambda B x$
   *  or $A x = \lambda x$ */
  typedef MatrixA matrix_a_type;

  /** The type of the matrix to be diagonalised */
  typedef MatrixDiag matrix_diag_type;
  ///@}

  // Assert that all matrices have the same underlying stored
  // matrix type:
  static_assert(
        std::is_same<stored_matrix_type, typename StoredTypeOf<MatrixDiag>::type>::value,
        "The Operator type and the MatrixA type need to have the "
        "same underlying stored matrix type.");

  // Assert that all matrix types are matrices
  static_assert(IsMatrix<MatrixA>::value, "MatrixA needs to be a matrix type");
  static_assert(IsMatrix<MatrixDiag>::value, "MatrixDiag needs to be a matrix type");

  /** May the eigensolver assume the operator to diagonalise is
   *  Hermitian? */
  static constexpr bool hermitian = isHermitian;

  /** May the eigensolver assume that all scalar types are real? */
  static constexpr bool real = !krims::IsComplexNumber<scalar_type>::value;

  /** The operator A of the eigenproblem $Ax = \lambda Bx$. */
  const MatrixA& A() const { return *m_A_ptr; }

  /** The matrix to diagonalise.
    *
    * This object may be different from A_ptr if some spectral
    * transformation (like a shift-and-invert) is employed.
    * E.g. we might want to apply the operation $(A - \sigma I)^{-1}$
    * if we look for eigenvalues of $A$ around $\sigma$.
    */
  const MatrixDiag& Diag() const { return *m_diag_ptr; }

  /** The number of eigenpairs to compute */
  size_type n_ep() const { return m_n_ep; }

  /** Return the dimensionality of the eigenproblem */
  size_type dim() const { return A().n_cols(); }

  /** Constructor
   *
   * If n_ep > dim()  the constructor automatically sets n_ep = dim()
   * */
  EigenproblemBase(const MatrixA& A, size_type n_ep, const MatrixDiag& diag)
        : m_A_ptr("EigenproblemA", A),
          m_diag_ptr("EigenproblemDiag", diag),
          m_n_ep(std::min(n_ep, dim())) {

    // Diag and A need to have the same number of columns
    assert_size(A.n_cols(), diag.n_cols());

#if DEBUG
    if (hermitian) {
      // Allow a tiny error when checking for hermiticity
      const real_type tolerance = 100 * Constants<real_type>::default_tolerance;
      assert_dbg(m_A_ptr->is_hermitian(tolerance), ExcMatrixNotHermitian());
      assert_dbg(m_diag_ptr->is_hermitian(tolerance), ExcMatrixNotHermitian());
    }
#endif
  }

 private:
  /** The matrix within $A x = \lambda x$ or $A x = \lambda B x$ */
  krims::SubscriptionPointer<const MatrixA> m_A_ptr;

  /** The matrix to diagonalise */
  krims::SubscriptionPointer<const MatrixDiag> m_diag_ptr;

  /** The number of eigenpairs to compute*/
  size_type m_n_ep;
};

//@{
/* \brief struct describing a generalised Eigenproblem
 *
 * \tparam isHermitian   Is the problem an Hermitian problem
 * \tparam MatrixA       The type of the operator A in the eigenproblem
 *                       $Ax = \lambda Bx$.
 * \tparam MatrixB       The type of the operator B in the eigenproblem
 *                       $Ax = \lambda Bx$ or void if it is a regular
 *                       Eigenproblem $Ax = \lambda x$.
 * \tparam MatrixDiag      The type of the operator to diagonalise
 */
template <bool isHermitian, typename MatrixA, typename MatrixB = void,
          typename MatrixDiag = MatrixA>
class Eigenproblem : public EigenproblemBase<isHermitian, MatrixA, MatrixDiag> {
 public:
  /** \name Type definitions */
  ///@{
  typedef EigenproblemBase<isHermitian, MatrixA, MatrixDiag> base_type;
  typedef typename base_type::stored_matrix_type stored_matrix_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::real_type real_type;

  /** The type of the matrix B in $A x = \lambda B x$ */
  typedef MatrixB matrix_b_type;
  ///@}

  static_assert(
        std::is_same<stored_matrix_type, typename StoredTypeOf<MatrixB>::type>::value,
        "For a general eigenproblem the MatrixB type and the MatrixA type "
        "need to have the same underlying stored matrix type.");
  static_assert(IsMatrix<MatrixB>::value, "MatrixB needs to be a matrix type");

  /** Is this a generalised eigenproblem.
   *
   * If this is false, than the B() function does not exist */
  static constexpr bool generalised = true;

  /** The operator B of the eigenproblem $Ax = \lambda Bx$. */
  const MatrixB& B() const { return *m_B_ptr; }

  /** Constructors */
  ///@{
  /** \brief Construct an eigenproblem $A x = \lambda B x$.
   *
   * This version allows to specify an actual matrix to be diagonalised.
   * This is may be used for spectral transformations.
   *
   *  \param A     The matrix A
   *  \param B     The matrix B
   *  \param n_ep  The number of eigenvalues to compute
   *               A value equal to IterationConstants<size_type>::all
   *               implies that all eigenpairs should be computed
   *  \param Diag    The operator to diagonalise.
   */
  Eigenproblem(const MatrixA& A, const MatrixB& B, size_type n_ep, const MatrixDiag& Diag)
        : base_type(A, n_ep, Diag), m_B_ptr("Eigenproblem", B) {
#if DEBUG
    // A and B need to have the same dimensionality.
    assert_size(A.n_cols(), B.n_cols());
    // B should be hermitian (with a slight error)
    const real_type tolerance = 100 * Constants<real_type>::default_tolerance;
    assert_dbg(B.is_hermitian(tolerance), ExcMatrixNotHermitian());
#endif
    //
    // TODO B should be positive definite as well (but we need to do an
    //      eigendecomposition to check for that ...)
  }

  /** \brief Construct an eigenproblem $A x = \lambda B x$
   *
   * \note This constructor is only enabled if MatrixA and MatrixDiag are the
   * same type.
   *
   *  \param A     The matrix A
   *  \param B     The matrix B
   *  \param n_ep  The number of eigenvalues to compute
   *               A value equal to IterationConstants<size_type>::all
   *               implies that all eigenpairs should be computed
   **/
  template <typename Matrix,
            typename = krims::enable_if_cond_same_t<
                  std::is_same<MatrixA, MatrixDiag>::value, Matrix, MatrixA>>
  Eigenproblem(const Matrix& A, const MatrixB& B,
               size_type n_ep = Constants<size_type>::all)
        : Eigenproblem(A, B, n_ep, A) {}
  ///@}

 private:
  /** The operator B of the eigenproblem $Ax = \lambda Bx$. */
  krims::SubscriptionPointer<const MatrixB> m_B_ptr;
};

template <bool isHermitian, typename MatrixA, typename MatrixDiag>
struct Eigenproblem<isHermitian, MatrixA, void, MatrixDiag>
      : public EigenproblemBase<isHermitian, MatrixA, MatrixDiag> {
  /** \name Type definitions */
  ///@{
  typedef EigenproblemBase<isHermitian, MatrixA, MatrixDiag> base_type;
  typedef typename base_type::size_type size_type;

  /** The type of the matrix B in $A x = \lambda B x$ */
  typedef void matrix_b_type;
  ///@}

  /** Is this a generalised eigenproblem.
   *
   * If this is false, than the B() member function does not exist */
  static constexpr bool generalised = false;

  /** Never call this routine. */
  const MatrixA& B() const {
    assert_dbg(false, krims::ExcDisabled("The method B() for a non-generalised "
                                         "eigenproblem should not be called."));
    return base_type::A();
  }

  /** \brief Construct an eigenproblem $A x = \lambda x$
  *
  * This version allows to specify an actual matrix to be diagonalised.
  * This is may be used for spectral transformations.
  *
  *  \param A     The matrix A
  *  \param n_ep  The number of eigenvalues to compute
  *               A value equal to IterationConstants<size_type>::all
  *               implies that all eigenpairs should be computed
  *  \param Diag    The operator to diagonalise.
  **/
  Eigenproblem(const MatrixA& A, size_type n_ep, const MatrixDiag& Diag)
        : base_type(A, n_ep, Diag) {}

  /** \brief Construct an eigenproblem $A x = \lambda x$
   *
   * \note This constructor is only enabled if MatrixA and MatrixDiag are the
   * same type.
   *
   *  \param A     The matrix A
   *  \param n_ep  The number of eigenvalues to compute
   *               A value equal to IterationConstants<size_type>::all
   *               implies that all eigenpairs should be computed
   **/
  template <typename Matrix,
            typename = krims::enable_if_cond_same_t<
                  std::is_same<MatrixA, MatrixDiag>::value, Matrix, MatrixA>>
  Eigenproblem(const Matrix& A, size_type n_ep = Constants<size_type>::all)
        : Eigenproblem(A, n_ep, A) {}
};
//@}

}  // namespace linalgwrap
