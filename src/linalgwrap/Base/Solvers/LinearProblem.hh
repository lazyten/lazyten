//
// Copyright (C) 2017 by the linalgwrap authors
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
#include "linalgwrap/Base/Interfaces/OperatorProperties.hh"
#include "linalgwrap/Matrix_i.hh"
#include "linalgwrap/MultiVector.hh"
#include "linalgwrap/TypeUtils.hh"

namespace linalgwrap {

/** \brief Class for linear problems.
 *
 * \tparam Matrix        The type of the operator A in the problem $Ax = b$.
 * \tparam Vector        The type of the vectors $x$,$b$ in the problem $Ax = b$.
 */
template <typename Matrix, typename Vector>
class LinearProblem {
 public:
  /** \name Type definitions */
  ///@{
  /** The matrix type used within the linear problem */
  typedef typename StoredTypeOf<Matrix>::type stored_matrix_type;

  /** The size type used in the matrices*/
  typedef typename stored_matrix_type::size_type size_type;

  /** The scalar type used in the matrices*/
  typedef typename stored_matrix_type::scalar_type scalar_type;

  /** The real type used in the matrices */
  typedef typename stored_matrix_type::real_type real_type;

  /** The multivector type used for the vectors x in
   *  $A x = b$ */
  typedef MultiVector<Vector> multivector_type;

  /** The multivector type used for the vectors b in
   *  $A x = b */
  typedef const MultiVector<const Vector> const_multivector_type;

  /** The type of the matrix A in $A x = \lambda x$ */
  typedef Matrix matrix_type;
  ///@}

  // Assert that all objects have the same scalar type
  static_assert(std::is_same<typename Vector::scalar_type, scalar_type>::value,
                "The Vector and the Matrix type need to have the same underlying stored "
                "matrix types.");

  // Assert that Matrices are matrices and vectors are vectors
  static_assert(IsMatrix<Matrix>::value, "Matrix needs to be a matrix type");
  static_assert(IsMutableVector<Vector>::value,
                "Vector needs to be a mutalbe vector type.");

  /** The operator A of the linear problem $Ax = b$. */
  const Matrix& A() const { return *m_A_ptr; }

  /** The RHS of the linear problem $Ax = b$, i.e. the vector $b$. */
  const_multivector_type& rhs() const { return *m_rhs_ptr; }

  /** Return the number of linear systems to solve
   *
   * Is equivalent to the number of RHS vectors
   * */
  size_type n_systems() const { return rhs().n_vectors(); }

  /** Return the dimensionality of the problem
   *
   * Equivalent to the number of columns of the problem matrix.
   */
  size_type dim() const { return A().n_cols(); }

  /** Constructor for the linear problem $Ax = b$
   *
   * \param A  The problem matrix
   * \param rhs  The rhs vector $b$
   * */
  LinearProblem(const Matrix& A, const_multivector_type& rhs)
        : m_A_ptr("LinearProblem", A), m_rhs_ptr("LinearProblem", rhs) {
    assert_size(A.n_rows(), rhs.n_elem());
  }

 private:
  /** The matrix within $A x = \lambda x$ or $A x = \lambda B x$ */
  krims::SubscriptionPointer<const Matrix> m_A_ptr;

  /** The rhs of the problem */
  krims::SubscriptionPointer<const_multivector_type> m_rhs_ptr;
};

}  // namespace linalgwrap
