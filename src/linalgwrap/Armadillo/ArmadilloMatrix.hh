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
#ifdef LINALGWRAP_HAVE_ARMADILLO

#include "ArmadilloTypes.hh"
#include "ArmadilloVector.hh"
#include "detail.hh"
#include "linalgwrap/Base/Interfaces/MutableMemoryVector_i.hh"
#include "linalgwrap/Base/Interfaces/Transposed.hh"
#include "linalgwrap/Constants.hh"
#include "linalgwrap/Exceptions.hh"
#include "linalgwrap/StoredMatrix_i.hh"
#include "linalgwrap/TypeUtils/mat_vec_apply_enabled_t.hh"
#include "linalgwrap/detail/scale_or_set.hh"
#include <initializer_list>
#include <memory>
#include <type_traits>

namespace linalgwrap {

// Forward-declare the interface class
template <typename Scalar>
class StoredMatrix_i;

/** A class for a dense stored matrix, currently implemented using armadillo.
 *
 * \note Armadillo is storing matrix data in column-major format
 *       (Fortran-style), but the default in our library is to
 *       use row-major storage (C-style) so in fact all data is
 *       stored in a transposed armadillo matrix such that memory
 *       access in a row-major fashion is contiguous.
 * */
template <typename Scalar>
class ArmadilloMatrix : public StoredMatrix_i<Scalar> {
  static_assert(std::is_same<double, Scalar>::value ||
                      std::is_same<float, Scalar>::value ||
                      std::is_same<std::complex<float>, Scalar>::value ||
                      std::is_same<std::complex<double>, Scalar>::value ||
                      std::is_same<short, Scalar>::value ||
                      std::is_same<int, Scalar>::value ||
                      std::is_same<long, Scalar>::value ||
                      std::is_same<unsigned short, Scalar>::value ||
                      std::is_same<unsigned int, Scalar>::value ||
                      std::is_same<unsigned long, Scalar>::value,
                "ArmadilloMatrix<Scalar> is currently only available for Scalar "
                "being one of double, float, complex<double>, "
                "complex<float>,  short, int, long, unsigned short, unsigned "
                "int, unsigned long");

 public:
  typedef StoredMatrix_i<Scalar> base_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::size_type size_type;

  /** The corresponding family of armadillo interface types */
  typedef ArmadilloTypes type_family;

  /** The corresponding vector type */
  typedef type_family::vector<scalar_type> vector_type;

  /** The type of the storage object used to store the data
   *  of the ArmadilloMatrix */
  typedef arma::Mat<Scalar> storage_type;

  // Swapping:
  template <typename S>
  friend void swap(ArmadilloMatrix<S>& first, ArmadilloMatrix<S>& second);

  /** \name Constructors
   */
  ///@{
  /** Construct a matrix of fixed size and optionally set the entries to
   * zero
   */
  ArmadilloMatrix(size_type n_rows, size_type n_cols, bool fill_zero = true);

  /** Construct a small matrix and copy all entries from ``mat`` which are
   *  not below the tolerance threshold.
   */
  ArmadilloMatrix(const ArmadilloMatrix& mat, scalar_type tolerance);

  /** \brief Construct from a nested initialiser list of scalars.
   *
   * The outermost layer gives the number of Columns, the innermost layer the
   * elements in each columns. An example would be
   * ```
   * ArmadilloMatrix mat{{1.0,2.0,0.5},{1.5,4.5,6.}})
   * ```
   * which produces a 2x3 matrix.
   */
  ArmadilloMatrix(
        std::initializer_list<std::initializer_list<scalar_type>> list_of_lists);

  /** Construct from transposed armadillo matrix
   *
   * \note The armadillo matrix has to be the transpose of the matrix
   * which is to be constructed. This is because internally we use a
   * transpose armadillo matrix in order to store the data.
   *
   * This is due to the fact that armadillo matrices are column-major,
   * but we are row-major.
   */
  explicit ArmadilloMatrix(storage_type inner) : m_arma(std::move(inner)) {}
  ///@}

  /** \name Matrix operations */
  ///@{
  /** Scale matrix by a scalar value */
  ArmadilloMatrix& operator*=(scalar_type s) {
    assert_finite(s);
    m_arma *= s;
    return *this;
  }

  /** Divide all matrix entries by a scalar value */
  ArmadilloMatrix& operator/=(scalar_type s) {
    assert_dbg(s != 0, krims::ExcDevideByZero());
    assert_finite(s);
    m_arma /= s;
    return *this;
  }

  /* Add a small matrix to this one */
  ArmadilloMatrix& operator+=(const ArmadilloMatrix& other) {
    assert_size(n_cols(), other.n_cols());
    assert_size(n_rows(), other.n_rows());
    m_arma += other.m_arma;
    return *this;
  }

  /* Subtract a small matrix from this one */
  ArmadilloMatrix& operator-=(const ArmadilloMatrix& other) {
    assert_size(n_cols(), other.n_cols());
    assert_size(n_rows(), other.n_rows());
    m_arma -= other.m_arma;
    return *this;
  }
  ///@}

  //
  // Relational operatiors
  //
  bool operator==(const ArmadilloMatrix& other) const {
    // does not work for some crazy arma reason
    // return (m_arma == other.m_arma);

    for (size_type i = 0; i < n_rows() * n_cols(); ++i) {
      if (m_arma[i] != other.m_arma[i]) return false;
    }
    return true;
  }

  bool operator!=(const ArmadilloMatrix& other) const { return !operator==(other); }

  //
  // matrix_i interface
  //
  /** \brief Number of rows of the matrix i
   */
  size_type n_rows() const override {
    // Note comment above the class definition why it
    // has to be this way round
    return m_arma.n_cols;
  }

  /** \brief Number of columns of the matrix
   */
  size_type n_cols() const override {
    // Note comment above the class definition why it
    // has to be this way round
    return m_arma.n_rows;
  }

  scalar_type operator()(size_type row, size_type col) const override {
    assert_greater(row, n_rows());
    assert_greater(col, n_cols());
    // Note comment above the class definition why it
    // has to be this way round
    return m_arma.at(col, row);
  }

  scalar_type operator[](size_type i) const override {
    assert_greater(i, n_cols() * n_rows());
    // Note that armadillo storage in column-major, but since
    // we store the transpose of what we represent internally
    // this give the correct interface (row-major access)
    return m_arma[i];
  }

  // We have a transpose operation mode available
  bool has_transpose_operation_mode() const override { return true; }

  /** Does this Matrix have an implemented inverse apply method?
   *
   * The idea is to allow special matrices to offer a method to take
   * advantage of their internal structure when one wants to apply
   * the inverse of such a matrix to a multivector. Possible examples
   * are diagonal and tridiagonal matrices.
   *
   * \note The inverse_apply should never be iterative.
   **/
  virtual bool has_apply_inverse() const override { return false; }

  /** \name Matrix application and matrix products
   */
  ///@{
  /** \brief Compute the Matrix-MultiVector application
   * For details see LazyMatrixExpression
   */
  template <typename VectorIn, typename VectorOut,
            mat_vec_apply_enabled_t<ArmadilloMatrix, VectorIn, VectorOut>...>
  void apply(const MultiVector<VectorIn>& x, MultiVector<VectorOut>& y,
             const Transposed mode = Transposed::None,
             const scalar_type c_this = Constants<scalar_type>::one,
             const scalar_type c_y = Constants<scalar_type>::zero) const;

  /** \brief Compute the application of the inverse of the matrix
   *  (or the inverse of the transpose of the matrix) to a MultiVector
   * For details see LazyMatrixExpression
   */
  template <typename VectorIn, typename VectorOut,
            mat_vec_apply_enabled_t<ArmadilloMatrix, VectorIn, VectorOut>...>
  void apply_inverse(const MultiVector<VectorIn>& /*x*/, MultiVector<VectorOut>& /*y*/,
                     const Transposed /*mode */ = Transposed::None,
                     const scalar_type /*c_this */ = 1,
                     const scalar_type /*c_y*/ = 0) const {
    // In general there is no easy way to do an inverse:
    assert_throw(false, krims::ExcDisabled("The apply_inverse function is in general "
                                           "very expensive and is only implemented in "
                                           "some cases. Use the function "
                                           "has_apply_inverse() to check when."));
  }

  /** Perform the Matrix-Vector product */
  template <typename Vector,
            typename = typename std::enable_if<IsStoredVector<Vector>::value>::type>
  MultiVector<typename std::remove_const<Vector>::type> operator*(
        const MultiVector<Vector>& v) const;

  /** Perform the Matrix-MultiVector product */
  template <typename Vector,
            typename = typename std::enable_if<IsStoredVector<Vector>::value>::type>
  Vector operator*(const Vector& v) const;

  /** \brief Compute the Matrix-Matrix product
   *
   * For more details see the docstring of the corresponding method
   * in LazyMatrixExpression */
  void mmult(const ArmadilloMatrix& in, ArmadilloMatrix& out,
             const Transposed mode = Transposed::None,
             const scalar_type c_this = Constants<scalar_type>::one,
             const scalar_type c_out = Constants<scalar_type>::zero) const;

  /** \brief Multiplication with a stored matrix */
  ArmadilloMatrix operator*(const ArmadilloMatrix& in) const {
    ArmadilloMatrix out(n_rows(), in.n_cols(), false);
    mmult(in, out);
    return out;
  }
  ///@}

  //
  // StoredMatrix_i interface
  //
  /** Set all elements to zero */
  void set_zero() override { m_arma.zeros(); }

  scalar_type& operator()(size_type row, size_type col) override {
    assert_greater(row, n_rows());
    assert_greater(col, n_cols());
    return m_arma.at(col, row);
  }

  /* Extract a block of the present matrix and copy it into a scaled
   * version of M.
   *
   * For more details see the docstring of the corresponding method
   * in LazyMatrixExpression */
  void extract_block(ArmadilloMatrix& M, const size_type start_row,
                     const size_type start_col, const Transposed mode = Transposed::None,
                     const scalar_type c_this = Constants<scalar_type>::one,
                     const scalar_type c_M = Constants<scalar_type>::zero) const;

  scalar_type& operator[](size_type i) override {
    assert_greater(i, n_cols() * n_rows());
    // Note that armadillo storage in column-major, but since
    // we store the transpose of what we represent internally
    // this give the correct interface (row-major access)
    return m_arma[i];
  }

  /** Read-only access to the inner storage
   *
   * \note Internally we store the data in transposed form, so
   * this will return an armadillo matrix representing the transpose
   * of the matrix this class represents.
   * */
  const storage_type& data() const { return m_arma; }

 private:
  storage_type m_arma;

  /** Performs
   *    y = c_y * y + c_this * A * x
   * If c_y is zero, then the values of y are never used
   */
  template <typename VectorIn, typename VectorOut>
  void apply_normal(const VectorIn& x, VectorOut& y, const scalar_type c_this,
                    const scalar_type c_y) const;

  /** Performs
   *    y = c_y * y + c_this * A^T * x
   * If c_y is zero, then the values of y are never used
   */
  template <typename VectorIn, typename VectorOut>
  void apply_transpose(const VectorIn& x, VectorOut& y, const scalar_type c_this,
                       const scalar_type c_y) const;

  /** Performs
   *    y = c_y * y + c_this * {A^{T})^* * x
   * If c_y is zero, then the values of y are never used
   */
  template <typename VectorIn, typename VectorOut>
  void apply_conjtranspose(const VectorIn& x, VectorOut& y, const scalar_type c_this,
                           const scalar_type c_y) const;
};

//
// Multiply by Scalar
//
template <typename Scalar>
ArmadilloMatrix<Scalar> operator*(Scalar s, ArmadilloMatrix<Scalar> m) {
  m *= s;
  return m;
}

template <typename Scalar>
ArmadilloMatrix<Scalar> operator*(ArmadilloMatrix<Scalar> m, Scalar s) {
  return s * m;
}

template <typename Scalar>
ArmadilloMatrix<Scalar> operator/(ArmadilloMatrix<Scalar> m, Scalar s) {
  m /= s;
  return m;
}

template <typename Scalar>
ArmadilloMatrix<Scalar> operator-(ArmadilloMatrix<Scalar> mat) {
  return -Constants<Scalar>::one * mat;
}

//
// Add and subtract small matrices
//
template <typename Scalar>
ArmadilloMatrix<Scalar> operator-(ArmadilloMatrix<Scalar> lhs,
                                  const ArmadilloMatrix<Scalar>& rhs) {
  lhs -= rhs;
  return lhs;
}

template <typename Scalar>
ArmadilloMatrix<Scalar> operator+(ArmadilloMatrix<Scalar> lhs,
                                  const ArmadilloMatrix<Scalar>& rhs) {
  lhs += rhs;
  return lhs;
}

//
// Specialisation of operations:
//
/** Compute the trace of the matrix
 *
 * \note only sensible for square matrices
 */
template <typename Scalar>
Scalar trace(const ArmadilloMatrix<Scalar>& m) {
  return arma::trace(m.data());
}

/** Accumulate all matrix values */
template <typename Scalar>
Scalar accumulate(const ArmadilloMatrix<Scalar>& m) {
  return arma::accu(m.data());
}

/** Compute the l1 norm (maximum of the sums over columns) */
template <typename Scalar>
Scalar norm_l1(const ArmadilloMatrix<Scalar>& m);

/** Calculate the linf norm (maximum of the sums over rows) */
template <typename Scalar>
Scalar norm_linf(const ArmadilloMatrix<Scalar>& m);

/** Calculate the Frobenius norm (sqrt of all matrix elements
 * squared
 *
 * \note This norm is not the matrix norm compatible to the l2 norm!
 */
template <typename Scalar>
Scalar norm_frobenius(const ArmadilloMatrix<Scalar>& m) {
  return arma::norm(m.data(), "fro");
}

/** Calculate the Frobenius norm squared
 *
 * \note This norm is not the matrix norm compatible to the l2 norm!
 */
template <typename Scalar>
Scalar norm_frobenius_squared(const ArmadilloMatrix<Scalar>& m) {
  // dot in arma is an elementwise dot product.
  return arma::cdot(m.data(), m.data());
}

//
// ---------------------------------------------------
//

template <typename Scalar>
void swap(ArmadilloMatrix<Scalar>& first, ArmadilloMatrix<Scalar>& second) {
  using std::swap;
  typedef typename ArmadilloMatrix<Scalar>::base_type base_type;
  swap(static_cast<base_type&>(first), static_cast<base_type&>(second));
  first.m_arma.swap(second.m_arma);
}

//
// Armadillo matrix
//
template <typename Scalar>
ArmadilloMatrix<Scalar>::ArmadilloMatrix(size_type n_rows, size_type n_cols,
                                         bool fill_zero)
      : m_arma(n_cols, n_rows, arma::fill::none) {
  // Note that the armadillo matrix stores the entries in transposed
  // form
  if (fill_zero) {
    // set all elements to zero
    m_arma.zeros();
  }
}

template <typename Scalar>
ArmadilloMatrix<Scalar>::ArmadilloMatrix(const ArmadilloMatrix& mat,
                                         scalar_type tolerance)
      : ArmadilloMatrix(mat.n_rows(), mat.n_cols(), false) {
  for (const auto elem : mat) {
    if (std::fabs(*elem) < tolerance) {
      (*this)(elem.row(), elem.col()) = Constants<scalar_type>::zero();
    } else {
      (*this)(elem.row(), elem.col()) = *elem;
    }
  }
}

template <typename Scalar>
ArmadilloMatrix<Scalar>::ArmadilloMatrix(
      std::initializer_list<std::initializer_list<scalar_type>> list_of_lists)
      : ArmadilloMatrix(list_of_lists.size(),
                        list_of_lists.size() > 0 ? list_of_lists.begin()->size() : 0,
                        false) {
#ifdef DEBUG
  size_type n_rows = list_of_lists.size();
  size_type n_cols = n_rows > 0 ? list_of_lists.begin()->size() : 0;
#endif

  // Assert all columns have equal length.
  assert_element_sizes(list_of_lists, n_cols);

  size_type i = 0;
  for (auto row : list_of_lists) {
    size_type j = 0;
    for (scalar_type elem : row) {
      (*this)(i, j) = elem;
      ++j;
    }
    assert_dbg(j == n_cols, krims::ExcInternalError());
    ++i;
  }
  assert_dbg(i == n_rows, krims::ExcInternalError());
}

// Matrix-Vector multiplication

template <typename Scalar>
template <typename VectorIn, typename VectorOut>
void ArmadilloMatrix<Scalar>::apply_normal(const VectorIn& x, VectorOut& y,
                                           const scalar_type c_this,
                                           const scalar_type c_y) const {
  const bool copy_into_arma = false;    // Do not copy the memory
  const bool fixed_vector_size = true;  // No memory reallocation

  const arma::Row<Scalar> x_arma(const_cast<Scalar*>(x.memptr()), x.size(),
                                 copy_into_arma, fixed_vector_size);
  arma::Row<Scalar> y_arma(y.memptr(), y.size(), copy_into_arma, fixed_vector_size);

  // Since m_arma is the transpose of what we represent, we
  // multiply from the left with a row-vector
  if (c_y == Constants<scalar_type>::zero) {
    y_arma = c_this * x_arma * m_arma;
  } else {
    y_arma = c_y * y_arma + c_this * x_arma * m_arma;
  }
}

template <typename Scalar>
template <typename VectorIn, typename VectorOut>
void ArmadilloMatrix<Scalar>::apply_transpose(const VectorIn& x, VectorOut& y,
                                              const scalar_type c_this,
                                              const scalar_type c_y) const {
  const bool copy_into_arma = false;    // Do not copy the memory
  const bool fixed_vector_size = true;  // No memory reallocation

  const arma::Col<Scalar> x_arma(const_cast<Scalar*>(x.memptr()), x.size(),
                                 copy_into_arma, fixed_vector_size);
  arma::Col<Scalar> y_arma(y.memptr(), y.size(), copy_into_arma, fixed_vector_size);

  // Since m_arma is the transpose of what we represent, we
  // multiply from the right with a column-vector
  // (no need to take transpose)
  if (c_y == Constants<scalar_type>::zero) {
    y_arma = c_this * m_arma * x_arma;
  } else {
    y_arma = c_y * y_arma + c_this * m_arma * x_arma;
  }
}

template <typename Scalar>
template <typename VectorIn, typename VectorOut>
void ArmadilloMatrix<Scalar>::apply_conjtranspose(const VectorIn& x, VectorOut& y,
                                                  const scalar_type c_this,
                                                  const scalar_type c_y) const {
  const bool copy_into_arma = false;    // Do not copy the memory
  const bool fixed_vector_size = true;  // No memory reallocation

  const arma::Col<Scalar> x_arma(const_cast<Scalar*>(x.memptr()), x.size(),
                                 copy_into_arma, fixed_vector_size);
  arma::Col<Scalar> y_arma(y.memptr(), y.size(), copy_into_arma, fixed_vector_size);

  // Since m_arma is the transpose of what we represent, we
  // multiply from the right with a column-vector
  // (no need to take transpose)
  if (c_y == Constants<scalar_type>::zero) {
    y_arma = c_this * arma::conj(m_arma) * x_arma;
  } else {
    y_arma = c_y * y_arma + c_this * arma::conj(m_arma) * x_arma;
  }
}

template <typename Scalar>
template <typename VectorIn, typename VectorOut,
          mat_vec_apply_enabled_t<ArmadilloMatrix<Scalar>, VectorIn, VectorOut>...>
void ArmadilloMatrix<Scalar>::apply(const MultiVector<VectorIn>& x,
                                    MultiVector<VectorOut>& y, const Transposed mode,
                                    const scalar_type c_this,
                                    const scalar_type c_y) const {
  assert_finite(c_this);
  assert_finite(c_y);
  assert_size(x.n_vectors(), y.n_vectors());
  if (mode == Transposed::Trans || mode == Transposed::ConjTrans) {
    assert_size(x.n_elem(), n_rows());
    assert_size(y.n_elem(), n_cols());
  } else {
    assert_size(x.n_elem(), n_cols());
    assert_size(y.n_elem(), n_rows());
  }
  assert_sufficiently_tested(mode != Transposed::ConjTrans);

  if (c_this == Constants<scalar_type>::zero) {
    for (auto& vec : y) detail::scale_or_set(vec, c_y);
    return;
  }  // c_this == 0

  for (size_type i = 0; i < x.n_vectors(); ++i) {
    switch (mode) {
      case Transposed::None:
        apply_normal(x[i], y[i], c_this, c_y);
        break;
      case Transposed::Trans:
        apply_transpose(x[i], y[i], c_this, c_y);
        break;
      case Transposed::ConjTrans:
        apply_conjtranspose(x[i], y[i], c_this, c_y);
        break;
    }  // mode
  }    // for
}

template <typename Scalar>
template <typename Vector, typename>
Vector ArmadilloMatrix<Scalar>::operator*(const Vector& v) const {
  assert_size(v.size(), n_cols());
  Vector out(n_rows(), false);
  apply_normal(v, out, Constants<scalar_type>::one, Constants<scalar_type>::zero);
  return out;
}

template <typename Scalar>
template <typename Vector, typename>
MultiVector<typename std::remove_const<Vector>::type> ArmadilloMatrix<Scalar>::operator*(
      const MultiVector<Vector>& mv) const {
  assert_size(mv.n_elem(), n_cols());
  MultiVector<typename std::remove_const<Vector>::type> out(n_rows(), mv.n_vectors(),
                                                            false);
  apply(mv, out, Transposed::None, Constants<scalar_type>::one,
        Constants<scalar_type>::zero);
  return out;
}

// mmult

template <typename Scalar>
void ArmadilloMatrix<Scalar>::mmult(const ArmadilloMatrix& in, ArmadilloMatrix& out,
                                    const Transposed mode, const scalar_type c_this,
                                    const scalar_type c_out) const {
  assert_finite(c_this);
  assert_finite(c_out);
  assert_size(in.n_cols(), out.n_cols());
  if (mode == Transposed::Trans || mode == Transposed::ConjTrans) {
    assert_size(n_rows(), in.n_rows());
    assert_size(n_cols(), out.n_rows());
  } else {
    assert_size(n_cols(), in.n_rows());
    assert_size(n_rows(), out.n_rows());
  }
  assert_sufficiently_tested(mode != Transposed::ConjTrans);

  if (c_this == Constants<scalar_type>::zero) {
    detail::scale_or_set(out, c_out);
    return;
  }  // c_this == 0

  // If c_out is zero we are not allowed to read the memory from out
  // since it could be uninitialised or nan
  const bool assign_out = c_out == Constants<scalar_type>::zero;

  // Transposing wrapper class
  typedef detail::ArmaTranspose<Scalar> AT;
  switch (mode) {
    case Transposed::None:
      // Since matrices are stored in a transposed sense
      // and B^T A^T = (AB)^T, i.e. what we need to store
      // in the product class
      if (assign_out) {
        out.m_arma = c_this * in.m_arma * m_arma;
      } else {
        out.m_arma = c_out * out.m_arma + c_this * in.m_arma * m_arma;
      }
      break;
    //
    case Transposed::Trans:
      if (assign_out) {
        out.m_arma = c_this * in.m_arma * AT::trans(m_arma);
      } else {
        out.m_arma = c_out * out.m_arma + c_this * in.m_arma * AT::trans(m_arma);
      }
      break;
    //
    case Transposed::ConjTrans:
      if (assign_out) {
        out.m_arma = c_this * in.m_arma * AT::conjtrans(m_arma);
      } else {
        out.m_arma = c_out * out.m_arma + c_this * in.m_arma * AT::conjtrans(m_arma);
      }
      break;
  }  // mode
}
template <typename Scalar>
void ArmadilloMatrix<Scalar>::extract_block(
      ArmadilloMatrix<Scalar>& M, const size_type start_row, const size_type start_col,
      const Transposed mode, const scalar_type c_this, const scalar_type c_M) const {
  assert_finite(c_this);
  assert_finite(c_M);
  // check that we do not overshoot the indices
  if (mode == Transposed::Trans || mode == Transposed::ConjTrans) {
    assert_greater_equal(start_row + M.n_rows(), n_cols());
    assert_greater_equal(start_col + M.n_cols(), n_rows());
  } else {
    assert_greater_equal(start_row + M.n_rows(), n_rows());
    assert_greater_equal(start_col + M.n_cols(), n_cols());
  }
  assert_sufficiently_tested(mode != Transposed::ConjTrans);

  // For empty matrices there is nothing to do
  if (M.n_rows() == 0 || M.n_cols() == 0) return;

  if (c_this == Constants<scalar_type>::zero) {
    detail::scale_or_set(M, c_M);
    return;
  }  // c_this == 0

  // If c_M is zero we are not allowed to read the memory from M
  // since it could be uninitialised or nan
  const bool assign_M = c_M == Constants<scalar_type>::zero;

  // Translate ranges of interesting rows and columns
  // to armadillo spans (which are closed intervals)
  arma::span rows(start_row, start_row + M.n_rows() - 1);
  arma::span cols(start_col, start_col + M.n_cols() - 1);

  typedef detail::ArmaTranspose<Scalar> AT;
  switch (mode) {
    case Transposed::None:
      // Note that the order of rows and cols has to be
      // interchanged, see class documentation for details
      if (assign_M) {
        M.m_arma = c_this * m_arma(cols, rows);
      } else {
        M.m_arma = c_M * M.m_arma + c_this * m_arma(cols, rows);
      }
      break;
    //
    case Transposed::Trans:
      if (assign_M) {
        M.m_arma = c_this * AT::trans(m_arma(rows, cols));
      } else {
        M.m_arma = c_M * M.m_arma + c_this * AT::trans(m_arma(rows, cols));
      }
      break;
    //
    case Transposed::ConjTrans:
      if (assign_M) {
        M.m_arma = c_this * AT::conjtrans(m_arma(rows, cols));
      } else {
        M.m_arma = c_M * M.m_arma + c_this * AT::conjtrans(m_arma(rows, cols));
      }
      break;
  }  // mode
}

//
// Out of class
//
template <typename Scalar>
Scalar norm_l1(const ArmadilloMatrix<Scalar>& m) {
  // l1 is the maximum of the sums over columns
  // linf is the maximum of the sums over rows
  //
  // Since m.data() is stored as the transpose of the
  // thing it actually represents, we can use linf
  // on the m.data() here --- in theory
  //
  if (m.data().n_rows == 1) {
    // Arma wants to be clever and treats matrices of one row
    // exactly like vectors (so here we do not need to transpose
    // and use the inf norm)
    return norm(m.data(), 1);
  } else {
    return norm(m.data(), "inf");
  }
}

template <typename Scalar>
Scalar norm_linf(const ArmadilloMatrix<Scalar>& m) {
  // See norm_l1 for reasons why we need the condition here
  // Note that m.data() is the *transpose* of the actual matrix
  // we represent. Therefore we calculate the l1 norm of it here.
  //
  if (m.data().n_rows == 1) {
    return norm(m.data(), "inf");
  } else {
    return norm(m.data(), 1);
  }
}

}  // namespace linalgwrap
#endif  // LINALGWRAP_HAVE_ARMADILLO
