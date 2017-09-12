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
#include "lazyten/config.hh"
#ifdef LAZYTEN_HAVE_BOHRIUM

#include "BohriumTypes.hh"
#include "BohriumVector.hh"
#include "common.hh"
#include "lazyten/Base/Interfaces/Transposed.hh"
#include "lazyten/StoredMatrix_i.hh"
#include "lazyten/TypeUtils/mat_vec_apply_enabled_t.hh"
#include "lazyten/detail/scale_or_set.hh"
#include <bhxx/bhxx.hpp>

namespace lazyten {

// Forward-declare the interface class
template <typename Scalar>
class StoredMatrix_i;

/** A class for Bohrium matrices */
template <typename Scalar>
class BohriumMatrix : public StoredMatrix_i<Scalar> {
  static_assert(std::is_same<double, Scalar>::value ||
                      std::is_same<float, Scalar>::value ||
                      std::is_same<std::complex<float>, Scalar>::value ||
                      std::is_same<std::complex<double>, Scalar>::value,
                "BohriumMatrix<Scalar> is currently only available for Scalar "
                "being one of double, float, complex<double>, complex<float>");

 public:
  typedef StoredMatrix_i<Scalar> base_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::size_type size_type;

  /** The family of Bohrium interface types */
  typedef BohriumTypes type_family;

  /** The corresponding vector type */
  typedef type_family::vector<scalar_type> vector_type;

  /** TODO Custom iteratiors for BohriumMatrix objects
   *    typedef Scalar* iterator;
   *    typedef const Scalar* const_iterator;
   */

  /** \name Constructors
   */
  ///@{
  /** Construct a matrix of fixed size and optionally set the entries to
   * zero
   */
  BohriumMatrix(size_type n_rows, size_type n_cols, bool /* fill_zero */ = true)
        : m_array({n_rows, n_cols}) {
    assert_greater(0, n_rows * n_cols);  // Bohrium arrays may not be empty
    set_zero();  // Here we fill the vector with zeros in all cases
  }

  /** Construct a small matrix and copy all entries from ``mat`` which are
   *  not below the tolerance threshold.
   */
  BohriumMatrix(const BohriumMatrix& mat, scalar_type tolerance);

  /** \brief Construct from a nested initialiser list of scalars.
   *
   * The outermost layer gives the number of Columns, the innermost layer the
   * elements in each columns. An example would be
   * ```
   * BohriumMatrix mat{{1.0,2.0,0.5},{1.5,4.5,6.}})
   * ```
   * which produces a 2x3 matrix.
   */
  BohriumMatrix(std::initializer_list<std::initializer_list<scalar_type>> list_of_lists);

  /** \brief Construct from a two-dimensional multi-array */
  explicit BohriumMatrix(bhxx::BhArray<Scalar> inner) : m_array(inner.shape) {
    assert_size(inner.rank(), 2);
    bhxx::identity(m_array, inner);
  }

  BohriumMatrix(BohriumMatrix&& other) = default;
  BohriumMatrix& operator=(BohriumMatrix&& other) = default;
  BohriumMatrix(const BohriumMatrix& other) : m_array({other.n_rows(), other.n_cols()}) {
    bhxx::identity(m_array, other.m_array);
  }

  BohriumMatrix& operator=(const BohriumMatrix& other) {
    if (m_array.shape != other.m_array.shape) {
      // If the shape changes we need to allocate a new
      // object first and only then call the identity operation
      m_array = bhxx::BhArray<Scalar>(other.m_array.shape);
    }
    bhxx::identity(m_array, other.m_array);
    return *this;
  }

  /** \name Matrix operations */
  ///@{
  /** Scale matrix by a scalar value */
  BohriumMatrix& operator*=(scalar_type s) {
    assert_finite(s);
    bhxx::multiply(m_array, m_array, s);
    return *this;
  }

  /** Divide all matrix entries by a scalar value */
  BohriumMatrix& operator/=(scalar_type s) {
    assert_dbg(std::abs(s) != 0, krims::ExcDevideByZero());
    assert_finite(s);
    bhxx::divide(m_array, m_array, s);
    return *this;
  }

  /* Add a matrix to this one */
  BohriumMatrix& operator+=(const BohriumMatrix& other) {
    assert_size(n_cols(), other.n_cols());
    assert_size(n_rows(), other.n_rows());
    bhxx::add(m_array, m_array, other.m_array);
    return *this;
  }

  /* Subtract a matrix from this one */
  BohriumMatrix& operator-=(const BohriumMatrix& other) {
    assert_size(n_cols(), other.n_cols());
    assert_size(n_rows(), other.n_rows());
    bhxx::subtract(m_array, m_array, other.m_array);
    return *this;
  }
  ///@}

  //
  // Relational operators
  //
  bool operator==(const BohriumMatrix& other) const;
  bool operator!=(const BohriumMatrix& other) const;

  //
  // Matrix_i interface
  //
  /** \brief Number of rows of the matrix */
  size_type n_rows() const override { return m_array.shape[0]; }

  /** \brief Number of columns of the matrix */
  size_type n_cols() const override { return m_array.shape[1]; }

  /** \name Size of the vector */
  size_type n_elem() const override { return m_array.n_elem(); }

  /** \name Size of the vector */
  size_type size() const override { return m_array.n_elem(); }

  scalar_type operator()(size_type row, size_type col) const override {
    assert_greater(row, n_rows());
    assert_greater(col, n_cols());
    return memptr()[row * n_cols() + col];
  }

  scalar_type operator[](size_type i) const override {
    assert_greater(i, n_elem());
    return memptr()[i];
  }

  // We have a transpose operation mode available
  bool has_transpose_operation_mode() const override { return true; }

  /** Does this Matrix have an implemented inverse apply method? */
  bool has_apply_inverse() const override { return false; }

  /** \name Matrix application and matrix products
   */
  ///@{
  /** \brief Compute the Matrix-MultiVector application
   * For details see LazyMatrixExpression
   */
  template <typename VectorIn, typename VectorOut,
            mat_vec_apply_enabled_t<BohriumMatrix, VectorIn, VectorOut>...>
  void apply(const MultiVector<VectorIn>& x, MultiVector<VectorOut>& y,
             const Transposed mode = Transposed::None, const scalar_type c_this = 1,
             const scalar_type c_y = 0) const;

  /** \brief Compute the application of the inverse of the matrix
   *  (or the inverse of the transpose of the matrix) to a MultiVector
   * For details see LazyMatrixExpression
   */
  template <typename VectorIn, typename VectorOut,
            mat_vec_apply_enabled_t<BohriumMatrix, VectorIn, VectorOut>...>
  void apply_inverse(const MultiVector<VectorIn>& /*x*/,
                     MultiVector<VectorOut>&
                     /*y*/,
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
        const MultiVector<Vector>& mv) const {
    assert_size(mv.n_elem(), n_cols());
    MultiVector<typename std::remove_const<Vector>::type> out(n_rows(), mv.n_vectors(),
                                                              false);
    apply(mv, out, Transposed::None, 1, 0);
    return out;
  }

  /** Perform the Matrix-MultiVector product */
  template <typename Vector,
            typename = typename std::enable_if<IsStoredVector<Vector>::value>::type>
  Vector operator*(const Vector& v) const {
    Vector out(v.size(), false);
    apply(m_array, v, out, 1, 0);
  }

  /** \brief Compute the Matrix-Matrix product
   *
   * For more details see the docstring of the corresponding method
   * in LazyMatrixExpression */
  void mmult(const BohriumMatrix& in, BohriumMatrix& out,
             const Transposed mode = Transposed::None, const scalar_type c_this = 1,
             const scalar_type c_out = 0) const;

  /** \brief Multiplication with a stored matrix */
  BohriumMatrix operator*(const BohriumMatrix& in) const {
    BohriumMatrix out(n_rows(), in.n_cols(), false);
    mmult(in, out);
    return out;
  }
  ///@}

  //  /** \name Iterators
  //   */
  //  ///@{
  //  /** Return an iterator to the beginning */
  //  iterator begin() { return memptr(); }
  //
  //  /** Return a const_iterator to the beginning */
  //  const_iterator begin() const { return memptr(); }
  //
  //  /** Return a const_iterator to the beginning */
  //  const_iterator cbegin() const { return memptr(); }
  //
  //  /** Return an iterator to the end */
  //  iterator end() { return memptr() + n_elem(); }
  //
  //  /** Return a const_iterator to the end */
  //  const_iterator end() const { return memptr() + n_elem(); }
  //
  //  /** Return a const_iterator to the end */
  //  const_iterator cend() const { return memptr() + n_elem(); }
  //  ///@}

  // StoredMatrix_i interface
  //
  /** Set all elements to zero */
  void set_zero() override { bhxx::identity(m_array, 0); }

  scalar_type& operator()(size_type row, size_type col) override {
    assert_greater(row, n_rows());
    assert_greater(col, n_cols());
    return memptr()[row * n_cols() + col];
  }

  /* Extract a block of the present matrix and copy it into a scaled
   * version of M.
   *
   * For more details see the docstring of the corresponding method
   * in LazyMatrixExpression */
  void extract_block(BohriumMatrix& M, const size_type start_row,
                     const size_type start_col, const Transposed mode = Transposed::None,
                     const scalar_type c_this = 1, const scalar_type c_M = 0) const;

  scalar_type& operator[](size_type i) override {
    assert_greater(i, n_elem());
    return memptr()[i];
  }

  /** Read-only access to the raw memory */
  const scalar_type* memptr() const {
    prepare_array_memory_for_use();
    return m_array.data();
  }

  /** Access to the raw memory */
  scalar_type* memptr() {
    prepare_array_memory_for_use();
    return m_array.data();
  }

  //@{
  /** Access to the inner BhArray */
  const bhxx::BhArray<Scalar>& bh_array() const { return m_array; }
  bhxx::BhArray<Scalar>& bh_array() { return m_array; }
  //@}
 private:
  bhxx::BhArray<Scalar> m_array;

  /** Initialise the m_array object and/or make it contiguous */
  void prepare_array_memory_for_use() const {
    // Cast the const away ...
    auto& array = const_cast<bhxx::BhArray<Scalar>&>(m_array);

    if (!array.is_contiguous()) {
      // Get contiguous representation of the array:
      bhxx::BhArray<Scalar> contiguous{array.shape};
      identity(contiguous, array);
      array = std::move(contiguous);
    }

    sync(array);
    bhxx::Runtime::instance().flush();

    assert_internal(m_array.is_contiguous());
    assert_internal(m_array.data() != nullptr);
  }

  //@{
  /** Performs
   *    y = c_y * y + c_this * A * x
   * If c_y is zero, then the values of y are never used
   */
  template <typename VectorIn, typename VectorOut>
  void apply(const bhxx::BhArray<Scalar>& A, const VectorIn& x, VectorOut& y,
             const scalar_type c_this, const scalar_type c_y) const;

  void apply(const bhxx::BhArray<Scalar>& A, const BohriumVector<Scalar>& x,
             BohriumVector<Scalar>& y, const scalar_type c_this,
             const scalar_type c_y) const {
    apply(A, x.bh_array(), y.bh_array(), c_this, c_y);
  }

  void apply(const bhxx::BhArray<Scalar>& A, const bhxx::BhArray<Scalar>& x,
             bhxx::BhArray<Scalar>& y, const scalar_type c_this,
             const scalar_type c_y) const;
  //@}
};

//
// Multiply by Scalar
//
template <typename Scalar>
BohriumMatrix<Scalar> operator*(Scalar s, BohriumMatrix<Scalar> m) {
  m *= s;
  return m;
}

template <typename Scalar>
BohriumMatrix<Scalar> operator*(BohriumMatrix<Scalar> m, Scalar s) {
  return s * m;
}

template <typename Scalar>
BohriumMatrix<Scalar> operator/(BohriumMatrix<Scalar> m, Scalar s) {
  m /= s;
  return m;
}

template <typename Scalar>
BohriumMatrix<Scalar> operator-(BohriumMatrix<Scalar> mat) {
  return Scalar(-1) * mat;
}

//
// Add and subtract bohrium matrices
//
template <typename Scalar>
BohriumMatrix<Scalar> operator-(BohriumMatrix<Scalar> lhs,
                                const BohriumMatrix<Scalar>& rhs) {
  lhs -= rhs;
  return lhs;
}

template <typename Scalar>
BohriumMatrix<Scalar> operator+(BohriumMatrix<Scalar> lhs,
                                const BohriumMatrix<Scalar>& rhs) {
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
Scalar trace(const BohriumMatrix<Scalar>& m) {
  assert_size(m.n_rows(), m.n_cols());
  bhxx::BhArray<Scalar> diagonal(m.bh_array().base, {m.n_cols()},
                                 {static_cast<int64_t>(m.n_rows() + 1)});
  return bhxx::as_scalar(bhxx::accumulate(diagonal));
}

/** Accumulate all matrix values */
template <typename Scalar>
Scalar accumulate(const BohriumMatrix<Scalar>& m) {
  return bhxx::as_scalar(bhxx::accumulate(m.bh_array()));
}

/** Compute the l1 norm (maximum of the sums over columns) */
template <typename Scalar>
Scalar norm_l1(const BohriumMatrix<Scalar>& m);

/** Calculate the linf norm (maximum of the sums over rows) */
template <typename Scalar>
Scalar norm_linf(const BohriumMatrix<Scalar>& m);

/** Calculate the Frobenius norm (sqrt of all matrix elements
 * squared
 *
 * \note This norm is not the matrix norm compatible to the l2 norm!
 */
template <typename Scalar>
Scalar norm_frobenius(const BohriumMatrix<Scalar>& m) {
  return std::sqrt(norm_frobenius_squared(m));
}

/** Calculate the Frobenius norm squared
 *
 * \note This norm is not the matrix norm compatible to the l2 norm!
 */
template <typename Scalar>
Scalar norm_frobenius_squared(const BohriumMatrix<Scalar>& m) {
  assert_implemented(!krims::IsComplexNumber<Scalar>::value);
  // TODO Need the complex conjugate on the lhs, i.e. tr(M^H * M)
  return bhxx::as_scalar(bhxx::inner_product(m.bh_array(), m.bh_array()));
}

// TODO Specialisation of elementwise functions like
//        sqrt, exp, ...
//        min and max

//
// ---------------------------------------------------
//

template <typename Scalar>
template <typename VectorIn, typename VectorOut>
void BohriumMatrix<Scalar>::apply(const bhxx::BhArray<Scalar>& A, const VectorIn& x,
                                  VectorOut& y, const scalar_type c_this,
                                  const scalar_type c_y) const {
  using bhxx::make_base_ptr;

  // Make Bohrium versions of x and y
  Scalar* x_ptr = const_cast<Scalar*>(x.memptr());
  bhxx::BhArray<Scalar> x_bh(make_base_ptr(x.size(), x_ptr), {x.size()});
  bhxx::BhArray<Scalar> y_bh(make_base_ptr(y.size(), y.memptr()), {y.size()});

  // Perform the application
  apply(A, x_bh, y_bh, c_this, c_y);

  // Sync and flush the output array
  bhxx::sync(y_bh);
  bhxx::Runtime::instance().flush();
}

template <typename Scalar>
template <typename VectorIn, typename VectorOut,
          mat_vec_apply_enabled_t<BohriumMatrix<Scalar>, VectorIn, VectorOut>...>
void BohriumMatrix<Scalar>::apply(const MultiVector<VectorIn>& x,
                                  MultiVector<VectorOut>& y, const Transposed mode,
                                  const scalar_type c_this, const scalar_type c_y) const {
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
  assert_implemented(mode != Transposed::ConjTrans);

  if (c_this == 0) {
    for (auto& vec : y) detail::scale_or_set(vec, c_y);
    return;
  }  // c_this == 0

  bhxx::BhArray<Scalar> trans = transpose(m_array);
  for (size_type i = 0; i < x.n_vectors(); ++i) {
    switch (mode) {
      case Transposed::None:
        apply(m_array, x[i], y[i], c_this, c_y);
        break;
      case Transposed::Trans:
        apply(trans, x[i], y[i], c_this, c_y);
        break;
      case Transposed::ConjTrans:
        assert_implemented(false);
        break;
    }  // mode
  }    // for
}

//
// Out of class
//
template <typename Scalar>
Scalar norm_l1(const BohriumMatrix<Scalar>& m) {
  // l1 is the maximum of the sums over columns
  bhxx::BhArray<Scalar> mabs(m.bh_array().shape);
  bhxx::absolute(mabs, m.bh_array());

  bhxx::BhArray<Scalar> colsum({mabs.shape[1]});
  bhxx::add_reduce(colsum, mabs, 0);

  bhxx::BhArray<Scalar> result({1});
  bhxx::maximum_reduce(result, colsum, 0);
  return bhxx::as_scalar(result);
}

template <typename Scalar>
Scalar norm_linf(const BohriumMatrix<Scalar>& m) {
  // linf is the maximum of the sums over rows
  bhxx::BhArray<Scalar> mabs(m.bh_array().shape);
  bhxx::absolute(mabs, m.bh_array());

  bhxx::BhArray<Scalar> rowsum({mabs.shape[0]});
  bhxx::add_reduce(rowsum, mabs, 1);

  bhxx::BhArray<Scalar> result({1});
  bhxx::maximum_reduce(result, rowsum, 0);
  return bhxx::as_scalar(result);
}

}  // namespace lazyten
#endif  // LAZYTEN_HAVE_BOHRIUM
