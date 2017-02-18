
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
#include "detail/scale_or_set.hh"
#include "linalgwrap/Constants.hh"
#include "linalgwrap/LazyMatrixExpression.hh"
#include "linalgwrap/LazyMatrixProduct.hh"
#include "linalgwrap/MultiVector.hh"
#include <algorithm>
#include <iterator>
#include <krims/GenMap.hh>
#include <krims/SubscriptionPointer.hh>
#include <vector>

namespace linalgwrap {

// forward declaration
template <typename StoredMatrix>
class LazyMatrixProduct;

template <typename StoredMatrix>
class LazyMatrixExpression;

// TODO: If sparsity should be incorporated, build up the sparsity pattern
//       as more and more terms are added to the sum object.
//
// TODO: Maybe have a scaling coefficient for the whole object as well
// this should improve numerical accuracy and maybe allows to easily
// handle sums of lazy objects without wrapping each of them inside a product

/** Macro to aid defining both in-place matrix-matrix addition and in-place
 * matrix-matrix subtraction */
#define LazyMatrixSum_inplace_addsub_op(othertype)        \
  LazyMatrixSum& operator+=(othertype other) {            \
    this->push_term(other);                               \
    return *this;                                         \
  }                                                       \
  LazyMatrixSum& operator-=(othertype other) {            \
    this->push_term(other, -Constants<scalar_type>::one); \
    return *this;                                         \
  }

/** Class to represent the sum of different MatrixProducts
 * It may include stored terms as well.
 */
template <typename StoredMatrix>
class LazyMatrixSum : public LazyMatrixExpression<StoredMatrix> {
 public:
  typedef LazyMatrixExpression<StoredMatrix> base_type;
  typedef typename base_type::stored_matrix_type stored_matrix_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::lazy_matrix_expression_ptr_type
        lazy_matrix_expression_ptr_type;

  /** The type of the lazy terms contained in this sum */
  typedef LazyMatrixProduct<stored_matrix_type> lazy_term_type;

 private:
  /** Type used for terms which are stored matrices
   *
   * Stores a reference to a stored matrix (by the means of a
   * subscription pointer) and a coefficient.
   *
   * @note: This implies that changing the stored matrix after
   *        enwrapping it in this term structure changes the
   *        value of this term as well!
   */
  class StoredTerm {
   public:
    StoredTerm(const stored_matrix_type& matrix, scalar_type coefficient)
          : m_coefficient{coefficient},
            m_matrix_ptr{make_subscription(matrix, "LazyMatrixSum")} {}

    //! get the coefficient
    scalar_type coefficient() const { return m_coefficient; }

    //! get the matrix as a reference
    const stored_matrix_type& matrix() const { return *m_matrix_ptr; }

    //! adjust the coefficient
    StoredTerm& operator*=(scalar_type c) {
      m_coefficient *= c;
      return *this;
    }

   private:
    scalar_type m_coefficient;
    krims::SubscriptionPointer<const stored_matrix_type> m_matrix_ptr;
    // TODO maybe use RCPWrapper here and allow to move matrix inside object
    //      as well?
  };

  /** Have an alias for the StoredTerm class */
  typedef StoredTerm stored_term_type;

 public:
  /** A swap function for LazyMatrixSums */
  friend void swap(LazyMatrixSum& first, LazyMatrixSum& second) {
    using std::swap;  // enable ADL

    swap(first.m_n_rows, second.m_n_rows);
    swap(first.m_n_cols, second.m_n_cols);
    swap(first.m_lazy_terms, second.m_lazy_terms);
    swap(first.m_stored_terms, second.m_stored_terms);
    swap(static_cast<base_type&>(first), static_cast<base_type&>(second));
  }

 public:
  /** \name Constructors for lazy sums
   */
  ///@{
  /** \brief Create an empty matrix sum object
   *
   * An empty matrix sum object has size 0 times 0.
   * It behaves as the additive identity (0), i.e. when
   * retrieving elements, they will all be zero.
   *
   * When the first summand is added it inherits its size.
   */
  explicit LazyMatrixSum() : m_n_rows(0), m_n_cols(0) {}

  /** \brief Create a matrix sum object
   *
   * @param term   The first matrix expression term
   */
  explicit LazyMatrixSum(lazy_term_type term)
        : m_n_rows(term.n_rows()), m_n_cols(term.n_cols()) {
    m_lazy_terms.push_back(std::move(term));
  }

  /** \brief Create a matrix sum object
   *
   * @param expr    The first matrix expression term
   * @param factor  The factor to multiply the matrix expression with.
   */
  explicit LazyMatrixSum(const LazyMatrixExpression<StoredMatrix>& expr,
                         scalar_type factor = Constants<scalar_type>::one)
        : LazyMatrixSum(lazy_term_type(expr, factor)) {}

  /** \brief Create a matrix sum object
   *
   * In this case a subscription to the term is made, meaning that
   * the term object should not be destroyed before the LazyMatrixSum
   * object. Otherwise an exception will throw in DEBUG mode.
   *
   * @param mat   The stored matrix to use as first term.
   * @param factor The factor to multiply the matrix with
   */
  explicit LazyMatrixSum(const stored_matrix_type& mat,
                         scalar_type factor = Constants<scalar_type>::one)
        : m_n_rows{mat.n_rows()}, m_n_cols{mat.n_cols()} {

    // Construct a term and subscribe to the reference
    stored_term_type term(mat, factor);

    // Push term onto the stored stack
    m_stored_terms.push_back(std::move(term));
  }

  //@{
  /** Defaults for the big five */
  ~LazyMatrixSum() = default;
  LazyMatrixSum(const LazyMatrixSum&) = default;
  LazyMatrixSum(LazyMatrixSum&&) = default;
  LazyMatrixSum& operator=(const LazyMatrixSum&) = default;
  LazyMatrixSum& operator=(LazyMatrixSum&&) = default;
  //@}
  ///@}

  //
  // Push back further terms
  //
  /** Push back a further lazy matrix term */
  void push_term(lazy_term_type term) {
    if (empty()) {
      m_n_rows = term.n_rows();
      m_n_cols = term.n_cols();
    }
    assert_size(n_cols(), term.n_cols());
    assert_size(n_rows(), term.n_rows());

    // Push back
    m_lazy_terms.push_back(std::move(term));
  }

  /** Push back a further lazy matrix expression */
  void push_term(const LazyMatrixExpression<StoredMatrix>& term,
                 scalar_type factor = Constants<scalar_type>::one) {
    // Reduce to version above
    push_term(lazy_term_type(term, factor));
  }

  /** \brief Push back a further stored matrix
   *
   * In this case a subscription to the term is made, meaning that
   * the term object should not be destroyed before the LazyMatrixSum
   * object. Otherwise an exception will throw in DEBUG mode.
   * */
  void push_term(const stored_matrix_type& mat,
                 scalar_type factor = Constants<scalar_type>::one) {
    if (empty()) {
      m_n_rows = mat.n_rows();
      m_n_cols = mat.n_cols();
    }
    assert_size(n_cols(), mat.n_cols());
    assert_size(n_rows(), mat.n_rows());

    // Construct a term and subscribe to the reference
    stored_term_type term(mat, factor);

    m_stored_terms.push_back(std::move(term));
  }

  void push_term(LazyMatrixSum sum) {
    if (empty()) {
      m_n_rows = sum.n_rows();
      m_n_cols = sum.n_cols();
    }
    assert_size(n_cols(), sum.n_cols());
    assert_size(n_rows(), sum.n_rows());

    // Move all lazy terms of sum to the end of this object
    std::move(std::begin(sum.m_lazy_terms), std::end(sum.m_lazy_terms),
              back_inserter(m_lazy_terms));

    // Move all stored terms of sum to the end of this object
    std::move(std::begin(sum.m_stored_terms), std::end(sum.m_stored_terms),
              back_inserter(m_stored_terms));
  }

  //
  // Other modifications of the sum
  //
  /** \brief scale the whole sum by the scalar \p c */
  void scale(const scalar_type c) {
    assert_finite(c);

    for (auto& term : m_stored_terms) {
      term *= c;
    }

    for (auto& term : m_lazy_terms) {
      term *= c;
    }
  }

  //
  // Matrix_i interface
  //

  /** \brief Number of rows of the matrix      */
  size_type n_rows() const override { return m_n_rows; }

  /** \brief Number of columns of the matrix     */
  size_type n_cols() const override { return m_n_cols; }

  /** \brief return an element of the matrix    */
  scalar_type operator()(size_type row, size_type col) const override;

  //
  // LazyMatrixExpression interface
  //
  /** Are operation modes Transposed::Trans and Transposed::ConjTrans
   *  supported for this matrix type.
   **/
  bool has_transpose_operation_mode() const override;

  /** Extract a block of a matrix and (optionally) add it to
   * a different matrix.
   *
   *  Loosely speaking we perform
   *  \[ M = c_M \cdot M + (A^{mode})_{rowrange,colrange} \]
   *  where
   *    - rowrange = [start_row, start_row+in.n_rows() ) and
   *    - colrange = [start_col, start_col+in.n_cols() )
   *
   * More details can be found in the same function in
   * LazyMatrixExpression
   */
  void extract_block(stored_matrix_type& M, const size_type start_row,
                     const size_type start_col, const Transposed mode = Transposed::None,
                     const scalar_type c_this = Constants<scalar_type>::one,
                     const scalar_type c_M = Constants<scalar_type>::zero) const override;

  /** \brief Compute the Matrix-Multivector application -- generic version
   *
   * Loosely speaking we perform
   * \[ y = c_this \cdot A^\text{mode} \cdot x + c_y \cdot y. \]
   *
   * See LazyMatrixExpression for more details
   *
   * \note Whenever the virtual apply method is overwritten, this method
   * should be implemented as well as it assures that conversion to
   * MultiVector<MutableMemoryVector_i<scalar_type>> can actually occur
   * automatically.
   */
  template <typename VectorIn, typename VectorOut,
            mat_vec_apply_enabled_t<LazyMatrixSum, VectorIn, VectorOut>...>
  void apply(const MultiVector<VectorIn>& x, MultiVector<VectorOut>& y,
             const Transposed mode = Transposed::None,
             const scalar_type c_this = Constants<scalar_type>::one,
             const scalar_type c_y = Constants<scalar_type>::zero) const {
    MultiVector<const MutableMemoryVector_i<scalar_type>> x_wrapped(x);
    MultiVector<MutableMemoryVector_i<scalar_type>> y_wrapped(y);
    apply(x_wrapped, y_wrapped, mode, c_this, c_y);
  }

  /** \brief Compute the Matrix-Multivector application
   *
   * Loosely speaking we perform
   * \[ y = c_this \cdot A^\text{mode} \cdot x + c_y \cdot y. \]
   *
   * See LazyMatrixExpression for more details
   */
  void apply(const MultiVector<const MutableMemoryVector_i<scalar_type>>& x,
             MultiVector<MutableMemoryVector_i<scalar_type>>& y,
             const Transposed mode = Transposed::None,
             const scalar_type c_this = Constants<scalar_type>::one,
             const scalar_type c_y = Constants<scalar_type>::zero) const override;

  /** Perform a matrix-matrix product.
   *
   * Loosely performs the operation
   * \[ out = c_this \cdot A^\text{mode} \cdot in + c_out \cdot out. \]
   *
   * See LazyMatrixExpression for more details
   */
  void mmult(const stored_matrix_type& in, stored_matrix_type& out,
             const Transposed mode = Transposed::None,
             const scalar_type c_this = Constants<scalar_type>::one,
             const scalar_type c_out = Constants<scalar_type>::zero) const override;

  /** \brief Update the internal data of all objects in this expression
   *         given the GenMap
   * */
  void update(const krims::GenMap& map) override {
    // Pass the call onto all factors:
    for (auto& expression : m_lazy_terms) {
      expression.update(map);
    }
  }

  /** \brief Clone the expression */
  lazy_matrix_expression_ptr_type clone() const override {
    // return a copy enwrapped in the pointer type
    return lazy_matrix_expression_ptr_type(new LazyMatrixSum(*this));
  }

  /** \brief Is this object empty? */
  bool empty() const { return m_stored_terms.size() == 0 && m_lazy_terms.size() == 0; }

  //
  // In-place scaling operators
  //
  /** \brief scale a matrix using a scalar
   */
  LazyMatrixSum& operator*=(scalar_type c) {
    assert_finite(c);
    this->scale(c);
    return *this;
  }

  /** \brief scale a matrix using a scalar
   */
  LazyMatrixSum& operator/=(scalar_type c) {
    assert_finite(c);
    assert_dbg(c != 0, krims::ExcDevideByZero());
    this->scale(Constants<scalar_type>::one / c);
    return *this;
  }

  //
  // Operators with += or -=
  //
  LazyMatrixSum_inplace_addsub_op(lazy_term_type);
  LazyMatrixSum_inplace_addsub_op(const LazyMatrixExpression<StoredMatrix>&);
  LazyMatrixSum_inplace_addsub_op(const stored_matrix_type&);
  LazyMatrixSum_inplace_addsub_op(LazyMatrixSum);

 private:
  //! The cached number of rows
  size_type m_n_rows;

  //! The cached number of columns
  size_type m_n_cols;

  //! Collection of all the terms which are subject to delayed evaluation
  std::vector<lazy_term_type> m_lazy_terms;

  /** Collection of stored matrix terms to which we reference */
  std::vector<stored_term_type> m_stored_terms;
};

//
// Operators which perform scaling
//
/** \brief Scale a lazy matrix sum */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator*(LazyMatrixSum<StoredMatrix> lhs,
                                      typename StoredMatrix::scalar_type rhs) {
  lhs *= rhs;
  return lhs;
}

/** \brief Scale a lazy matrix sum */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator*(typename StoredMatrix::scalar_type lhs,
                                      LazyMatrixSum<StoredMatrix> rhs) {
  return rhs * lhs;
}

/** \brief Scale a lazy matrix sum */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator/(LazyMatrixSum<StoredMatrix> lhs,
                                      typename StoredMatrix::scalar_type rhs) {
  lhs /= rhs;
  return lhs;
}

//
// Unary operator-
//
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator-(LazyMatrixSum<StoredMatrix> mat) {
  typedef typename StoredMatrix::scalar_type scalar_type;
  return -Constants<scalar_type>::one * mat;
}

//
// Macro definition for adding and subtracting things to/from LazyMatrixSums
//
#define LazyMatrixSum_outofplace_addsub_op(othertype)                      \
  template <typename StoredMatrix>                                         \
  LazyMatrixSum<StoredMatrix> operator+(LazyMatrixSum<StoredMatrix> sum,   \
                                        othertype other) {                 \
    sum += other;                                                          \
    return sum;                                                            \
  }                                                                        \
  template <typename StoredMatrix>                                         \
  LazyMatrixSum<StoredMatrix> operator-(LazyMatrixSum<StoredMatrix> sum,   \
                                        othertype other) {                 \
    sum -= other;                                                          \
    return sum;                                                            \
  }                                                                        \
  template <typename StoredMatrix>                                         \
  LazyMatrixSum<StoredMatrix> operator+(othertype other,                   \
                                        LazyMatrixSum<StoredMatrix> sum) { \
    return sum + other;                                                    \
  }                                                                        \
  template <typename StoredMatrix>                                         \
  LazyMatrixSum<StoredMatrix> operator-(othertype other,                   \
                                        LazyMatrixSum<StoredMatrix> sum) { \
    return (-sum) + other;                                                 \
  }

//
// Macro definition for adding and subtracting things to/from LazyMatrixProducts
//
#define LazyMatrixProduct_outofplace_addsub_op(othertype)                       \
  template <typename StoredMatrix>                                              \
  LazyMatrixSum<StoredMatrix> operator+(LazyMatrixProduct<StoredMatrix> prod,   \
                                        othertype other) {                      \
    LazyMatrixSum<StoredMatrix> sum(std::move(prod));                           \
    return sum + other;                                                         \
  }                                                                             \
  template <typename StoredMatrix>                                              \
  LazyMatrixSum<StoredMatrix> operator-(LazyMatrixProduct<StoredMatrix> prod,   \
                                        othertype other) {                      \
    LazyMatrixSum<StoredMatrix> sum(std::move(prod));                           \
    return sum - other;                                                         \
  }                                                                             \
  template <typename StoredMatrix>                                              \
  LazyMatrixSum<StoredMatrix> operator+(othertype other,                        \
                                        LazyMatrixProduct<StoredMatrix> prod) { \
    return prod + other;                                                        \
  }                                                                             \
  template <typename StoredMatrix>                                              \
  LazyMatrixSum<StoredMatrix> operator-(othertype other,                        \
                                        LazyMatrixProduct<StoredMatrix> prod) { \
    return (-prod) + other;                                                     \
  }

//
// Further addition operators
//
/* clang-format off */
LazyMatrixSum_outofplace_addsub_op(LazyMatrixProduct<StoredMatrix>)
LazyMatrixSum_outofplace_addsub_op(const StoredMatrix&)
LazyMatrixSum_outofplace_addsub_op(const LazyMatrixExpression<StoredMatrix>&)

LazyMatrixProduct_outofplace_addsub_op(const StoredMatrix&)
LazyMatrixProduct_outofplace_addsub_op(
      const LazyMatrixExpression<StoredMatrix>&)

template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator+(LazyMatrixSum<StoredMatrix> lhs,
                                      LazyMatrixSum<StoredMatrix> rhs) {
  /* clang-format on */
  lhs += rhs;
  return lhs;
}

template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator-(LazyMatrixSum<StoredMatrix> lhs,
                                      LazyMatrixSum<StoredMatrix> rhs) {
  lhs -= rhs;
  return lhs;
}

template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator+(LazyMatrixProduct<StoredMatrix> lhs,
                                      LazyMatrixProduct<StoredMatrix> rhs) {
  LazyMatrixSum<StoredMatrix> sum(std::move(lhs));
  sum += rhs;
  return sum;
}

template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator-(LazyMatrixProduct<StoredMatrix> lhs,
                                      LazyMatrixProduct<StoredMatrix> rhs) {
  LazyMatrixSum<StoredMatrix> sum(std::move(lhs));
  sum -= rhs;
  return sum;
}

//
// ------------------------------------------------------------
//
//
// LazyMatrixSum
//
template <typename StoredMatrix>
typename LazyMatrixSum<StoredMatrix>::scalar_type LazyMatrixSum<StoredMatrix>::operator()(
      size_type row, size_type col) const {
  if (empty()) return Constants<scalar_type>::zero;

  assert_range(0u, row, n_rows());
  assert_range(0u, col, n_cols());
  stored_matrix_type block(1, 1, false);
  extract_block(block, row, col);
  return block(0, 0);
}

template <typename StoredMatrix>
bool LazyMatrixSum<StoredMatrix>::has_transpose_operation_mode() const {
  // If all have it, than we have it.
  const bool lazy = std::all_of(
        std::begin(m_lazy_terms), std::end(m_lazy_terms),
        [](const lazy_term_type& t) { return t.has_transpose_operation_mode(); });
  const bool stored = std::all_of(std::begin(m_stored_terms), std::end(m_stored_terms),
                                  [](const stored_term_type& t) {
                                    return t.matrix().has_transpose_operation_mode();
                                  });
  return lazy && stored;
}

template <typename StoredMatrix>
void LazyMatrixSum<StoredMatrix>::extract_block(
      stored_matrix_type& M, const size_type start_row, const size_type start_col,
      const Transposed mode, const scalar_type c_this, const scalar_type c_M) const {
  assert_dbg(mode == Transposed::None || has_transpose_operation_mode(),
             ExcUnsupportedOperationMode(mode));
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

  if (c_this == Constants<scalar_type>::zero || empty()) {
    detail::scale_or_set(M, c_M);
    return;
  }  // c_this == 0

  // We need to use the c_M of the input only once,
  // so store a local copy and set this to 1 after first run
  scalar_type our_cM = c_M;

  // Extract all the terms in turn and
  // accumulate results in M matrix
  for (const auto& stored_term : m_stored_terms) {
    const scalar_type coeff = stored_term.coefficient();
    const stored_matrix_type& mat = stored_term.matrix();
    mat.extract_block(M, start_row, start_col, mode, c_this * coeff, our_cM);

    // All future multiplication results will just be
    // accumulated in out
    our_cM = Constants<scalar_type>::one;
  }

  for (const auto& lazy_term : m_lazy_terms) {
    lazy_term.extract_block(M, start_row, start_col, mode, c_this, our_cM);
    our_cM = Constants<scalar_type>::one;
  }
}

template <typename StoredMatrix>
void LazyMatrixSum<StoredMatrix>::apply(
      const MultiVector<const MutableMemoryVector_i<scalar_type>>& x,
      MultiVector<MutableMemoryVector_i<scalar_type>>& y, const Transposed mode,
      const scalar_type c_this, const scalar_type c_y) const {
  assert_dbg(mode == Transposed::None || has_transpose_operation_mode(),
             ExcUnsupportedOperationMode(mode));
  assert_finite(c_this);
  assert_finite(c_y);
  assert_dbg(!empty(), krims::ExcInvalidState("LazyMatrixSum is empty."));
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

  // Local, modifiable copy of c_y
  scalar_type our_cy = c_y;

  // Apply the terms to the input and
  // accumulate results in y multivector
  for (const auto& stored_term : m_stored_terms) {
    const scalar_type coeff = stored_term.coefficient();
    const stored_matrix_type& mat = stored_term.matrix();
    mat.apply(x, y, mode, c_this * coeff, our_cy);

    // All future multiplication results will just be
    // accumulated in y
    our_cy = Constants<scalar_type>::one;
  }

  for (const auto& lazy_term : m_lazy_terms) {
    lazy_term.apply(x, y, mode, c_this, our_cy);
    our_cy = Constants<scalar_type>::one;
  }
}

template <typename StoredMatrix>
void LazyMatrixSum<StoredMatrix>::mmult(const stored_matrix_type& in,
                                        stored_matrix_type& out, const Transposed mode,
                                        const scalar_type c_this,
                                        const scalar_type c_out) const {
  assert_dbg(mode == Transposed::None || has_transpose_operation_mode(),
             ExcUnsupportedOperationMode(mode));
  assert_finite(c_this);
  assert_finite(c_out);
  assert_dbg(!empty(), krims::ExcInvalidState("LazyMatrixSum is empty."));
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
  }

  // Local, modifiable copy of c_out
  scalar_type our_cout = c_out;

  // Multiply the terms with the input and
  // accumulate results in out matrix
  for (const auto& stored_term : m_stored_terms) {
    const scalar_type coeff = stored_term.coefficient();
    const stored_matrix_type& mat = stored_term.matrix();
    mat.mmult(in, out, mode, c_this * coeff, our_cout);

    // All future multiplication results will just be
    // accumulated in out
    our_cout = Constants<scalar_type>::one;
  }

  for (const auto& lazy_term : m_lazy_terms) {
    lazy_term.mmult(in, out, mode, c_this, our_cout);
    our_cout = Constants<scalar_type>::one;
  }
}

}  // namespace linalgwrap
