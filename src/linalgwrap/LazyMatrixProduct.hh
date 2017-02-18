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
#include "detail/scale_or_set.hh"
#include "linalgwrap/Constants.hh"
#include "linalgwrap/LazyMatrixExpression.hh"
#include "linalgwrap/MultiVector.hh"
#include <algorithm>
#include <iterator>
#include <krims/GenMap.hh>
#include <memory>
#include <vector>

namespace linalgwrap {

// Forward declaration
template <typename StoredMatrix>
class LazyMatrixExpression;

/**
 * \brief Class to represent a product of matrices, which are scaled
 * by a common coefficient.
 */
template <typename StoredMatrix>
class LazyMatrixProduct : public LazyMatrixExpression<StoredMatrix> {
 public:
  typedef LazyMatrixExpression<StoredMatrix> base_type;
  typedef typename base_type::stored_matrix_type stored_matrix_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::vector_type vector_type;
  typedef typename base_type::lazy_matrix_expression_ptr_type
        lazy_matrix_expression_ptr_type;

  /** A swap function for LazyMatrixProducts */
  friend void swap(LazyMatrixProduct& first, LazyMatrixProduct& second) {
    using std::swap;  // enable ADL

    swap(first.m_coefficient, second.m_coefficient);
    swap(first.m_factors, second.m_factors);
    swap(static_cast<base_type&>(first), static_cast<base_type&>(second));
  }

  /** \name Constructors, desctructors and assignment
   */
  ///@{
  /** \brief Create a matrix product object:
   *
   * @param expr   The first matrix expression factor
   * @param factor The scalar factor to initialise the object with
   */
  explicit LazyMatrixProduct(const LazyMatrixExpression<StoredMatrix>& expr,
                             scalar_type factor = Constants<scalar_type>::one)
        : m_coefficient(factor) {

    // Push back the first factor:
    m_factors.push_back(std::move(expr.clone()));
  }

  /** \brief Copy and scale constructor constructor
   *
   * Copy the contents of another product and scale it altogether
   * by an extra factor.
   *
   *  @param factor   The factor to scale the product with
   * */
  explicit LazyMatrixProduct(LazyMatrixProduct prod, scalar_type factor) {
    swap(*this, prod);
    m_coefficient *= factor;
  }

  //@{
  /** Defaults for the big five */
  ~LazyMatrixProduct() = default;
  LazyMatrixProduct(const LazyMatrixProduct&) = default;
  LazyMatrixProduct(LazyMatrixProduct&&) = default;
  LazyMatrixProduct& operator=(const LazyMatrixProduct& other) = default;
  LazyMatrixProduct& operator=(LazyMatrixProduct&& other) = default;
  //@}
  ///@}

  //
  // Push back further factors
  //
  /** \brief Add a factor by copying or moving the matrix \p m into
   * the current object.
   */
  void push_factor(const LazyMatrixExpression<StoredMatrix>& e) {
    // Check fitting dimensionality of matrix expressions:
    assert_size(n_cols(), e.n_rows());

    // Place into m_factors by moving a copy there
    m_factors.push_back(std::move(e.clone()));
  }

  /** \brief Push back all factors of a product onto another product
   */
  void push_factor(LazyMatrixProduct prod) {
    // Check fitting dimensionality of matrix expressions:
    assert_size(n_cols(), prod.n_rows());

    // Move all factors of m to the very end:
    std::move(std::begin(prod.m_factors), std::end(prod.m_factors),
              back_inserter(m_factors));

    // Adjust the scaling:
    m_coefficient *= prod.m_coefficient;
  }

  //
  // Other modifications of the product
  //
  /** \brief Scale the whole matrix_factor by the scalar \p c
   */
  void scale(const scalar_type c) {
    assert_finite(c);
    m_coefficient *= c;
  }

  //
  // Matrix_i interface
  //
  /** \brief Number of rows of the matrix */
  size_type n_rows() const override {
    // The number of rows of the product
    // equals the number of rows of the first element in the product.
    return m_factors.front()->n_rows();
  }

  /** \brief Number of columns of the matrix  */
  size_type n_cols() const override {
    // The number of columns of the product
    // equals the nuber of columns of the last element in the product.
    return m_factors.back()->n_cols();
  }

  /** \brief return an element of the matrix    */
  scalar_type operator()(size_type row, size_type col) const override;

  //
  // LazyMatrixExpression interface
  //

  /** Are operation modes Transposed::Trans and Transposed::ConjTrans
   *  supported for this matrix type.
   **/
  bool has_transpose_operation_mode() const override {
    // If all have it, than we have it.
    return std::all_of(
          std::begin(m_factors), std::end(m_factors),
          [](const factor_ptr_type& p) { return p->has_transpose_operation_mode(); });
  }

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
            mat_vec_apply_enabled_t<LazyMatrixProduct, VectorIn, VectorOut>...>
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
    for (auto& expression : m_factors) {
      expression->update(map);
    }
  }

  /** \brief Clone the expression */
  lazy_matrix_expression_ptr_type clone() const override {
    // return a copy enwrapped in the pointer type
    return lazy_matrix_expression_ptr_type(new LazyMatrixProduct(*this));
  }

  //
  // In-place scalar operators:
  //
  /** \brief scale a matrix using a scalar
   */
  LazyMatrixProduct& operator*=(scalar_type c) {
    assert_finite(c);
    this->scale(c);
    return *this;
  }

  LazyMatrixProduct& operator/=(scalar_type c) {
    assert_finite(c);
    assert_dbg(c != 0, krims::ExcDevideByZero());
    this->scale(Constants<scalar_type>::one / c);
    return *this;
  }

 private:
  /** \name Inner functions level 2
   * Where the assertions have been checked, trivial cases have been dealt
   * with
   * and we just need to run over the factors and recursively deal with the
   * problem
   * at hand (for normal mode we usually parse m_factors in reverse order,
   * else in normal order and with the mode passed on
   */
  ///@{

  /** Run over iterator range from begin to end and extract the appropriate
   * block of the product.
   */
  template <typename BidirectIterator>
  void extract_block_inner(BidirectIterator begin, BidirectIterator end,
                           stored_matrix_type& M, const size_type start_row,
                           const size_type start_col, const Transposed mode,
                           const scalar_type c_this, const scalar_type c_M) const;

  /** Run over iterator range from begin to end and multiply the
   * encountered factors recursively with in, i.e. forms the product
   *
   * out = c_out * out + c_this * fac1^{mode} * fac2^{mode} * fac3^{mode} * in
   *
   * For mode == Transposed == None we usually iterate over m_factors in the
   * reverse order using a reverse iterator(BidirectIterator == reverse
   * iterator),
   * else in the normal order (since A0^T A1^T A2^T = (A2 A1 A0)^T )
   */
  template <typename BidirectIterator>
  void mmult_inner(BidirectIterator begin, BidirectIterator end,
                   const stored_matrix_type& in, stored_matrix_type& out,
                   const Transposed mode, const scalar_type c_this,
                   const scalar_type c_out) const;

  /** Run over iterator range from begin to end and multiply the
   * encountered factors recursively with x, i.e. forms the product
   *
   * y = c_y * y + c_this * fac1^{mode} * fac2^{mode} * fac3^{mode} * x
   *
   * For mode == Transposed == None we usually iterate over m_factors in the
   * reverse order using a reverse iterator(BidirectIterator == reverse
   * iterator),
   * else in the normal order (since A0^T A1^T A2^T = (A2 A1 A0)^T )
   */
  template <typename BidirectIterator>
  void apply_inner(BidirectIterator begin, BidirectIterator end,
                   const MultiVector<const MutableMemoryVector_i<scalar_type>>& x,
                   MultiVector<MutableMemoryVector_i<scalar_type>>& y,
                   const Transposed mode, const scalar_type c_this,
                   const scalar_type c_y) const;
  ///@}

  /** \name Apply and multiply level 1
   *
   * Just recursively multiply or apply a chain of factors (represented
   * by the iterator range) in-place.
   */
  ///@{
  /** \brief mmult a range of factors given by the iterator range.
   *
   * The function is equivalent to
   * m = (fac1)^{mode} * (fac2)^{mode} * (fac3)^{mode} * ... * m
   */
  template <typename BidirectIterator>
  static void multiply_in_place(BidirectIterator begin, BidirectIterator end,
                                stored_matrix_type& m, Transposed mode);

  /** \brief apply a range of factors given by the iterator range.
   *
   * The function is equivalent to
   * v = (fac1)^{mode} * (fac2)^{mode} * (fac3)^{mode} * ... * v
   */
  template <typename BidirectIterator>
  static void apply_in_place(BidirectIterator begin, BidirectIterator end,
                             MultiVector<vector_type>& mv, Transposed mode);
  ///@}

  //! The type of the lazy matrix expression pointers inside the vector
  typedef std::shared_ptr<LazyMatrixExpression<StoredMatrix>> factor_ptr_type;

  //! The vector containing the lazy matrix expression factors
  std::vector<factor_ptr_type> m_factors;

  //! The global scaling coefficient of all factors
  scalar_type m_coefficient;
};

/** \brief Multiply two lazy matrix products */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator*(LazyMatrixProduct<StoredMatrix> lhs,
                                          LazyMatrixProduct<StoredMatrix> rhs) {
  lhs.push_factor(std::move(rhs));
  return lhs;
}

/** \brief Multiply a lazy matrix product and an expression */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator*(LazyMatrixProduct<StoredMatrix> lhs,
                                          const LazyMatrixExpression<StoredMatrix>& rhs) {
  lhs.push_factor(rhs);
  return lhs;
}

/** \brief Multiply a lazy matrix product and an expression */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator*(const LazyMatrixExpression<StoredMatrix>& rhs,
                                          LazyMatrixProduct<StoredMatrix> lhs) {
  return lhs * rhs;
}

//
// Operator*
//
/** \brief Scale a lazy matrix product  */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator*(LazyMatrixProduct<StoredMatrix> m,
                                          typename StoredMatrix::scalar_type s) {
  m *= s;
  return m;
}

/** \brief Scale a lazy matrix product  */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator*(typename StoredMatrix::scalar_type s,
                                          LazyMatrixProduct<StoredMatrix> m) {
  return m * s;
}

/** \brief Scale a lazy matrix product  */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator/(LazyMatrixProduct<StoredMatrix> lhs,
                                          typename StoredMatrix::scalar_type rhs) {
  lhs /= rhs;
  return lhs;
}

//
// Unary operator-
//
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator-(LazyMatrixProduct<StoredMatrix> mat) {
  typedef typename StoredMatrix::scalar_type scalar_type;
  return -Constants<scalar_type>::one * mat;
}

//
// ------------------------------------------------------------
//
//
// LazyMatrixProduct
//
template <typename StoredMatrix>
typename LazyMatrixProduct<StoredMatrix>::scalar_type LazyMatrixProduct<StoredMatrix>::
operator()(size_type row, size_type col) const {
  assert_range(0u, row, n_rows());
  assert_range(0u, col, n_cols());
  stored_matrix_type block(1, 1, false);
  extract_block(block, row, col);
  return block(0, 0);
}

template <typename StoredMatrix>
void LazyMatrixProduct<StoredMatrix>::extract_block(
      stored_matrix_type& M, const size_type start_row, const size_type start_col,
      const Transposed mode, const scalar_type c_this, const scalar_type c_M) const {
  // For empty matrices there is nothing to do
  if (M.n_rows() == 0 || M.n_cols() == 0) return;

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

  if (c_this == Constants<scalar_type>::zero) {
    detail::scale_or_set(M, c_M);
    return;
  }  // c_this == 0

  // TODO Here we assume that m_factors.size() > 0
  assert_dbg(m_factors.size() > 0, krims::ExcInternalError());

  if (m_factors.size() == 1) {
    // Forward to the one factor
    m_factors.front()->extract_block(M, start_row, start_col, mode,
                                     c_this * m_coefficient, c_M);
    return;
  }

  if (mode == Transposed::None) {
    // Go about factors in reverse order (i.e. extract the required
    // block on the right factor and multiply recursively until
    // we reach the first factor)
    extract_block_inner(m_factors.rbegin(), m_factors.rend(), M, start_row, start_col,
                        mode, c_this, c_M);
  } else {
    // Here we exploit that
    // (A0 A1 A2 ... An)^T M = An^T ... A2^T A1^T A0^T M
    extract_block_inner(m_factors.begin(), m_factors.end(), M, start_row, start_col, mode,
                        c_this, c_M);
  }
}

template <typename StoredMatrix>
template <typename BidirectIterator>
void LazyMatrixProduct<StoredMatrix>::extract_block_inner(
      BidirectIterator begin, BidirectIterator end, stored_matrix_type& M,
      const size_type start_row, const size_type start_col, const Transposed mode,
      const scalar_type c_this, const scalar_type c_M) const {
  assert_dbg(end != begin, krims::ExcInternalError());
  assert_dbg((end - 1) != begin, krims::ExcInternalError());
  assert_dbg(end != (begin + 1), krims::ExcInternalError());

  // Do the first factor. On it we call extract_block
  // with the full number of rows
  const size_type rows =
        mode == Transposed::None ? (*begin)->n_rows() : (*begin)->n_cols();
  stored_matrix_type tmp(rows, M.n_cols(), false);
  (*begin)->extract_block(tmp, 0, start_col, mode);

  // In between just multiply
  multiply_in_place((begin + 1), (end - 1), tmp, mode);

  // Extract required values of the last factor
  auto last = (end - 1);
  const size_type cols = mode == Transposed::None ? (*last)->n_cols() : (*last)->n_rows();
  stored_matrix_type last_vals(M.n_rows(), cols, false);
  (*last)->extract_block(last_vals, start_row, 0, mode);

  assert_dbg(last_vals.n_rows() == M.n_rows(), krims::ExcInternalError());
  assert_dbg(last_vals.n_cols() == tmp.n_rows(), krims::ExcInternalError());
  assert_dbg(tmp.n_cols() == M.n_cols(), krims::ExcInternalError());

  // Perform final multiplication
  last_vals.mmult(tmp, M, Transposed::None, c_this * m_coefficient, c_M);

  // If extracting the required part of the first factor is much
  // more expansive that doing another lazy*stored multiplication
  // and binning the rows we do not require, than perhaps it is
  // in fact better to let the for loop run over all factors but
  // the last and use another extract_block on the cache to
  // extract the appropriate row elements.
  // (i.e. reverse the order of the last 2 steps)
}

template <typename StoredMatrix>
void LazyMatrixProduct<StoredMatrix>::apply(
      const MultiVector<const MutableMemoryVector_i<scalar_type>>& x,
      MultiVector<MutableMemoryVector_i<scalar_type>>& y, const Transposed mode,
      const scalar_type c_this, const scalar_type c_y) const {
  assert_dbg(mode == Transposed::None || has_transpose_operation_mode(),
             ExcUnsupportedOperationMode(mode));
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

  // TODO Here we assume that m_factors.size() > 0
  assert_dbg(m_factors.size() > 0, krims::ExcInternalError());

  if (m_factors.size() == 1) {
    // Forward to the one factor:
    m_factors.front()->apply(x, y, mode, c_this * m_coefficient, c_y);
    return;
  }

  if (mode == Transposed::None) {
    // Go about factors in reverse order (i.e. applying from right to left
    // to the supplied input recursively)
    apply_inner(m_factors.rbegin(), m_factors.rend(), x, y, mode, c_this, c_y);
  } else {
    // Here we exploit that
    // (A0 A1 A2 ... An)^T M = An^T ... A2^T A1^T A0^T M
    apply_inner(m_factors.begin(), m_factors.end(), x, y, mode, c_this, c_y);
  }
}

template <typename StoredMatrix>
template <typename BidirectIterator>
void LazyMatrixProduct<StoredMatrix>::apply_inner(
      BidirectIterator begin, BidirectIterator end,
      const MultiVector<const MutableMemoryVector_i<scalar_type>>& x,
      MultiVector<MutableMemoryVector_i<scalar_type>>& y, const Transposed mode,
      const scalar_type c_this, const scalar_type c_y) const {
  assert_dbg(end != begin, krims::ExcInternalError());
  assert_dbg((end - 1) != begin, krims::ExcInternalError());
  assert_dbg(end != (begin + 1), krims::ExcInternalError());

  // Deal with first factor:
  const size_type rows =
        mode == Transposed::None ? (*begin)->n_rows() : (*begin)->n_cols();
  MultiVector<vector_type> tmp(rows, x.n_vectors(), false);
  (*begin)->apply(x, tmp, mode);

  // Deal with stuff in the middle
  // (end-1) points to last but one
  apply_in_place((begin + 1), (end - 1), tmp, mode);

  // Deal with last factor:
  const auto last = end - 1;
  if (mode == Transposed::None) {
    assert_dbg(tmp.n_elem() == (*last)->n_cols(), krims::ExcInternalError());
  } else {
    assert_dbg(tmp.n_elem() == (*last)->n_rows(), krims::ExcInternalError());
  }
  (*last)->apply(tmp, y, mode, c_this * m_coefficient, c_y);
}

template <typename StoredMatrix>
void LazyMatrixProduct<StoredMatrix>::mmult(const stored_matrix_type& in,
                                            stored_matrix_type& out,
                                            const Transposed mode,
                                            const scalar_type c_this,
                                            const scalar_type c_out) const {
  assert_dbg(mode == Transposed::None || has_transpose_operation_mode(),
             ExcUnsupportedOperationMode(mode));
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
  }

  // TODO Here we assume that m_factors.size() > 0
  assert_dbg(m_factors.size() > 0, krims::ExcInternalError());

  if (m_factors.size() == 1) {
    // Forward to the one factor:
    m_factors.front()->mmult(in, out, mode, c_this * m_coefficient, c_out);
    return;
  }

  if (mode == Transposed::None) {
    // Go about factors in reverse order (i.e. applying from right to left
    // to the supplied input recursively)
    mmult_inner(m_factors.rbegin(), m_factors.rend(), in, out, mode, c_this, c_out);
  } else {
    // Here we exploit that
    // (A0 A1 A2 ... An)^T M = An^T ... A2^T A1^T A0^T M
    mmult_inner(m_factors.begin(), m_factors.end(), in, out, mode, c_this, c_out);
  }
}

template <typename StoredMatrix>
template <typename BidirectIterator>
void LazyMatrixProduct<StoredMatrix>::mmult_inner(
      BidirectIterator begin, BidirectIterator end, const stored_matrix_type& in,
      stored_matrix_type& out, const Transposed mode, const scalar_type c_this,
      const scalar_type c_out) const {
  assert_dbg(end != begin, krims::ExcInternalError());
  assert_dbg((end - 1) != begin, krims::ExcInternalError());
  assert_dbg(end != (begin + 1), krims::ExcInternalError());

  // Deal with first factor:
  const size_type rows =
        mode == Transposed::None ? (*begin)->n_rows() : (*begin)->n_cols();
  stored_matrix_type tmp(rows, in.n_cols(), false);
  (*begin)->mmult(in, tmp, mode);

  // Deal with stuff in the middle
  // (end-1) points to last but one
  multiply_in_place((begin + 1), (end - 1), tmp, mode);

  // Deal with the last factor
  const auto last = end - 1;
  if (mode == Transposed::None) {
    assert_dbg(tmp.n_rows() == (*last)->n_cols(), krims::ExcInternalError());
  } else {
    assert_dbg(tmp.n_rows() == (*last)->n_rows(), krims::ExcInternalError());
  }
  (*last)->mmult(tmp, out, mode, c_this * m_coefficient, c_out);
}

template <typename StoredMatrix>
template <typename BidirectIterator>
void LazyMatrixProduct<StoredMatrix>::multiply_in_place(BidirectIterator begin,
                                                        BidirectIterator end,
                                                        stored_matrix_type& m,
                                                        Transposed mode) {
  for (; begin != end; ++begin) {
    const size_type rows =
          mode == Transposed::None ? (*begin)->n_rows() : (*begin)->n_cols();
    stored_matrix_type tmp(rows, m.n_cols(), false);
    (*begin)->mmult(m, tmp, mode);
    m = std::move(tmp);
  }

  /*
      if (begin == end) return;
      do {
          --end;
          const size_type rows =
                mode == Transposed::None ? (*end)->n_rows() :
     (*end)->n_cols();
          stored_matrix_type tmp(rows, m.n_cols(), false);
          (*end)->mmult(m, tmp, mode);  // TODO <----------------- Error here
          m = std::move(tmp);           // Move tmp into m and free m's
     storage
      } while (end != begin);
      */
}

template <typename StoredMatrix>
template <typename BidirectIterator>
void LazyMatrixProduct<StoredMatrix>::apply_in_place(BidirectIterator begin,
                                                     BidirectIterator end,
                                                     MultiVector<vector_type>& mv,
                                                     Transposed mode) {
  for (; begin != end; ++begin) {
    const size_type rows =
          mode == Transposed::None ? (*begin)->n_rows() : (*begin)->n_cols();
    MultiVector<vector_type> tmp(rows, mv.n_vectors(), false);
    (*begin)->apply(mv, tmp, mode);
    // Move tmp into mv and free mv's storage
    mv = std::move(tmp);
  }
}

}  // namespace linalgwrap
