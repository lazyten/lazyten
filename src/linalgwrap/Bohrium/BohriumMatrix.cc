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

#include "BohriumMatrix.hh"
#ifdef LINALGWRAP_HAVE_BOHRIUM

namespace linalgwrap {

//
// Constructors
//
template <typename Scalar>
BohriumMatrix<Scalar>::BohriumMatrix(const BohriumMatrix& mat, scalar_type tolerance)
      : BohriumMatrix(mat.n_rows(), mat.n_cols(), false) {
  // TODO something better?
  assert_implemented(false);

  /*
  for (const auto elem : mat) {
    if (std::fabs(*elem) < tolerance) {
      (*this)(elem.row(), elem.col()) = 0;
    } else {
      (*this)(elem.row(), elem.col()) = *elem;
    }
  }
  */
}

template <typename Scalar>
BohriumMatrix<Scalar>::BohriumMatrix(
      std::initializer_list<std::initializer_list<scalar_type>> list_of_lists)
      : BohriumMatrix(list_of_lists.size(),
                      list_of_lists.size() > 0 ? list_of_lists.begin()->size() : 0,
                      false) {
#ifdef DEBUG
  size_type n_rows = list_of_lists.size();
  size_type n_cols = n_rows > 0 ? list_of_lists.begin()->size() : 0;

  // Assert all columns have equal length.
  assert_element_sizes(list_of_lists, n_cols);
#endif

  size_type i = 0;
  Scalar* mem = memptr();
  for (auto row : list_of_lists) {
    std::copy(row.begin(), row.end(), mem + i);
    i += row.size();
  }
  assert_internal(i == n_rows * n_cols);
}

template <typename Scalar>
bool BohriumMatrix<Scalar>::operator==(const BohriumMatrix& other) const {
  if (n_rows() != other.n_rows()) return false;
  if (n_cols() != other.n_cols()) return false;

  // Use the inner_product function with equal as the elementwise operation
  // and logical_and_reduce as the accumulate operation
  return bhxx::as_scalar(bhxx::inner_product(
        m_array, other.m_array, bhxx::Equal<Scalar>{}, bhxx::LogicalAndReduce<bool>{}));
}

template <typename Scalar>
bool BohriumMatrix<Scalar>::operator!=(const BohriumMatrix& other) const {
  if (n_rows() != other.n_rows()) return true;
  if (n_cols() != other.n_cols()) return true;

  // Use the inner_product function with not_equal as the elementwise operation
  // and logical_or_reduce as the accumulate operation
  return bhxx::as_scalar(bhxx::inner_product(
        m_array, other.m_array, bhxx::NotEqual<Scalar>{}, bhxx::LogicalOrReduce<bool>{}));
}

//
// Application and multiplication
//

template <typename Scalar>
void BohriumMatrix<Scalar>::apply(const bhxx::BhArray<Scalar>& A,
                                  const bhxx::BhArray<Scalar>& x,
                                  bhxx::BhArray<Scalar>& y, const scalar_type c_this,
                                  const scalar_type c_y) const {
  bhxx::BhArray<Scalar> Ax = bhxx::matmul(A, x);
  if (c_this != 0) bhxx::multiply(Ax, Ax, c_this);

  if (c_y == 0) {
    y = std::move(Ax);
  } else {
    bhxx::BhArray<Scalar> scaled(y.shape);
    bhxx::multiply(scaled, y, c_y);
    bhxx::add(y, Ax, scaled);
  }
}

template <typename Scalar>
void BohriumMatrix<Scalar>::mmult(const BohriumMatrix& in, BohriumMatrix& out,
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
  assert_implemented(mode != Transposed::ConjTrans);

  if (c_this == 0) {
    detail::scale_or_set(out, c_out);
    return;
  }

  bhxx::BhArray<Scalar> A = bh_array();
  bhxx::BhArray<Scalar> B = in.bh_array();
  auto& C = out.bh_array();

  if (mode == Transposed::Trans) A = bhxx::transpose(std::move(A));
  if (c_this != 1) bhxx::multiply(A, A, c_this);

  bhxx::BhArray<Scalar> AB = bhxx::matmul(std::move(A), std::move(B));
  assert_internal(AB.rank() == 2);

  if (c_out != 0) {
    bhxx::multiply(C, C, c_out);
    assert_internal(C.shape == AB.shape);
    bhxx::add(AB, C, AB);
  }
  C = std::move(AB);
}

template <typename Scalar>
void BohriumMatrix<Scalar>::extract_block(BohriumMatrix<Scalar>& M,
                                          const size_type start_row,
                                          const size_type start_col,
                                          const Transposed mode, const scalar_type c_this,
                                          const scalar_type c_M) const {
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
  assert_implemented(mode != Transposed::ConjTrans);

  // For empty matrices there is nothing to do
  if (M.n_rows() == 0 || M.n_cols() == 0) return;

  if (c_this == 0) {
    detail::scale_or_set(M, c_M);
    return;
  }  // c_this == 0

  assert_internal(m_array.is_contiguous());

  // Construct the view of m_array, which we wish to extract
  size_t offset = start_row * n_cols() + start_col;
  bhxx::Stride view_stride = m_array.stride;
  if (mode == Transposed::Trans) {
    // For transposed the strides need to be swapped and the offset
    // adjusted keeping in mind that start_row and start_col now refer
    // to the starting row/column in the transposed object.
    offset = start_col * n_cols() + start_row;
    view_stride = {m_array.stride[1], m_array.stride[0]};
  }
  bhxx::BhArray<Scalar> Aview(m_array.base, {M.n_rows(), M.n_cols()}, view_stride,
                              offset);

  bhxx::BhArray<Scalar> Ascaled(Aview.shape);
  bhxx::multiply(Ascaled, Aview, c_this);
  if (c_M == 0) {
    bhxx::identity(M.bh_array(), Ascaled);
  } else {
    assert_internal(M.bh_array().is_data_initialised());
    bhxx::BhArray<Scalar> Mscaled(M.bh_array().shape);
    bhxx::multiply(Mscaled, M.bh_array(), c_M);
    bhxx::add(M.bh_array(), Mscaled, Ascaled);
  }
}

//
// Explicit instantiation
//
#define INSTANTIATE(TYPE) template class BohriumMatrix<TYPE>
INSTANTIATE(double);
// INSTANTIATE(float);
// INSTANTIATE(std::complex<float>);
// INSTANTIATE(std::complex<double>);
#undef INSTANTIATE

}  // namespace linalgwrap

#endif  // LINALGWRAP_HAVE_BOHRIUM
