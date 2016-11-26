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
#include "linalgwrap/TypeUtils.hh"
#include <krims/RCPWrapper.hh>

namespace linalgwrap {

namespace detail {

/** ProxyBase class for stored matrices */
template <typename Matrix, typename StoredMatrix = typename StoredTypeOf<Matrix>::type,
          bool isStoredMatrix = IsStoredMatrix<Matrix>::value,
          bool isReadonly = std::is_const<Matrix>::value>
class ProxyBase : public LazyMatrixExpression<StoredMatrix> {
  static_assert(IsMatrix<Matrix>::value, "Matrix needs to be a matrix type");

 public:
  //! Matrix type underlying this TransposeProxyBase class
  typedef Matrix matrix_type;
  typedef StoredMatrix stored_matrix_type;

  /** \brief Update the internal data
   *
   *  In this case does nothing.
   * */
  void update(const krims::ParameterMap&) override {}

  /** Does this class own the inner matrix
   *
   * In other words are we the one who is responsible for
   * managing the inner storage or not.V
   */
  bool owns_inner_matrix() const { return m_inner_ptr.is_shared_ptr(); }

  /** Access to the inner object */
  Matrix& inner_matrix() { return *m_inner_ptr; }

  /** Const access to the inner object */
  const Matrix& inner_matrix() const { return *m_inner_ptr; }

 protected:
  ProxyBase(Matrix& matrix)
        : m_inner_ptr(krims::make_subscription(matrix, "ProxyBase")) {}

  ProxyBase(Matrix&& matrix) : m_inner_ptr(std::make_shared<Matrix>(matrix)) {}

  krims::RCPWrapper<Matrix> m_inner_ptr;
};

/** ProxyBase class for const lazy matrix classes */
template <typename Matrix, typename StoredMatrix>
class ProxyBase<Matrix, StoredMatrix, /*isStored*/ false,
                /*isReadonly*/ true> : public LazyMatrixExpression<StoredMatrix> {
  static_assert(IsMatrix<Matrix>::value, "Matrix needs to be a matrix type");

 public:
  //! Matrix type underlying this ProxyBase class
  typedef Matrix matrix_type;
  typedef StoredMatrix stored_matrix_type;

  /** \brief Update the internal data
   *
   *  In this case does nothing, since the internal object is
   *  either a stored matrix or const.
   * */
  void update(const krims::ParameterMap&) override {
    assert_dbg(false, krims::ExcDisabled("Update is not possible for this Matrix "
                                         "expression, since the matrix inside the "
                                         "Proxy object is const."));
  }

  /** Does this class own the inner matrix
   *
   * In other words are we the one who is responsible for
   * managing the inner storage or not.V
   *
   * In this case it always owns it since it is a lazy matrix object.
   */
  constexpr bool owns_inner_matrix() const { return true; }

  /** Access to the inner object */
  Matrix& inner_matrix() { return m_inner; }

  /** Const access to the inner object */
  const Matrix& inner_matrix() const { return m_inner; }

 protected:
  ProxyBase(const Matrix& matrix) : m_inner(matrix) {}
  ProxyBase(Matrix&& matrix) : m_inner(std::move(matrix)) {}

  /** Matrix is lazy, so we can store a full copy */
  Matrix m_inner;
};

/** TransposeProxyBase class for mutable lazy matrix classes */
template <typename Matrix, typename StoredMatrix>
class ProxyBase<Matrix, StoredMatrix, /*isStored*/ false,
                /*isReadonly*/ false> : public LazyMatrixExpression<StoredMatrix> {
  static_assert(IsMatrix<Matrix>::value, "Matrix needs to be a matrix type");

 public:
  //! Matrix type underlying this TransposeProxyBase class
  typedef Matrix matrix_type;
  typedef StoredMatrix stored_matrix_type;

  /** \brief Update the internal data
   *
   *  In this case does nothing, since the internal object is
   *  either a stored matrix or const.
   * */
  void update(const krims::ParameterMap& map) override { m_inner.update(map); }

  /** Does this class own the inner matrix
   *
   * In other words are we the one who is responsible for
   * managing the inner storage or not.V
   *
   * In this case it always owns it since it is a lazy matrix object.
   */
  constexpr bool owns_inner_matrix() const { return true; }

  /** Access to the inner object */
  Matrix& inner_matrix() { return m_inner; }

  /** Const access to the inner object */
  const Matrix& inner_matrix() const { return m_inner; }

 protected:
  ProxyBase(const Matrix& matrix) : m_inner(matrix) {}
  ProxyBase(Matrix&& matrix) : m_inner(std::move(matrix)) {}

  /** Matrix is lazy, so we can store a full copy */
  Matrix m_inner;
};
}  // namespace detail
}  // namespace linalgwrap
