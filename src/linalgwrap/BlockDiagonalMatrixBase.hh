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
#include "detail/GenericFunctionals.hh"
#include "linalgwrap/BlockDiagonalMatrix.hh"
#include "linalgwrap/LazyMatrixExpression.hh"
#include "linalgwrap/StoredMatrix_i.hh"
#include <krims/TupleUtils.hh>

namespace linalgwrap {
// Forward definition of BlockDiagonalMatrix class
template <typename... Matrices>
class BlockDiagonalMatrix;

namespace detail {
/** Definition of type traits for all BlockDiagonalMatrix classes */
template <typename... Matrices>
class BlockDiagonalMatrixTypeTraits;
}  // namespace detail

/** \brief struct representing a type (std::true_type, std::false_type) which
 *  indicates whether T is a block diagonal matrix
 *
 *  \note This works, since by construction all BlockDiagonalMatrix classes
 *  are indirectly derived off BlockDiagonalMatrixTypeTraits<>.
 *  */
template <typename T>
struct IsBlockDiagonalMatrix
      : public std::is_base_of<detail::BlockDiagonalMatrixTypeTraits<>, T> {};

namespace detail {
/** Explicit specialisation of type traits class for a matrix object with
 *  no matrices
 *
 * \note This is by construction the common base class for all block
 * diagonal matrices */
template <>
class BlockDiagonalMatrixTypeTraits<> {};

/** Explicit specialisation of type traits class for a matrix object with
 *  exactly one matrix  */
template <typename Matrix>
class BlockDiagonalMatrixTypeTraits<Matrix>
      : public BlockDiagonalMatrixTypeTraits<> {
  public:
    typedef typename std::remove_reference<Matrix>::type matrix_type;
    typedef typename matrix_type::size_type size_type;
    typedef typename matrix_type::scalar_type scalar_type;

    // Check that the matrices are of a sensible type
    static_assert(IsStoredMatrix<Matrix>::value || IsLazyMatrix<Matrix>::value,
                  "Constituents of BlockDiagonalMatrix objects need to be "
                  "stored or lazy.");
    static_assert(!IsBlockDiagonalMatrix<Matrix>::value,
                  "Constituents of BlockDiagonalMatrix objects cannot be "
                  "BlockDiagonalMatrix objects themself.");

    //! Are all matrices in this diagonal matrix stored?
    static bool constexpr all_stored = IsStoredMatrix<Matrix>::value;
};

/** Explicit specialisation of type traits class for a matrix object with
 * arbitrarily many matrices */
template <typename Matrix, typename... Matrices>
class BlockDiagonalMatrixTypeTraits<Matrix, Matrices...>
      : public BlockDiagonalMatrixTypeTraits<Matrices...> {
  public:
    typedef BlockDiagonalMatrixTypeTraits<Matrices...> base_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;

    // Check that size and scalar types agree.
    static_assert(std::is_same<size_type, typename Matrix::size_type>::value,
                  "All constituent matrices of a block diagonal matrix need to "
                  "have the same size type");
    static_assert(
          std::is_same<scalar_type, typename Matrix::scalar_type>::value,
          "All constituent matrices of a block diagonal matrix need to have "
          "the same scalar type");

    // Check that the matrices are of a sensible type
    static_assert(
          IsMatrix<Matrix>::value,
          "Constituents of BlockDiagonalMatrix objects need to be matrices.");
    static_assert(!IsBlockDiagonalMatrix<Matrix>::value,
                  "Constituents of BlockDiagonalMatrix objects cannot be "
                  "BlockDiagonalMatrix objects themself.");

    //! Are all matrices in this diagonal matrix stored matrices?
    static bool constexpr all_stored =
          base_type::all_stored && IsStoredMatrix<Matrix>::value;
};
}  // namespace detail

/** \brief Class implementing the common interface and functionality of
 *  all BlockDiagonalMatrixBase classes
 *
 * This base class implements all common functionality of the base class
 *
 * \note If all constituent matrices are stored matrices, then this class
 *       derives off StoredMatrix_i,else it just represents a Matrix_i.
 */
// TODO better iterator for this type of matrices.
template <typename... Matrices>
class BlockDiagonalMatrixBase
      : detail::BlockDiagonalMatrixTypeTraits<Matrices...>,
        public Matrix_i<typename detail::BlockDiagonalMatrixTypeTraits<
              Matrices...>::scalar_type> {
    /* TODO Implementation of StoredMatrix_i interface
public std::conditional<
  BlockDiagonalMatrixTypeTraits<Matrices...>::all_stored,
  StoredMatrix_i<typename BlockDiagonalMatrixTypeTraits<
        Matrices...>::scalar_type>,
  Matrix_i<typename BlockDiagonalMatrixTypeTraits<
        Matrices...>::scalar_type>>::type {
        */
  public:
    typedef detail::BlockDiagonalMatrixTypeTraits<Matrices...> traits_type;
    typedef typename traits_type::size_type size_type;
    typedef typename traits_type::scalar_type scalar_type;

    typedef typename std::conditional<traits_type::all_stored,
                                      StoredMatrix_i<scalar_type>,
                                      Matrix_i<scalar_type>>::type base_type;
    static_assert(
          std::is_same<scalar_type, typename base_type::scalar_type>::value,
          "Scalar type of base type and of all constituent matrices does not "
          "agree.");
    static_assert(std::is_same<size_type, typename base_type::size_type>::value,
                  "Size type of base type and of all constituent matrices does "
                  "not agree.");

    /** \brief Const access to the matrix blocks */
    virtual std::tuple<const Matrices&...> blocks() const = 0;
};

/** \name BlockDiagonalMatrix operations */
///@{

//! Add two block diagonal matrices
template <typename... MatricesLHS, typename... MatricesRHS>
auto operator+(const BlockDiagonalMatrixBase<MatricesLHS...>& lhs,
               const BlockDiagonalMatrixBase<MatricesRHS...>& rhs);

//! Subtract two block diagonal matrices
template <typename... MatricesLHS, typename... MatricesRHS>
auto operator-(const BlockDiagonalMatrixBase<MatricesLHS...>& lhs,
               const BlockDiagonalMatrixBase<MatricesRHS...>& rhs);

//! Multiply two block diagonal matrices
template <typename... MatricesLHS, typename... MatricesRHS>
auto operator*(const BlockDiagonalMatrixBase<MatricesLHS...>& lhs,
               const BlockDiagonalMatrixBase<MatricesRHS...>& rhs);

//! Unary minus on a block diagonal matrix
template <typename... Matrices>
auto operator-(const BlockDiagonalMatrixBase<Matrices...>& bm);

//@{
//! Scale a block diagonal matrix
template <typename Matrix, typename... OtherMatrices>
auto operator*(const BlockDiagonalMatrixBase<Matrix, OtherMatrices...>& bm,
               typename Matrix::scalar_type s);

template <typename Matrix, typename... OtherMatrices>
auto operator*(typename Matrix::scalar_type s,
               const BlockDiagonalMatrixBase<Matrix, OtherMatrices...>& bm);

template <typename Matrix, typename... OtherMatrices>
auto operator/(const BlockDiagonalMatrixBase<Matrix, OtherMatrices...>& bm,
               typename Matrix::scalar_type s);
//@}
///@}

//
// --------------------------------------------------
//

// TODO: problem for C++11 since no explicit return type
template <typename... MatricesLHS, typename... MatricesRHS>
auto operator+(const BlockDiagonalMatrixBase<MatricesLHS...>& lhs,
               const BlockDiagonalMatrixBase<MatricesRHS...>& rhs) {
    // Generic lambdas are not part of c++11, so we cannot use them
    // auto plus = [](auto& x, auto& y) { return x + y; };

    // Use generic plus functor instead
    detail::PlusFctr plus;

    // Apply it to all elements of the tuple:
    auto result = tuple_map(std::move(plus), lhs.blocks(), rhs.blocks());

    // Return the result enwrapped inside a BlockDiagonalMatrix:
    return make_block_diagonal(std::move(result));
}

template <typename... MatricesLHS, typename... MatricesRHS>
auto operator-(const BlockDiagonalMatrixBase<MatricesLHS...>& lhs,
               const BlockDiagonalMatrixBase<MatricesRHS...>& rhs) {
    // Generic lambdas are not part of c++11, so we cannot use them
    // auto bminus = [](auto& x, auto& y) { return x - y; };

    // use generic minus functor:
    detail::MinusFctr bminus;

    // Apply it to all elements of the tuple:
    auto result = tuple_map(std::move(bminus), lhs.blocks(), rhs.blocks());

    // Return the result enwrapped inside a BlockDiagonalMatrix:
    return make_block_diagonal(std::move(result));
}

template <typename... MatricesLHS, typename... MatricesRHS>
auto operator*(const BlockDiagonalMatrixBase<MatricesLHS...>& lhs,
               const BlockDiagonalMatrixBase<MatricesRHS...>& rhs) {
    // Avoid generic lambda (C++14 and above only)
    // auto times = [](auto& x, auto& y) { return x * y; };
    detail::MultipliesFctr times;

    // Apply it to all elements of the tuple:
    auto result = tuple_map(std::move(times), lhs.blocks(), rhs.blocks());

    // Return the result enwrapped inside a BlockDiagonalMatrix:
    return make_block_diagonal(std::move(result));
}

template <typename... Matrices>
auto operator-(const BlockDiagonalMatrixBase<Matrices...>& bm) {
    // Avoid generic lambda (C++14 and above only)
    // auto uminus = [](auto& x) { return -x; };
    detail::NegateFctr uminus;

    // Apply it to all elements of the tuple:
    auto result = tuple_map(std::move(uminus), bm.blocks());

    // Return the result enwrapped inside a BlockDiagonalMatrix:
    return make_block_diagonal(std::move(result));
}

template <typename Matrix, typename... OtherMatrices>
auto operator*(const BlockDiagonalMatrixBase<Matrix, OtherMatrices...>& bm,
               typename Matrix::scalar_type s) {
    // Avoid generic lambda (C++14 and above only)
    // auto scalemult = [&](auto& x) { return s * x; };
    detail::ScaleByFctr<typename Matrix::scalar_type> scalemult{s};

    // Apply it to all elements of the tuple:
    auto result = tuple_map(std::move(scalemult), bm.blocks());

    // Return the result enwrapped inside a BlockDiagonalMatrix:
    return make_block_diagonal(std::move(result));
}

template <typename Matrix, typename... OtherMatrices>
auto operator*(typename Matrix::scalar_type s,
               const BlockDiagonalMatrixBase<Matrix, OtherMatrices...>& bm) {
    return bm * s;
}

template <typename Matrix, typename... OtherMatrices>
auto operator/(const BlockDiagonalMatrixBase<Matrix, OtherMatrices...>& bm,
               typename Matrix::scalar_type s) {
    // Avoid generic lambda (C++14 and above only)
    // auto scalediv = [&](auto& x) { return x / s; };
    detail::DivideByFctr<typename Matrix::scalar_type> scalediv{s};

    // Apply it to all elements of the tuple:
    auto result = tuple_map(std::move(scalediv), bm.blocks());

    // Return the result enwrapped inside a BlockDiagonalMatrix:
    return make_block_diagonal(std::move(result));
}

}  // namespace linalgwrap
