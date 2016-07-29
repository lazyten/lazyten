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
#include "linalgwrap/Range.hh"
#include "linalgwrap/VectorOfSpecialise.hh"
#include <initializer_list>
#include <vector>

namespace linalgwrap {

/** Class to represent a row vector specialisation of an arbitrary matrix */
template <typename MatrixType>
class VectorOf : public detail::VectorOfSpecialise<MatrixType> {
  public:
    typedef detail::VectorOfSpecialise<MatrixType> base_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::size_type size_type;

    /** \name Constructors
     */
    ///@{
    /** Construct a matrix of fixed size and optionally set the entries to
     * zero
     */
    explicit VectorOf(size_type size, bool fill_zero = true);

    /** \brief Construct from an initialiser list.
     */
    explicit VectorOf(std::initializer_list<scalar_type> list);

    /** \brief Construct from std::vector */
    explicit VectorOf(const std::vector<scalar_type>& v);

    /** \brief Construct from Matrix type
     *
     * The matrix is asserted to be a 1xn matrix.
     */
    explicit VectorOf(const MatrixType& m);

    /** \brief Construct from iterators */
    template <class InputIterator>
    VectorOf(InputIterator first, InputIterator last);

    ///@}

    /** \name Vector operations
     *  \note For the outer Product of a Row and a column vector we take the
     * default implementation in the \t MatrixType class.
     * The inner product is implemented in VectorOfBase or VectorOfSpecialise
     * by the \t dot_with() function. Additionally a specialisation for the
     * operation TransposeView<VectorOf> * VectorOf exists and internally calls
     * the dot_with() function of this class.
     * */
    ///@{
    /** Scale vector by a scalar value */
    VectorOf& operator*=(scalar_type s);

    /** Divide all vector entries by a scalar value */
    VectorOf& operator/=(scalar_type s);

    /* Add a vector to this one */
    VectorOf& operator+=(const VectorOf& other);

    /* Subtract a small matrix from this one */
    VectorOf& operator-=(const VectorOf& other);
    ///@}
};

//
// Multiply by Scalar
//
template <typename MatrixType>
VectorOf<MatrixType> operator*(typename MatrixType::scalar_type s,
                               VectorOf<MatrixType> m) {
    m *= s;
    return m;
}

template <typename MatrixType>
VectorOf<MatrixType> operator*(VectorOf<MatrixType> m,
                               typename MatrixType::scalar_type s) {
    return s * m;
}

template <typename MatrixType>
VectorOf<MatrixType> operator/(VectorOf<MatrixType> m,
                               typename MatrixType::scalar_type s) {
    m /= s;
    return m;
}

//
// Add and subtract small matrices
//
template <typename MatrixType>
VectorOf<MatrixType> operator-(VectorOf<MatrixType> lhs,
                               const VectorOf<MatrixType>& rhs) {
    lhs -= rhs;
    return lhs;
}

template <typename MatrixType>
VectorOf<MatrixType> operator+(VectorOf<MatrixType> lhs,
                               const VectorOf<MatrixType>& rhs) {
    lhs += rhs;
    return lhs;
}

//
// ------------------------------------------
//

//
// VectorOf
//
template <typename MatrixType>
VectorOf<MatrixType>::VectorOf(size_type size, bool fill_zero)
      : base_type{size, fill_zero} {}

template <typename MatrixType>
template <class InputIterator>
VectorOf<MatrixType>::VectorOf(InputIterator first, InputIterator last)
      : VectorOf(std::distance(first, last), false) {
    assert_greater(0, std::distance(first, last));
    std::copy(first, last, base_type::begin());
}

template <typename MatrixType>
VectorOf<MatrixType>::VectorOf(std::initializer_list<scalar_type> list)
      : VectorOf(list.size(), false) {
    auto itbase = base_type::begin();
    for (auto it = list.begin(); it != list.end(); ++it, ++itbase) {
        *itbase = *it;
    }
}

template <typename MatrixType>
VectorOf<MatrixType>::VectorOf(const std::vector<scalar_type>& v)
      : VectorOf(v.begin(), v.end()) {}

template <typename MatrixType>
VectorOf<MatrixType>::VectorOf(const MatrixType& m)
      : VectorOf(m.n_rows(), false) {
    // Assert that the matrix size is correct.
    assert_size(m.n_cols(), 1);

    // If the distance between begin and end is not equal to the number
    // of rows the matrix contained unexpected sparsity.
    // TODO use some kind of "has sparsity" check here once it exists.
    assert_dbg(std::distance(std::begin(m), std::end(m)) == m.n_rows(),
               ExcMatrixNotDense());

    std::copy(std::begin(m), std::end(m), base_type::begin());
}

template <typename MatrixType>
VectorOf<MatrixType>& VectorOf<MatrixType>::operator*=(scalar_type s) {
    assert_finite(s);
    base_type::operator*=(s);
    return (*this);
}

template <typename MatrixType>
VectorOf<MatrixType>& VectorOf<MatrixType>::operator/=(scalar_type s) {
    assert_dbg(s != 0, ExcDevideByZero());
    assert_finite(s);
    base_type::operator/=(s);
    return (*this);
}

template <typename MatrixType>
VectorOf<MatrixType>& VectorOf<MatrixType>::operator+=(const VectorOf& other) {
    assert_size(this->size(), other.size());
    base_type::operator+=(other);
    return *this;
}

template <typename MatrixType>
VectorOf<MatrixType>& VectorOf<MatrixType>::operator-=(const VectorOf& other) {
    assert_size(this->size(), other.size());
    base_type::operator-=(other);
    return *this;
}

}  // end namespace linalgwrap
