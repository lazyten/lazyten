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
#include "linalgwrap/Constants.hh"
#include <algorithm>

namespace linalgwrap {
namespace detail {

/** Base class for representing a row vector specialisation of an arbitrary
 * matrix
 *
 * Implements all common functionality and gives defaults.
 * */
template <typename MatrixType>
class VectorOfBase : public MatrixType {
  public:
    typedef MatrixType base_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::size_type size_type;

    /** \name Constructors
     */
    ///@{
    /** Construct a matrix of fixed size and optionally set the entries to
     * zero
     */
    explicit VectorOfBase(size_type size, bool fill_zero);
    ///@}

    /** \name Vector information */
    ///@{
    /** Size of the vector */
    size_type size() const;

    /** Number of columns */
    size_type constexpr n_cols() const override;
    ///@}

    //@{
    /** Vector element access (identical to [])*/
    scalar_type operator()(size_type i) const;
    scalar_type& operator()(size_type i);
    //@}

    /** \name Scalar product and norms */
    ///@{
    /** Calculate the scalar product
     *
     * \note TransposeView<VectorOf> * VectorOf calls this routine,
     *       so there is no difference between these two.
     * */
    scalar_type dot(const VectorOfBase& other) const;

    /** Calculate the l2 norm squared */
    scalar_type norm_squared() const;

    /** Calculate the l2 norm */
    scalar_type l2_norm() const;

    /** Calculate the l1 norm */
    scalar_type l1_norm() const;

    /** Calculate the linf norm */
    scalar_type linf_norm() const;
    ///@}
};

//
// ------------------------------------------------------
//

template <typename MatrixType>
VectorOfBase<MatrixType>::VectorOfBase(size_type size, bool fill_zero)
      : base_type{size, 1, fill_zero} {}

template <typename MatrixType>
typename MatrixType::size_type VectorOfBase<MatrixType>::size() const {
    return base_type::n_rows();
}

template <typename MatrixType>
typename MatrixType::scalar_type VectorOfBase<MatrixType>::operator()(
      size_type i) const {
    return base_type::operator[](i);
}

template <typename MatrixType>
typename MatrixType::scalar_type& VectorOfBase<MatrixType>::operator()(
      size_type i) {
    return base_type::operator[](i);
}

template <typename MatrixType>
constexpr typename MatrixType::size_type VectorOfBase<MatrixType>::n_cols()
      const {
    return 1;
}

template <typename MatrixType>
typename MatrixType::scalar_type VectorOfBase<MatrixType>::dot(
      const VectorOfBase& other) const {
    assert_size(other.size(), size());

    scalar_type res = Constants<scalar_type>::zero;
    for (size_type i = 0; i < size(); ++i) {
        res += (*this)[i] * other[i];
    }
    return res;
}

template <typename MatrixType>
typename MatrixType::scalar_type VectorOfBase<MatrixType>::norm_squared()
      const {
    return dot(*this);
}

template <typename MatrixType>
typename MatrixType::scalar_type VectorOfBase<MatrixType>::l2_norm() const {
    return sqrt(norm_squared());
}

template <typename MatrixType>
typename MatrixType::scalar_type VectorOfBase<MatrixType>::l1_norm() const {
    scalar_type res = Constants<scalar_type>::zero;
    for (size_type i = 0; i < size(); ++i) {
        res += std::fabs((*this)[i]);
    }
    return res;
}

template <typename MatrixType>
typename MatrixType::scalar_type VectorOfBase<MatrixType>::linf_norm() const {
    scalar_type res = Constants<scalar_type>::zero;
    for (size_type i = 0; i < size(); ++i) {
        res = std::max(res, std::fabs((*this)[i]));
    }
    return res;
}

}  // namespace detail
}  // namespace linalgwrap
