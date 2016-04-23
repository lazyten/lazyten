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
#include "linalgwrap/ArmadilloMatrix.hh"
#include "linalgwrap/VectorOfBase.hh"

namespace linalgwrap {
namespace detail {

/** Class implementing the type-specific VectorOf functionality. */
template <typename MatrixType>
class VectorOfSpecialise : public VectorOfBase<MatrixType> {
  public:
    typedef VectorOfBase<MatrixType> base_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::size_type size_type;

    /** Construct a VectorOfSpecialise, pass all data to base_type */
    explicit VectorOfSpecialise(size_type size, bool fill_zero);
};

//
// ArmadialloMatrix
//
#ifdef LINALGWRAP_HAVE_ARMADILLO
/** \brief Class implementing the armadillo specific VectorOf functionality.
 *
 * \tparam Scalar   The scalar type to use
 * */
template <typename Scalar>
class VectorOfSpecialise<ArmadilloMatrix<Scalar>>
      : public VectorOfBase<ArmadilloMatrix<Scalar>> {
    static_assert(!std::is_same<std::complex<float>, Scalar>::value &&
                        !std::is_same<std::complex<double>, Scalar>::value,
                  "VectorOfSpecialise<ArmadilloMatrix<Scalar>> is currently "
                  "not well-tested for complex-valued data types. Hence "
                  "complex<double> and complex<float> cannot be used");

  public:
    typedef VectorOfBase<ArmadilloMatrix<Scalar>> base_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::size_type size_type;

    /** Construct a VectorOfSpecialise, pass all data to base_type */
    explicit VectorOfSpecialise(size_type size, bool fill_zero);

    /** \name Scalar product and norms */
    ///@{
    /** Calculate the scalar product
     *
     * \note TransposeView<VectorOf> * VectorOf calls this routine,
     *       so there is no difference between these two.
     * */
    scalar_type dot(const VectorOfBase<ArmadilloMatrix<Scalar>>& other) const;

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
#endif

//
// -------------------------------------------------------
//

template <typename MatrixType>
VectorOfSpecialise<MatrixType>::VectorOfSpecialise(size_type size,
                                                   bool fill_zero)
      : base_type{size, fill_zero} {}

//
// ArmadialloMatrix
//
#ifdef LINALGWRAP_HAVE_ARMADILLO
template <typename Scalar>
VectorOfSpecialise<ArmadilloMatrix<Scalar>>::VectorOfSpecialise(size_type size,
                                                                bool fill_zero)
      : base_type{size, fill_zero} {}

template <typename Scalar>
typename VectorOfSpecialise<ArmadilloMatrix<Scalar>>::scalar_type
VectorOfSpecialise<ArmadilloMatrix<Scalar>>::dot(
      const VectorOfBase<ArmadilloMatrix<Scalar>>& other) const {
    typedef typename ArmadilloMatrix<Scalar>::storage_type storage_type;
    const storage_type& thismat = this->data();
    const storage_type& othermat = other.data();
    // TODO if thismat is complex, we might need cdot here!
    return arma::dot(thismat, othermat);
}

template <typename Scalar>
typename VectorOfSpecialise<ArmadilloMatrix<Scalar>>::scalar_type
VectorOfSpecialise<ArmadilloMatrix<Scalar>>::norm_squared() const {
    return dot(*this);
}

template <typename Scalar>
typename VectorOfSpecialise<ArmadilloMatrix<Scalar>>::scalar_type
VectorOfSpecialise<ArmadilloMatrix<Scalar>>::l2_norm() const {
    return arma::norm(this->data());
}

template <typename Scalar>
typename VectorOfSpecialise<ArmadilloMatrix<Scalar>>::scalar_type
VectorOfSpecialise<ArmadilloMatrix<Scalar>>::l1_norm() const {
    return arma::norm(this->data(), 1);
}

template <typename Scalar>
typename VectorOfSpecialise<ArmadilloMatrix<Scalar>>::scalar_type
VectorOfSpecialise<ArmadilloMatrix<Scalar>>::linf_norm() const {
    return arma::norm(this->data(), "inf");
}

#endif

}  // namespace detail
}  // namespace linalgwrap
