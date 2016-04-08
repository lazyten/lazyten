#pragma once
#include "Range.hh"
#include "VectorOfSpecialise.hh"
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
    ///@}

    /** \name Vector operations
     *  \note For the outer Product of a Row and a column vector we take the
     * default implementation in the \t MatrixType class. i
     * The inner product is implemented in VectorOfBase or VectorOfSpecialise
     * by the \t dot() function. Additionally a specialisation for the operation
     * TransposeView<VectorOf> * VectorOf exists and internally calls the dot()
     * function of this class.
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
VectorOf<MatrixType>::VectorOf(std::initializer_list<scalar_type> list)
      : VectorOf(list.size(), false) {
    auto itbase = base_type::begin();
    for (auto it = list.begin(); it != list.end(); ++it, ++itbase) {
        *itbase = *it;
    }
}

template <typename MatrixType>
VectorOf<MatrixType>::VectorOf(const std::vector<scalar_type>& v)
      : VectorOf(v.size(), false) {
    std::copy(v.begin(), v.end(), base_type::begin());
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
