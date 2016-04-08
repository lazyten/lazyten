#pragma once
#include "Range.hh"
#include <initializer_list>
#include <vector>

namespace linalgwrap {

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

/** Class to represent a row vector specialisation of an arbitrary matrix */
template <typename MatrixType>
class VectorOf : public VectorOfSpecialise<MatrixType> {
  public:
    typedef VectorOfSpecialise<MatrixType> base_type;
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
// VectorOfBase
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

//
// VectorOfSpecialise
//
template <typename MatrixType>
VectorOfSpecialise<MatrixType>::VectorOfSpecialise(size_type size,
                                                   bool fill_zero)
      : base_type{size, fill_zero} {}

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

// Include the specialisations for VectorOf.hh
#include "VectorOf.special.hh"
