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
#ifdef LINALGWRAP_HAVE_ARMADILLO

#include "ArmadilloMatrix.hh"
#include "linalgwrap/BaseInterfaces.hh"
#include <armadillo>
#include <iterator>

/** A class for an armadillo column vector with the interface of linalgwrap */
namespace linalgwrap {

// Forward-declare
template <typename Scalar>
class ArmadilloMatrix;

template <typename Scalar>
class ArmadilloVector : public StoredVector_i<Scalar> {
    static_assert(
          std::is_same<double, Scalar>::value ||
                std::is_same<float, Scalar>::value ||
                std::is_same<std::complex<float>, Scalar>::value ||
                std::is_same<std::complex<double>, Scalar>::value ||
                std::is_same<short, Scalar>::value ||
                std::is_same<int, Scalar>::value ||
                std::is_same<long, Scalar>::value ||
                std::is_same<unsigned short, Scalar>::value ||
                std::is_same<unsigned int, Scalar>::value ||
                std::is_same<unsigned long, Scalar>::value,
          "ArmadilloVector<Scalar> is currently only available for Scalar "
          "being one of double, float, complex<double>, "
          "complex<float>,  short, int, long, unsigned short, unsigned "
          "int, unsigned long");

  public:
    typedef StoredVector_i<Scalar> base_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::real_type real_type;

    /** The corresponding matrix type */
    typedef ArmadilloMatrix<Scalar> matrix_type;

    /** The type of the storage object used to store the data
     *  of the ArmadilloVector */
    typedef arma::Col<Scalar> storage_type;
    typedef const storage_type const_storage_type;

    //! The iterator types from armadillo
    typedef typename storage_type::iterator iterator;
    typedef typename storage_type::const_iterator const_iterator;

    /** \name Constructors */
    ///@{
    /** \brief Construct vector of fixed size and optionally fill with zeros.
     *
     * \param fill_zero   If true all entries are set to zero
     */
    explicit ArmadilloVector(size_type size, bool fill_zero = true)
          : m_arma(size, arma::fill::none) {
        if (fill_zero) m_arma.zeros();
    }

    /** \brief Construct from initialiser list */
    explicit ArmadilloVector(std::initializer_list<Scalar> list)
          : m_arma(list) {}

    /** \brief Construct from std::vector */
    explicit ArmadilloVector(std::vector<scalar_type> v)
          : m_arma(std::move(v)) {}

    /** \brief Construct from input iterator */
    template <class InputIterator>
    ArmadilloVector(InputIterator first, InputIterator last)
          : ArmadilloVector(std::distance(first, last), false) {
        std::copy(first, last, m_arma.begin());
    }

    /** \brief Construct from Arbitrary Indexable_i */
    template <typename Indexable, typename = typename std::enable_if<
                                        IsIndexable<Indexable>::value>::type>
    explicit ArmadilloVector(Indexable i) : ArmadilloVector(i.n_elem(), false) {
        std::move(i.begin(), i.end(), m_arma.begin());
    }

    /** \brief Construct from inner storage object */
    explicit ArmadilloVector(storage_type inner) : m_arma(std::move(inner)) {}
    ///@}

    /** \name Size of the vector */
    size_type n_elem() const override { return m_arma.n_elem; }

    /** \name Size of the vector */
    size_type size() const override { return m_arma.n_elem; }

    /** \name Data access
     */
    ///@{
    /** \brief return an element of the vector    */
    scalar_type operator()(size_type i) const override {
        assert_range(0, i, n_elem());
        return m_arma[i];
    }

    /** \brief return an element of the vector */
    scalar_type operator[](size_type i) const override {
        assert_range(0, i, n_elem());
        return m_arma[i];
    }

    /** \brief return an element of the vector    */
    scalar_type& operator()(size_type i) override {
        assert_range(0, i, n_elem());
        return m_arma[i];
    }

    /** \brief return an element of the vector */
    scalar_type& operator[](size_type i) override {
        assert_range(0, i, n_elem());
        return m_arma[i];
    }
    ///@}

    /** Set all elements of the vector to zero */
    void set_zero() override { m_arma.zeros(); }

    /** \name Iterators
     */
    ///@{
    /** Return an iterator to the beginning */
    iterator begin() { return m_arma.begin(); }

    /** Return a const_iterator to the beginning */
    const_iterator begin() const { return m_arma.begin(); }

    /** Return a const_iterator to the beginning */
    const_iterator cbegin() const { return m_arma.begin(); }

    /** Return an iterator to the end */
    iterator end() { return m_arma.end(); }

    /** Return a const_iterator to the end */
    const_iterator end() const { return m_arma.end(); }

    /** Return a const_iterator to the end */
    const_iterator cend() const { return m_arma.end(); }
    ///@}

    /** \name Vector operations */
    ///@{
    /** Scale vector by a scalar value */
    ArmadilloVector& operator*=(scalar_type s) {
        assert_finite(s);
        m_arma *= s;
        return *this;
    }

    /** Divide all vector entries by a scalar value */
    ArmadilloVector& operator/=(scalar_type s) {
        assert_nonzero(s);
        assert_finite(s);
        m_arma /= s;
        return *this;
    }

    /* Add a vector to this one */
    ArmadilloVector& operator+=(const ArmadilloVector& other) {
        assert_size(n_elem(), other.n_elem());
        m_arma += other.m_arma;
        return *this;
    }

    /* Add a vector to this one */
    ArmadilloVector& operator-=(const ArmadilloVector& other) {
        assert_size(n_elem(), other.n_elem());
        m_arma -= other.m_arma;
        return *this;
    }

    bool operator==(const ArmadilloVector& other) const {
        if (m_arma.n_elem != other.m_arma.n_elem) return false;
        return std::equal(m_arma.begin(), m_arma.end(), other.m_arma.begin());
    }

    // Note != taken from default implementation

    ///@}

    /** Read-only access to the inner storage */
    const_storage_type& data() const { return m_arma; }

    /** Read-write access to the inner storage (use with caution) */
    storage_type& data() { return m_arma; }

  protected:
    storage_type m_arma;
};

//
// Standard operations
//
template <typename Scalar>
ArmadilloVector<Scalar> operator*(Scalar s, ArmadilloVector<Scalar> m) {
    m *= s;
    return m;
}

template <typename Scalar>
ArmadilloVector<Scalar> operator*(ArmadilloVector<Scalar> m, Scalar s) {
    return s * m;
}

template <typename Scalar>
ArmadilloVector<Scalar> operator/(ArmadilloVector<Scalar> m, Scalar s) {
    m /= s;
    return m;
}

template <typename Scalar>
ArmadilloVector<Scalar> operator-(ArmadilloVector<Scalar> mat) {
    return -Constants<Scalar>::one * mat;
}

//
// Add and subtract small matrices
//
template <typename Scalar>
ArmadilloVector<Scalar> operator-(ArmadilloVector<Scalar> lhs,
                                  const ArmadilloVector<Scalar>& rhs) {
    lhs -= rhs;
    return lhs;
}

template <typename Scalar>
ArmadilloVector<Scalar> operator+(ArmadilloVector<Scalar> lhs,
                                  const ArmadilloVector<Scalar>& rhs) {
    lhs += rhs;
    return lhs;
}

//
// Specialisations of fallback operations
//

// -- dot
template <typename Scalar1, typename Scalar2>
typename std::common_type<Scalar1, Scalar2>::type dot(
      const ArmadilloVector<Scalar1>& A, const ArmadilloVector<Scalar2>& B) {
    return arma::dot(A.data(), B.data());
}

template <typename Scalar1, typename Scalar2>
typename std::common_type<Scalar1, Scalar2>::type cdot(
      const ArmadilloVector<Scalar1>& A, const ArmadilloVector<Scalar2>& B) {
    // arma does the complex conjugate in the first argument as well.
    return arma::cdot(A.data(), B.data());
}

// -- accumulate
template <typename Scalar>
Scalar accumulate(const ArmadilloVector<Scalar>& A) {
    return arma::accu(A.data());
}

// -- minmax
template <typename Scalar>
Scalar min(const ArmadilloVector<Scalar>& A) {
    return arma::min(A.data());
}

template <typename Scalar>
Scalar max(const ArmadilloVector<Scalar>& A) {
    return arma::max(A.data());
}

// -- norms
template <typename Scalar>
typename ArmadilloVector<Scalar>::real_type norm_l1(
      const ArmadilloVector<Scalar>& A) {
    return arma::norm(A.data(), 1);
}

template <typename Scalar>
typename ArmadilloVector<Scalar>::real_type norm_linf(
      const ArmadilloVector<Scalar>& A) {
    return arma::norm(A.data(), "inf");
}

// for norm_l2_squared use default implementation.

template <typename Scalar>
typename ArmadilloVector<Scalar>::real_type norm_l2(
      const ArmadilloVector<Scalar>& A) {
    return arma::norm(A.data(), 2);
}

// -- elementwise
// (TODO this implementation introduces temporaries ... maybe not what we want)
template <typename Scalar>
ArmadilloVector<typename krims::RealTypeOf<Scalar>::type> abs(
      const ArmadilloVector<Scalar>& v) {
    typedef typename krims::RealTypeOf<Scalar>::type real_type;
    return ArmadilloVector<real_type>(arma::abs(v.data()));
}

template <typename Scalar>
ArmadilloVector<Scalar> conj(const ArmadilloVector<Scalar>& v) {
    return ArmadilloVector<Scalar>(arma::conj(v.data()));
}

template <typename Scalar>
ArmadilloVector<Scalar> sqrt(const ArmadilloVector<Scalar>& v) {
    return ArmadilloVector<Scalar>(arma::sqrt(v.data()));
}

template <typename Scalar>
ArmadilloVector<Scalar> square(const ArmadilloVector<Scalar>& v) {
    return ArmadilloVector<Scalar>(arma::square(v.data()));
}

}  // end namespace linalgwrap
#endif
