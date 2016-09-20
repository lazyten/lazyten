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

#include "linalgwrap/BaseInterfaces.hh"
#include <iterator>
#include <memory>

namespace linalgwrap {

/** A very basic class for a stored vector */
template <typename Scalar>
class BuiltinVector : public Stored_i, public MutableVector_i<Scalar> {
  public:
    typedef StoredVector_i<Scalar> base_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::real_type real_type;

    // TODO currently there is none
    /** The corresponding matrix type */
    typedef void matrix_type;

    /** The type of the storage object used to store the data */
    typedef std::unique_ptr<scalar_type[]> storage_type;
    typedef const std::unique_ptr<const scalar_type[]> const_storage_type;

    //! The iterator types
    typedef scalar_type* iterator;
    typedef const scalar_type* const_iterator;

    /** \name Constructors */
    ///@{
    /** \brief Construct vector of fixed size and optionally fill with zeros.
     *
     * \param fill_zero   If true all entries are set to zero
     */
    explicit BuiltinVector(size_type size, bool fill_zero = true)
          : m_data(new scalar_type[size]), m_size(size) {
        if (fill_zero) set_zero();
    }

    /** \brief Construct from initialiser list */
    BuiltinVector(std::initializer_list<Scalar> list)
          : BuiltinVector(list.size(), false) {
        std::move(std::begin(list), std::end(list), begin());
    }

    /** \brief Construct from std::vector */
    explicit BuiltinVector(std::vector<scalar_type> v)
          : BuiltinVector(v.size(), false) {
        std::move(std::begin(v), std::end(v), begin());
    }

    /** \brief Construct from input iterator */
    template <class InputIterator>
    BuiltinVector(InputIterator first, InputIterator last)
          : BuiltinVector(std::distance(first, last), false) {
        std::move(first, last, begin());
    }

    /** \brief Construct from Arbitrary Indexable_i */
    template <typename Indexable, typename = typename std::enable_if<
                                        IsIndexable<Indexable>::value>::type>
    explicit BuiltinVector(Indexable i) : BuiltinVector(i.n_elem(), false) {
        std::move(i.begin(), i.end(), begin());
    }
    ///@}

    ~BuiltinVector() = default;
    BuiltinVector& operator=(BuiltinVector&&) = default;
    BuiltinVector(BuiltinVector&&) = default;

    /** Copy assignment operator */
    BuiltinVector& operator=(const BuiltinVector& other) {
        // TODO test this function!
        //
        if (m_size != other.m_size) {
            // size is different: Reallocate memory:
            m_size = other.m_size;
            m_data.reset(new scalar_type[other.m_size]);
        }

        std::copy(other.begin(), other.end(), begin());
    }

    /** Copy constructor */
    BuiltinVector(const BuiltinVector& other)
          : BuiltinVector(other.m_size, false) {
        std::copy(other.begin(), other.end(), begin());
    }

    /** \name Size of the vector */
    size_type n_elem() const override { return m_size; }

    /** \name Size of the vector */
    size_type size() const override { return m_size; }

    /** \name Data access
     */
    ///@{
    /** \brief return an element of the vector    */
    scalar_type operator()(size_type i) const override {
        assert_range(0, i, n_elem());
        return m_data[i];
    }

    /** \brief return an element of the vector */
    scalar_type operator[](size_type i) const override {
        assert_range(0, i, n_elem());
        return m_data[i];
    }

    /** \brief return an element of the vector    */
    scalar_type& operator()(size_type i) override {
        assert_range(0, i, n_elem());
        return m_data[i];
    }

    /** \brief return an element of the vector */
    scalar_type& operator[](size_type i) override {
        assert_range(0, i, n_elem());
        return m_data[i];
    }
    ///@}

    /** Set all elements of the vector to zero */
    virtual void set_zero() override {
        std::fill(begin(), end(), Constants<scalar_type>::zero);
    }

    /** \name Iterators
     */
    ///@{
    /** Return an iterator to the beginning */
    iterator begin() { return m_data.get(); }

    /** Return a const_iterator to the beginning */
    const_iterator begin() const { return m_data.get(); }

    /** Return a const_iterator to the beginning */
    const_iterator cbegin() const { return m_data.get(); }

    /** Return an iterator to the end */
    iterator end() { return m_data.get() + m_size; }

    /** Return a const_iterator to the end */
    const_iterator end() const { return m_data.get() + m_size; }

    /** Return a const_iterator to the end */
    const_iterator cend() const { return m_data.get() + m_size; }
    ///@}

    /** \name Vector operations */
    ///@{
    /** Scale vector by a scalar value */
    BuiltinVector& operator*=(scalar_type s) {
        assert_finite(s);
        std::transform(begin(), end(), begin(),
                       [&](scalar_type& v) { return v * s; });
        return *this;
    }

    /** Divide all vector entries by a scalar value */
    BuiltinVector& operator/=(scalar_type s) {
        assert_nonzero(s);
        assert_finite(s);
        std::transform(begin(), end(), begin(),
                       [&](scalar_type& v) { return v / s; });
        return *this;
    }

    /* Add a vector to this one */
    BuiltinVector& operator+=(const BuiltinVector& other) {
        assert_size(n_elem(), other.n_elem());
        std::transform(
              begin(), end(), std::begin(other), begin(),
              [](scalar_type& v1, scalar_type& v2) { return v1 + v2; });
        return *this;
    }

    /* Add a vector to this one */
    BuiltinVector& operator-=(const BuiltinVector& other) {
        assert_size(n_elem(), other.n_elem());
        std::transform(
              begin(), end(), std::begin(other), begin(),
              [](scalar_type& v1, scalar_type& v2) { return v1 - v2; });
        return *this;
    }

    bool operator==(const BuiltinVector& other) const {
        if (m_size != other.m_size) return false;
        return std::equal(begin(), end(), std::begin(other));
    }

    // Note != taken from default implementation

    ///@}

    /** Read-only access to the inner storage */
    const_storage_type& data() const { return m_data; }

    /** Read-write access to the inner storage (use with caution) */
    storage_type& data() { return m_data; }

  protected:
    storage_type m_data;
    size_type m_size;
};

//
// Standard operations
//
template <typename Scalar>
BuiltinVector<Scalar> operator*(Scalar s, BuiltinVector<Scalar> m) {
    m *= s;
    return m;
}

template <typename Scalar>
BuiltinVector<Scalar> operator*(BuiltinVector<Scalar> m, Scalar s) {
    return s * m;
}

template <typename Scalar>
BuiltinVector<Scalar> operator/(BuiltinVector<Scalar> m, Scalar s) {
    m /= s;
    return m;
}

template <typename Scalar>
BuiltinVector<Scalar> operator-(BuiltinVector<Scalar> mat) {
    return -Constants<Scalar>::one * mat;
}

//
// Add and subtract small matrices
//
template <typename Scalar>
BuiltinVector<Scalar> operator-(BuiltinVector<Scalar> lhs,
                                const BuiltinVector<Scalar>& rhs) {
    lhs -= rhs;
    return lhs;
}

template <typename Scalar>
BuiltinVector<Scalar> operator+(BuiltinVector<Scalar> lhs,
                                const BuiltinVector<Scalar>& rhs) {
    lhs += rhs;
    return lhs;
}

}  // end namespace linalgwrap
