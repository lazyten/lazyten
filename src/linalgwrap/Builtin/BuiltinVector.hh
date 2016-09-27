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

#include "BuiltinTypes.hh"
#include "linalgwrap/BaseInterfaces.hh"
#include "linalgwrap/VectorMemoryWrapper.hh"
#include <iterator>
#include <memory>

namespace linalgwrap {

/** \brief A very basic class for a stored vector
 *
 * \note This class is not intended to be fast. It just is a fallback.
 * */
template <typename Scalar>
class BuiltinVector : public VectorMemoryWrapper<Scalar*>, public Stored_i {
  public:
    typedef VectorMemoryWrapper<Scalar*> base_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::real_type real_type;

    /** The corresponding family of builtin linear algebra types */
    typedef BuiltinTypes type_family;

    //@{
    /** The type of the storage object used to store the data */
    typedef std::unique_ptr<scalar_type[]> storage_type;
    typedef const std::unique_ptr<const scalar_type[]> const_storage_type;
    //@}

    /** \name Constructors */
    ///@{
    /** \brief Construct vector of fixed size and optionally fill with zeros.
     *
     * \param fill_zero   If true all entries are set to zero
     */
    explicit BuiltinVector(size_type size, bool fill_zero = true)
          : base_type(), m_data(new scalar_type[size]) {
        base_type::initialise(m_data.get(), size);
        if (fill_zero) base_type::set_zero();
    }

    /** \brief Construct from initialiser list */
    BuiltinVector(std::initializer_list<Scalar> list)
          : BuiltinVector(list.size(), false) {
        std::move(std::begin(list), std::end(list), base_type::begin());
    }

    /** \brief Construct from std::vector */
    explicit BuiltinVector(std::vector<scalar_type> v)
          : BuiltinVector(v.size(), false) {
        std::move(std::begin(v), std::end(v), base_type::begin());
    }

    /** \brief Construct from input iterator */
    template <class InputIterator>
    BuiltinVector(InputIterator first, InputIterator last)
          : BuiltinVector(std::distance(first, last), false) {
        std::move(first, last, base_type::begin());
    }

    /** \brief Construct from Arbitrary Indexable_i */
    template <typename Indexable, typename = typename std::enable_if<
                                        IsIndexable<Indexable>::value>::type>
    explicit BuiltinVector(Indexable i) : BuiltinVector(i.n_elem(), false) {
        std::move(i.begin(), i.end(), base_type::begin());
    }
    ///@}

    ~BuiltinVector() = default;
    BuiltinVector& operator=(BuiltinVector&&) = default;
    BuiltinVector(BuiltinVector&&) = default;

    /** Copy assignment operator */
    BuiltinVector& operator=(const BuiltinVector& other) {
        if (base_type::size() != other.size()) {
            // size is different: Reallocate memory:
            base_type::m_size = other.m_size;
            m_data.reset(new scalar_type[other.m_size]);
        }
        std::copy(other.begin(), other.end(), base_type::begin());
        return *this;
    }

    /** Copy constructor */
    BuiltinVector(const BuiltinVector& other)
          : BuiltinVector(other.m_size, false) {
        std::copy(other.begin(), other.end(), base_type::begin());
    }

    /** \name Vector operations */
    ///@{
    /** Scale vector by a scalar value */
    BuiltinVector& operator*=(scalar_type s) {
        base_type::operator*=(s);
        return *this;
    }

    /** Divide all vector entries by a scalar value */
    BuiltinVector& operator/=(scalar_type s) {
        base_type::operator/=(s);
        return *this;
    }

    /* Add a vector to this one */
    BuiltinVector& operator+=(const BuiltinVector& other) {
        base_type::operator+=(other);
        return *this;
    }

    /* Add a vector to this one */
    BuiltinVector& operator-=(const BuiltinVector& other) {
        base_type::operator-=(other);
        return *this;
    }
    ///@}

    /** Read-only access to the inner storage */
    const_storage_type& data() const { return m_data; }

    /** Read-write access to the inner storage (use with caution) */
    storage_type& data() { return m_data; }

  protected:
    storage_type m_data;
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
