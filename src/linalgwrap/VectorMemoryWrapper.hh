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

namespace linalgwrap {
namespace detail {
template <typename T>
struct MakeConstIterator {
    typedef std::iterator_traits<T> tt;

    struct type : public T {
        typedef typename tt::iterator_category iterator_category;
        typedef
              typename std::add_const<typename tt::value_type>::type value_type;
        typedef typename tt::difference_type difference_type;
        typedef const typename tt::pointer pointer;
        typedef const typename tt::reference reference;

        type(T t) : T(t) {}
        reference operator*() const { return T::operator*(); }
        pointer operator->() const { return T::operator->(); }
    };
};

template <typename T>
struct MakeConstIterator<T*> {
    typedef const T* type;
};

template <typename T>
struct MakeConstIterator<const T*> {
    typedef const T* type;
};
}  // namespace detail

// TODO VectorMemoryWrapper is maybe not the best name for this guy
// Think about renaming it to something else some day.

/** \brief Simple class to offer a vector-like view into an arbitrary stride
 * of memory.
 *
 * The memory is either specified by an iterator range or
 * by an initial position and a size.
 *
 * \tparam Iterator   the iterator which represents the range of memory
 * \tparam ConstIterator   the type of the equivalent const iterator.
 *                         If missing it will either be deduced or if not
 *                         possible a constant access view will be employed.
 */
template <typename Iterator,
          typename ConstIterator =
                typename detail::MakeConstIterator<Iterator>::type>
class VectorMemoryWrapper
      : public MutableVector_i<typename std::decay<
              typename std::iterator_traits<Iterator>::value_type>::type> {
  public:
    //@{
    //! The iterator types
    typedef Iterator iterator;
    typedef ConstIterator const_iterator;
    //@}

    typedef typename std::decay<typename std::iterator_traits<
          iterator>::value_type>::type scalar_type;
    typedef MutableVector_i<scalar_type> base_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::real_type real_type;

    /** \name Constructors */
    ///@{
    /** \brief Construct a VectorMemoryWrapper over an iterator range.  */
    VectorMemoryWrapper(iterator begin, iterator end) {
        initialise(std::forward<iterator>(begin), std::forward<iterator>(end));
    }

    /** \brief Construct a VectorMemoryWrapper by providing a start of the range
     * and a size
     *
     * \param size Number of elements this container holds, including the one
     * pointed to by begin
     * */
    VectorMemoryWrapper(iterator begin, size_type size) {
        initialise(std::forward<iterator>(begin), size);
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
        // assertion done in operator[]
        return (*this)[i];
    }

    /** \brief return an element of the vector */
    scalar_type operator[](size_type i) const override;

    /** \brief return an element of the vector    */
    scalar_type& operator()(size_type i) override {
        // assertion done in operator[]
        return (*this)[i];
    }

    /** \brief return an element of the vector */
    scalar_type& operator[](size_type i) override;
    ///@}

    /** Set all elements of the vector to zero */
    virtual void set_zero() override {
        std::fill(begin(), end(), Constants<scalar_type>::zero);
    }

    /** \name Iterators
     */
    ///@{
    /** Return an iterator to the beginning */
    iterator begin() { return m_begin; }

    /** Return a const_iterator to the beginning */
    const_iterator begin() const { return cbegin(); }

    /** Return a const_iterator to the beginning */
    const_iterator cbegin() const { return const_iterator(m_begin); }

    /** Return an iterator to the end */
    iterator end() { return m_end; }

    /** Return a const_iterator to the end */
    const_iterator end() const { return cend(); }

    /** Return a const_iterator to the end */
    const_iterator cend() const { return const_iterator(m_end); }
    ///@}

    /** \name Vector operations */
    ///@{
    /** Scale vector by a scalar value */
    VectorMemoryWrapper& operator*=(scalar_type s) {
        assert_finite(s);
        std::transform(begin(), end(), begin(),
                       [&](scalar_type& v) { return v * s; });
        return *this;
    }

    /** Divide all vector entries by a scalar value */
    VectorMemoryWrapper& operator/=(scalar_type s) {
        assert_nonzero(s);
        assert_finite(s);
        std::transform(begin(), end(), begin(),
                       [&](scalar_type& v) { return v / s; });
        return *this;
    }

    /* Add a vector to this one */
    VectorMemoryWrapper& operator+=(const VectorMemoryWrapper& other) {
        assert_size(n_elem(), other.n_elem());
        std::transform(
              begin(), end(), std::begin(other), begin(),
              [](scalar_type& v1, scalar_type& v2) { return v1 + v2; });
        return *this;
    }

    /* Add a vector to this one */
    VectorMemoryWrapper& operator-=(const VectorMemoryWrapper& other) {
        assert_size(n_elem(), other.n_elem());
        std::transform(
              begin(), end(), std::begin(other), begin(),
              [](scalar_type& v1, scalar_type& v2) { return v1 - v2; });
        return *this;
    }

    bool operator==(const VectorMemoryWrapper& other) const {
        if (m_size != other.m_size) return false;
        return std::equal(begin(), end(), std::begin(other));
    }

  protected:
    /** Allow deferred initialisation construction
     *
     * Before this class is used you *must* call initialise.
     **/
    VectorMemoryWrapper() {}

    /** Initialise from range */
    void initialise(iterator begin, iterator end);

    /** Initialise from begin and size */
    void initialise(iterator begin, size_type size);

    iterator m_begin;
    iterator m_end;
    size_type m_size;
};

/** Convenience function to construct a vector memory wrapper class */
template <typename Iterator>
VectorMemoryWrapper<Iterator> make_vector_mem_wrap(Iterator begin,
                                                   Iterator end) {
    return VectorMemoryWrapper<Iterator>(begin, end);
}

/** Convenience function to construct a vector memory wrapper class */
template <typename Iterator, typename Size>
VectorMemoryWrapper<Iterator> make_vector_mem_wrap(Iterator begin, Size size) {
    return VectorMemoryWrapper<Iterator>(begin, size);
}

//
// Standard operations
//
template <typename Iterator, typename ConstIterator>
VectorMemoryWrapper<Iterator, ConstIterator> operator*(
      typename VectorMemoryWrapper<Iterator, ConstIterator>::scalar_type s,
      VectorMemoryWrapper<Iterator, ConstIterator> m) {
    m *= s;
    return m;
}

template <typename Iterator, typename ConstIterator>
VectorMemoryWrapper<Iterator, ConstIterator> operator*(
      VectorMemoryWrapper<Iterator, ConstIterator> m,
      typename VectorMemoryWrapper<Iterator, ConstIterator>::scalar_type s) {
    return s * m;
}

template <typename Iterator, typename ConstIterator>
VectorMemoryWrapper<Iterator, ConstIterator> operator/(
      VectorMemoryWrapper<Iterator, ConstIterator> m,
      typename VectorMemoryWrapper<Iterator, ConstIterator>::scalar_type s) {
    m /= s;
    return m;
}

template <typename Iterator, typename ConstIterator>
VectorMemoryWrapper<Iterator, ConstIterator> operator-(
      VectorMemoryWrapper<Iterator, ConstIterator> mat) {
    typedef typename VectorMemoryWrapper<Iterator, ConstIterator>::scalar_type
          scalar_type;
    return -Constants<scalar_type>::one * mat;
}

//
// Add and subtract
//
template <typename Iterator, typename ConstIterator>
VectorMemoryWrapper<Iterator, ConstIterator> operator-(
      VectorMemoryWrapper<Iterator, ConstIterator> lhs,
      const VectorMemoryWrapper<Iterator, ConstIterator>& rhs) {
    lhs -= rhs;
    return lhs;
}

template <typename Iterator, typename ConstIterator>
VectorMemoryWrapper<Iterator, ConstIterator> operator+(
      VectorMemoryWrapper<Iterator, ConstIterator> lhs,
      const VectorMemoryWrapper<Iterator, ConstIterator>& rhs) {
    lhs += rhs;
    return lhs;
}

//
// -------------------------------------------------------------
//

template <typename Iterator, typename ConstIterator>
typename VectorMemoryWrapper<Iterator, ConstIterator>::scalar_type
      VectorMemoryWrapper<Iterator, ConstIterator>::operator[](
            size_type i) const {
    assert_range(0, i, n_elem());
    iterator cpy(m_begin);
    std::advance(cpy, i);
    return *cpy;
}

template <typename Iterator, typename ConstIterator>
typename VectorMemoryWrapper<Iterator, ConstIterator>::scalar_type&
      VectorMemoryWrapper<Iterator, ConstIterator>::operator[](size_type i) {
    assert_range(0, i, n_elem());
    iterator cpy(m_begin);
    std::advance(cpy, i);
    return *cpy;
}

template <typename Iterator, typename ConstIterator>
void VectorMemoryWrapper<Iterator, ConstIterator>::initialise(iterator begin,
                                                              iterator end) {
    m_begin = std::move(begin);
    m_end = std::move(end);
    m_size = std::distance(m_begin, m_end);
}

template <typename Iterator, typename ConstIterator>
void VectorMemoryWrapper<Iterator, ConstIterator>::initialise(iterator begin,
                                                              size_type size) {
    m_begin = std::move(begin);
    m_end = m_begin;
    m_size = size;
    std::advance(m_end, m_size);
}

}  // end namespace linalgwrap
