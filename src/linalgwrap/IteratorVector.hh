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
#include "linalgwrap/Base/Interfaces.hh"

// TODO IteratorVector is maybe not the best name for this guy
// Think about renaming it to something else some day.

namespace linalgwrap {
namespace detail {
template <typename T>
struct MakeConstIterator {
  typedef std::iterator_traits<T> tt;

  struct type : public T {
    typedef typename tt::iterator_category iterator_category;
    typedef typename std::add_const<typename tt::value_type>::type value_type;
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

//
// -------------------------------------------
//

/** \brief Simple class to offer a vector-like view into an arbitrary stride
 * of memory, which we walk over using an iterator range.
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
          typename ConstIterator = typename detail::MakeConstIterator<Iterator>::type>
class IteratorVector : public MutableVector_i<typename std::decay<
                             typename std::iterator_traits<Iterator>::value_type>::type> {

 public:
  //@{
  //! The iterator types
  typedef Iterator iterator;
  typedef ConstIterator const_iterator;
  //@}

  typedef typename std::decay<typename std::iterator_traits<iterator>::value_type>::type
        scalar_type;
  typedef MutableVector_i<scalar_type> base_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::real_type real_type;

  // TODO lift this at some point?
  static_assert(!std::is_const<scalar_type>::value,
                "The scalar_type underlying the IteratorVector may not be constant.");

  /** \name Constructors */
  ///@{
  /** \brief Construct a IteratorVector over an iterator range.  */
  IteratorVector(iterator begin, iterator end) {
    initialise(std::forward<iterator>(begin), std::forward<iterator>(end));
  }

  /** \brief Construct a IteratorVector by providing a start of the range
   * and a size
   *
   * \param size Number of elements this container holds, including the one
   * pointed to by begin
   * */
  IteratorVector(iterator begin, size_type size) {
    initialise(std::forward<iterator>(begin), size);
  }
  ///@}

  /** \name Size of the vector */
  size_type n_elem() const override { return m_size; }

  /** \name Size of the vector */
  size_type size() const override { return m_size; }

  /** \name Data access
   */
  ///@{
  /** \brief return an element of the vector    */
  scalar_type operator()(size_type i) const override { return (*this)[i]; }

  /** \brief return an element of the vector */
  scalar_type operator[](size_type i) const override;

  /** \brief return an element of the vector    */
  scalar_type& operator()(size_type i) override { return (*this)[i]; }

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
  bool operator==(const IteratorVector& other) const {
    if (m_size != other.m_size) return false;
    return std::equal(begin(), end(), std::begin(other));
  }

  /** Scale vector by a scalar value */
  IteratorVector& operator*=(scalar_type s);

  /** Divide all vector entries by a scalar value */
  IteratorVector& operator/=(scalar_type s);

  /** Add a vector to this one */
  IteratorVector& operator+=(const IteratorVector& other);

  /** Add a vector to this one */
  IteratorVector& operator-=(const IteratorVector& other);
  ///@}

 protected:
  /** Allow deferred initialisation construction
   *
   * Before this class is used you *must* call initialise.
   **/
  IteratorVector() {}

  /** Initialise from range */
  void initialise(iterator begin, iterator end);

  /** Initialise from begin and size */
  void initialise(iterator begin, size_type size);

  iterator m_begin;
  iterator m_end;
  size_type m_size;
};

//
// -------------------------------------------------------------
//

/** \brief Simple class to offer a vector-like view into an arbitrary stride
 * of memory, which we walk over using an iterator range.
 *
 * Specialisation for pointers. This class also provides access to the
 * underlying pointer range via the memptr() function and hence
 * classifies as a full MutableMemoryVector.
 */
template <typename T, typename CT>
class IteratorVector<T*, CT> : public MutableMemoryVector_i<T> {
 public:
  typedef MutableMemoryVector_i<T> base_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;

  // TODO lift this at some point?
  static_assert(!std::is_const<scalar_type>::value,
                "The scalar_type underlying the IteratorVector may not be constant.");

  /** \name Constructors */
  ///@{
  /** \brief Construct a IteratorVector over an iterator range.  */
  IteratorVector(T* begin, T* end) : m_begin(begin), m_size(end - begin) {
    assert_greater(0, end - begin);
  }

  /** \brief Construct a IteratorVector by providing a start of the range
   * and a size
   *
   * \param size Number of elements this container holds, including the one
   * pointed to by begin
   * */
  IteratorVector(T* begin, size_type size) : m_begin(begin), m_size(size) {}
  ///@}

  /** \name Size of the vector */
  size_type n_elem() const override { return m_size; }

  /** \name Size of the vector */
  size_type size() const override { return m_size; }

  /** \name Vector operations */
  ///@{
  /** Scale vector by a scalar value */
  IteratorVector& operator*=(scalar_type s);

  /** Divide all vector entries by a scalar value */
  IteratorVector& operator/=(scalar_type s);

  /** Add a vector to this one */
  IteratorVector& operator+=(const IteratorVector& other);

  /** Add a vector to this one */
  IteratorVector& operator-=(const IteratorVector& other);

  bool operator==(const IteratorVector& other) const {
    if (m_size != other.m_size) return false;
    return std::equal(m_begin, m_begin + m_size, other.m_begin);
  }

  bool operator!=(const IteratorVector& other) const { return !(*this == other); }
  ///@}

  /** Read-only access to the raw memory */
  const T* memptr() const override { return m_begin; }

  /** Access to the raw memory */
  T* memptr() override { return m_begin; }

 protected:
  /** Allow deferred initialisation construction
   * Before this class is used you *must* call initialise.
   **/
  IteratorVector() : base_type() {}

  /** Initialise from range */
  void initialise(T* begin, T* end) {
    m_begin = begin;
    m_size = (end - begin);
  }

  /** Initialise from begin and size */
  void initialise(T* begin, size_type size) {
    m_size = size;
    m_begin = begin;
  }

  scalar_type* m_begin;
  size_type m_size;
};

//
// Standard operations
//
template <typename Iterator, typename ConstIterator>
IteratorVector<Iterator, ConstIterator> operator*(
      typename IteratorVector<Iterator, ConstIterator>::scalar_type s,
      IteratorVector<Iterator, ConstIterator> m) {
  m *= s;
  return m;
}

template <typename Iterator, typename ConstIterator>
IteratorVector<Iterator, ConstIterator> operator*(
      IteratorVector<Iterator, ConstIterator> m,
      typename IteratorVector<Iterator, ConstIterator>::scalar_type s) {
  return s * m;
}

template <typename Iterator, typename ConstIterator>
IteratorVector<Iterator, ConstIterator> operator/(
      IteratorVector<Iterator, ConstIterator> m,
      typename IteratorVector<Iterator, ConstIterator>::scalar_type s) {
  m /= s;
  return m;
}

template <typename Iterator, typename ConstIterator>
IteratorVector<Iterator, ConstIterator> operator-(
      IteratorVector<Iterator, ConstIterator> mat) {
  typedef typename IteratorVector<Iterator, ConstIterator>::scalar_type scalar_type;
  return -Constants<scalar_type>::one * mat;
}

//
// Add and subtract
//
template <typename Iterator, typename ConstIterator>
IteratorVector<Iterator, ConstIterator> operator-(
      IteratorVector<Iterator, ConstIterator> lhs,
      const IteratorVector<Iterator, ConstIterator>& rhs) {
  lhs -= rhs;
  return lhs;
}

template <typename Iterator, typename ConstIterator>
IteratorVector<Iterator, ConstIterator> operator+(
      IteratorVector<Iterator, ConstIterator> lhs,
      const IteratorVector<Iterator, ConstIterator>& rhs) {
  lhs += rhs;
  return lhs;
}

//
// --------------------------------------------------------
//

template <typename Iterator, typename ConstIterator>
typename IteratorVector<Iterator, ConstIterator>::scalar_type
      IteratorVector<Iterator, ConstIterator>::operator[](size_type i) const {
  assert_range(0, i, n_elem());
  iterator cpy(m_begin);
  std::advance(cpy, i);
  return *cpy;
}

template <typename Iterator, typename ConstIterator>
typename IteratorVector<Iterator, ConstIterator>::scalar_type&
      IteratorVector<Iterator, ConstIterator>::operator[](size_type i) {
  assert_range(0, i, n_elem());
  iterator cpy(m_begin);
  std::advance(cpy, i);
  return *cpy;
}

template <typename Iterator, typename ConstIterator>
void IteratorVector<Iterator, ConstIterator>::initialise(iterator begin, iterator end) {
  m_begin = std::move(begin);
  m_end = std::move(end);
  m_size = std::distance(m_begin, m_end);
}

template <typename Iterator, typename ConstIterator>
void IteratorVector<Iterator, ConstIterator>::initialise(iterator begin, size_type size) {
  m_begin = std::move(begin);
  m_end = m_begin;
  m_size = size;
  std::advance(m_end, m_size);
}

template <typename Iterator, typename ConstIterator>
IteratorVector<Iterator, ConstIterator>& IteratorVector<Iterator, ConstIterator>::
operator*=(scalar_type s) {
  assert_finite(s);
  std::transform(base_type::begin(), base_type::end(), base_type::begin(),
                 [&](scalar_type& v) { return v * s; });
  return *this;
}

template <typename Iterator, typename ConstIterator>
IteratorVector<Iterator, ConstIterator>& IteratorVector<Iterator, ConstIterator>::
operator/=(scalar_type s) {
  assert_nonzero(s);
  assert_finite(s);
  std::transform(base_type::begin(), base_type::end(), base_type::begin(),
                 [&](scalar_type& v) { return v / s; });
  return *this;
}

template <typename Iterator, typename ConstIterator>
IteratorVector<Iterator, ConstIterator>& IteratorVector<Iterator, ConstIterator>::
operator+=(const IteratorVector<Iterator, ConstIterator>& other) {
  assert_size(this->n_elem(), other.n_elem());
  std::transform(base_type::begin(), base_type::end(), std::begin(other),
                 base_type::begin(),
                 [](scalar_type& v1, scalar_type& v2) { return v1 + v2; });
  return *this;
}

template <typename Iterator, typename ConstIterator>
IteratorVector<Iterator, ConstIterator>& IteratorVector<Iterator, ConstIterator>::
operator-=(const IteratorVector<Iterator, ConstIterator>& other) {
  assert_size(this->n_elem(), other.n_elem());
  std::transform(base_type::begin(), base_type::end(), std::begin(other),
                 base_type::begin(),
                 [](scalar_type& v1, scalar_type& v2) { return v1 - v2; });
  return *this;
}

//
// --------------------------------------------------------
//

template <typename T, typename CT>
IteratorVector<T*, CT>& IteratorVector<T*, CT>::operator*=(scalar_type s) {
  assert_finite(s);
  std::transform(m_begin, m_begin + m_size, m_begin,
                 [&](scalar_type& v) { return v * s; });
  return *this;
}

template <typename T, typename CT>
IteratorVector<T*, CT>& IteratorVector<T*, CT>::operator/=(scalar_type s) {
  assert_nonzero(s);
  assert_finite(s);
  std::transform(m_begin, m_begin + m_size, m_begin,
                 [&](scalar_type& v) { return v / s; });
  return *this;
}

template <typename T, typename CT>
IteratorVector<T*, CT>& IteratorVector<T*, CT>::operator+=(const IteratorVector& other) {
  assert_size(this->n_elem(), other.n_elem());
  std::transform(m_begin, m_begin + m_size, other.m_begin, m_begin,
                 [](scalar_type& v1, scalar_type& v2) { return v1 + v2; });
  return *this;
}

template <typename T, typename CT>
IteratorVector<T*, CT>& IteratorVector<T*, CT>::operator-=(const IteratorVector& other) {
  assert_size(this->n_elem(), other.n_elem());
  std::transform(m_begin, m_begin + m_size, other.m_begin, m_begin,
                 [](scalar_type& v1, scalar_type& v2) { return v1 - v2; });
  return *this;
}
}  // end namespace linalgwrap
