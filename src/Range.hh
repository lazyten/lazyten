#pragma once
#include "Constants.hh"
#include "Exceptions.hh"
#include <cstddef>
#include <iostream>
#include <type_traits>
#include <utility>

namespace linalgwrap {

template <typename T>
class RangeIterator;

/** A range of integral values
 *
 * \note Empty ranges are allowed, but calling any function in order to
 * access elements of the range(``first()``,``last()``, ``operator[]``)
 * leads to undefined behaviour. In this case ``begin()`` is furthermore
 * equivalent to ``end()``.
 * */
template <typename T>
class Range {
  public:
    static_assert(std::is_integral<T>::value,
                  "T needs to be an integral data type");

    typedef T value_type;
    typedef size_t size_type;

    DefExceptionMsg(
          ExcEmptyRange,
          "The range object you attempted to use represents an empty range"
          " and hence cannot be used in this way.");

    /** \brief Construct a range
     *
     * Note that the interval is half-open, i.e. the range
     * is including the first element but not the last.
     * */
    Range(value_type first, value_type last);
    /** \brief Construct from pair
     *
     * The first element is inclusive, the last exclusive
     * (half-open interval)
     */
    explicit Range(std::pair<value_type, value_type> first_last);

    /** \brief Return the effective number of elements in the range
     *
     * i.e. if this class represents [3,5) it returns 2.
     * */
    size_type length() const;

    /** \brief Return the effective number of elements in the range
     *
     * alias to length()
     */
    size_type size() const;

    /** Get the first element, which is inclusive */
    value_type first() const;

    /** Get the last element, which is exclusive */
    value_type last() const;

    /** Is this range empty */
    bool is_empty() const;

    /** Return the ith value of the range */
    value_type operator[](size_type i) const;

    /** Return an iterator to the first element of the range */
    RangeIterator<T> begin() const;

    /** Return an iterator to the last element of the range */
    RangeIterator<T> end() const;

  private:
    value_type m_first;  // inclusive
    value_type m_last;   // exclusive
};

/** Output operator for ranges */
template <typename T>
std::ostream& operator<<(std::ostream& o, const Range<T>& r);

/** Iterator for ranges */
template <typename T>
class RangeIterator : public std::iterator<std::input_iterator_tag, T> {
  public:
    typedef std::iterator<std::input_iterator_tag, T> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::reference reference;
    typedef typename base_type::pointer pointer;

    /** Constructs an iterator-past-the-end */
    RangeIterator();

    /** Construct an iterator which currently points at
     * current and runs until last-1 (i.e. last is *not*
     * included any more.
     * */
    RangeIterator(value_type current, value_type last);

    /** Prefix increment to the next value */
    RangeIterator& operator++();

    /** Postfix increment to the next value */
    RangeIterator operator++(int);

    /** Return the value of the element we point to. */
    value_type operator*() const;

    /** Access the members of the element we point to. */
    const pointer operator->() const;

    /** Check if two iterators are equal */
    bool operator==(const RangeIterator& other) const;

    /** Check whether two iterators are unequal */
    bool operator!=(const RangeIterator& other) const;

  private:
    /** Does this data structure represent an
     *  iterator-past-the-end */
    bool is_past_the_end() const;

    /** Assert that the internal state is valid */
    void assert_valid_state() const;

    value_type m_current;
    value_type m_last;
};

//
// Helper functions for ranges:
//
/** Return a range interval from 0 to \p t, i.e. 0 is included, but \p t not. */
template <typename T>
Range<T> range(const T& t);

/** Return a range interval from \p t1 to \p t2, where \p t1 is included,
 * but \p t2 not. */
template <typename T>
Range<T> range(const T& t1, const T& t2);

//
// ---------------------------------------------------
//

template <typename T>
Range<T>::Range(value_type first, value_type last)
      : m_first(first), m_last(last) {
    assert_greater_equal(m_first, m_last);
};

template <typename T>
Range<T>::Range(std::pair<value_type, value_type> first_last)
      : m_first(first_last.first), m_last(first_last.second) {
    assert_greater_equal(m_first, m_last);
};

template <typename T>
typename Range<T>::value_type Range<T>::first() const {
    assert_dbg(!is_empty(), ExcEmptyRange());
    return m_first;
}

template <typename T>
typename Range<T>::value_type Range<T>::last() const {
    assert_dbg(!is_empty(), ExcEmptyRange());
    return m_last;
}

template <typename T>
typename Range<T>::size_type Range<T>::length() const {
    return m_last - m_first;
}

template <typename T>
typename Range<T>::size_type Range<T>::size() const {
    return length();
}

template <typename T>
bool Range<T>::is_empty() const {
    return length() <= 0;
}

template <typename T>
typename Range<T>::value_type Range<T>::operator[](size_type i) const {
    assert_dbg(!is_empty(), ExcEmptyRange());
    assert_range(0, i, length());
    return m_first + i;
}

template <typename T>
RangeIterator<T> Range<T>::begin() const {
    if (is_empty()) return end();
    return RangeIterator<T>{m_first, m_last};
}

template <typename T>
RangeIterator<T> Range<T>::end() const {
    return RangeIterator<T>{};
}

template <typename T>
std::ostream& operator<<(std::ostream& o, const Range<T>& r) {
    if (r.is_empty()) {
        o << "[0,0)";
    } else {
        o << "[" << r.first() << "," << r.last() << ")";
    }
    return o;
}

//
// ------------------------------------------------
//

template <typename T>
void RangeIterator<T>::assert_valid_state() const {
    assert_dbg(!is_past_the_end(), ExcIteratorPastEnd());
}

template <typename T>
bool RangeIterator<T>::is_past_the_end() const {
    return m_current >= m_last;
}

template <typename T>
RangeIterator<T>::RangeIterator()
      : m_current{Constants<T>::invalid}, m_last{Constants<T>::invalid} {}

template <typename T>
RangeIterator<T>::RangeIterator(value_type current, value_type last)
      : m_current{current}, m_last{last} {
    assert_valid_state();
}

template <typename T>
RangeIterator<T>& RangeIterator<T>::operator++() {
    assert_valid_state();
    m_current++;
    return (*this);
}

template <typename T>
RangeIterator<T> RangeIterator<T>::operator++(int) {
    assert_valid_state();
    RangeIterator<T> copy{*this};
    ++(*this);
    return copy;
}

template <typename T>
typename RangeIterator<T>::value_type RangeIterator<T>::operator*() const {
    assert_valid_state();
    return m_current;
}

template <typename T>
const typename RangeIterator<T>::pointer RangeIterator<T>::operator->() const {
    assert_valid_state();
    return &m_current;
}

template <typename T>
bool RangeIterator<T>::operator==(const RangeIterator& other) const {
    // The iterators are equal if they are either both past the end
    // or their m_current and their m_last agrees

    const bool both_past_the_end = is_past_the_end() && other.is_past_the_end();
    const bool identical_values =
          (m_current == other.m_current && m_last == other.m_last);
    return both_past_the_end || identical_values;
}

template <typename T>
bool RangeIterator<T>::operator!=(const RangeIterator& other) const {
    return !(*this == other);
}

//
// Range helper functions.
//

/** Return a range interval from 0 to \p t, i.e. 0 is included, but \p t not. */
template <typename T>
Range<T> range(const T& t) {
    return Range<T>{0, t};
}

/** Return a range interval from \p t1 to \p t2, where \p t1 is included,
 * but \p t2 not. */
template <typename T>
Range<T> range(const T& t1, const T& t2) {
    return Range<T>{t1, t2};
}

}  // namespace linalgwrap
