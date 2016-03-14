#pragma once
#include <cstddef>
#include "Exceptions.hh"
#include <utility>
#include <type_traits>

namespace linalgwrap {

template <typename T>
class RangeIterator;

/** A range of integral values
 *
 * \note Empty ranges are allowed, but calling any function except
 * length(), size() or is_empty() on them leads to undefined behaviour.
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
    Range(std::pair<value_type, value_type> first_last);

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
    void assert_valid_state() const;

    value_type m_current;
    value_type m_last;
};

//
// ---------------------------------------------------
//

template <typename T>
Range<T>::Range(value_type first, value_type last)
      : m_first(first), m_last(last) {
    assert_lower_bound(first, last);
};

template <typename T>
Range<T>::Range(std::pair<value_type, value_type> first_last)
      : m_first(first_last.first), m_last(first_last.second){};

template <typename T>
typename Range<T>::value_type Range<T>::first() const {
    assert_dbg(length() > 0, ExcEmptyRange());
    return m_first;
}

template <typename T>
typename Range<T>::value_type Range<T>::last() const {
    assert_dbg(length() > 0, ExcEmptyRange());
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
    return length() == 0;
}

template <typename T>
typename Range<T>::value_type Range<T>::operator[](size_type i) const {
    assert_dbg(length() > 0, ExcEmptyRange());
    assert_range(0, i, length());
    return m_first + i;
}

template <typename T>
RangeIterator<T> Range<T>::begin() const {
    assert_dbg(length() > 0, ExcEmptyRange());
    return RangeIterator<T>{m_first, m_last};
}

template <typename T>
RangeIterator<T> Range<T>::end() const {
    assert_dbg(length() > 0, ExcEmptyRange());
    return RangeIterator<T>{};
}

//
// ------------------------------------------------
//

template <typename T>
void RangeIterator<T>::assert_valid_state() const {
    assert_dbg(m_current < m_last, ExcIteratorPastEnd());
}

template <typename T>
RangeIterator<T>::RangeIterator()
      : m_current{0}, m_last{0} {}

template <typename T>
RangeIterator<T>::RangeIterator(value_type current, value_type last)
      : m_current{current}, m_last{last} {
    assert_valid_state();
}

template <typename T>
RangeIterator<T>& RangeIterator<T>::operator++() {
    assert_valid_state();
    m_current++;
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
    return m_current == other.m_current;
}

template <typename T>
bool RangeIterator<T>::operator!=(const RangeIterator& other) const {
    return !(*this == other);
}

}  // namespace linalgwrap
