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
#include "MutableVector_i.hh"

namespace linalgwrap {

/** \brief Abstract vector interface class to represent a vector,
 * from which we can both read/write entries and get access to the
 * underlying memory as well via the memptr function.
 *
 * Unlike a StoredVector_i, we cannot necessarily allocate memory
 * for this object, since that might be managed elsewhere.
 *
 * One should consider to override the following methods:
 *  - ``begin()``, ``cbegin()``  Return an iterator or a constant iterator
 *    to the beginning of the stride of data.
 *  - ``end()``, ``cend()``   Return an iterator/constant iterator to the
 *    position past-the-end of the stride of data.
 *  - The types iterator and const_iterator should also be altered accordingly
 *
 * The following operators should be implemented between vectors of the
 * implementing type:
 *   - Addition(+), Subtraction(-), Scalar multiplication and division
 *   - In-place addition(+=) and subtraction(-=)
 *   - In-place scalar multiplication (*=)
 */
template <typename Scalar>
class MutableMemoryVector_i : public MutableVector_i<Scalar> {
 public:
  typedef MutableVector_i<Scalar> base_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::size_type size_type;

  //! The default iterator types for MutableMemoryVectors
  //@{
  typedef scalar_type* iterator;
  typedef const scalar_type* const_iterator;
  //@}

  /** \name Data access
   *        Access to vector data
   */
  ///@{
  /** Access to the underlying memory */
  virtual scalar_type* memptr() = 0;

  /** Const access to the underlying memory */
  virtual const scalar_type* memptr() const = 0;

  /** \brief return an element of the vector   */
  virtual scalar_type operator()(size_type i) const override { return (*this)[i]; }

  /** \brief return an element of the vector   */
  virtual scalar_type operator[](size_type i) const override {
    assert_range(0, i, this->n_elem());
    return *(memptr() + i);
  }

  /** \brief return an element of the vector    */
  virtual scalar_type& operator()(size_type i) override { return (*this)[i]; }

  /** \brief return an element of the vector   */
  virtual scalar_type& operator[](size_type i) override {
    assert_range(0, i, this->n_elem());
    return *(memptr() + i);
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
  iterator begin() { return memptr(); }

  /** Return a const_iterator to the beginning */
  const_iterator begin() const { return cbegin(); }

  /** Return a const_iterator to the beginning */
  const_iterator cbegin() const { return memptr(); }

  /** Return an iterator to the end */
  iterator end() { return (memptr() + this->n_elem()); }

  /** Return a const_iterator to the end */
  const_iterator end() const { return cend(); }

  /** Return a const_iterator to the end */
  const_iterator cend() const { return (memptr() + this->n_elem()); }
  ///@}
};

//@{
/** \brief struct representing a type (std::true_type, std::false_type) which
 *  indicates whether T is a MutableMemoryVector.
 *  */
template <typename T, typename = void>
struct IsMutableMemoryVector : public std::false_type {};

template <typename T>
struct IsMutableMemoryVector<T, krims::VoidType<typename T::scalar_type>>
      : public std::is_base_of<MutableMemoryVector_i<typename T::scalar_type>, T> {};
//@}
}  // namespace linalgwrap
