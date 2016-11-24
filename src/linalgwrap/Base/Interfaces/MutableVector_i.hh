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
#include "Vector_i.hh"

namespace linalgwrap {
/** \brief Abstract vector interface class to represent a vector,
 * from which we can both read entries and modify them as well
 * (if the object is not const)
 *
 * The idea is that unlike the Vector_i entries may not only be read, but
 * also modified, i.e. written. Unlike the StoredVector_i, however, the
 * storage is handled externally, such that it is not possible to construct
 * such an object by the constructors laid out in StoredVector_i.
 *
 * The following types should also be defined:
 *  - iterator  The type returned by begin()
 *  - const_iterator The type returned by cbegin()
 *
 * The following methods should be implemented:
 *  - ``begin()``, ``cbegin()``  Return an iterator or a constant iterator
 *    to the beginning of the stride of data.
 *  - ``end()``, ``cend()``   Return an iterator/constant iterator to the
 *    position past-the-end of the stride of data.
 *
 * The following operators should be implemented between vectors of the
 * implementing type:
 *   - Addition(+), Subtraction(-), Scalar multiplication and division
 *   - In-place addition(+=) and subtraction(-=)
 *   - In-place scalar multiplication (*=)
 */
template <typename Scalar>
class MutableVector_i : public Vector_i<Scalar> {
 public:
  typedef Vector_i<Scalar> base_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;

  /** \name Data access
   *        Access to vector data
   */
  ///@{
  /** \brief return an element of the vector    */
  virtual scalar_type& operator()(size_type i) = 0;

  /** \brief return an element of the vector   */
  virtual scalar_type& operator[](size_type i) = 0;
  ///@}

  /** Set all elements of the vector to zero */
  virtual void set_zero() = 0;
};

//@{
/** \brief struct representing a type (std::true_type, std::false_type) which
 *  indicates whether T is a mutable vector.
 *  */
template <typename T, typename = void>
struct IsMutableVector : public std::false_type {};

template <typename T>
struct IsMutableVector<T, krims::VoidType<typename T::scalar_type>>
      : public std::is_base_of<MutableVector_i<typename T::scalar_type>, T> {};
//@}

}  // namespace linalgwrap
