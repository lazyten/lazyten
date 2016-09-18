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

/** \brief Abstract vector interface class to represent a stored column vector
 *
 * We expect any implementing class to also provide the following constructors:
 * - Construct vector of fixed size and optionally fill with zeros or leave
 *   memory unassigned:
 *   ```
 *   Vector_i(size, fill_zero);
 *   ```
 *
 * - Construct from an initialiser list
 *   ```
 *   Vector_i(std::initializer_list<Scalar> list);
 *   ```
 *
 * - Construct from a range provided as iterators:
 *   ```
 *   template <class InputIterator>
 *   Vector_i(InputIterator first, InputIterator last);
 *   ```
 *
 * - Construct from a std::vector<Scalar>:
 *   ```
 *   Vector_i(const std::vector<scalar_type>& v);
 *   ```
 *
 * - Construct from an arbitrary indexable:
 *   ```
 *   template <typename Indexable, typename = typename std::enable_if<
 *                                      IsIndexable<Indexable>::value>::type>
 *   explicit ArmadilloVector(const Indexable& i);
 *   ```
 *
 * Furthermore an implementation of the default vector operations
 * +, -, +=, -= as well as multiplication and division by a scalar
 * both in-place (*=) and out-of-place (*) should be provided
 *
 * The following types should also be defined:
 *  - iterator  The type returned by begin()
 *  - const_iterator The type returned by cbegin()
 *  - matrix_type  The type of the corresponding matrix.
 *
 * The following methods should be implemented:
 *  - ``begin()``, ``cbegin()``  Return an iterator or a constant iterator
 *    to the beginning of the stride of memory.
 *  - ``end()``, ``cend()``   Return an iterator/constant iterator to the
 *    position past-the-end of the stride of memory.
 */
template <typename Scalar>
class StoredVector_i : public Vector_i<Scalar> {
  public:
    typedef Vector_i<Scalar> base_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;

    /** \name Data access
     *        Access to vector data
     */
    ///@{
    /** \brief return an element of the vector   */
    virtual scalar_type operator()(size_type i) const = 0;

    /** \brief return an element of the vector    */
    virtual scalar_type& operator()(size_type i) = 0;

    /** \brief return an element of the vector   */
    virtual scalar_type& operator[](size_type i) = 0;
    ///@}
};

//@{
/** \brief struct representing a type (std::true_type, std::false_type) which
 *  indicates whether T is a stored vector.
 *  */
template <typename T, typename = void>
struct IsStoredVector : public std::false_type {};

template <typename T>
struct IsStoredVector<T, krims::VoidType<typename T::scalar_type>>
      : public std::is_base_of<StoredVector_i<typename T::scalar_type>, T> {};
//@}

}  // namespace linalgwrap
