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
#include "Stored_i.hh"

namespace linalgwrap {

/** \brief struct representing a type (std::true_type, std::false_type) which
 *  indicates whether T is a stored vector.
 *
 *  A stored vector is a mutable vector, which also inherits off the
 *  Stored_i marker class.
 *
 *  The expectation is that it satisfies the following:
 *
 * If the object is not const both read and write access to the data is
 * given.
 *
 * We expect any stored vector class to also provide the following constructors:
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
 *   Vector_i(std::vector<scalar_type> v);
 *   ```
 *
 * - Construct from an arbitrary indexable:
 *   ```
 *   template <typename Indexable, typename = typename std::enable_if<
 *                                      IsIndexable<Indexable>::value>::type>
 *   explicit ArmadilloVector(Indexable i);
 *   ```
 * The default vector operations (see MutalbeVector_i for details) should
 * also be implemented between stored vectors of the same kind.
 *
 * The following types should also be defined:
 *  - iterator  The type returned by begin()
 *  - const_iterator The type returned by cbegin()
 *  - type_family  The struct containing the corresponding vector and matrix
 *                 types which are compatible to this vector type
 *                 (see Armadillo folder for an example)
 *
 * The following methods should be implemented:
 *  - ``begin()``, ``cbegin()``  Return an iterator or a constant iterator
 *    to the beginning of the stride of memory.
 *  - ``end()``, ``cend()``   Return an iterator/constant iterator to the
 *    position past-the-end of the stride of memory.
 */

template <typename T>
struct IsStoredVector : public std::integral_constant<
                              bool, IsMutableVector<T>::value &&
                                          std::is_base_of<Stored_i, T>::value> {
};

}  // namespace linalgwrap
