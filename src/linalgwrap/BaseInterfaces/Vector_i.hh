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
#include "Indexable_i.hh"

namespace linalgwrap {

/** \brief Abstract vector interface class to represent a vector,
 * from which we can only read data.
 *
 * The following types should also be defined:
 *  - iterator  The type returned by begin()
 *  - const_iterator The type returned by cbegin()
 *
 * The following methods should be implemented:
 *  - ``begin()``, ``cbegin()``  Return an iterator or a constant iterator
 *    to the beginning of the stride of memory.
 *  - ``end()``, ``cend()``   Return an iterator/constant iterator to the
 *    position past-the-end of the stride of memory.
 */
template <typename Scalar>
class Vector_i : public Indexable_i<Scalar> {
  public:
    typedef Indexable_i<Scalar> base_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;

    /** The value type of the vectors (for compatibility with STL vectors */
    typedef scalar_type value_type;

    /** Size of the vector */
    virtual size_type size() const = 0;

    /** \name Data access
     *        Access to vector data
     */
    ///@{
    /** \brief return an element of the vector   */
    virtual scalar_type operator()(size_type i) const = 0;
    ///@}
};

//@{
/** \brief struct representing a type (std::true_type, std::false_type) which
 *  indicates whether T is a stored vector
 *
 * The definition is done using SFINAE, such that even for types not having a
 * typedef scalar_type this expression is valid.
 *  */
template <typename T, typename = void>
struct IsVector : public std::false_type {};

template <typename T>
struct IsVector<T, krims::VoidType<typename T::scalar_type>>
      : public std::is_base_of<Vector_i<typename T::scalar_type>, T> {};
//@}

/** Type alias to obtain the scalar type of a valid Vector object */
template <typename T>
using ValidVectorScalarT =
      typename std::enable_if<IsVector<T>::value,
                              typename T::scalar_type>::type;

/** \brief Simple output operator, that plainly shows all entries of
 *  the Matrix one by one.
 *
 *  Rows are separated by a newline and entries by spaces.
 *  The last row is not terminated by a newline character.
 *  */
template <typename Scalar>
std::ostream& operator<<(std::ostream& o, const Vector_i<Scalar>& v) {
    for (size_t i = 0; i < v.size(); ++i) {
        o << v[i] << " ";
    }

    // TODO extend
    // assert_dbg(false, krims::ExcNotImplemented());
    // io::MatrixPrinter().print(m, o);
    return o;
}

}  // namespace linalgwrap
