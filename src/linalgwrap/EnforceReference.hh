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

namespace linalgwrap {
namespace detail {

/** \brief Class to enforce a reference to be returned
 *
 * The class requires prior knowledge about the situation
 * Only in the case that we get the data by-value it really
 * does anything (it copies the data to internal storage
 * and returns a reference to it)
 *
 * \tparam T The data type
 * \tparam FromValue Do we get the data by-value(true) or by-reference(false)
 * */
template <typename T, bool FromValue>
struct EnforceReference {};

/** \brief Class to enforce a reference to be returned
 *
 * This version is an identity operation, just returning
 * the reference it got.
 */
template <typename T>
struct EnforceReference<T, false> {
    typedef T& result_type;
    typedef T& argument_type;
    T& operator()(T& t) const;
};

/** \brief Class to enforce a reference to be returned
 *
 * This version copies the value to internal storage and
 * returns a reference to this internal value.
*/
template <typename T>
struct EnforceReference<T, true> {
    typedef T& result_type;
    typedef T argument_type;
    T& operator()(T t) const;

  private:
    mutable T dummy;
};

//
// -----------------------------------------------
//

template <typename T>
T& EnforceReference<T, false>::operator()(T& t) const {
    return t;
}

template <typename T>
T& EnforceReference<T, true>::operator()(T t) const {
    dummy = t;
    return dummy;
}

}  // namespace detail
}  // namespace linalgwrap
