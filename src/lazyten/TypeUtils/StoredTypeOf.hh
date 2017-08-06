//
// Copyright (C) 2016-17 by the lazyten authors
//
// This file is part of lazyten.
//
// lazyten is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// lazyten is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with lazyten. If not, see <http://www.gnu.org/licenses/>.
//

#pragma once
#include "lazyten/LazyMatrixExpression.hh"
#include <krims/TypeUtils.hh>

namespace lazyten {

//@{
/** \brief helper struct to extract the underlying stored matrix
 *         type from a potentially lazy matrix.
 *
 * References or array extents are stripped off via std::decay
 *
 * The stored matrix type is accessible via the member type "type".
 *  */
template <typename Matrix, typename = void>
struct StoredTypeOf {
  // Neither lazy nor stored -> don't have any type
};

template <typename Matrix>
struct StoredTypeOf<Matrix, typename std::enable_if<IsStoredMatrix<
                                  typename std::decay<Matrix>::type>::value>::type> {
  typedef typename std::decay<Matrix>::type type;
};

template <typename Matrix>
struct StoredTypeOf<Matrix, typename std::enable_if<IsLazyMatrix<
                                  typename std::decay<Matrix>::type>::value>::type> {
  typedef typename std::decay<Matrix>::type::stored_matrix_type type;
};
//@}

}  // namespace lazyten
