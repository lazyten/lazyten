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

// Forward declarations:
template <typename IteratorCore>
class MatrixIterator;

namespace detail {
template <typename Matrix, bool Constness>
class MatrixIteratorDefaultCore;
}  // namespace detail

//! The default matrix iterator
template <typename Matrix>
using DefaultMatrixIterator =
      MatrixIterator<detail::MatrixIteratorDefaultCore<Matrix, false>>;

//! The default matrix const iterator
template <typename Matrix>
using DefaultMatrixConstIterator =
      MatrixIterator<detail::MatrixIteratorDefaultCore<Matrix, true>>;

}  // namespace linalgwrap

// This is intentional, in order to make sure that the cyclic dependency
// between the iterators and the Matrix_i class are correctly resolved.
#include "linalgwrap/MatrixIterator.hh"
