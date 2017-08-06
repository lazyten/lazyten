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
#include <krims/ExceptionSystem.hh>

namespace lazyten {

/** Should the matrix be applied in a transposed sense? */
enum class Transposed {
  /** The matrix is applied/extracted/multiplied directly */
  None,
  /** The transpose of the matrix is applied/extracted/multiplied */
  Trans,
  /** The conjugate transpose (Hermitian conjugate) of the matrix
   * is applied/extracted/multiplied */
  ConjTrans,
};

// Print a transposed object
std::ostream& operator<<(std::ostream& o, Transposed mode);

/** A matrix that was expected to be Hermetian is not */
DefException1(ExcUnsupportedOperationMode, Transposed,
              << "Operation mode " << arg1 << " not valid for this matrix.");

}  // namespace lazyten
