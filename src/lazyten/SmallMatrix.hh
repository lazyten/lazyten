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
#include "lazyten/config.hh"

#include "lazyten/Armadillo/ArmadilloMatrix.hh"

namespace lazyten {
#if defined LAZYTEN_HAVE_ARMADILLO
template <typename Scalar>
class ArmadilloMatrix;

template <typename Scalar>
using SmallMatrix = ArmadilloMatrix<Scalar>;
#else
template <typename Scalar>
class SmallMatrix {
  static_assert(false, "No default implementation for SmallMatrix.");
};
#endif

}  // namespace lazyten
