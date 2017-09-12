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

// Macro to insert the required code
#define SMALL_MATRIX_CODE(MATRIX)     \
  namespace lazyten {                 \
  template <typename Scalar>          \
  class MATRIX;                       \
  template <typename Scalar>          \
  using SmallMatrix = MATRIX<Scalar>; \
  }  // end namespace lazyten

#if defined LAZYTEN_SMALL_MATRIX_ARMADILLO
#include "lazyten/Armadillo/ArmadilloMatrix.hh"
SMALL_MATRIX_CODE(ArmadilloMatrix)

//
#elif defined LAZYTEN_SMALL_MATRIX_BOHRIUM
#include "lazyten/Bohrium/BohriumMatrix.hh"
SMALL_MATRIX_CODE(BohriumMatrix)

//
#else
#error "No alias for SmallMatrix defined"
#endif

#undef SMALL_MATRIX_CODE
