//
// Copyright (C) 2017 by the lazyten authors
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
#ifdef LAZYTEN_HAVE_BOHRIUM
namespace lazyten {

// Forward-declare
template <typename Scalar>
class BohriumMatrix;

// Forward-declare
template <typename Scalar>
class BohriumVector;

struct BohriumTypes {
  template <typename Scalar>
  using vector = BohriumVector<Scalar>;

  template <typename Scalar>
  using matrix = BohriumMatrix<Scalar>;
};

}  // namespace lazyten
#endif  // LAZYTEN_HAVE_BOHRIUM
