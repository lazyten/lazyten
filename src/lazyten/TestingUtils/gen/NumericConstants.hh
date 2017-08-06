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
#include <cstddef>

namespace lazyten {
namespace gen {

/* \brief The maximum an arbitrary numeric value may be. */
static constexpr long double max_value = 1e5;

/** \brief The minimum an arbitrary numeric value may be */
static constexpr long double min_value = 1e-5;

/** \brief The maximum l2 norm all values in a numeric container/tensor may have
 *         (objects with larger norms will have their values scaled down)
 * */
static constexpr long double max_norm = max_value;

}  // namespace gen
}  // namespace lazyten
