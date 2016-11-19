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
#include "Numeric.hh"
#include "NumericSize.hh"

namespace linalgwrap {
namespace gen {

/** \brief Generator for a container filled with numeric values
 *
 * Generate ``count`` numeric values and place them into the container.
 */
template <typename Container>
rc::Gen<Container> numeric_container(size_t count) {
    return rc::gen::container<Container>(
          count, numeric<typename Container::value_type>());
}

/** \brief Generator for a container filled with numeric values
 *
 * The number of values the container contains is determined
 * using numeric_size()
 */
template <typename Container>
rc::Gen<Container> numeric_container() {
    return rc::gen::exec(
          [] { return *numeric_container<Container>(*numeric_size<1>()); });
}
}  // gen
}  // linalgwrap
