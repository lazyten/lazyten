//
// Copyright (C) 2016-17 by the linalgwrap authors
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
#include <krims/Functionals.hh>

namespace linalgwrap {
namespace gen {

/** \brief Generator for a container filled with numeric values
 *
 * Generate ``count`` numeric values and place them into the container.
 */
template <typename Container>
rc::Gen<Container> numeric_container(size_t count) {
  auto gen_container =
        rc::gen::container<Container>(count, numeric<typename Container::value_type>());
  auto fix_container_norm = [](Container c) {
    long double norm = 0.0;
    krims::ConjFctr conj;
    for (const auto& elem : c) norm += std::real(conj(elem) * elem);

    if (norm > max_norm * max_norm) {
      for (auto& elem : c) elem /= std::sqrt(norm);
    }

    return c;
  };

  return rc::gen::map(gen_container, std::move(fix_container_norm));
}

/** \brief Generator for a container filled with numeric values
 *
 * The number of values the container contains is determined
 * using numeric_size()
 */
template <typename Container>
rc::Gen<Container> numeric_container() {
  return rc::gen::mapcat(numeric_size<1>(), [](size_t count) {
    return numeric_container<Container>(count);
  });
}
}  // gen
}  // linalgwrap
