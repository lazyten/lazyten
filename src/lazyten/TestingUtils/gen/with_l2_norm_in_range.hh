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
#include "NumericTensor.hh"

namespace lazyten {
namespace gen {

/** Create a numeric tensor which has a norm within a provided range (both sides are
 * inclusive)
 *  The raw tensors are drawn from the generator gen.
 */
template <typename Tensor>
rc::Gen<Tensor> with_l2_norm_in_range(typename Tensor::real_type min_norm,
                                      typename Tensor::real_type max_norm,
                                      rc::Gen<Tensor> gen) {
  auto map_to_norm_range = [min_norm, max_norm](Tensor t) {
    auto t_norm = norm_l2(t);

    // Discard cases with too small or too large norms
    // => cannot scale without great numerical error
    RC_PRE(t_norm >= min_norm / 1e6);
    RC_PRE(t_norm <= 1e6 * max_norm);

    // Bring to min_norm
    if (t_norm < min_norm) t *= min_norm / t_norm;

    // Bring to max_norm
    if (t_norm > max_norm) t *= max_norm / t_norm;

    return t;
  };
  return rc::gen::map(gen, std::move(map_to_norm_range));
}

/** Create a numeric tensor which has a norm within a provided range (both sides are
 * inclusive)
 *  The raw tensors are drawn from the generator numeric_tensor
 */
template <typename Tensor>
rc::Gen<Tensor> with_l2_norm_in_range(typename Tensor::real_type min_norm,
                                      typename Tensor::real_type max_norm) {
  return with_l2_norm_in_range(min_norm, max_norm, numeric_tensor<Tensor>());
}

}  // gen
}  // namespace lazyten
