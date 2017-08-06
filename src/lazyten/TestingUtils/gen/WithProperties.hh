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
#include "lazyten/Base/Interfaces/OperatorProperties.hh"

namespace lazyten {
namespace gen {

namespace detail {
template <typename Matrix>
struct WithProperties {
  static_assert(IsStoredMatrix<Matrix>::value,
                "Matrix needs to be a stored matrix type.");

  static rc::Gen<Matrix> with_properties(rc::Gen<Matrix> gen, OperatorProperties prop) {
    switch (prop) {
      case OperatorProperties::None:
        return gen;
      case OperatorProperties::RealSymmetric: {
        auto make_symmetric = [](Matrix m) {
          if (m.n_rows() != m.n_cols()) {
            RC_DISCARD("Matrix not quadratic: Could not make matrix symmetric");
          }

          m.symmetrise();
          if (!m.check_properties_satisfied(OperatorProperties::RealSymmetric)) {
            RC_DISCARD("Could not make matrix symmetric");
          }
          m.add_properties(OperatorProperties::RealSymmetric);
          return m;
        };
        return rc::gen::map(gen, make_symmetric);
      }
      default:
        RC_DISCARD("Property not yet available");
        return gen;
    }
  }
};
}  // namespace detail

/** Generate a matrix with the given matrix properties using the provided generator
 *  to generate the raw matrix */
template <typename Matrix>
rc::Gen<Matrix> with_properties(rc::Gen<Matrix> gen, OperatorProperties prop) {
  return detail::WithProperties<Matrix>::with_properties(std::move(gen), prop);
}

/** Generate a matrix with the given matrix properties using the numeric tensor
 *  generator as the raw generator */
template <typename Matrix>
rc::Gen<Matrix> with_properties(OperatorProperties prop) {
  return with_properties(gen::numeric_tensor<Matrix>(), prop);
}

}  // namespace gen
}  // namespace lazyten
