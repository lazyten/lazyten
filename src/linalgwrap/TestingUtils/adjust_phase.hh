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
#include "krims_NumComp.hh"
#include <krims/TypeUtils.hh>

namespace linalgwrap {

/** Read element by element from object1 and and Object2
 * and if the phase differs, adjust it such that it
 * hopefully agrees.
 *
 * The idea behind this function is that eigensolvers can only compute
 * the eigenvector up to the sign and this ambiguity has to be fixed before
 * comparing the solver results and a reference
 */
template <typename Object1, typename Object2>
void adjust_phase(const Object1& from, Object2& to) {
  using namespace krims;

  typedef typename Object1::scalar_type scalar_type;
  static_assert(std::is_same<scalar_type, typename Object2::scalar_type>::value,
                "Both submitted objects should have the same scalar type.");

  auto itfrom = std::begin(from);
  auto itto = std::begin(to);
  for (; itfrom != std::end(from); ++itfrom, ++itto) {
    // Skip if the element is numerically zero
    if (numcomp(*itfrom).failure_action(NumCompActionType::Return) ==
        Constants<scalar_type>::zero) {
      continue;
    }

    // TODO The case of complex eigenvalues has not been properly thought
    // through when it comes to sign normalisation: It seems that one could
    // perform any rotation of real and imaginary part on the unit circle,
    // but I am not sure whether this is correct ... mfh
    assert_implemented(!IsComplexNumber<scalar_type>::value);

    // The sign of the first important element is different:
    if ((std::real(*itfrom) < 0. && std::real(*itto) > 0.) ||
        (std::real(*itfrom) > 0. && std::real(*itto) < 0.)) {
      to *= -1;
      return;
    }
  }
}

}  // namespace linalgwrap
