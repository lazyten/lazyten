//
// Copyright (C) 2017 by the linalgwrap authors
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

#include "OperatorProperties.hh"

namespace linalgwrap {

#define CASENAME(prop)           \
  case OperatorProperties::prop: \
    o << #prop;                  \
    return o;

std::ostream& operator<<(std::ostream& o, OperatorProperties prop) {
  switch (prop) {
    CASENAME(PositiveDefinite);
    CASENAME(PositiveSemiDefinite);
    CASENAME(RealSymmetric);
    CASENAME(Real);
    CASENAME(Hermitian);
    CASENAME(AntiHermitian);
    CASENAME(None);

    default:
      o << "<Unknown property>";
      return o;
  }
}

}  // namespace linalgwrap
