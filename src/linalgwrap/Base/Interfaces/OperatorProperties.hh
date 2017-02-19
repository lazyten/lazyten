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

#pragma once
#include <ostream>
#include <type_traits>

namespace linalgwrap {

/** Enum flag which describes the kind of operator at hand
 *
 * \note Do not trust the meaning of the bits in the enum class.
 *       This may change from time to time.
 * */
enum class OperatorProperties : uint32_t {
  /** operator problem, i.e. no implicit assumptions about
   *  scalar type or structure of the operator can be made. */
  None = 0x00,

  /** The operator is Hermitian (or Symmetric if it is real as well) */
  Hermitian = 0x01,

  /** The operator is real valued */
  Real = 0x02,

  /** The operator is real and symmetric */
  RealSymmetric = (Real + Hermitian),

  /** The operator is real, symmetric and positive semi-definite
   *
   * \note This flag implies Hermitian and Real.
   **/
  PositiveSemiDefinite = (0x04 + RealSymmetric),

  /** The operator is real, symmetric and positive definite
   *
   * \note This flag implies Hermitian and PositiveSemiDefinite.
   * */
  PositiveDefinite = (PositiveSemiDefinite + 0x08),

  /** The operator is anti-Hermitian (or anti-Symmetric) */
  AntiHermitian = 0x80,
};

inline OperatorProperties operator&(OperatorProperties l, OperatorProperties r) {
  using type = typename std::underlying_type<OperatorProperties>::type;
  return static_cast<OperatorProperties>(static_cast<type>(l) & static_cast<type>(r));
}

inline OperatorProperties operator|(OperatorProperties l, OperatorProperties r) {
  using type = typename std::underlying_type<OperatorProperties>::type;
  return static_cast<OperatorProperties>(static_cast<type>(l) | static_cast<type>(r));
}

inline OperatorProperties& operator&=(OperatorProperties& l, OperatorProperties r) {
  l = l & r;
  return l;
}

inline OperatorProperties& operator|=(OperatorProperties& l, OperatorProperties r) {
  l = l | r;
  return l;
}

/** Return true if the first set of properties is contained in the second
 *
 * ``other`` may contain further flags, but all flags in ``mask`` have to be set.
 * */
inline bool props_contained_in(OperatorProperties mask, OperatorProperties other) {
  return (mask & other) == mask;
}

/** Stream output operator for OperatorProperties objects */
std::ostream& operator<<(std::ostream& o, OperatorProperties prop);

}  // namespace linalgwrap
