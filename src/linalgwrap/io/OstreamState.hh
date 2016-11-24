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
#include <ios>
#include <ostream>

namespace linalgwrap {
namespace io {

/** Class to store a current internal state of an ostream
 * and restore it upon destruction of the class.
 */
class OstreamState {
 public:
  //! Construct from an ostream, taking its current state
  OstreamState(std::ostream& o);

  //! Destruct class, set state back to the ostream.
  ~OstreamState() { restore_to(m_o); }

  //! Default copy constructor
  OstreamState(const OstreamState&) = default;

  //! Default move constructor
  OstreamState(OstreamState&&) = default;

  //! Default copy assignment
  OstreamState& operator=(const OstreamState&) = default;

  //! Default move assignment
  OstreamState& operator=(OstreamState&&) = default;

  /** Set original state to this ostream as well
   *
   *  \note This class automatically resets the state of the ostream passed
   *  upon construction as soon as it is destructed.
   **/
  void restore_to(std::ostream& o) const;

 private:
  //! The ostream we manage
  std::ostream& m_o;

  //! The original format flags it had:
  std::ios::fmtflags m_orig_flags;

  //! The original stream width it had:
  std::streamsize m_orig_width;

  //! The original precision it had:
  std::streamsize m_orig_precision;

  //! The original fill character it had:
  char m_orig_fill;
};

//
// --------------------------------------
//

inline OstreamState::OstreamState(std::ostream& o)
      : m_o(o),
        m_orig_flags{o.flags()},
        m_orig_width{o.width()},
        m_orig_precision{o.precision()},
        m_orig_fill{o.fill()} {}

inline void OstreamState::restore_to(std::ostream& o) const {
  o.flags(m_orig_flags);
  o.width(m_orig_width);
  o.precision(m_orig_precision);
  o.fill(m_orig_fill);
}

}  // namespace io
}  // namespace linalgwrap
