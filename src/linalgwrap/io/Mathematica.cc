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

#include "Mathematica.hh"

namespace linalgwrap {
namespace io {

Mathematica::Mathematica() : m_thresh{1e-16}, m_check_for_thresh{false} {};

Mathematica::Mathematica(double thresh)
      : m_thresh{thresh}, m_check_for_thresh{true} {};

const std::vector<std::string> Mathematica::extensions{"m"};

}  // namespace io
}  // namespace linalgwrap
