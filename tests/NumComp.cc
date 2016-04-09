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

#include "NumComp.hh"
#include <sstream>

namespace linalgwrap {
namespace tests {

NumCompException::NumCompException(double lhs_, double rhs_, double error_,
                                   double tolerance_,
                                   const std::string description_) noexcept
      : lhs{lhs_},
        rhs{rhs_},
        error{error_},
        tolerance{tolerance_},
        description{description_},
        what_str("") {}

const char* NumCompException::what() const noexcept {
    try {
        if (what_str == "") {
            std::stringstream ss;

            ss << std::scientific << std::setprecision(15) << lhs
               << " == " << rhs << " returned false. Error(" << error
               << "), Tolerance(" << tolerance
               << "), Description: " << description;

            what_str = ss.str();
        }
        return what_str.c_str();
    } catch (...) {
        return "Error building what_str";
    }
}
}
}
