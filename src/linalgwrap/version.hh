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
#include <string>

#define VERSION_MAJOR 0
#define VERSION_MINOR 1
#define VERSION_REVISION 0

struct version {
    static int constexpr major{VERSION_MAJOR};
    static int constexpr minor{VERSION_MINOR};
    static int constexpr revision{VERSION_REVISION};

    // Return the version as a string
    static std::string version_string();
};
