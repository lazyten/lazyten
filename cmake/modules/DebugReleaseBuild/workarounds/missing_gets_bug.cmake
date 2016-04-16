## ---------------------------------------------------------------------
##
## Copyright (C) 2016 by the linalgwrap authors
##
## This file is part of linalgwrap.
##
## linalgwrap is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## linalgwrap is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with linalgwrap. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

option(PROVIDE_DUMMY_GETS "Do we have the missing gets bug?" OFF)

#
# This bug is only in glibc 2.16 in combination with the libstdc++ of gcc 4.8, 
# which are both shipped with Ubuntu LTS 14.04.
#
# See http://gcc.gnu.org/bugzilla/show_bug.cgi?id=51785 for some details.
#

if (PROVIDE_DUMMY_GETS) 
	FILE(WRITE "${${CMAKE_PROJECT_NAME}_BINARY_DIR}/gets_dummy.hh"
"#pragma once
extern \"C\" char* gets (char* __s);")

	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -include ${${CMAKE_PROJECT_NAME}_BINARY_DIR}/gets_dummy.hh")
endif()
