## ---------------------------------------------------------------------
##
## Copyright (C) 2016 by the layzten authors
##
## This file is part of layzten.
##
## layzten is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## layzten is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with layzten. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

# adds entries to these things
#
#       LAZYTEN_DEPENDENCIES			everyone needs these libraries
#       LAZYTEN_DEPENDENCIES_DEBUG		debug mode needs these extras
#       LAZYTEN_DEPENDENCIES_RELEASE		release mode needs these extras
#       LAZYTEN_DEPENDENCIES_TEST		tests need these extra libraries
#

####################
#-- C++ standard --#
####################
if (NOT CMAKE_CXX_STANDARD VERSION_LESS 14)
	message(STATUS "Detected C++14 support: Setting LAZYTEN_HAVE_CXX14")
	set(LAZYTEN_HAVE_CXX14 ON)
endif()
if (NOT CMAKE_CXX_STANDARD VERSION_LESS 17)
	message(STATUS "Detected C++17 support: Setting LAZYTEN_HAVE_CXX17")
	set(LAZYTEN_HAVE_CXX17 ON)
endif()

################
#--  ARPACK  --#
################
SET(ARPACK_DIR "" CACHE PATH "An optional hint to an ARPACK installation")
disable_feature(arpack)

find_library(
	ARPACK_LIBRARY
	NAMES arpack
	HINTS ${ARPACK_DIR}
	PATH_SUFFIXES lib64 lib
)
if(NOT "${ARPACK_LIBRARY}" MATCHES "-NOTFOUND")
	message(STATUS "Found ARPACK at ${ARPACK_LIBRARY}")
	set(LAZYTEN_DEPENDENCIES ${LAZYTEN_DEPENDENCIES} ${ARPACK_LIBRARY})
	enable_feature(arpack)
endif()

#################
#--  Bohrium  --#
#################
option(LAZYTEN_ENABLE_BOHRIUM "Enable bohrium as a computational backend." OFF)
disable_feature(bohrium)

# TODO Provide way to search for Bohrium and enable it only if found
if (LAZYTEN_ENABLE_BOHRIUM)
	# Forward parameters to included module
	set(BOHRIUM_VERSION 0.8.3) # We need at least this version
	include(cmake/findBohrium.cmake)
	unset(BOHRIUM_VERSION)

	set(LAZYTEN_DEPENDENCIES ${LAZYTEN_DEPENDENCIES} ${bohrium_TARGET})
	enable_feature(bohrium)
endif()
