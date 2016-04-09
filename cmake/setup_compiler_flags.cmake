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

# Determine the compiler type and load corresponding flags.
# 
# sets the following variables
#       CMAKE_CXX_FLAGS				Flags for the c++ compiler, all builds
#	CMAKE_CXX_FLAGS_DEBUG			Extra flags for the debug build
#	CMAKE_CXX_FLAGS_RELEASE			Extra flags for the release build
#
#	CMAKE_SHARED_LINKER_FLAGS		Flags for the shared linker, all builds
#	CMAKE_SHARED_LINKER_FLAGS_DEBUG		Extra flags for the debug build
#	CMAKE_SHARED_LINKER_FLAGS_RELEASE	Extra flags for the release build
#
#	CMAKE_STATIC_LINKER_FLAGS		Flags for the static linker, all builds
#	CMAKE_STATIC_LINKER_FLAGS_DEBUG		Extra flags for the debug build
#	CMAKE_STATIC_LINKER_FLAGS_RELEASE	Extra flags for the release build
#
#	CMAKE_EXE_LINKER_FLAGS			Flags for linking executables, all builds
#	CMAKE_EXE_LINKER_FLAGS_DEBUG		Extra flags for the debug build
#	CMAKE_EXE_LINKER_FLAGS_RELEASE		Extra flags for the release build
#
#	LINALGWRAP_DEFINITIONS			Extra definitions for all builds
#	LINALGWRAP_DEFINITIONS_DEBUG		Extra definitions for debug builds
# 	LINALGWRAP_DEFINITIONS_RELEASE		Extra definitions for release builds

#################################
#--     Modules and macros    --#
#################################
include(CheckCXXCompilerFlag)

# macro to enforce c++ standard:
macro(USE_CXX_STANDARD STANDARD)
	if (CMAKE_VERSION VERSION_LESS "3.1")
		set (FLAG "c++${STANDARD}")

		if (STANDARD EQUAL "14")
			if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.5")
				# Clang <3.5 needs special treatment here, 
				# since c++14 flag is not yet know there
				set(FLAG "c++1y")
			endif()
		endif()
		set (CMAKE_CXX_FLAGS "-std=${FLAG} ${CMAKE_CXX_FLAGS}")
		unset(FLAG)
	else()
		# set the standard and enforce it:
		set(CMAKE_CXX_STANDARD ${STANDARD})
		set(CMAKE_CXX_STANDARD_REQUIRED ON)

		# make sure we no not use the gnu++ extensions:
		# this adds c++n instead of gnu++n to the command line
		set(CMAKE_CXX_EXTENSIONS OFF)
	endif()
endmacro(USE_CXX_STANDARD)

macro(enable_if_compiles VARIABLE FLAG)
	string(REGEX REPLACE "[^a-zA-Z0-9]" "" FLAG_CLEAN "${FLAG}")
	CHECK_CXX_COMPILER_FLAG("-Werror ${FLAG}" LINALGWRAP_HAVE_FLAG_${FLAG_CLEAN})
	if (LINALGWRAP_HAVE_FLAG_${FLAG_CLEAN})
		set(${VARIABLE} "${${VARIABLE}} ${FLAG}")
	endif()
	unset(FLAG_CLEAN)
endmacro(enable_if_compiles)


# TODO have something similar for the linker

########################################################
#--     Set the compiler flags and build options     --#
########################################################

#
# extra flags for the specific compiler types:
#
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	set(LINALGWRAP_KNOWN_COMPILER TRUE)
	include(cmake/setup_compiler_flags_standard.cmake)
else()
	set(LINALGWRAP_KNOWN_COMPILER FALSE)
	message(WARNING "Untested compiler: ${CMAKE_CXX_COMPILER_ID}, you are on your own.")
	message(WARNING "Currently we only support clang and gnu compilers.")
	include(cmake/setup_compiler_flags_standard.cmake)
endif()
