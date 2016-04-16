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

macro(drb_init)
	# Setup a DebugReleaseBuild environment. Does the following:
	# 	1. It sets up the cache variable CMAKE_BUILD_TYPE
	#	   (Which can now take the values Debug, Release and DebugRelease
	#
	#	2. It checks wheather this variable has an appropriate value.
	#
	#	3. It hides all those CMAKE options that we won't use in our builds
	#	   from the user.
	#
	#	4. Sets the variable DRB_BUILD_TYPES to the list of build types to perform
	#	   (i.e the list contains DEBUG and/or RELEASE)
	#	   The idea is that you can use a FOREACH loop to setup your targets.
	#

	# Setup CMAKE_BUILD_TYPE and check for value
	set(CMAKE_BUILD_TYPE "DebugRelease"
		CACHE STRING "The type of build, options are Debug, Release and DebugRelease.")

	if( NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Release" AND
			NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" AND
			NOT "${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease" )
		message(FATAL_ERROR "CMAKE_BUILD_TYPE must be one of Release, Debug or DebugRelease.")
	endif()

	# Hide options we explicitly do not use
	set(DRB_HIDDEN_OPTIONS
		CMAKE_CXX_FLAGS_MINSIZEREL
		CMAKE_CXX_FLAGS_RELWITHDEBINFO
		CMAKE_C_FLAGS_MINSIZEREL
		CMAKE_C_FLAGS_RELWITHDEBINFO
		CMAKE_Fortran_FLAGS_MINSIZEREL
		CMAKE_Fortran_FLAGS_RELWITHDEBINFO
		CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL
		CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO
		CMAKE_STATIC_LINKER_FLAGS_MINSIZEREL
		CMAKE_STATIC_LINKER_FLAGS_RELWITHDEBINFO
		CMAKE_EXE_LINKER_FLAGS_MINSIZEREL
		CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO
	)

	foreach(_flag ${DRB_HIDDEN_OPTIONS})
		set(${_flag} ${${_flag}} CACHE INTERNAL "" FORCE)
		set(${_flag} "")
	endforeach()

	# Parse the build type and generate a list of build types
	# to perform
	if (CMAKE_BUILD_TYPE MATCHES "Debug")
		LIST(APPEND DRB_BUILD_TYPES "DEBUG")
	endif()
	if (CMAKE_BUILD_TYPE MATCHES "Release")
		LIST(APPEND DRB_BUILD_TYPES "RELEASE")
	endif()

	set(DRB_INITIALISED ON)
endmacro(drb_init)
