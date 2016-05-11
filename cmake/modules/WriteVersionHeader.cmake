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

include(CMakeParseArguments)

function(WRITE_VERSION_HEADER _filename)
	# Write a header file which defines the macro contstants
	# ${NAME}_VERSION_MAJOR, ${NAME}_VERSION_MINOR and
	# ${NAME}_VERSION_PATCH according to the current 
	# project's version and name or by a name and version
	# supplied via the VERSION and NAME flags.


	set(options )
	set(oneValueArgs VERSION NAME)
	set(multiValueArgs )

	cmake_parse_arguments(WVH "${options}" "${oneValueArgs}" "${multiValueArgs}"  ${ARGN})

	if(WVH_UNPARSED_ARGUMENTS)
		message(FATAL_ERROR "Unknown keywords given to WRITE_VERSION_HEADER: \"${WVH_UNPARSED_ARGUMENTS}\"")
	endif()

	if("${WVH_VERSION}" STREQUAL "")
		if ("${PROJECT_VERSION}" STREQUAL "")
			message(FATAL_ERROR "No VERSION specified for WRITE_VERSION_HEADER()")
		else()
			set(WVH_VERSION "${PROJECT_VERSION}")
		endif()
	endif()

	if("${WVH_NAME}" STREQUAL "")
		if ("${PROJECT_NAME}" STREQUAL "")
			message(FATAL_ERROR "No NAME specified for WRITE_VERSION_HEADER()")
		else()
			set(WVH_NAME "${PROJECT_NAME}")
		endif()
	endif()

	# Split version into individual vars:
	string(REPLACE "." ";" VERSION_LIST ${WVH_VERSION})
	list(GET VERSION_LIST 0 WVH_MAJOR)
	list(GET VERSION_LIST 1 WVH_MINOR)
	list(GET VERSION_LIST 2 WVH_PATCH)

	# Now dump the header:

	file(WRITE "${_filename}"
"#pragma once
#define ${WVH_NAME}_VERSION_MAJOR ${WVH_MAJOR}
#define ${WVH_NAME}_VERSION_MINOR ${WVH_MINOR}
#define ${WVH_NAME}_VERSION_PATCH ${WVH_PATCH}")
endfunction()
