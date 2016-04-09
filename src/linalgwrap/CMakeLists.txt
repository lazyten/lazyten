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

#
# set the source files.
#
set(LINALGWRAP_SOURCES
	ParameterMap.cc
	ExceptionBase.cc
	ExceptionSystem.cc
	version.cc
)	

#
# add the Release and Debug targets and set their properties
#
FOREACH(build ${LINALGWRAP_BUILD_TYPES})
	# add the library target
	add_library(${LINALGWRAP_${build}_TARGET} ${LINALGWRAP_SOURCES})

	if (BUILD_SHARED_LIBS)
		set(LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS_${build}}")
	else()
		set(LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} ${CMAKE_STATIC_LINKER_FLAGS_${build}}")
	endif()

	# set the compilation properties
	set_target_properties(
		"${LINALGWRAP_${build}_TARGET}"
		PROPERTIES
		# version for the so library:
		SOVERSION "${PROJECT_VERSION}"
		LINK_FLAGS "${LINKER_FLAGS}"
		LINKER_LANGUAGE "CXX"
		COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${build}}"
		COMPILE_DEFINITIONS "${LINALGWRAP_DEFINITIONS};${LINALGWRAP_DEFINITIONS_${build}}"
	)

	# link in the dependencies:
	if(NOT "${LINALGWRAP_DEPENDENCIES_${build}} ${LINALGWRAP_DEPENDENCIES}" STREQUAL " ")
		# link, but only if there is anything to link
		target_link_libraries(
			${LINALGWRAP_${build}_TARGET}
			${LINALGWRAP_DEPENDENCIES_${build}}
			${LINALGWRAP_DEPENDENCIES}
		)
	endif()

# 	TODO this export stuff does not work
#	FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/cmake)	
#	export(TARGETS ${LINALGWRAP_${build}_TARGET} 
#		FILE ${CMAKE_BINARY_DIR}/cmake/${CMAKE_PROJECT_NAME}_libraries.cmake
#		APPEND
#	)


	# set the installation property for this target:
	install(TARGETS ${LINALGWRAP_${build}_TARGET}
		DESTINATION lib 
#		EXPORT ${CMAKE_PROJECT_NAME}_libraries
		COMPONENT library
	)
ENDFOREACH()

# install the includes:
install(DIRECTORY . DESTINATION "include/${CMAKE_PROJECT_NAME}"
	FILES_MATCHING PATTERN "*.hh")

## install the cmake configuration:
#install(EXPORT ${CMAKE_PROJECT_NAME}_libraries
#	DESTINATION "include/${CMAKE_PROJECT_NAME}/cmake" 
#	COMPONENT library)
