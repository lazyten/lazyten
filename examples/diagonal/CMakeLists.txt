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
# The sources for this example executable
#
set(EXAMPLE_DIAGONAL_SOURCES
	diagonal.cc)

# Compile and link:
add_executable(diagonal ${EXAMPLE_DIAGONAL_SOURCES})

set_target_properties(
	diagonal
	PROPERTIES
	LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS_${EXAMPLES_BUILD_TYPE}}" #${CMAKE_EXE_LINKER_FLAGS} are added
	LINKER_LANGUAGE "CXX"
	COMPILE_FLAGS "${CMAKE_CXX_FLAGS_${EXAMPLES_BUILD_TYPE}}" #${CMAKE_CXX_FLAGS} are added automaticaly
)

# setup definitions.
target_compile_definitions(diagonal PRIVATE "${LINALGWRAP_DEFINITIONS};${LINALGWRAP_DEFINITIONS_${EXAMPLES_BUILD_TYPE}}")

# link it with the appropriate linalgwrap library target
target_link_libraries(diagonal 
	${linalgwrap_${EXAMPLES_BUILD_TYPE}_TARGET}
	${LINALGWRAP_DEPENDENCIES}
	${LINALGWRAP_DEPENDENCIES_${EXAMPLES_BUILD_TYPE}}
)
