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

cmake_minimum_required(VERSION 3.0.0)

# Setup options, has to be done BEFORE project definition.
include(cmake/setup_options.cmake)

# enable the new version feature introduced in cmake 3.0.0
cmake_policy(SET CMP0048 NEW)
project(linalgwrap VERSION 0.1.0)

# enable testing of this project
enable_testing()

##########################################################################
# Setup compiler and dependencies:

# default default compiler flags:
include(cmake/setup_compiler_flags.cmake)

# setup hard dependencies
include(cmake/setup_dependencies.cmake)

# setup optional dependencies and features.
include(cmake/setup_optional.cmake)

##########################################################################
# Library and binary dir for this project.

# Setup the main include directory of this project.
set(LINALGWRAP_INCLUDE_DIR "${linalgwrap_SOURCE_DIR}/src")
include_directories(${LINALGWRAP_INCLUDE_DIR})

# Setup the main binary directory
set(LINALGWRAP_LIBRARY_DIR "${linalgwrap_BINARY_DIR}/src")
link_directories(${LINALGWRAP_LIBRARY_DIR})

# The name of the library targets in Debug and Release form
set(LINALGWRAP_DEBUG_TARGET "${CMAKE_PROJECT_NAME}.g")
set(LINALGWRAP_RELEASE_TARGET "${CMAKE_PROJECT_NAME}")

##########################################################################
# Pass on to subdirectories:

# The src directory does most of the work: Here the library is build and 
# installation is configured.
add_subdirectory(src)

# The external dependencies via submodules or similar
add_subdirectory(external)

# Add subdirectories for the test.
if(LINALGWRAP_ENABLE_TESTS) 
	add_subdirectory(tests)
endif()

# Add subdirectories for the examples.
if(LINALGWRAP_ENABLE_EXAMPLES)
	add_subdirectory(examples)
endif()
