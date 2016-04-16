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

# Set the module path
list(APPEND CMAKE_MODULE_PATH  "${CMAKE_CURRENT_LIST_DIR}/cmake/modules")

# Load the DebugReleaseBuild module for a simultaneous Debug and Release
# and release build and initialise it.
include(DebugReleaseBuild)
drb_init()

##########################################################################
# Global options

# Path where cmake wants the modules and config packages:
set(PackageModuleLocation share/cmake)

# setup global options (needs to be done before call to project())
option(LINALGWRAP_ENABLE_TESTS "Build linalgwrap tests" ON)
option(LINALGWRAP_ENABLE_EXAMPLES "Build linalgwrap examples" ON)

##########################################################################
# Project setup

# enable the new version feature introduced in cmake 3.0.0
cmake_policy(SET CMP0048 NEW)
project(linalgwrap VERSION 0.1.0)
drb_setup_project(linalgwrap)

# Global include directories for this project:
include_directories("${linalgwrap_SOURCE_DIR}/src")

# enable testing of this project
enable_testing()

##########################################################################
# Setup hard and optional dependencies and find components

# setup hard dependencies
include(cmake/setup_dependencies.cmake)

# setup optional dependencies and features.
include(cmake/setup_optional.cmake)

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

##########################################################################
# Finally install modules and packages:
include(cmake/install_modules_packages.cmake)
