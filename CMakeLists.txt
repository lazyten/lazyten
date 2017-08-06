## ---------------------------------------------------------------------
##
## Copyright (C) 2016-17 by the lazyten authors
##
## This file is part of lazyten.
##
## lazyten is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## lazyten is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with lazyten. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

cmake_minimum_required(VERSION 3.0.0)
cmake_policy(VERSION 3.0.0)

# Load the DebugReleaseBuild module for a simultaneous Debug and Release
# and release build and initialise it.
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/IncludeKrimsCmakeModule.cmake")
include_krims_cmake_module(DebugReleaseBuild)
drb_init()

##########################################################################
# Project setup

project(lazyten VERSION 0.3.0)
message("\n# Configuring ${PROJECT_NAME}\n#")

# Global include directories for this project:
include_directories("${PROJECT_SOURCE_DIR}/src")  # for default sources
include_directories("${PROJECT_BINARY_DIR}/src")  # for generated sources

# enable testing of this project
enable_testing()

##########################################################################
# Global options

# Path where cmake wants the modules and config packages:
set(PackageModuleLocation share/cmake)

# setup global options
include_krims_cmake_module(GlobalOption)
global_option(ENABLE_DOCUMENTATION "Build documentation"         OFF )
global_option(ENABLE_EXAMPLES      "Build example executables"   ON  )
global_option(ENABLE_TESTS         "Build unit test executables" ON  )

##########################################################################
# Setup hard and optional dependencies and find components

# Determine and setup compiler flags, enforcing C++11
drb_setup_compiler_flags(11)

# Export compile commands if cmake supports it.
if (CMAKE_VERSION VERSION_GREATER 3.5.0)
	set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

# setup hard dependencies
include(cmake/setup_dependencies.cmake)

# setup optional dependencies and features.
include(cmake/setup_optional.cmake)

##########################################################################
# Pass on to subdirectories:

# The src directory does most of the work: Here the library is build and
# installation is configured.
add_subdirectory(src)
set(LAW_SUBDIRS src)

# Add subdirectories for the testing.
if(LAZYTEN_ENABLE_TESTS)
	add_subdirectory(tests)
	set(LAW_SUBDIRS ${LAW_SUBDIRS} tests)
endif()

# Add subdirectories for the examples.
if(LAZYTEN_ENABLE_EXAMPLES)
	add_subdirectory(examples)
	set(LAW_SUBDIRS ${LAW_SUBDIRS} examples)
endif()

if(LAZYTEN_ENABLE_DOCUMENTATION)
	add_subdirectory(doc)
endif()

##########################################################################
# Setup extra targets for clang tools
include_krims_cmake_module(SetupClangTargets)
add_available_clang_targets_for(${PROJECT_NAME}
	DIRECTORIES ${LAW_SUBDIRS})

# Finally install package:
include(cmake/install_package.cmake)

# and print a closing message
message("#\n# Configuring ${PROJECT_NAME} done\n")
