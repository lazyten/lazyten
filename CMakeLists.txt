cmake_minimum_required(VERSION 2.8.8)

if (CMAKE_VERSION VERSION_LESS "3.0.0")
	# work around the new version feature introduced
	# for cmake 3.0.0
	project(linalgwrap)
	set(PROJECT_VERSION "0.0.0")
else()
	# enable the new feature:
	cmake_policy(SET CMP0048 NEW)
	project(linalgwrap VERSION 0.0.0)
endif()

# enable testing of this project
enable_testing()

# Setup options
include(cmake/setup_options.cmake)

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
