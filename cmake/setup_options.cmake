# this file is greatly inspired by cmake/setup_cached_variables.cmake of deal.ii

#
# Mess with build type
#
set(CMAKE_BUILD_TYPE "DebugRelease"
	CACHE STRING "The type of build, options are Debug, Release and DebugRelease.")

if( NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Release" AND
		NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" AND
		NOT "${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease" )
	message(FATAL_ERROR "CMAKE_BUILD_TYPE must be one of Release, Debug or DebugRelease.")
endif()

#
# Parse the build type and generate a list of build types
# to perform
#
if (CMAKE_BUILD_TYPE MATCHES "Debug")
	LIST(APPEND LINALGWRAP_BUILD_TYPES "DEBUG")
endif()
if (CMAKE_BUILD_TYPE MATCHES "Release")
	LIST(APPEND LINALGWRAP_BUILD_TYPES "RELEASE")
endif()

#
# Hide options we explicitly do not use
#
set(LINALGWRAP_HIDDEN_OPTIONS
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

foreach(_flag ${LINALGWRAP_HIDDEN_OPTIONS})
	set(${_flag} ${${_flag}} CACHE INTERNAL "" FORCE)
	set(${_flag} "")
endforeach()

#
# Some default options
# TODO: Promote them to cache
#
set(LINALGWRAP_DEFINITIONS "")	# empty list
LIST(APPEND LINALGWRAP_DEFINITIONS_DEBUG "DEBUG")
LIST(APPEND LINALGWRAP_DEFINITIONS_RELEASE "RELEASE")

#
# Setup build type for examples
#
if (CMAKE_BUILD_TYPE MATCHES "Release")
	set(EXAMPLES_BUILD_TYPE "RELEASE")
else()
	set(EXAMPLES_BUILD_TYPE "DEBUG")
endif()

#
# Extra options
#
option(LINALGWRAP_ENABLE_TESTS "Build linalgwrap tests" ON)
option(LINALGWRAP_ENABLE_EXAMPLES "Build linalgwrap examples" ON)
option(LINALGWRAP_COMPILE_PEDANTIC "Use extremely pedantic compiler flags" OFF)
option(LINALGWRAP_WITH_EXTERNAL_RAPIDCHECK "Use an external rapidcheck library for the tests." OFF)
option(LINALGWRAP_WITH_EXTERNAL_CATCH "Use an external catch library for the tests." OFF)
