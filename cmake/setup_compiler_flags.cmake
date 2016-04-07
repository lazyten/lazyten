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
	CHECK_CXX_COMPILER_FLAG(${FLAG} LINALGWRAP_HAVE_FLAG_${FLAG})
	if (LINALGWRAP_HAVE_FLAG_${FLAG})
		set(${VARIABLE} "${${VARIABLE}} ${FLAG}")
	endif()
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
