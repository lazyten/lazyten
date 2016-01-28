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
#	LINALGWRAP_DEFINITIONS			Extra definitions for all builds
#	LINALGWRAP_DEFINITIONS_DEBUG		Extra definitions for debug builds
# 	LINALGWRAP_DEFINITIONS_RELEASE		Extra definitions for release builds

######################
#--     Macros     --#
######################


message(AUTHOR_WARNING "TODO: Have an enable_if kind of adding mechanism")

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
