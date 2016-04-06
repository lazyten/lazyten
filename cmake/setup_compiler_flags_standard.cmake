# Compiler flags for the default compilers, i.e. gnu or clang

#######################
#--  Version check  --#
#######################

set(GCC_MIN_VERSION "4.9")
set(CLANG_MIN_VERSION "3.5")

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
	if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS GCC_MIN_VERSION)
		message(FATAL_ERROR "gcc version ${GCC_MIN_VERSION} or higher is required for compilation (C++14 support)")
	endif()
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS CLANG_MIN_VERSION)
		message(FATAL_ERROR "clang version ${CLANG_MIN_VERSION} or higher is required for compilation (C++14 support)")
	endif()
endif() 

###############
#--  Common --#
###############
# We definitely need c++14
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

#
# Warning policy 
#
# The warning policy we use here is based on a very valuable
# post by a clang developer on stack exchange
# https://programmers.stackexchange.com/questions/122608#124574

# Show high confidence warning
enable_if_compiles(CMAKE_CXX_FLAGS  "-Wall")

# Show valuable extra warnings
enable_if_compiles(CMAKE_CXX_FLAGS "-Wextra")

# Turn on warnings about language extensions
enable_if_compiles(CMAKE_CXX_FLAGS "-pedantic")

# But silence some rather annoying warnings
enable_if_compiles(CMAKE_CXX_FLAGS "-Wno-unused-macros")
enable_if_compiles(CMAKE_CXX_FLAGS "-Wno-unused-parameter")

# Make warnings errors, such that we cannot ignore them
enable_if_compiles(CMAKE_CXX_FLAGS "-Werror")

##############
#--  Debug --#
##############
#
# Extra stuff for debug:
#
if (CMAKE_BUILD_TYPE MATCHES "Debug")
	enable_if_compiles(CMAKE_CXX_FLAGS_DEBUG "-O0")
	enable_if_compiles(CMAKE_CXX_FLAGS_DEBUG "-Og")

	# Common linker flags for all of debug:
	set(COMMON_LINKER_FLAGS_DEBUG "${COMMON_LINKER_FLAGS_DEBUG} -g")
	set(COMMON_LINKER_FLAGS_DEBUG "${COMMON_LINKER_FLAGS_DEBUG} -ggdb")

	set(CMAKE_STATIC_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${COMMON_LINKER_FLAGS_DEBUG}")
	set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${COMMON_LINKER_FLAGS_DEBUG}")
	set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${COMMON_LINKER_FLAGS_DEBUG}")

	unset(COMMON_LINKER_FLAGS_DEBUG)
endif()

################
#--  Release --#
################
#
# Extra stuff for release:
#
if (CMAKE_BUILD_TYPE MATCHES "Release")
	# nothing atm
endif()
