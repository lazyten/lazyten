# Compiler flags for the default compilers, i.e. gnu or clang


###############
#--  Common --#
###############
# We definitely need c++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Show high confidence warning
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

# Make them errors such that we cannot ignore them
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")

# Show valuable extra warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")

# Turn on warnings about language extensions
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic")

# But silence some rather annoying warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-macros")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")

#
# Warning policy 
#
# The warning policy we use here is based on a very valuable
# post by a clang developer on stack exchange
# https://programmers.stackexchange.com/questions/122608#124574

if (LINALGWRAP_COMPILE_PEDANTIC)
else()
endif()


##############
#--  Debug --#
##############
#
# Extra stuff for debug:
#
if (CMAKE_BUILD_TYPE MATCHES "Debug")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

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
