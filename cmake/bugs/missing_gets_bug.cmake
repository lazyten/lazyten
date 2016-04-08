option(PROVIDE_DUMMY_GETS "Do we have the missing gets bug?" OFF)

#
# This bug is only in glibc 2.16 in combination with the libstdc++ of gcc 4.8, 
# which are both shipped with Ubuntu LTS 14.04.
#
# See http://gcc.gnu.org/bugzilla/show_bug.cgi?id=51785 for some details.
#

if (PROVIDE_DUMMY_GETS) 
FILE(WRITE "${linalgwrap_BINARY_DIR}/gets_dummy.hh"
"#pragma once
extern \"C\" char* gets (char* __s);")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -include ${linalgwrap_BINARY_DIR}/gets_dummy.hh")
endif()
