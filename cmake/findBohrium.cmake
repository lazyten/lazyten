## ---------------------------------------------------------------------
##
## Copyright (C) 2017 by the linalgwrap authors
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

# Finds and sets up Bohrium under the target names stored in
#      bohrium_TARGET
# such that just linking against it as a dependency does everything
# automatically.
#
# In case the krims library is not found and AUTOCHECKOUT_MISSING_LIBS is set to
# on, krims is automatically checked out and built.
# Otherwise a fatal error is produced.
#

#
# Options and properties required
#
#option(AUTOCHECKOUT_MISSING_REPOS "Automatically checkout missing repositories" OFF)

# TODO This file is pretty much only a stub
# 
# TODO Improve search and use external repos
#      Check version
# TODO mark_as_advanced on some variables

# TODO Checkout via ./get_bohrium.sh but configure and build via external project.


#
# -------
#

if (TARGET "${bohrium_TARGET}")
	message(STATUS "Found bohrium targets, assume bohrium already configured for build.")
	return()
endif()

find_library(
	BOHRIUM_LIBRARY
	NAMES bh
	PATHS
	$ENV{BOHRIUM_DIR}
	$ENV{HOME}/.local
	${BOHRIUM_DIR}
	PATH_SUFFIXES lib64 lib
	DOC "Bohrium library directory"
)
if(${BOHRIUM_LIBRARY} MATCHES "-NOTFOUND")
	message(FATAL_ERROR "Bohrium library not found. Try environment variable BOHRIUM_DIR.")
endif()
get_filename_component(BOHRIUM_LIBDIR "${BOHRIUM_LIBRARY}" DIRECTORY)

find_path(
	BOHRIUM_INCLUDE bohrium/bh_memory.h
	HINTS ${BOHRIUM_LIBDIR}/../include
	DOC "Bohrium include directory"
)
if(${BOHRIUM_INCLUDE} MATCHES "-NOTFOUND")
	message(FATAL_ERROR "Bohrium include not found. Try environment variable BOHRIUM_DIR.")
endif()
message(STATUS "Found bohrium at ${BOHRIUM_LIBRARY}")

set(bohrium_TARGET bohrium)
add_library(${bohrium_TARGET} INTERFACE IMPORTED GLOBAL)
set_target_properties(${bohrium_TARGET} PROPERTIES
	INTERFACE_LINK_LIBRARIES "${BOHRIUM_LIBRARY};${BOHRIUM_LIBDIR}/libbhxx.so"
	INTERFACE_INCLUDE_DIRECTORIES "${BOHRIUM_INCLUDE}/bohrium;${BOHRIUM_INCLUDE}"
)
MESSAGE(WARNING "This is a bit of a hack.")

# TODO Dummy for now ... better do this via a module
file(WRITE "${PROJECT_BINARY_DIR}/ycm_extra_includes.yaml"
"---
- \"${BOHRIUM_INCLUDE}\"
- \"${BOHRIUM_INCLUDE}/bohrium\"
...")


