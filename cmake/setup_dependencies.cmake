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

# sets these things
#
# 	LINALGWRAP_DEPENDENCIES			everyone needs these libraries
# 	LINALGWRAP_DEPENDENCIES_DEBUG		debug mode needs these extras
# 	LINALGWRAP_DEPENDENCIES_RELEASE		release mode needs these extras
# 	LINALGWRAP_DEPENDENCIES_TEST		tests need these extra libraries
#
#       LINALGWRAP_DEFINITIONS			definitions for all compilation
#       LINALGWRAP_DEFINITIONS_DEBUG		definitions for debug mode
#       LINALGWRAP_DEFINITIONS_RELEASE		definitions for release mode
#       

####################
#-- Empty it all --#
####################
set(LINALGWRAP_DEPENDENCIES "")
set(LINALGWRAP_DEPENDENCIES_DEBUG "")
set(LINALGWRAP_DEPENDENCIES_RELEASE "")
set(LINALGWRAP_DEPENDENCIES_TEST "")
set(LINALGWRAP_DEFINITIONS "")
set(LINALGWRAP_DEFINITIONS_DEBUG "")
set(LINALGWRAP_DEFINITIONS_RELEASE "")

##################
#-- rapidcheck --#
##################

if (LINALGWRAP_WITH_EXTERNAL_RAPIDCHECK AND LINALGWRAP_ENABLE_TESTS)
	# We need to compile rapidcheck further down in the externals
	# add it to the dependencies right here:
	set(LINALGWRAP_DEPENDENCIES_TEST ${LINALGWRAP_DEPENDENCIES_TEST} rapidcheck)
elseif(LINALGWRAP_ENABLE_TESTS)
	message(FATAL_ERROR "Currently only the build with rapidcheck as a submodule is supported")
	#TODO change this
endif()

##############
#-- catch  --#
##############

if (LINALGWRAP_WITH_EXTERNAL_CATCH AND LINALGWRAP_ENABLE_TESTS)
	# We need to setup catch further down in the externals
	# so add it to the dependencies
	set(LINALGWRAP_DEPENDENCIES_TEST ${LINALGWRAP_DEPENDENCIES_TEST} law_catch)
elseif(LINALGWRAP_ENABLE_TESTS)
	# Note: We need the target to have the name law_catch and not just "catch"
	#       since rapidcheck uses this target name internally.
	add_library(law_catch INTERFACE)

	find_path(catch_INCLUDE_DIR catch.hpp
		PATHS
		$ENV{catch_INCLUDE_DIR}
		DOC "catch include directory"
	)

	if ("${catch_INCLUDE_DIR}" STREQUAL "catch_INCLUDE_DIR-NOTFOUND")
		message(FATAL_ERROR "Could not find catch include directory. 
Either disable LINALGWRAP_WITH_EXTERNAL_CATCH or provide a hint where \
the catch include file can be found via the environment variable \
catch_INCLUDE_DIR.")
	endif()

	message(STATUS "Found system-provided catch at ${catch_INCLUDE_DIR}/catch.hpp")

	target_include_directories(law_catch INTERFACE ${catch_INCLUDE_DIR})
	set(GSCF_DEPENDENCIES_TEST  ${LINALGWRAP_DEPENDENCIES_TEST} law_catch)
endif()



#################
#-- armadillo --#
#################
include(FindArmadillo)
find_package(Armadillo REQUIRED)

# We need at least 4.000
if(ARMADILLO_VERSION_STRING VERSION_LESS "4.000")
	message(FATAL_ERROR "Armadillo version is too old. Expect at least version ${ARMADILLO_MIN_VERSION}.")
endif()

# add to general dependencies and include string
set(LINALGWRAP_DEPENDENCIES ${LINALGWRAP_DEPENDENCIES} ${ARMADILLO_LIBRARIES})
include_directories(${ARMADILLO_INCLUDE_DIRS})

# enable armadillo-dependant code:
LIST(APPEND LINALGWRAP_DEFINITIONS "LINALGWRAP_HAVE_ARMADILLO")
