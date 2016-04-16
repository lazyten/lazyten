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

macro(drb_get_matching_build_types WHICH MATCHLIST)
	# Macro which parses the DRB_BUILD_TYPES list and stores only those
	# entries inside MATCHLIST, which match the WHICH statement.
	#
	# if WHICH==ALL all entries of DRB_BUILD_TYPES are just copied to
	# MATCHLIST.
	#
	# if WHICH==DEBUG, MATCHLIST may only contain debug, but may also be
	# empty (if DRB_BUILD_TYPES=="RELEASE")
	#
	# similarly for WHICH==RELEASE

	if(NOT DRB_INITIALISED)
		message(FATAL_ERROR "You have to call drb_init before any other DebugReleaseBuild function.")
	endif()

	# Set MATCHLIST to be a copy.
	set(MATCHLIST "${DRB_BUILD_TYPES}")

	if (${WHICH} STREQUAL "DEBUG")
		# remove RELEASE:
		list(REMOVE_ITEM MATCHLIST "RELEASE")
	elseif(${WHICH} STREQUAL "RELEASE")
		# remove DEBUG:
		list(REMOVE_ITEM MATCHLIST "DEBUG")
	elseif(${WHICH} STREQUAL "ALL")
		# remove nothing
	else()
		message(FATAL_ERROR "Unknown value to WHICH: ${WHICH}. Should be ALL, DEBUG or RELEASE")
	endif()
endmacro(drb_get_matching_build_types)
