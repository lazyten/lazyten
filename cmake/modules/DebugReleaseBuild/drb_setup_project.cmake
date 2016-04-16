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

macro(DRB_SETUP_PROJECT PROJECTNAME)
	# This macro sets up a project for a debug-and-release build.
	#
	# Sets the following variables:
	#       PROJECTNAME_DEFINITIONS  Default definitions for the project
	#                                in both debug and release network
	#       PROJECTNAME_DEFINITIONS_DEBUG  Variable for definitions 
	#                                      which are used only in DEBUG
	#                                      mode on top of the ones above
	#       PROJECTNAME_DEFINITIONS_RELEASE Variable for definitions
	#                                       which are used only in RELEASE
	#
	# then calls drb_setup_compiler_flags() to setup the compiler flags
	# CMAKE_CXX_FLAGS, CMAKE_SHARED_LINKER_FLAGS, CMAKE_CXX_FLAGS_DEBUG ...
	# see drb_setup_compiler_flags.cmake for details.
	#

	# check that drb has been initialised
	if(NOT DRB_INITIALISED)
		message(FATAL_ERROR "You have to call drb_init before any other DebugReleaseBuild function (here: drb_setup_project)")
	endif()

	#
	# Some default options
	# TODO: Promote them to cache
	#
	set(${PROJECTNAME}_DEFINITIONS "")	# empty list
	LIST(APPEND ${PROJECTNAME}_DEFINITIONS_DEBUG "DEBUG")
	LIST(APPEND ${PROJECTNAME}_DEFINITIONS_RELEASE "RELEASE")

	# setup the compilers
	drb_setup_compiler_flags()

	# note that the project has been set up for drb
	set(${PROJECTNAME}_HAS_DRB_SETUP TRUE)
endmacro(DRB_SETUP_PROJECT)
