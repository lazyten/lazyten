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

# Check we have at least the required version:
cmake_minimum_required(VERSION 3.0.0)

# Set basedir for this module
set(DRB_DIR "${CMAKE_CURRENT_LIST_DIR}/DebugReleaseBuild")

# include other files:
include("${DRB_DIR}/drb_init.cmake")
include("${DRB_DIR}/drb_utils.cmake")
include("${DRB_DIR}/drb_setup_compiler_flags.cmake")
include("${DRB_DIR}/drb_setup_project.cmake")
include("${DRB_DIR}/drb_targets.cmake")
