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

# Installs the cmake modules and cmake package information this project
# provides
#
# Requires the variable PackageModuleLocation to be set.

# Installing DebugReleaseBuild cmake module
install(DIRECTORY "${linalgwrap_SOURCE_DIR}/cmake/modules"
	DESTINATION ${PackageModuleLocation}
	COMPONENT devel
	FILES_MATCHING PATTERN "*.cmake"
)

# Write a basic version file for linalgwrap
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
	"${CMAKE_CURRENT_BINARY_DIR}/linalgwrapConfigVersion.cmake"
	COMPATIBILITY AnyNewerVersion
)

# Adjust a configure file
configure_file(cmake/linalgwrapConfig.cmake.in
	"${CMAKE_CURRENT_BINARY_DIR}/linalgwrapConfig.cmake"
	COPYONLY
)

# Set an export location:
install(FILES
	"${CMAKE_CURRENT_BINARY_DIR}/linalgwrapConfig.cmake"
	"${CMAKE_CURRENT_BINARY_DIR}/linalgwrapConfigVersion.cmake"
	DESTINATION "${PackageModuleLocation}/linalgwrap"
	COMPONENT devel
)

