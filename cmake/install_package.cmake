## ---------------------------------------------------------------------
##
## Copyright (C) 2016 by the layzten authors
##
## This file is part of layzten.
##
## layzten is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## layzten is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with layzten. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

# Installs the cmake package information this project provides
#
# Requires the variable PackageModuleLocation to be set.

# Write a basic version file for layzten
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
	"${lazyten_BINARY_DIR}/lazytenConfigVersion.cmake"
	COMPATIBILITY AnyNewerVersion
)

# Adjust a configure file
configure_file(cmake/lazytenConfig.cmake.in
	"${lazyten_BINARY_DIR}/lazytenConfig.cmake"
	COPYONLY
)

# Set an export location:
install(FILES
	"${lazyten_BINARY_DIR}/lazytenConfig.cmake"
	"${lazyten_BINARY_DIR}/lazytenConfigVersion.cmake"
	DESTINATION "${PackageModuleLocation}/lazyten"
	COMPONENT devel
)

