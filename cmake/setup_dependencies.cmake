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

# sets these things
#
#       LAZYTEN_DEPENDENCIES			everyone needs these libraries
#       LAZYTEN_DEPENDENCIES_DEBUG		debug mode needs these extras
#       LAZYTEN_DEPENDENCIES_RELEASE		release mode needs these extras
#       LAZYTEN_DEPENDENCIES_TEST		tests need these extra libraries
#

####################
#-- Empty it all --#
####################
set(LAZYTEN_DEPENDENCIES "")
set(LAZYTEN_DEPENDENCIES_DEBUG "")
set(LAZYTEN_DEPENDENCIES_RELEASE "")
set(LAZYTEN_DEPENDENCIES_TEST "")
include_krims_cmake_module(ProjectFeatures)

############################
#-- rapidcheck and catch --#
############################
if (LAZYTEN_ENABLE_TESTS)
	# We need to setup rapidcheck and catch for the tests:
	include(cmake/findRapidcheck.cmake)
	set(LAZYTEN_DEPENDENCIES_TEST ${LAZYTEN_DEPENDENCIES_TEST} ${rapidcheck_TARGET})

	include(cmake/findCatch.cmake)
	set(LAZYTEN_DEPENDENCIES_TEST ${LAZYTEN_DEPENDENCIES_TEST} ${catch_TARGET})
endif()

#############
#-- krims --#
#############
# Find at least version 0.0.0
set(KRIMS_VERSION 0.0.0)
include(cmake/findKrims.cmake)

foreach (build ${DRB_BUILD_TYPES})
	set(LAZYTEN_DEPENDENCIES_${build} ${LAZYTEN_DEPENDENCIES_${build}} ${krims_${build}_TARGET})
endforeach()

#########################
#--  LAPACK and BLAS  --#
#########################
set(BLAS_VENDOR "All" CACHE STRING "The BLAS and LAPACK vendor to use \
(e.g. Intel, ATLAS, OpenBLAS, ACML, Apple)")
set(BLA_VENDOR ${BLAS_VENDOR})

find_package(LAPACK REQUIRED)
set(LAZYTEN_DEPENDENCIES ${LAZYTEN_DEPENDENCIES} ${LAPACK_LIBRARIES})
enable_feature(lapack)

unset(BLA_VENDOR)

#################
#-- armadillo --#
#################
find_package(Armadillo 4.000 REQUIRED)
set(LAZYTEN_DEPENDENCIES ${LAZYTEN_DEPENDENCIES} ${ARMADILLO_LIBRARIES})
include_directories(${ARMADILLO_INCLUDE_DIRS})
enable_feature(armadillo)
