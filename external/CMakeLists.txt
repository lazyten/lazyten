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

#
# rapidcheck
#
if (LINALGWRAP_ENABLE_TESTS AND NOT LINALGWRAP_WITH_EXTERNAL_RAPIDCHECK)
	#
	# Add subdirectory
	#
	
	# Check wether rapidcheck is already there.
	if(NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/rapidcheck/CMakeLists.txt")
		message(FATAL_ERROR "Could not include rapidcheck from submodule. Try \"git submodule update --init --recursive\"")
	endif()

	# Set flags such that both rapidcheck tests and catch are built.
	set(RC_ENABLE_TESTS ON CACHE BOOL "Build RapidCheck tests")
	set(RC_ENABLE_CATCH ON CACHE BOOL "Build Catch.hpp support")

	# Add the rapidcheck subdirectory
	add_subdirectory(rapidcheck)

	#
	# Adjust compiler flags
	#

	# Set some compiler options to make rapidcheck compile
	enable_if_compiles(LINALGWRAP_RAPIDCHECK_EXTRA_FLAGS "-Wno-gnu-zero-variadic-macro-arguments")
	#enable_if_compiles(LINALGWRAP_RAPIDCHECK_EXTRA_FLAGS "-Wno-error=pedantic")
	enable_if_compiles(LINALGWRAP_RAPIDCHECK_EXTRA_FLAGS "-Wno-pedantic")
	set(LINALGWRAP_RAPIDCHECK_EXTRA_FLAGS "${LINALGWRAP_RAPIDCHECK_EXTRA_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")

	# Activate extra flags for all potential rapidcheck targets.
	foreach(_potential_target
			rapidcheck
			rapidcheck_tests
			rapidcheck_test_utils
			rapidcheck_boost_tests
			rapidcheck_gmock_tests
			boost_test_integration
			database
			counter
			gtest_integration
			classify
			mapparser
	)
		if(TARGET ${_potential_target})
			set_target_properties(${_potential_target} PROPERTIES COMPILE_FLAGS "${LINALGWRAP_RAPIDCHECK_EXTRA_FLAGS}")
		endif()
	endforeach()

	# unset the extra flags
	unset(LINALGWRAP_RAPIDCHECK_EXTRA_FLAGS)

	#
	# Install rapidcheck
	#
	if(BUILD_SHARED_LIBS)
		install(FILES "${rapidcheck_BINARY_DIR}/librapidcheck.so"
			DESTINATION lib
			COMPONENT library-rapidcheck
		)
	else()
		install(FILES "${rapidcheck_BINARY_DIR}/librapidcheck.a"
			DESTINATION lib
			COMPONENT library-rapidcheck
		)
	endif()

	# Export definitions:
	export(TARGETS rapidcheck
		NAMESPACE Upstream::
		FILE "${CMAKE_CURRENT_BINARY_DIR}/rapidcheckTargets.cmake"
	)

	# TODO this does not work, since we cannot define the export group
	#      for the rapidcheck target.
	# and install them appropriately:
	#install(EXPORT rapidcheck
	#	FILE "rapidcheckTargets.cmake"
	#	NAMESPACE Upstream::
	#	DESTINATION "${PackageModuleLocation}/rapidcheck"
	#)

	# Install catch
	install(FILES ${rapidcheck_SOURCE_DIR}/ext/catch/single_include/catch.hpp DESTINATION include)
endif()

if(LINALGWRAP_ENABLE_TESTS AND NOT LINALGWRAP_WITH_EXTERNAL_CATCH)
	# TODO download catch somewhere sensible and use that one
endif()
