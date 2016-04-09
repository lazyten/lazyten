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
	# compile the rapidcheck we have in the submodule.

	# TODO Doing it this way is a hell of a mess ...
	#      Is there no better way to do this?

	# check it out if it is not already there:
	if(NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/rapidcheck/CMakeLists.txt")
		message(FATAL_ERROR "Could not include rapidcheck from submodule. Try \"git submodule update --init --recursive\"")
	endif()

	# make sure rapidcheck gets checked as well ...
	set(RC_ENABLE_TESTS "ON")
	set(RC_ENABLE_CATCH "ON")

	# Add the rapidcheck stuff
	add_subdirectory(rapidcheck)

	# Set some compiler options to make rapidcheck compile
	set(LINALGWRAP_RAPIDCHECK_EXTRA_FLAGS "${LINALGWRAP_RAPIDCHECK_EXTRA_FLAGS} -Wno-gnu-zero-variadic-macro-arguments")
	set(LINALGWRAP_RAPIDCHECK_EXTRA_FLAGS "${LINALGWRAP_RAPIDCHECK_EXTRA_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")

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

	# This is a hell of a hack to install the compiled library
	if(BUILD_SHARED_LIBS)
		install(FILES "${rapidcheck_BINARY_DIR}/librapidcheck.so" DESTINATION lib)
	else()
		install(FILES "${rapidcheck_BINARY_DIR}/librapidcheck.a" DESTINATION lib)
	endif()

	# Install catch
	install(FILES ${rapidcheck_SOURCE_DIR}/ext/catch/single_include/catch.hpp DESTINATION include/${CMAKE_PROJECT_NAME})
endif()

if(LINALGWRAP_ENABLE_TESTS AND NOT LINALGWRAP_WITH_EXTERNAL_CATCH)
	# TODO download catch somewhere sensible and use that one
endif()
