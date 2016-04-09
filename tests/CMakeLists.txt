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
# The sources for the test executable
#
set(LINALGWRAP_TESTS_SOURCES
	NumComp.cc
	ExceptionTests.cc
	RangeTests.cc
	SubscriptionTests.cc
	ArmadilloMatrixTests.cc
	ArmadilloVectorTests.cc
	ParameterMapTests.cc
	LazyMatrixWrapperTests.cc
	LazyMatrixProductTests.cc
	LazyMatrixSumTests.cc
	LazyMatrix_i_Tests.cc
	ScaleViewTests.cc
	TransposeViewTests.cc
	BlockViewTests.cc
	BlockDiagonalMatrixTests.cc
	main.cc
)

# TODO have a benchmark test!

# 
# For each build type compile seperate tests.
#
FOREACH(build ${LINALGWRAP_BUILD_TYPES})
	# add the test executable
	add_executable(linalgwrap_tests_${build} ${LINALGWRAP_TESTS_SOURCES})

	set_target_properties(
		linalgwrap_tests_${build}
		PROPERTIES
		LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_${build}}"
		LINKER_LANGUAGE "CXX"
		COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${build}}"
		COMPILE_DEFINITIONS "${LINALGWRAP_DEFINITIONS};${LINALGWRAP_DEFINITIONS_${build}}"
	)

	# link it with the appropriate linalgwrap library target
	target_link_libraries(linalgwrap_tests_${build}
		${LINALGWRAP_${build}_TARGET}
	       	${LINALGWRAP_DEPENDENCIES}
	       	${LINALGWRAP_${build}_DEPENDENCIES}
	       	${LINALGWRAP_DEPENDENCIES_TEST}
	)

	# add it to list of test cases
	add_test(NAME linalgwrap_tests_${build} COMMAND linalgwrap_tests_${build})
ENDFOREACH()
		
