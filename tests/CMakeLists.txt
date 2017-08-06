## ---------------------------------------------------------------------
##
## Copyright (C) 2016-17 by the lazyten authors
##
## This file is part of lazyten.
##
## lazyten is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published
## by the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## lazyten is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with lazyten. If not, see <http://www.gnu.org/licenses/>.
##
## ---------------------------------------------------------------------

#
# The sources for the test executable
#
set(LAZYTEN_TESTS_SOURCES
	# Low-level interface classes
	OperatorPropertiesTests.cc
	versionTest.cc

	# Vector and Matrix LA interfaces
# TODO	IteratorVectorTests.cc
	MultiVectorTests.cc
	BuiltinVectorTests.cc

	# Armadillo vector/matrix
	ArmadilloMatrixTests.cc
	ArmadilloVectorTests.cc

	# Helper classes
	TypeUtilsTests.cc
	RandomTests.cc
	orthoTests.cc

	# Lazy matrices
	LazyMatrix_i_Tests.cc
	LazyMatrixWrapperTests.cc
	LazyMatrixProductTests.cc
	LazyMatrixSumTests.cc

	# Proxies and other matrix functionality
	TransposeProxyTests.cc
	transTests.cc
	inverseTests.cc

	# Lazy matrix implementations
	BlockDiagonalMatrixTests.cc
	DiagonalMatrixTests.cc

	# Eigensolver
	ArpackEigensolverTests.cc
	ArmadilloEigensolverTests.cc
	LapackEigensolverTests.cc
	eigensystemTests.cc

	# linear solver
	solveTests.cc

	# main of the test suite
	main.cc
)

# TODO have a benchmark test!

drb_add_executable(lazyten_tests
	DBGSUFFIX  "_DEBUG"
	RELSUFFIX  "_RELEASE" 
	FILES ${LAZYTEN_TESTS_SOURCES}
)

drb_target_compile_definitions(ALL     lazyten_tests
	PRIVATE ${LAZYTEN_DEFINITIONS})
drb_target_compile_definitions(DEBUG   lazyten_tests
	PRIVATE ${LAZYTEN_DEFINITIONS_DEBUG})
drb_target_compile_definitions(RELEASE lazyten_tests
	PRIVATE ${LAZYTEN_DEFINITIONS_RELEASE})

drb_target_link_libraries(ALL     lazyten_tests
	${LAZYTEN_DEPENDENCIES}   ${LAZYTEN_DEPENDENCIES_TEST})
drb_target_link_libraries(DEBUG   lazyten_tests
	${lazyten_DEBUG_TARGET}   ${LAZYTEN_DEPENDENCIES_DEBUG})
drb_target_link_libraries(RELEASE lazyten_tests
	${lazyten_RELEASE_TARGET} ${LAZYTEN_DEPENDENCIES_RELEASE})

# Add the test executables
FOREACH(target ${lazyten_tests_TARGETS})
	add_test(NAME ${target} COMMAND ${target})
ENDFOREACH()

