## ---------------------------------------------------------------------
##
## Copyright (C) 2016-17 by the linalgwrap authors
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
	# Low-level interface classes
	OperatorPropertiesTests.cc

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

drb_add_executable(linalgwrap_tests
	DBGSUFFIX  "_DEBUG"
	RELSUFFIX  "_RELEASE" 
	FILES ${LINALGWRAP_TESTS_SOURCES}
)

drb_target_compile_definitions(ALL     linalgwrap_tests 
	PRIVATE ${LINALGWRAP_DEFINITIONS})
drb_target_compile_definitions(DEBUG   linalgwrap_tests 
	PRIVATE ${LINALGWRAP_DEFINITIONS_DEBUG})
drb_target_compile_definitions(RELEASE linalgwrap_tests 
	PRIVATE ${LINALGWRAP_DEFINITIONS_RELEASE})

drb_target_link_libraries(ALL     linalgwrap_tests 
	${LINALGWRAP_DEPENDENCIES} ${LINALGWRAP_DEPENDENCIES_TEST})
drb_target_link_libraries(DEBUG   linalgwrap_tests 
	${linalgwrap_DEBUG_TARGET} ${LINALGWRAP_DEPENDENCIES_DEBUG})
drb_target_link_libraries(RELEASE linalgwrap_tests 
	${linalgwrap_RELEASE_TARGET} ${LINALGWRAP_DEPENDENCIES_RELEASE})

# Add the test executables
FOREACH(target ${linalgwrap_tests_TARGETS})
	add_test(NAME ${target} COMMAND ${target})
ENDFOREACH()

