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
		
