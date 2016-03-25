#
# The sources for this example executable
#
set(EXAMPLE_DIAGONAL_SOURCES
	diagonal.cc)

# Compile and link:
add_executable(example_diagonal ${EXAMPLE_DIAGONAL_SOURCES})

set_target_properties(
	example_diagonal
	PROPERTIES
	LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_${EXAMPLES_BUILD_TYPE}}"
	LINKER_LANGUAGE "CXX"
	COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${EXAMPLES_BUILD_TYPE}}"
	COMPILE_DEFINITIONS "${LINALGWRAP_DEFINITIONS};${LINALGWRAP_DEFINITIONS_${EXAMPLES_BUILD_TYPE}}"
)

# link it with the appropriate linalgwrap library target
target_link_libraries(example_diagonal 
	${LINALGWRAP_${EXAMPLES_BUILD_TYPE}_TARGET}
	${LINALGWRAP_DEPENDENCIES}
	${LINALGWRAP_${EXAMPLES_BUILD_TYPE}_DEPENDENCIES}
)
