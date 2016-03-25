# sets these things
#
# 	LINALGWRAP_DEPENDENCIES			everyone needs these libraries
# 	LINALGWRAP_DEPENDENCIES_DEBUG		debug mode needs these extras
# 	LINALGWRAP_DEPENDENCIES_RELEASE		release mode needs these extras
# 	LINALGWRAP_DEPENDENCIES_TEST		tests need these extra libraries

##############
#-- catch  --#
##############

if (LINALGWRAP_WITH_EXTERNAL_CATCH)
	message(FATAL_ERROR "Currently only the build with the catch from rapidcheck as a submodule is supported")
	#TODO change this
else()
	# We need to compile catch further down in the externals
	# so add it to the dependencies
	set(LINALGWRAP_DEPENDENCIES_TEST ${LINALGWRAP_DEPENDENCIES_TEST} catch)
endif()


##################
#-- rapidcheck --#
##################

if (LINALGWRAP_WITH_EXTERNAL_RAPIDCHECK)
	message(FATAL_ERROR "Currently only the build with rapidcheck as a submodule is supported")
	#TODO change this
else()
	# We need to compile rapidcheck further down in the externals
	# add it to the dependencies right here:
	set(LINALGWRAP_DEPENDENCIES_TEST ${LINALGWRAP_DEPENDENCIES_TEST} rapidcheck catch)
endif()


#################
#-- armadillo --#
#################

SET(ARMADILLO_DIR "$ENV{ARPACK_DIR}" CACHE PATH "Optional hint to find and armadillo installation")

find_library(ARMADILLO_LIBRARY armadillo
	DOC "Path to the armadillo library."
	HINTS ${ARMADILLO_DIR} ${ARMADILLO_DIR}/lib
)
if ("${ARMADILLO_LIBRARY}" STREQUAL "ARMADILLO_LIBRARY-NOTFOUND")
	message(FATAL_ERROR "Could not find armadillo library. Specify ARMADILLO_DIR for a hint")
endif()

find_path(ARMADILLO_INCLUDE armadillo
	DOC "Path to the armadillo include directory."
	HINTS ${ARMADILLO_DIR} ${ARMADILLO_DIR}/include
)
if ("${ARMADILLO_INCLUDE}" STREQUAL "ARMADILLO_INCLUDE-NOTFOUND")
	message(FATAL_ERROR "Could not find armadillo include files. Specify ARMADILLO_DIR for a hint")
endif()

# add to general dependencies:
set(LINALGWRAP_DEPENDENCIES ${LINALGWRAP_DEPENDENCIES} ${ARMADILLO_LIBRARY})

# enable armadillo-dependant code:
LIST(APPEND LINALGWRAP_DEFINITIONS "LINALGWRAP_HAVE_ARMADILLO")
