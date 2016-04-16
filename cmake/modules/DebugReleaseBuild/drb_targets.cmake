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

include(CMakeParseArguments)

macro(drb_add_target TARGET KIND)
	# Add a library or executable to a DebugReleaseBuild
	# For internal usage only. Interface may change at any
	# time. Use drb_add_library or drb_add_executable instead.

	# Check that we have a drb setup
	if(NOT ${CMAKE_PROJECT_NAME}_HAS_DRB_SETUP)
		message(FATAL_ERROR "You have to setup the DebugReleaseBuild \
for this project(${CMAKE_PROJECT_NAME}) before setting up any targets via \
drb_add_target. Try calling drb_setup_project first.")
	endif()

	# Set the target variables
	set(${TARGET}_DEBUG_TARGET "${TARGET}.g" CACHE INTERNAL "Target name of debug version of ${TARGET}")
	set(${TARGET}_RELEASE_TARGET "${TARGET}" CACHE INTERNAL "Target name of release version of ${TARGET}")

	# Set the empty TARGETS list.
	set(${TARGET}_TARGETS "" CACHE INTERNAL "Name of target (Debug/Release) of ${TARGET} we actually build")

	# Convert ARGN into a real list.
	set(ARGN_list "${ARGN}")

	if ("${ARGV2}" STREQUAL "DBGSUFFIX" OR "${ARGV2}" STREQUAL "RELSUFFIX" OR "${ARGV2}" STREQUAL "FILES") 
		cmake_parse_arguments(ADD_TARGET "" "DBGSUFFIX;RELSUFFIX" "FILES" ${ARGN})

		if(NOT "${ADD_TARGET_DBGSUFFIX}" STREQUAL "")
			set(${TARGET}_DEBUG_TARGET "${TARGET}${ADD_TARGET_DBGSUFFIX}")
		endif()

		if(NOT "${ADD_TARGET_RELSUFFIX}" STREQUAL "")
			set(${TARGET}_RELEASE_TARGET "${TARGET}${ADD_TARGET_RELSUFFIX}")
		endif()

		if("${ADD_TARGET_FILES}" STREQUAL "")
			message(FATAL_ERROR "Did not supply a FILES block.")
		else()
			set(ARGN_list "${ADD_TARGET_FILES}")
		endif()

		unset(ADD_TARGET_DBGSUFFIX)
		unset(ADD_TARGET_RELSUFFIX)
		unset(ADD_TARGET_FILES)
	endif()

	foreach(build ${DRB_BUILD_TYPES})
		if ("${KIND}" STREQUAL "lib")
			# Add the actual library target:
			add_library("${${TARGET}_${build}_TARGET}" ${ARGN_list})

			# Extract linker flags
			if (BUILD_SHARED_LIBS)
				set(LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS_${build}}")
			else()
				set(LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} ${CMAKE_STATIC_LINKER_FLAGS_${build}}")
			endif()
		elseif("${KIND}" STREQUAL "exe")
			add_executable("${${TARGET}_${build}_TARGET}" ${ARGN_list})

			# Set linker flags:
			set(LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_${build}}")
		else()
			message(FATAL_ERROR "Unknown kind: ${KIND}")
		endif()

		# and append it to the list
		list(APPEND ${TARGET}_TARGETS "${${TARGET}_${build}_TARGET}")

		# set the compilation properties we are responsible for
		set_target_properties(
			"${${TARGET}_${build}_TARGET}"
			PROPERTIES
			LINK_FLAGS "${LINKER_FLAGS}"
			LINKER_LANGUAGE "CXX"
			COMPILE_FLAGS "${CMAKE_CXX_FLAGS_${build}}" #${CMAKE_CXX_FLAGS} is added by default
			COMPILE_DEFINITIONS "${build}"
		)
	endforeach()

	unset(ARGN_list)
	unset(LINKER_FLAGS)
endmacro(drb_add_target)

macro(drb_add_library TARGET)
	# Add a library to a DebugReleaseBuild.
	# This sets up a target for the debug build and a target for the 
	# release build. Compiler flags and linker flags are set appropriately.
	#
	# In the default setting the call signature is:
	#     drb_add_library(<TARGET> <files> ...)
	# and this sets up the targets <TARGET>.g for the debug build and 
	# the <TARGET> for the release build.
	#
	# For the setup purpose the function reads the following variables:
	#       CMAKE_SHARED_LINKER_FLAGS CMAKE_SHARED_LINKER_FLAGS_DEBUG
	#		CMAKE_SHARED_LINKER_FLAGS_RELEASE
	#       CMAKE_STATIC_LINKER_FLAGS CMAKE_STATIC_LINKER_FLAGS_DEBUG
	#		CMAKE_STATIC_LINKER_FLAGS_RELEASE
	#	CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG
	#		CMAKE_CXX_FLAGS_RELEASE
	#
	# If one wants to change the target name of debug and release build one
	# can use the alternative signature:
	#    drb_add_library(<TARGET>
	#                    DBGSUFFIX <debug_suffix>
	#                    RELSUFFIX <release_suffix>
	#                    FILES <files> ...)
	# Then the target for the debug build is <TARGET><debug_suffix> and the
	# target for the release build is <TARGET><release_suffix>.
	#	
	# The macro sets the variables
	#       ${TARGET}_DEBUG_TARGET		full name of the Debug target	
	#	${TARGET}_RELEASE_TARGET	full name of the Release target
	#	${TARGET}_TARGETS
	#		list of the above targets which have been configured
	#		using this macro. I.e. if we only build DEBUG 
	#		(DRB_BUILD_TYPES=="DEBUG"), then it contains only 
	#		the name of the debug target. Helpful when looping over
	#		all targets configured by this macro
	#

	# Check that we have a drb setup
	if(NOT ${CMAKE_PROJECT_NAME}_HAS_DRB_SETUP)
		message(FATAL_ERROR "You have to setup the DebugReleaseBuild \
for this project(${CMAKE_PROJECT_NAME}) before setting up any targets via \
drb_add_library. Try calling drb_setup_project first.")
	endif()

	# Make ARGN string a list, which we then can properly pass on
	string(REPLACE " " ";" VALUES "${ARGN}")

	# Call the common function
	drb_add_target(${TARGET} "lib" ${VALUES})

	unset(VALUES)
endmacro(drb_add_library)

#
# ------------------------------
#

macro(drb_add_executable TARGET)
	# Add an executable to a DebugReleaseBuild.
	# This sets up a target for the debug build and a target for the 
	# release build. Compiler flags and linker flags are set appropriately.
	#
	# In the default setting the call signature is:
	#     drb_add_executable(<TARGET> <files> ...)
	# and this sets up the targets <TARGET>.g for the debug build and 
	# the <TARGET> for the release build.
	#
	# For the setup purpose the function reads the following variables:
	#	CMAKE_EXE_LINKER_FLAGS CMAKE_EXE_LINKER_FLAGS_DEBUG
	#		CMAKE_EXE_LINKER_FLAGS_RELEASE
	#	CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG
	#		CMAKE_CXX_FLAGS_RELEASE
	#
	# If one wants to change the target name of debug and release build one
	# can use the alternative signature:
	#    drb_add_executable(<TARGET>
	#                       DBGSUFFIX <debug_suffix>
	#                       RELSUFFIX <release_suffix>
	#                       FILES <files> ...)
	# Then the target for the debug build is <TARGET><debug_suffix> and the
	# target for the release build is <TARGET><release_suffix>.
	#	
	# The macro sets the variables
	#       ${TARGET}_DEBUG_TARGET		full name of the Debug target	
	#	${TARGET}_RELEASE_TARGET	full name of the Release target
	#	${TARGET}_TARGETS
	#		list of the above targets which have been configured
	#		using this macro. I.e. if we only build DEBUG 
	#		(DRB_BUILD_TYPES=="DEBUG"), then it contains only 
	#		the name of the debug target. Helpful when looping over
	#		all targets configured by this macro
	#
	#	
	# The macro sets the variables
	#       ${TARGET}_DEBUG_TARGET		full name of the Debug target	
	#	${TARGET}_RELEASE_TARGET	full name of the Release target
	#	${TARGET}_TARGETS
	#		list of the above targets which have been configured
	#		using this macro. I.e. if we only build DEBUG 
	#		(DRB_BUILD_TYPES=="DEBUG"), then it contains only 
	#		the name of the debug target. Helpful when looping over
	#		all targets configured by this macro
	#

	# check that drb has been initialised
	if(NOT ${CMAKE_PROJECT_NAME}_HAS_DRB_SETUP)
		message(FATAL_ERROR "You have to setup the DebugReleaseBuild \
for this project(${CMAKE_PROJECT_NAME}) before setting up any targets via \
drb_add_executable. Try calling drb_setup_project first.")
	endif()

	# Make ARGN string a list, which we then can properly pass on
	string(REPLACE " " ";" VALUES "${ARGN}")

	# Call the common function
	drb_add_target(${TARGET} "exe" ${VALUES})

	unset(VALUES)
endmacro(drb_add_executable)

#
# ------------------------------
#

macro(drb_set_target_properties WHICH TARGET)
	# Call target_link_libraries for the proper target names 
	# depinding on WHICH.
	#
	# if WHICH == ALL      calls it for all configured configurations
	# if WHICH == DEBUG    calls it for the DEBUG configurations
	#                      (if this configuration is enabled)
	# if WHICH == RELEASE  calls it for the RELEASE configurations
	#                      (if this configuration is enabled)
	#
	# TARGET should be exactly the target base name given to 
	# drb_add_library and drb_add_executable
	#

	# check that drb has been initialised
	if("${${TARGET}_TARGETS}" STREQUAL "")
		message(FATAL_ERROR "TARGET(${TARGET}) has not been setup with DebugReleaseBuild.")
	endif()

	# Check that we do not want to override the COMPILE_DEFINITIONS, since we need them
	# the way they are to distinguish debug and release builds.
	if("${ARGN}" MATCHES "COMPILE_DEFINITIONS")
		message(FATAL_ERROR "Setting COMPILE_DEFINITIONS using \
drb_set_target_properties is not supported.")
	endif()

	if (NOT "${ARGN}" STREQUAL "")
		# Extract the build types we want to change:
		# This effectively filters DRB_BUILD_TYPES according
		# to WHICH
		drb_get_matching_build_types(${WHICH} MATCHLIST)

		# Make ARGN string a list, which we then can properly pass to 
		# target_compile_definitions.
		string(REPLACE " " ";" VALUES "${ARGN}")

		foreach(build ${MATCHLIST})
			set_target_properties(${${TARGET}_${build}_TARGET} ${VALUES})
		endforeach()

		unset(MATCHLIST)
		unset(VALUES)
	endif()
endmacro(drb_set_target_properties)

#
# ------------------------------
#

macro(drb_target_link_libraries WHICH TARGET)
	# Call target_link_libraries for the proper target names 
	# depinding on WHICH.
	#
	# if WHICH == ALL      calls it for all configured configurations
	# if WHICH == DEBUG    calls it for the DEBUG configurations
	#                      (if this configuration is enabled)
	# if WHICH == RELEASE  calls it for the RELEASE configurations
	#                      (if this configuration is enabled)
	#
	# TARGET should be exactly the target base name given to 
	# drb_add_library and drb_add_executable
	#

	# check that drb has been initialised
	if("${${TARGET}_TARGETS}" STREQUAL "")
		message(FATAL_ERROR "TARGET(${TARGET}) has not been setup with DebugReleaseBuild.")
	endif()

	if (NOT "${ARGN}" STREQUAL "")
		# Extract the build types we want to change:
		# This effectively filters DRB_BUILD_TYPES according
		# to WHICH
		drb_get_matching_build_types(${WHICH} MATCHLIST)

		# Make ARGN string a list, which we then can properly pass to 
		# target_compile_definitions.
		string(REPLACE " " ";" VALUES "${ARGN}")

		foreach(build ${MATCHLIST})
			target_link_libraries(${${TARGET}_${build}_TARGET} ${VALUES})
		endforeach()

		unset(MATCHLIST)
		unset(VALUES)
	endif()
endmacro(drb_target_link_libraries)

#
# ------------------------------
#

macro(drb_target_compile_definitions WHICH TARGET)
	# Call target_compile_definitions for the proper target names 
	# depinding on WHICH.
	#
	# if WHICH == ALL      calls it for all configured configurations
	# if WHICH == DEBUG    calls it for the DEBUG configurations
	#                      (if this configuration is enabled)
	# if WHICH == RELEASE  calls it for the RELEASE configurations
	#                      (if this configuration is enabled)
	#
	# TARGET should be exactly the target base name given to 
	# drb_add_library and drb_add_executable
	#

	# check that drb has been initialised
	if("${${TARGET}_TARGETS}" STREQUAL "")
		message(FATAL_ERROR "TARGET(${TARGET}) has not been setup with DebugReleaseBuild.")
	endif()

	if (NOT "${ARGN}" STREQUAL "")
		# Extract the build types we want to change:
		# This effectively filters DRB_BUILD_TYPES according
		# to WHICH
		drb_get_matching_build_types(${WHICH} MATCHLIST)

		# Make ARGN string a list, which we then can properly pass to 
		# target_compile_definitions.
		string(REPLACE " " ";" VALUES "${ARGN}")

		foreach(build ${MATCHLIST})
			target_compile_definitions(${${TARGET}_${build}_TARGET} ${VALUES})
		endforeach()

		unset(MATCHLIST)
		unset(VALUES)
	endif()
endmacro(drb_target_compile_definitions)
