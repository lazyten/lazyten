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
# set the source files.
#
set(LAZYTEN_SOURCES
	version.cc
	io/Mathematica.cc
	Base/Interfaces/OperatorProperties.cc
	Base/Interfaces/Transposed.cc
	Base/Solvers/EigensolverBaseKeys.cc
	Base/Solvers/IterativeWrapperKeys.cc
	Base/Solvers/LinearSolverBaseKeys.cc
	Base/Solvers/SolverStateBase.cc
	Arpack/ArpackEigensolver.cc
	Lapack/LapackEigensolver.cc
	Lapack/detail/lapack.cc
	LinearSolver.cc
	EigensystemSolver.cc
	rescue.cc
)

# Write the config file config.hh and the version header (into the binary dir)
configure_file("config.hh.in"  "config.hh")
feature_list_cxx(PROJECT_FEATURE_LIST)
configure_file("version.hh.in" "version.hh")

# add the libraries for Debug and Release:
drb_add_library(lazyten
	DBGSUFFIX  ".g"
	FILES ${LAZYTEN_SOURCES}
)

drb_set_target_properties(ALL lazyten
	PROPERTIES
	VERSION "${PROJECT_VERSION}"
)

drb_target_compile_definitions(ALL     lazyten
	PUBLIC ${LAZYTEN_DEFINITIONS})
drb_target_compile_definitions(DEBUG   lazyten
	PUBLIC ${LAZYTEN_DEFINITIONS_DEBUG})
drb_target_compile_definitions(RELEASE lazyten
	PUBLIC ${LAZYTEN_DEFINITIONS_RELEASE})

drb_target_link_libraries(ALL     lazyten ${LAZYTEN_DEPENDENCIES})
drb_target_link_libraries(DEBUG   lazyten ${LAZYTEN_DEPENDENCIES_DEBUG})
drb_target_link_libraries(RELEASE lazyten ${LAZYTEN_DEPENDENCIES_RELEASE})

# set the installation property for these targets:
install(TARGETS ${lazyten_TARGETS}
	EXPORT LazytenTargets
	COMPONENT library
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	INCLUDES DESTINATION include
)

# install the includes:
install(DIRECTORY .
	DESTINATION "include/lazyten"
	COMPONENT devel
	FILES_MATCHING PATTERN "*.hh"
)
install(FILES
	"${lazyten_BINARY_DIR}/src/lazyten/config.hh"
	"${lazyten_BINARY_DIR}/src/lazyten/version.hh"
	DESTINATION "include/lazyten"
	COMPONENT devel
)

# Export the target specifications for lazyten.
export(EXPORT LazytenTargets
	FILE "${lazyten_BINARY_DIR}/lazytenTargets.cmake"
	NAMESPACE Upstream::
)

# and install them appropriately:
install(EXPORT LazytenTargets
	FILE "lazytenTargets.cmake"
	NAMESPACE Upstream::
	DESTINATION "${PackageModuleLocation}/lazyten"
)

