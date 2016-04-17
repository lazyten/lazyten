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

# adds entries to these things
#
# 	LINALGWRAP_DEPENDENCIES			everyone needs these libraries
# 	LINALGWRAP_DEPENDENCIES_DEBUG		debug mode needs these extras
# 	LINALGWRAP_DEPENDENCIES_RELEASE		release mode needs these extras
# 	LINALGWRAP_DEPENDENCIES_TEST		tests need these extra libraries
#
#       LINALGWRAP_DEFINITIONS			definitions for all compilation
#       LINALGWRAP_DEFINITIONS_DEBUG		definitions for debug mode
#       LINALGWRAP_DEFINITIONS_RELEASE		definitions for release mode
#       


#################################
#--     Modules and macros    --#
#################################
include(CheckCXXSourceCompiles)

########################
#-- Exception system --#
########################

#
# Check whether stacktrace information is available for the exception
# system. We expect the interface of glibc. 
# If it is than we also add the -rdynamic flag since this is required
# if one wants meaningful backtraces.
#
CHECK_CXX_SOURCE_COMPILES(
	"
	#include <execinfo.h>
	#include <stdlib.h>
	const int asize = 25;
	void * array[asize];
	int nBT = backtrace(array, asize);
	char ** bt_raw = backtrace_symbols(array, nBT);
	int main() { free(bt_raw); return 0; }
	"
	LINALGWRAP_HAVE_GLIBC_STACKTRACE)

if(LINALGWRAP_HAVE_GLIBC_STACKTRACE)
	# needed for meaningful stacktraces
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -rdynamic")
	set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -rdynamic")
	LIST(APPEND LINALGWRAP_DEFINITIONS "LINALGWRAP_HAVE_GLIBC_STACKTRACE")
endif()

#
# Check for demangling symbols from within the program.
# We expect the interface from the libstdc++
#
# The test example is taken from 
#   http://gcc.gnu.org/onlinedocs/libstdc++/18_support/howto.html#6
#
if(LINALGWRAP_HAVE_GLIBC_STACKTRACE)
	CHECK_CXX_SOURCE_COMPILES(
		"
		#include <exception>
		#include <iostream>
		#include <cxxabi.h>
		#include <cstdlib>

		struct empty { };

		template <typename T, int N>
		struct bar { };

		int     status;
		char   *realname;

		int main()
		{
		// exception classes not in <stdexcept>, thrown by the implementation
		// instead of the user
		std::bad_exception  e;
		realname = abi::__cxa_demangle(e.what(), 0, 0, &status);
		free(realname);


		// typeid
		bar<empty,17>          u;
		const std::type_info  &ti = typeid(u);

		realname = abi::__cxa_demangle(ti.name(), 0, 0, &status);
		free(realname);

		return 0;
		}
		"
		LINALGWRAP_HAVE_LIBSTDCXX_DEMANGLER)

	if(LINALGWRAP_HAVE_LIBSTDCXX_DEMANGLER)
		LIST(APPEND LINALGWRAP_DEFINITIONS "LINALGWRAP_HAVE_LIBSTDCXX_DEMANGLER")
	endif()
endif()

