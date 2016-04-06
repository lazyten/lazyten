# Setup optional dependencies and features

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
	set(CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} -rdynamic")
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

