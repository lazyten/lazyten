# Setup optional dependencies and features

########################
#-- Exception system --#
########################

#
# Check whether stacktrace information is available for the exception
# system. We expect the interface of glibc. 
# If it is than we also add the -rdynamic flag since this is required
# if one wants meaningful backtraces.
#
#CHECK_CXX_SOURCE_COMPILES(
#	"
#	#include <execinfo.h>
#	#include <stdlib.h>
#	const int asize = 25;
#	void * array[asize];
#	int nBT = backtrace(array, asize);
#	char ** bt_raw = backtrace_symbols(array, nBT);
#	int main() { free(symbols); return 0; }
#	"
#	LINALGWRAP_HAVE_GLIBC_STACKTRACE)
message(AUTHOR_WARNING "TODO Test for glibc stacktrace")
set(LINALGWRAP_HAVE_GLIBC_STACKTRACE)

if(LINALGWRAP_HAVE_GLIBC_STACKTRACE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -rdynamic")
endif()

#
# Check for demangling symbols from within the program.
# We expect the interface from the libstdc++
#
# The test example is taken from 
#   http://gcc.gnu.org/onlinedocs/libstdc++/18_support/howto.html#6
#
if(LINALGWRAP_HAVE_GLIBC_STACKTRACE)
#	CHECK_CXX_SOURCE_COMPILES(
#		"
#		#include <exception>
#		#include <iostream>
#		#include <cxxabi.h>
#		#include <cstdlib>
#
#		struct empty { };
#
#		template <typename T, int N>
#		struct bar { };
#
#		int     status;
#		char   *realname;
#
#		int main()
#		{
#		// exception classes not in <stdexcept>, thrown by the implementation
#		// instead of the user
#		std::bad_exception  e;
#		realname = abi::__cxa_demangle(e.what(), 0, 0, &status);
#		free(realname);
#
#
#		// typeid
#		bar<empty,17>          u;
#		const std::type_info  &ti = typeid(u);
#
#		realname = abi::__cxa_demangle(ti.name(), 0, 0, &status);
#		free(realname);
#
#		return 0;
#		}
#		"
#		LINALGWRAP_HAVE_LIBSTDCXX_DEMANGLER)
	message(AUTHOR_WARNING "TODO Test for libc++ demangler")
	set(LINALGWRAP_HAVE_LIBSTDCXX_DEMANGLER)
endif()

