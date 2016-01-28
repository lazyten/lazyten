#include "ExceptionBase.hh"

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <sstream>

#ifdef LINALGWRAP_HAVE_GLIBC_STACKTRACE
#include <execinfo.h>
#endif

#ifdef HAVE_LIBSTDCXX_DEMANGLER
#include <cxxabi.h>
#endif

namespace linalgwrap {

ExceptionBase::ExceptionBase()
      : m_name(""),
        m_file(""),
        m_line(0),
        m_function(""),
        m_failed_condition(""),
        m_n_stacktrace_frames(0),
        m_what_str("") {}

void ExceptionBase::add_exc_data(const char* file, int line,
                                 const char* function,
                                 const char* failed_condition,
                                 const char* exception_name) {
    m_name = exception_name;
    m_file = file;
    m_line = line;
    m_function = function;
    m_failed_condition = failed_condition;

#ifdef LINALGWRAP_HAVE_GLIBC_STACKTRACE
    // If the system supports it, get a stacktrace how we got here
    // We defer the symbol lookup via backtrace_symbols() since
    // this loads external libraries which can take up to seconds
    // on some machines.
    // See generate_stacktrace() for the place where this is done.
    m_n_stacktrace_frames = backtrace(m_raw_stacktrace, 25);
#endif
}

const char* ExceptionBase::what() const noexcept {
    // if string is empty: build it.
    if (m_what_str == "") m_what_str = generate_message();

    // TODO: Is this noexcept???
    return m_what_str.c_str();
}

const char* ExceptionBase::name() const { return m_name; }

void ExceptionBase::print_extra(std::ostream& out) const noexcept {
    out << "(none)";
}

void ExceptionBase::print_exc_data(std::ostream& out) const noexcept {
    out << "The assertion" << std::endl << "   " << m_failed_condition
        << std::endl << "failed in line " << m_line << " of file \"" << m_file
        << "\" while executing the function" << std::endl << "   " << m_function
        << std::endl << "This raised the exception" << std::endl << "   "
        << m_name << std::endl;
}

void ExceptionBase::print_stacktrace(std::ostream& out) const {
    // return if no stracktrace frames
    if (m_n_stacktrace_frames <= 0) return;

#ifdef LINALGWRAP_HAVE_GLIBC_STACKTRACE
    // We have deferred the symbol lookup to this point to avoid costly
    // runtime penalties due to linkage of external libraries by
    // backtrace_symbols.
    char** stacktrace =
          backtrace_symbols(m_raw_stacktrace, m_n_stacktrace_frames);

    // if there is a stackframe stored, print it
    out << std::endl;
    out << "Backtrace:" << std::endl;
    out << "----------" << std::endl;
    out << "#   file  functionname  @  address" << std::endl;
    out << "----------------------------------" << std::endl << std::endl;

    // Skip the frames we are not interested in:
    int initframe = 0;
    for (; initframe < m_n_stacktrace_frames; ++initframe) {
        if (std::strstr(stacktrace[initframe], "linalgwrap") &&
            std::strstr(stacktrace[initframe], "exception") &&
            std::strstr(stacktrace[initframe], "raise")) {
            // The current frame contains linalgwrap, exceptions and raise, ie.
            // it is the one
            // corresponding to
            // linalgwrap::exceptions::raise
            // so the next is the first frame we are interested in
            ++initframe;
            break;
        }
    }

    for (int frame = initframe; frame < m_n_stacktrace_frames; ++frame) {
        out << "#" << frame - initframe << "  ";

        const std::string stacktrace_entry(stacktrace[frame]);
        // The stacktrace frames are in the format
        //     filename(functionname+offset) [address]
        // Try to extract the functionname
        const size_t pos_bracket = stacktrace_entry.find("(");
        const size_t pos_bracketcl = stacktrace_entry.find(")");
        const size_t pos_plus = stacktrace_entry.find("+");

        // If functionname+offset is the empty string, then the user omitted the
        // flag -rdynamic on compilation. Warn him and continue:
        if (pos_bracketcl == pos_bracket + 1) {
            out << stacktrace_entry << "   (add flag \"-rdynamic\" on linking "
                                       "to improve stacktrace)" << std::endl;
            continue;
        }

        // Split up the stacktrace string:
        const std::string functionname = stacktrace_entry.substr(
              pos_bracket + 1, pos_plus - pos_bracket - 1);
        const std::string file = stacktrace_entry.substr(0, pos_bracket);
        const std::string offset = stacktrace_entry.substr(
              pos_plus + 1, pos_bracketcl - pos_plus - 1);

        // print the file name:
        out << file;

#ifdef LINALGWRAP_HAVE_LIBSTDCXX_DEMANGLER
        // try to demangle the function name:
        int status;
        char* p = abi::__cxa_demangle(functionname.c_str(), 0, 0, &status);

        if (status == 0) {
            // all well
            out << "  " << p;
        } else {
            out << "  " << functionname;
        }

        // Free the allocated char*
        free(p);
#else
        // output using the original function name:
        out << "  " << functionname;
#endif

        // print the raw address:
        out << " @ " << m_raw_stacktrace[frame];

        // finish the line:
        out << std::endl;

        // stop once we are in main()
        if (functionname == "main") break;
    }

    // TODO: have this hint?
    out << std::endl << "Hint: Use \"addr2line -e <file> <address>\" to get "
                        "line number in BT." << std::endl;

    // Free the memory for the above.
    free(stacktrace);
#endif
}

std::string ExceptionBase::generate_message() const noexcept {
    // Build c_string inside a try block since this function is noexcept
    try {
        std::ostringstream converter;

        converter << std::endl
                  << "--------------------------------------------------------"
                  << std::endl;

        // Print first the general data we hold:
        print_exc_data(converter);

        // Now print the extra information:
        converter << std::endl << "Extra information:" << std::endl;
        print_extra(converter);

        // Finally print a stacktrace:
        print_stacktrace(converter);

        converter << std::endl
                  << "--------------------------------------------------------"
                  << std::endl;

        return converter.str();
    } catch (...) {
        // Deal with error by printing some generic message
        std::string message(
              "Failed to generate the exception message in "
              "ExceptionBase::generate_message()");
        return message;
    }
}

}  // linalgwrap
