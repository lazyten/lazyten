#include "version.hh"
#include <sstream>

std::string version::version_string() {
    std::stringstream ss;
    ss << major << "." << minor << "." << revision;
    return ss.str();
}
