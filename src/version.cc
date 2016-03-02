#include "version.hh"
#include <string>
#include <sstream>

struct version {
    int major{VERSION_MAJOR};
    int minor{VERSION_MINOR};
    int revision{VERSION_REVISION};

    std::string version_string() {
        std::stringstream ss;
        ss << major << "." << minor << "." << revision;
        return ss.str();
    }
};
