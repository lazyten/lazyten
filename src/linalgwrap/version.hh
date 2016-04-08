#pragma once
#include <string>

#define VERSION_MAJOR 0
#define VERSION_MINOR 1
#define VERSION_REVISION 0

struct version {
    static int constexpr major{VERSION_MAJOR};
    static int constexpr minor{VERSION_MINOR};
    static int constexpr revision{VERSION_REVISION};

    // Return the version as a string
    static std::string version_string();
};
