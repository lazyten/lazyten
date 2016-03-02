#include "NumComp.hh"
#include <sstream>

namespace linalgwrap {
namespace tests {

NumCompException::NumCompException(double lhs_, double rhs_, double error_,
                                   double tolerance_,
                                   const std::string description_) noexcept
      : lhs{lhs_},
        rhs{rhs_},
        error{error_},
        tolerance{tolerance_},
        description{description_},
        what_str("") {}

const char* NumCompException::what() const noexcept {
    try {
        if (what_str == "") {
            std::stringstream ss;

            ss << std::scientific << std::setprecision(15) << lhs
               << " == " << rhs << " returned false. Error(" << error
               << "), Tolerance(" << tolerance
               << "), Description: " << description;

            what_str = ss.str();
        }
        return what_str.c_str();
    } catch (...) {
        return "Error building what_str";
    }
}
}
}
