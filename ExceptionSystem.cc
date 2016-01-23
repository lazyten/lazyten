#include "ExceptionSystem.hh"
#include <iostream>

namespace linalgwrap {

namespace exceptions {

void handle(const ExceptionBase& ex, ExceptionEffect effect) {
    if (effect == ExceptionEffect::PRINT) {
        // Just print exception
        std::cerr << "Exception: " << ex.name() << std::endl;
    } else if (effect == ExceptionEffect::THROW) {
        // Throw it
        throw ex;
    } else {
        // print what is going on and abort.
        std::cerr << ex.what() << std::endl;
        std::abort();
    }
}

}  // exceptions
}  // linalgwrap
