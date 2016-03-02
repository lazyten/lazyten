#pragma once
#include <rapidcheck.h>
#include <rapidcheck/state.h>

namespace linalgwrap {
namespace tests {
using namespace rc;

/** Namespace for utility components for
 * testing stuff with rapidcheck
 **/
namespace rapidcheck_utils {

// TODO use until something better exists
template <typename Model, typename Sut, typename GenFunc>
void state_check_scaled(const Model &initialState, Sut &sut, double scale,
                        GenFunc &&generationFunc) {
    // Undo the scaling on the generationFunc
    const auto generationFuncUnscaled = [=](const Model &state) {
        auto gen = generationFunc(state);
        return rc::gen::scale(1. / scale, gen);
    };

    // Define the commands generator:
    const auto commandsGen =
          rc::state::gen::commands<rc::state::Command<Model, Sut>>(
                initialState, generationFuncUnscaled);

    // Scale the generator:
    const auto commandsGenScaled = *rc::gen::scale(scale, commandsGen);

    // Run the thing:
    runAll(commandsGenScaled, initialState, sut);
}

}  // matrix_test_utils
}  // tests
}  // linalgwrap
