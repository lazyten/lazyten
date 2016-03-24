#pragma once
#include <rapidcheck.h>
#include <rapidcheck/state.h>

namespace rc {
namespace state {
namespace gen {

// TODO implement into rapidcheck in some sensible way
template <typename Cmd, typename GenerationFunc>
Gen<Commands<Cmd>> commandsScaledLength(const typename Cmd::Model &initialState,
                                        double scale,
                                        GenerationFunc &&genFunc) {

    /// Generate a sequence of commands where the commands themself
    /// get passed the size ``size``.
    auto commands_with_size = [=](int size) {
        return commands<Cmd>(initialState,
                             [=](const typename Cmd::Model &state) {
            return rc::gen::resize(size, genFunc(state));
        });
    };

    return rc::gen::withSize([=](int size) {
        return rc::gen::scale(scale, commands_with_size(size));
    });
}

}  // namespace gen
}  // namespace state
}  // namespace rc

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

    const auto commandsGenScaled =
          *rc::state::gen::commandsScaledLength<rc::state::Command<Model, Sut>>(
                initialState, scale, generationFunc);

    // Run the thing:
    runAll(commandsGenScaled, initialState, sut);
}

}  // matrix_test_utils
}  // tests
}  // linalgwrap
