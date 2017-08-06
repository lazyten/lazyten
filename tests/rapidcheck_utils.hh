//
// Copyright (C) 2016-17 by the lazyten authors
//
// This file is part of lazyten.
//
// lazyten is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// lazyten is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with lazyten. If not, see <http://www.gnu.org/licenses/>.
//

#pragma once
#include <rapidcheck.h>
#include <rapidcheck/state.h>

/** Make a rapidcheck assertion with captures disabled
 *
 * This is required for the krims::numcomp function to work
 * properly --- which I personally value more)
 * */
#define RC_ASSERT_NC(expression) RC_ASSERT((expression));

namespace rc {
namespace state {
namespace gen {

// TODO implement into rapidcheck in some sensible way
template <typename Model, typename GenerationFunc>
auto commandsScaledLength(const Model& initialState, double scale,
                          GenerationFunc&& genFunc)
      -> decltype(commands<Model>(initialState, genFunc)) {

  /// Generate a sequence of commands where the commands themself
  /// get passed the size ``size``.
  auto commands_with_size = [=](int size) {
    return commands<Model>(initialState, [=](const Model& state) {
      return rc::gen::resize(size, genFunc(state));
    });
  };

  return rc::gen::withSize(
        [=](int size) { return rc::gen::scale(scale, commands_with_size(size)); });
}

}  // namespace gen
}  // namespace state
}  // namespace rc

namespace lazyten {
namespace tests {
using namespace rc;

/** Namespace for utility components for
 * testing stuff with rapidcheck
 **/
namespace rapidcheck_utils {

// TODO use until something better exists
template <typename Model, typename Sut, typename GenFunc>
void state_check_scaled(const Model& initialState, Sut& sut, double scale,
                        GenFunc&& generationFunc) {

  const auto commandsGenScaled =
        *rc::state::gen::commandsScaledLength(initialState, scale, generationFunc);

  // Run the thing:
  runAll(commandsGenScaled, initialState, sut);
}

}  // namespace matrix_test_utils
}  // namespace tests
}  // namespace lazyten
