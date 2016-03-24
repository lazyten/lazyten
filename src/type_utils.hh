#pragma once
#include <type_traits>

namespace linalgwrap {

/* Using statement we need for some SFINAE
 * (substitution failure is not an error) */
template <typename... Ts>
using void_t = void;
}  // namespace linalgwrap
