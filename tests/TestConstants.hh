#pragma once
#include <cstddef>
#include <limits>

namespace linalgwrap {
namespace tests {

struct TestConstants {
    /** The maximal matrix size arbitrary matrices have */
    static constexpr std::size_t max_matrix_size = 10;

    /** The default numeric tolerance when comparing
     * computed values for equality. Usually comparison
     * operations will interpret this as tolerance towards
     * the *relative* error */
    static constexpr double default_num_tol = 1e-12;
};
}
}
