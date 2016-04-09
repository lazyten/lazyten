//
// Copyright (C) 2016 by the linalgwrap authors
//
// This file is part of linalgwrap.
//
// linalgwrap is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// linalgwrap is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with linalgwrap. If not, see <http://www.gnu.org/licenses/>.
//

#pragma once
#include <cstddef>
#include <limits>

namespace linalgwrap {
namespace tests {

struct TestConstants {
    //
    // Random matrix properties:
    //
    /* \brief The maximal matrix size arbitrary matrices have
     *
     * This means that both columns and rows, so a matrix can have
     * a maximum of 100 entries if this is chosen to be 10.
     */
    static constexpr std::size_t max_matrix_size = 10;

    /* \brief The maximum magnitude an arbitrary matrix entry may have
     *
     * If this is set to 6, then all entries will have an absolute
     * value smaller than $10^6$.
     */
    static constexpr int max_matrix_entry_magnitude = 6;
    static constexpr double max_matrix_entry = 1e6;

    /** \brief The minimal magnitude an arbitrary matrix entry may have
     *
     * If this is set to -6 that all entries will either be zero or have
     * an absolute value larger than $10^{-6}$.
     */
    static constexpr int min_matrix_entry_magnitude = -6;
    static constexpr double min_matrix_entry = 1e-6;

    //
    // Numerical tolerance for numerical comparison:
    //
    /** The default numeric tolerance when comparing computed values
     * for equality. Usually comparison operations will interpret
     * this as tolerance towards the *relative* error
     */
    static constexpr double default_num_tol = 1e-12;
};
}
}
