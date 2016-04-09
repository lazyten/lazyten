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

#include "../TestConstants.hh"
#include "MatrixElement.hh"
#include <linalgwrap/ArmadilloMatrix.hh>
#include <rapidcheck.h>

namespace rc {
template <typename Scalar>
struct Arbitrary<::linalgwrap::ArmadilloMatrix<Scalar>> {
    typedef ::linalgwrap::ArmadilloMatrix<Scalar> small_matrix_type;
    typedef typename small_matrix_type::size_type size_type;
    typedef typename small_matrix_type::scalar_type scalar_type;

    static Gen<::linalgwrap::ArmadilloMatrix<Scalar>> arbitrary() {
        // Define a lambda that returns a new matrix object.
        auto callable = [] {
            const size_type maxsize =
                  ::linalgwrap::tests::TestConstants::max_matrix_size;

            // rhs of inRange is exclusive
            const auto n_rows = *gen::inRange<size_type>(1, maxsize + 1);
            const auto n_cols = *gen::inRange<size_type>(1, maxsize + 1);

            small_matrix_type matrix(n_rows, n_cols, false);

            // set to arbitrary values
            for (size_type i = 0; i < matrix.n_rows() * matrix.n_cols(); ++i) {
                matrix[i] = *gen::matrix_element<scalar_type>();
            }

            return matrix;
        };

        // Return the callable wrapped in gen::exec.
        return gen::exec(callable);
    };
};
}
