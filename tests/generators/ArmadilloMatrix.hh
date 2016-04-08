#pragma once

#include "../TestConstants.hh"
#include "MatrixElement.hh"
#include <ArmadilloMatrix.hh>
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
