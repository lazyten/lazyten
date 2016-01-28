#pragma once

#include <SmallMatrix.hh>
#include <rapidcheck.h>
#include "../tests_configuration.hh"

namespace rc {
template <typename Scalar>
struct Arbitrary<::linalgwrap::SmallMatrix<Scalar>> {
    static Gen<::linalgwrap::SmallMatrix<Scalar>> arbitrary() {
        // Define a lambda that returns a new matrix object.
        auto callable = [] {
            typedef ::linalgwrap::SmallMatrix<Scalar> small_matrix_type;
            typedef typename small_matrix_type::size_type size_type;
            typedef typename small_matrix_type::scalar_type scalar_type;

            const scalar_type maxsize =
                  ::linalgwrap::tests::tests_constants::max_matrix_size + 1;

            const auto n_rows = *gen::inRange<size_type>(1, maxsize);
            const auto n_cols = *gen::inRange<size_type>(1, maxsize);

            small_matrix_type matrix(n_rows, n_cols, false);

            // set to arbitrary values
            for (size_type i = 0; i < matrix.n_rows() * matrix.n_cols(); ++i) {
                matrix[i] = *rc::gen::arbitrary<scalar_type>();
            }

            return matrix;
        };

        // Return the callable wrapped in gen::exec.
        return gen::exec(callable);
    };
};
}
