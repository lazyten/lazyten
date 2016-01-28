#pragma once
#include <rapidcheck.h>
#include <StoredMatrix_i.hh>

namespace rc {

template <typename Matrix>
struct FixedSize {
    static_assert(
          std::is_base_of<
                ::linalgwrap::StoredMatrix_i<typename Matrix::scalar_type>,
                Matrix>::value,
          "Matrix must be a child class of StoredMatrix_i");

    typedef typename Matrix::size_type size_type;
    typedef typename Matrix::scalar_type scalar_type;

    static Gen<Matrix> fixed_size(size_type n_rows, size_type n_cols) {
        auto callable = [=] {
            // allocate memory, but don't initialise
            Matrix m(n_rows, n_cols, false);

            // set to arbitrary values
            for (size_type i = 0; i < m.n_rows() * m.n_cols(); ++i) {
                m[i] = *rc::gen::arbitrary<scalar_type>();
            }

            return m;
        };
        return gen::exec(callable);
    }
};
}
