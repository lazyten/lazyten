#pragma once
#include "MatrixElement.hh"
#include <StoredMatrix_i.hh>
#include <rapidcheck.h>

namespace rc {

template <typename Matrix>
struct FixedSizeMatrix {
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
                m[i] = *gen::matrix_element<scalar_type>();
            }

            return m;
        };
        return gen::exec(callable);
    }
};

namespace gen {
template <typename Matrix>
Gen<Matrix> fixed_size(typename Matrix::size_type n_rows,
                       typename Matrix::size_type n_cols) {
    return FixedSizeMatrix<Matrix>::fixed_size(n_rows, n_cols);
}

}  // namespace gen
}  // namespace rc
