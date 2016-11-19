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
#include <linalgwrap/LazyMatrixWrapper.hh>
#include <linalgwrap/TestingUtils.hh>
#include <rapidcheck.h>

namespace rc {
// TODO deprecate in favour for gen::numeric_tensor and gen::construct

template <typename Matrix>
struct FixedSizeMatrix {
    typedef typename Matrix::size_type size_type;
    static Gen<Matrix> fixed_size(size_type n_rows, size_type n_cols) {
        return linalgwrap::gen::numeric_tensor<Matrix>(n_rows, n_cols);
    }
};

template <typename StoredMatrix>
struct FixedSizeMatrix<::linalgwrap::LazyMatrixWrapper<StoredMatrix>> {
    typedef typename ::linalgwrap::LazyMatrixWrapper<StoredMatrix>
          generated_type;
    typedef typename generated_type::size_type size_type;

    static Gen<generated_type> fixed_size(size_type n_rows, size_type n_cols) {
        return gen::map(
              linalgwrap::gen::numeric_tensor<StoredMatrix>(n_rows, n_cols),
              [](StoredMatrix m) { return generated_type(std::move(m)); });
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
