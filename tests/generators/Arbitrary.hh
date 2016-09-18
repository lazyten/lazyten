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
#include <linalgwrap/ArmadilloMatrix.hh>
#include <linalgwrap/LazyMatrixWrapper.hh>
#include <linalgwrap/TestingUtils.hh>
#include <rapidcheck.h>

namespace rc {
template <typename Scalar>
struct Arbitrary<::linalgwrap::ArmadilloMatrix<Scalar>> {
    typedef ::linalgwrap::ArmadilloMatrix<Scalar> small_matrix_type;
    static Gen<small_matrix_type> arbitrary() {
        return linalgwrap::gen::numeric_tensor<small_matrix_type>();
    };
};

template <typename StoredMatrix, typename InnerMatrix>
struct Arbitrary<::linalgwrap::LazyMatrixWrapper<StoredMatrix, InnerMatrix>> {
    typedef typename ::linalgwrap::LazyMatrixWrapper<StoredMatrix, InnerMatrix>
          generated_type;
    typedef typename generated_type::size_type size_type;
    typedef typename generated_type::scalar_type scalar_type;

    static Gen<generated_type> arbitrary() {
        // Define a lambda that returns a new matrix object.
        auto callable = [] {
            // Generate an arbitary inner matrix
            auto inner = *gen::arbitrary<InnerMatrix>();

            // Enwrap and return
            return ::linalgwrap::LazyMatrixWrapper<StoredMatrix, InnerMatrix>(
                  std::move(inner));
        };

        // Return the callable wrapped in gen::exec.
        return gen::exec(callable);
    }
};
}
