//
// Copyright (C) 2016-17 by the lazyten authors
//
// This file is part of lazyten.
//
// lazyten is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// lazyten is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with lazyten. If not, see <http://www.gnu.org/licenses/>.
//

#pragma once
#include <lazyten/Armadillo/ArmadilloMatrix.hh>
#include <lazyten/LazyMatrixWrapper.hh>
#include <lazyten/TestingUtils.hh>
#include <rapidcheck.h>

namespace rc {
// TODO deprecate these functions some day in favour for gen::numeric_tensor and
// gen::construct

template <typename Scalar>
struct Arbitrary<::lazyten::ArmadilloMatrix<Scalar>> {
  typedef ::lazyten::ArmadilloMatrix<Scalar> small_matrix_type;
  static Gen<small_matrix_type> arbitrary() {
    return lazyten::gen::numeric_tensor<small_matrix_type>();
  };
};

template <typename StoredMatrix>
struct Arbitrary<::lazyten::LazyMatrixWrapper<StoredMatrix>> {
  typedef typename ::lazyten::LazyMatrixWrapper<StoredMatrix> generated_type;
  static Gen<generated_type> arbitrary() {
    return gen::map(lazyten::gen::numeric_tensor<StoredMatrix>(),
                    [](StoredMatrix m) { return generated_type(std::move(m)); });
  }
};
}
