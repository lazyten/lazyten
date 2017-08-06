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

#include <lazyten/LazyMatrixWrapper.hh>
#include <lazyten/SmallMatrix.hh>
#include <lazyten/SmallVector.hh>
#include <lazyten/TypeUtils.hh>

namespace lazyten {
namespace tests {
namespace type_utils_tests {

template <typename Scalar>
class Tests {
  typedef SmallMatrix<Scalar> mat_type;
  typedef SmallVector<Scalar> vec_type;
  typedef LazyMatrixWrapper<mat_type> lazy_mat_type;

  typedef krims::RealTypeOf<Scalar> real_type;
  static constexpr bool is_complex = krims::IsComplexNumber<Scalar>::value;

  static_assert(std::is_same<mat_type, typename StoredTypeOf<mat_type>::type>::value,
                "StoredTypeOf with stored matrix failed.");
  static_assert(std::is_same<mat_type, typename StoredTypeOf<lazy_mat_type>::type>::value,
                "StoredTypeOf with lazy matrix failed.");
  static_assert(HasComplexScalar<mat_type>::value || !is_complex,
                "HasComplexScalar with stored matrix");
  static_assert(HasComplexScalar<vec_type>::value || !is_complex,
                "HasComplexScalar with stored vector");
  static_assert(HasComplexScalar<lazy_mat_type>::value || !is_complex,
                "HasComplexScalar with lazy matrix");

  // References and const
  static_assert(std::is_same<mat_type, typename StoredTypeOf<mat_type&>::type>::value,
                "StoredTypeOf with stored matrix lvalue failed.");
  static_assert(
        std::is_same<mat_type, typename StoredTypeOf<lazy_mat_type&&>::type>::value,
        "StoredTypeOf with lazy matrix rvalue failed.");
  static_assert(HasComplexScalar<const mat_type>::value || !is_complex,
                "HasComplexScalar with const stored matrix");
  static_assert(HasComplexScalar<const vec_type&&>::value || !is_complex,
                "HasComplexScalar with const stored vector reference");
  static_assert(HasComplexScalar<lazy_mat_type&>::value || !is_complex,
                "HasComplexScalar with lazy matrix reference");
};

// Explicit instantiation of test cases:
template class Tests<double>;
template class Tests<float>;
template class Tests<std::complex<float>>;
template class Tests<std::complex<double>>;

}  // namespace type_utils_tests
}  // namespace tests
}  // namespace lazyten
