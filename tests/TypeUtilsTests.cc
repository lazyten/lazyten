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

#include <linalgwrap/LazyMatrixWrapper.hh>
#include <linalgwrap/SmallMatrix.hh>
#include <linalgwrap/SmallVector.hh>
#include <linalgwrap/TypeUtils.hh>

namespace linalgwrap {
namespace tests {
namespace type_utils_tests {

template <typename Scalar>
class Tests {
  typedef SmallMatrix<Scalar> mat_type;
  typedef SmallVector<Scalar> vec_type;
  typedef LazyMatrixWrapper<mat_type> lazy_mat_type;

  typedef krims::RealTypeOf<Scalar> real_type;
  static constexpr bool isComplex = krims::IsComplexNumber<Scalar>::value;

  static_assert(std::is_same<mat_type, typename StoredTypeOf<mat_type>::type>::value,
                "StoredTypeOf with stored matrix failed.");
  static_assert(std::is_same<mat_type, typename StoredTypeOf<lazy_mat_type>::type>::value,
                "StoredTypeOf with lazy matrix failed.");
  static_assert(HasComplexScalar<mat_type>::value || !isComplex,
                "HasComplexScalar with stored matrix");
  static_assert(HasComplexScalar<vec_type>::value || !isComplex,
                "HasComplexScalar with stored vector");
  static_assert(HasComplexScalar<lazy_mat_type>::value || !isComplex,
                "HasComplexScalar with lazy matrix");

  // References and const
  static_assert(std::is_same<mat_type, typename StoredTypeOf<mat_type&>::type>::value,
                "StoredTypeOf with stored matrix lvalue failed.");
  static_assert(
        std::is_same<mat_type, typename StoredTypeOf<lazy_mat_type&&>::type>::value,
        "StoredTypeOf with lazy matrix rvalue failed.");
  static_assert(HasComplexScalar<const mat_type>::value || !isComplex,
                "HasComplexScalar with const stored matrix");
  static_assert(HasComplexScalar<const vec_type&&>::value || !isComplex,
                "HasComplexScalar with const stored vector reference");
  static_assert(HasComplexScalar<lazy_mat_type&>::value || !isComplex,
                "HasComplexScalar with lazy matrix reference");
};

// Explicit instantiation of test cases:
template class Tests<double>;
template class Tests<float>;
template class Tests<std::complex<float>>;
template class Tests<std::complex<double>>;

}  // type_utils_tests
}  // tests
}  // linalgwrap
