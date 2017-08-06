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
#include <krims/NumComp.hh>
#include <lazyten/Base/Interfaces.hh>
#include <lazyten/Matrix_i.hh>
#include <lazyten/MultiVector.hh>

// This file provides specialisations for the NumComp system to work with
// lazyten matrices out of the box.
// Include it in tests making use of krims/NumComp.hh in which you want to be
// able to numerically compare lazyten matrices as well.

// Note since this specialises the NumEqual class in krims/NumComp.hh this code
// is deliberately put into the krims namespace.
namespace krims {

/** \brief Functor to check that two indexables are numerically equal
 */
template <typename Ible1, typename Ible2>
struct NumEqual<
      Ible1, Ible2,
      typename std::enable_if<
            lazyten::IsIndexable<Ible1>::value && lazyten::IsIndexable<Ible2>::value &&
            !lazyten::IsMatrix<Ible1>::value && !lazyten::IsMatrix<Ible2>::value &&
            !lazyten::IsVector<Ible1>::value && !lazyten::IsVector<Ible2>::value &&
            std::is_same<typename Ible1::scalar_type,
                         typename Ible2::scalar_type>::value>::type>
      : private NumEqualContainerBase<Ible1, Ible2> {
  typedef const Ible1& first_argument_type;
  typedef const Ible2& second_argument_type;
  typedef bool result_type;

  typedef NumEqualContainerBase<Ible1, Ible2> base_type;
  typedef typename base_type::real_type real_type;

  NumEqual(const real_type tolerance, const NumCompActionType failure_action)
        : base_type{tolerance, failure_action} {}

  bool operator()(const Ible1& lhs, const Ible2& rhs) const {
    return base_type::number_elem_match(lhs, rhs, "indexables") &&
           base_type::element_values_match(lhs, rhs, "indexables");
  }
};

/** \brief Functor to check that two matrices (of possibly different type)
 * are numerically equal
 */
template <typename Mat1, typename Mat2>
struct NumEqual<Mat1, Mat2,
                typename std::enable_if<
                      lazyten::IsMatrix<Mat1>::value && lazyten::IsMatrix<Mat2>::value &&
                      std::is_same<typename Mat1::scalar_type,
                                   typename Mat2::scalar_type>::value>::type>
      : private NumEqualContainerBase<Mat1, Mat2> {
  typedef const Mat1& first_argument_type;
  typedef const Mat2& second_argument_type;
  typedef bool result_type;

  typedef NumEqualContainerBase<Mat1, Mat2> base_type;
  typedef typename base_type::real_type real_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::value_type scalar_type;

  NumEqual(const real_type tolerance, const NumCompActionType failure_action)
        : base_type{tolerance, failure_action} {}

  bool operator()(const Mat1& lhs, const Mat2& rhs) const;
};

/** \brief Functor to check that two vectors (of possibly different type)
 * are numerically equal
 */
template <typename Vec1, typename Vec2>
struct NumEqual<Vec1, Vec2,
                typename std::enable_if<
                      lazyten::IsVector<Vec1>::value && lazyten::IsVector<Vec2>::value &&
                      std::is_same<typename Vec1::scalar_type,
                                   typename Vec2::scalar_type>::value>::type>
      : private NumEqualContainerBase<Vec1, Vec2> {
  typedef const Vec1& first_argument_type;
  typedef const Vec2& second_argument_type;
  typedef bool result_type;

  typedef NumEqualContainerBase<Vec1, Vec2> base_type;
  typedef typename base_type::real_type real_type;

  NumEqual(const real_type tolerance, const NumCompActionType failure_action)
        : base_type(tolerance, failure_action) {}

  bool operator()(const Vec1& lhs, const Vec2& rhs) const {
    return base_type::number_elem_match(lhs, rhs, "vectors") &&
           base_type::element_values_match(lhs, rhs, "vectors");
  }
};

/** \brief Functor to check that two multivectors
 * (with vectors of possibly different type)
 * are numerically equal
 */
template <typename Vec1, typename Vec2>
struct NumEqual<lazyten::MultiVector<Vec1>, lazyten::MultiVector<Vec2>,
                typename std::enable_if<
                      lazyten::IsVector<Vec1>::value && lazyten::IsVector<Vec2>::value &&
                      std::is_same<typename Vec1::scalar_type,
                                   typename Vec2::scalar_type>::value>::type>
      : private NumEqualContainerBase<Vec1, Vec2> {
  typedef const lazyten::MultiVector<Vec1>& first_argument_type;
  typedef const lazyten::MultiVector<Vec2>& second_argument_type;
  typedef bool result_type;

  typedef NumEqualContainerBase<Vec1, Vec2> base_type;
  typedef typename base_type::real_type real_type;
  typedef typename base_type::size_type size_type;

  NumEqual(const real_type tolerance, const NumCompActionType failure_action)
        : base_type(tolerance, failure_action) {}

  bool operator()(const lazyten::MultiVector<Vec1>& lhs,
                  const lazyten::MultiVector<Vec2>& rhs) const;
};

//
// -----------------------------------------------------------------------
//

template <typename Mat1, typename Mat2>
bool NumEqual<Mat1, Mat2,
              typename std::enable_if<
                    lazyten::IsMatrix<Mat1>::value && lazyten::IsMatrix<Mat2>::value &&
                    std::is_same<typename Mat1::scalar_type,
                                 typename Mat2::scalar_type>::value>::type>::
operator()(const Mat1& lhs, const Mat2& rhs) const {
  // First compare the sizes. If we get an exception just pass it upwards
  // appending the matrix values in case the user wants detailed throw
  // information.
  try {
    if (!base_type::sizes_match(lhs.n_rows(), rhs.n_rows(), "rows") ||
        !base_type::sizes_match(rhs.n_cols(), lhs.n_cols(), "columns")) {
      return false;
    }
  } catch (NumCompExceptionBase& e) {
    if (base_type::failure_action == NumCompActionType::ThrowVerbose) {
      std::stringstream ss;
      ss << " when comparing matrices " << std::endl
         << lhs << std::endl
         << "and" << std::endl
         << rhs << "." << std::endl;
      e.append_extra(ss.str());
    }
    throw;
  }

  // Now compare elements.
  // If one is not equal return false or catch the exception and amend
  // the data we are interested in before rethrowing.
  NumEqual<scalar_type, scalar_type> is_equal{base_type::tolerance,
                                              base_type::failure_action};
  for (size_type i = 0; i < lhs.n_rows(); ++i) {
    for (size_type j = 0; j < lhs.n_cols(); ++j) {
      try {
        if (!is_equal(lhs(i, j), rhs(i, j))) {
          return false;
        }
      } catch (NumCompExceptionBase& e) {
        std::stringstream ss;
        ss << "Matrix entry (" << i << "," << j << ") not equal";
        if (base_type::failure_action == NumCompActionType::ThrowVerbose) {
          ss << " when comparing matrices " << std::endl
             << lhs << std::endl
             << "and" << std::endl
             << rhs << "." << std::endl;
        } else {
          ss << ".";
        }
        e.append_extra(ss.str());
        throw;
      }
    }
  }
  return true;
}

template <typename Vec1, typename Vec2>
bool NumEqual<lazyten::MultiVector<Vec1>, lazyten::MultiVector<Vec2>,
              typename std::enable_if<
                    lazyten::IsVector<Vec1>::value && lazyten::IsVector<Vec2>::value &&
                    std::is_same<typename Vec1::scalar_type,
                                 typename Vec2::scalar_type>::value>::type>::
operator()(const lazyten::MultiVector<Vec1>& lhs,
           const lazyten::MultiVector<Vec2>& rhs) const {
  // First compare the sizes. If we get an exception just pass it upwards
  // appending the matrix values in case the user wants detailed throw
  // information.
  try {
    if (!base_type::sizes_match(lhs.n_vectors(), rhs.n_vectors(), "vectors")) {
      return false;
    }
  } catch (NumCompExceptionBase& e) {
    if (base_type::failure_action == NumCompActionType::ThrowVerbose) {
      std::stringstream ss;
      ss << "when comparing MultiVectors " << std::endl
         << lhs << std::endl
         << "and" << std::endl
         << rhs << "." << std::endl;
      e.append_extra(ss.str());
    }
    throw;
  }

  // Now compare vectors
  // If one is not equal return false or catch the exception and amend
  // the data we are interested in before rethrowing.
  NumEqual<Vec1, Vec2> is_equal{base_type::tolerance, base_type::failure_action};
  for (size_type i = 0; i < lhs.n_vectors(); ++i) {
    try {
      if (!is_equal(lhs[i], rhs[i])) {
        return false;
      }
    } catch (NumCompExceptionBase& e) {
      std::stringstream ss;
      ss << "Vector (" << i << ") not equal";
      if (base_type::failure_action == NumCompActionType::ThrowVerbose) {
        ss << "when comparing MultiVectors " << std::endl
           << lhs << std::endl
           << "and" << std::endl
           << rhs << "." << std::endl;
      } else {
        ss << ".";
      }
      e.append_extra(ss.str());
      throw;
    }
  }
  return true;
}

}  // namespace krims
