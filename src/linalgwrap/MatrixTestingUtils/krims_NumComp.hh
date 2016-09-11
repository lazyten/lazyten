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
#include <krims/NumComp.hh>
#include <linalgwrap/Matrix_i.hh>

// This file provides specialisations for the NumComp system to work with
// linalgwrap matrices out of the box.
// Include it in tests making use of krims/NumComp.hh in which you want to be
// able to numerically compare linalgwrap matrices as well.

// Note since this specialises the NumEqual class in krims/NumComp.hh this code
// is deliberately put into the krims namespace.
namespace krims {

/** \brief Functor to check that two matrices (of possibly different type)
 * are numerically equal
 */
template <typename Mat1, typename Mat2>
struct NumEqual<Mat1, Mat2,
                typename std::enable_if<
                      linalgwrap::IsMatrix<Mat1>::value &&
                      linalgwrap::IsMatrix<Mat2>::value &&
                      std::is_same<typename Mat1::scalar_type,
                                   typename Mat2::scalar_type>::value>::type> {
    typedef const Mat1& first_argument_type;
    typedef const Mat2& second_argument_type;
    typedef bool result_type;
    typedef typename RealTypeOf<typename Mat1::scalar_type>::type real_type;

    NumEqual(const real_type tolerance, const NumCompActionType failure_action)
          : m_tolerance(tolerance), m_failure_action(failure_action) {}

    bool operator()(const Mat1& lhs, const Mat2& rhs) const;

  private:
    static_assert(std::is_same<typename Mat1::size_type,
                               typename Mat2::size_type>::value,
                  "The size types of Mat1 and Mat2 do not agree.");
    typedef typename Mat1::size_type size_type;

    bool sizes_match(size_type lhs_value, size_type rhs_value,
                     std::string what) const;

    const real_type m_tolerance;
    const NumCompActionType m_failure_action;
};

//
// -----------------------------------------------------------------------
//

template <typename Mat1, typename Mat2>
bool NumEqual<Mat1, Mat2,
              typename std::enable_if<
                    linalgwrap::IsMatrix<Mat1>::value &&
                    linalgwrap::IsMatrix<Mat2>::value &&
                    std::is_same<typename Mat1::scalar_type,
                                 typename Mat2::scalar_type>::value>::type>::
      sizes_match(size_type lhs_value, size_type rhs_value,
                  std::string what) const {
    if (lhs_value == rhs_value) {
        return true;
    } else if (m_failure_action == NumCompActionType::ThrowNormal ||
               m_failure_action == NumCompActionType::ThrowVerbose) {

        const size_type diff = lhs_value < rhs_value ? rhs_value - lhs_value
                                                     : lhs_value - rhs_value;
        NumCompException<size_type> e(lhs_value, rhs_value, diff, 0, "==",
                                      "Size mismatch in number of " + what);
        e.add_exc_data(__FILE__, __LINE__, __PRETTY_FUNCTION__);
        throw e;
    } else {
        return false;
    }
}

template <typename Mat1, typename Mat2>
bool NumEqual<Mat1, Mat2,
              typename std::enable_if<
                    linalgwrap::IsMatrix<Mat1>::value &&
                    linalgwrap::IsMatrix<Mat2>::value &&
                    std::is_same<typename Mat1::scalar_type,
                                 typename Mat2::scalar_type>::value>::type>::
operator()(const Mat1& lhs, const Mat2& rhs) const {
    // First compare the sizes. If we get an exception just pass it upwards
    // appending the matrix values in case the user wants detailed throw
    // information.
    try {
        if (!sizes_match(lhs.n_rows(), rhs.n_rows(), "rows") ||
            !sizes_match(rhs.n_cols(), lhs.n_cols(), "columns")) {
            return false;
        }
    } catch (NumCompException<size_type>& e) {
        if (m_failure_action == NumCompActionType::ThrowVerbose) {
            std::stringstream ss;
            ss << "when comparing matrices " << std::endl
               << lhs << std::endl
               << "and" << std::endl
               << rhs << std::endl;
            e.append(ss.str());
            throw;
        }
    }

    // Now compare elements.
    // If one is not equal return false or catch the exception and amend
    // the data we are interested in before rethrowing.
    NumEqual<real_type, real_type> is_equal{m_tolerance, m_failure_action};
    for (size_type i = 0; i < lhs.n_rows(); ++i) {
        for (size_type j = 0; j < lhs.n_cols(); ++j) {
            try {
                if (!is_equal(lhs(i, j), rhs(i, j))) {
                    return false;
                }
            } catch (NumCompException<real_type>& e) {
                std::stringstream ss;
                ss << "Matrix entry (" << i << "," << j << ") not equal";
                if (m_failure_action == NumCompActionType::ThrowVerbose) {
                    ss << "when comparing matrices " << std::endl
                       << lhs << std::endl
                       << "and" << std::endl
                       << rhs << std::endl;
                } else {
                    ss << ".";
                }

                e.append(ss.str());
                throw;
            }
        }
    }
    return true;
}
}  // namespace krims
