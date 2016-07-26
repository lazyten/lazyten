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
#include "TestConstants.hh"
#include <algorithm>
#include <cmath>
#include <linalgwrap/Matrix_i.hh>

#include <iostream>

namespace linalgwrap {
namespace tests {

/** Exception raised by the NumComp operations if the objects are not
 *  numerically equal. */
struct NumCompException : public std::exception {
  public:
    NumCompException(double lhs_, double rhs_, double error_, double tolerance_,
                     const std::string description_) noexcept;
    const char* what() const noexcept;

    //! The value of the lhs
    double lhs;

    //! The value of the rhs
    double rhs;

    //! The error that was obtained
    double error;

    //! The tolerance we applied
    double tolerance;

    //! The description that was addititonally supplied
    const std::string description;

  private:
    mutable std::string what_str;
};

/* Do numeric-error aware comparison of certain objects. */
struct NumComp {
    template <typename Scalar>
    class NumEqual {
      public:
        typedef Scalar first_argument_type;
        typedef Scalar second_argument_type;
        typedef bool result_type;

        /** \brief Functor to check that two values are numerically equal
         *
         * \param tolerance    The tolerance in absolute error
         *                     (if we are comparing against zero)
         *                     or relative error (else).
         *
         * \param do_throw    Should the comparison throw an exception
         *                    if its not valid or just return false.
         *                    (default is true)
         *
         */
        NumEqual(double tolerance = TestConstants::default_num_tol,
                 bool do_throw = true)
              : m_tolerance(tolerance), m_throw(do_throw) {}

        bool operator()(const Scalar& lhs, const Scalar& rhs) const {
            // TODO Make a constexpr in C++1y
            const Scalar minabs = std::min(std::fabs(lhs), std::fabs(rhs));

            // Do we compare against zero?
            const bool compare_zero = minabs < m_tolerance;

            // Absolute error between the two numbers
            const double abserror = std::fabs(lhs - rhs);

            // Error value used to qualify the agreement between the numbers
            // (absolute for comparing against zero, else relative)
            const double error = abserror / (compare_zero ? 1. : minabs);

            // Are the numbers equal (including numerical tolerance?
            bool equal = compare_zero ?
                                      // Check absolute error
                               abserror < m_tolerance
                                      :
                                      // Check relative error
                               abserror < (m_tolerance * minabs);

            if (!equal && m_throw) {
                throw NumCompException(lhs, rhs, error, m_tolerance, "");
            }
            return equal;
        }

      private:
        const double m_tolerance;
        const bool m_throw;
    };

    /* Specialisation of NumEqual for Matrices. */
    template <typename Scalar>
    class NumEqualMatrix {
        typedef Matrix_i<Scalar> first_argument_type;
        typedef Matrix_i<Scalar> second_argument_type;
        typedef bool result_type;

      public:
        /** Construct a NumEqualMatrix object
         *
         * \param tolerance   The tolerance to compare with
         * \param do_throw    Should the comparison throw an exception
         *                    if its not valid or just return false.
         *                    (default is true)
         * \param verbose_throw  Be extra verbose on a throw (e.g. print the
         *                       full matrices instead of just the problematic
         *                       value.
         *
         */
        NumEqualMatrix(double tolerance = TestConstants::default_num_tol,
                       bool do_throw = true, bool verbose_throw = false)
              : m_tolerance(tolerance),
                m_throw(do_throw),
                m_verbose_throw(verbose_throw) {}

        bool operator()(const Matrix_i<Scalar>& lhs,
                        const Matrix_i<Scalar>& rhs) const {
            // TODO Make a constexpr in C++1y

            typedef typename Matrix_i<Scalar>::size_type size_type;

            // Check sizes:
            if (lhs.n_rows() != rhs.n_rows()) return false;
            if (lhs.n_cols() != rhs.n_cols()) return false;

            // Check all entries:
            // TODO use matrix iterator
            for (size_type i = 0; i < lhs.n_rows(); ++i) {
                for (size_type j = 0; j < lhs.n_cols(); ++j) {
                    try {
                        if (!is_equal(lhs(i, j), rhs(i, j), m_tolerance,
                                      m_throw))
                            return false;
                    } catch (const NumCompException& e) {
                        std::stringstream ss;
                        ss << "Matrix entry (" << i << "," << j
                           << ") not equal";
                        if (m_verbose_throw) {
                            ss << "when comparing matrices " << std::endl
                               << lhs << std::endl
                               << "and" << std::endl
                               << rhs << std::endl;
                        } else {
                            ss << ".";
                        }

                        throw NumCompException(e.lhs, e.rhs, e.error,
                                               e.tolerance, ss.str());
                    }
                }
            }

            return true;
        }

      private:
        const double m_tolerance;
        // TODO unify into a "throw-type" of none, normal and verbose.
        const bool m_throw;
        const bool m_verbose_throw;
    };

    /** \brief Check that two values are numerically equal
     *
     * \param tolerance    The tolerance in absolute error
     *                     (if we are comparing floating point
     *                     entries against zero)
     *                     or relative error (else).
     * \param do_throw    Should the comparison throw an exception
     *                    if its not valid or just return false.
     *                    (default is true)
     */
    template <typename Scalar>
    static bool is_equal(
          const Scalar& lhs, const Scalar& rhs,
          const double tolerance = TestConstants::default_num_tol,
          const bool do_throw = true) {
        NumEqual<Scalar> eq_comp(tolerance, do_throw);
        return eq_comp(lhs, rhs);
    }

    /** \brief Check that two matrices are equivalent p to numeric
     * error.
     *
     * Uses is_equal above to perform the equality test on each value.
     */
    template <typename Scalar>
    static bool is_equal_matrix(
          const Matrix_i<Scalar>& lhs, const Matrix_i<Scalar>& rhs,
          const double tolerance = TestConstants::default_num_tol,
          const bool do_throw = true, const bool verbose_throw = false) {
        NumEqualMatrix<Scalar> eq_comp(tolerance, do_throw, verbose_throw);
        return eq_comp(lhs, rhs);
    }
};
}
}
