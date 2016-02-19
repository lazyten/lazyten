#pragma once
#include "TestConstants.hh"
#include <cmath>
#include <algorithm>
#include <Matrix_i.hh>

#include <iostream>

namespace linalgwrap {
namespace tests {

/* Do numeric-error aware comparison of certain objects. */
struct NumComp {
    /** \brief Functor to check that two values are numerically equal
     *
     * \param tolerance    The tolerance in absolute error
     *                     (if we are comparing against zero)
     *                     or relative error (else).
     */
    template <typename Scalar>
    class NumEqual {
        // TODO derive this off binary functor!

      public:
        NumEqual(double tolerance = TestConstants::default_num_tol)
              : m_tolerance(tolerance) {}

        bool operator()(const Scalar& lhs, const Scalar& rhs) const {
            // TODO Make a constexpr in C++1y
            auto min = std::min(lhs, rhs);
            auto minabs = std::fabs(min);

            if (minabs < m_tolerance) {
                // We are comparing against zero
                // Check absolute error
                return std::fabs(lhs - rhs) < m_tolerance;
            } else {
                // We are comparing two number:
                // Check relative error
                return (std::fabs(lhs - rhs) < m_tolerance * minabs);
            }
        }

      private:
        double m_tolerance;
    };

    /* Specialisation of NumEqual for Matrices. */
    template <typename Scalar>
    class NumEqualMatrix {
        // TODO derive off binary functor as well.

      public:
        NumEqualMatrix(double tolerance = TestConstants::default_num_tol)
              : m_tolerance(tolerance) {}

        bool operator()(const Matrix_i<Scalar>& lhs,
                        const Matrix_i<Scalar>& rhs) const {
            // TODO Make a constexpr in C++1y

            typedef typename Matrix_i<Scalar>::size_type size_type;

            // Check sizes:
            if (lhs.n_rows() != rhs.n_rows()) return false;
            if (lhs.n_cols() != rhs.n_cols()) return false;

            // Calculate number of entries:
            size_type entries = lhs.n_rows() * lhs.n_cols();

            // Check all entries:
            for (size_type i = 0; i < entries; ++i) {
                if (!is_equal(lhs[i], rhs[i], m_tolerance)) return false;
            }

            return true;
        }

      private:
        double m_tolerance;
    };

    /** \brief Check that two values are numerically equal
     *
     * \param tolerance    The tolerance in absolute error
     *                     (if we are comparing floating point
     *                     entries against zero)
     *                     or relative error (else).
     */
    template <typename Scalar>
    static bool is_equal(
          const Scalar& lhs, const Scalar& rhs,
          const double tolerance = TestConstants::default_num_tol) {
        NumEqual<Scalar> eq_comp(tolerance);
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
          const double tolerance = TestConstants::default_num_tol) {
        NumEqualMatrix<Scalar> eq_comp(tolerance);
        return eq_comp(lhs, rhs);
    }
};
}
}
