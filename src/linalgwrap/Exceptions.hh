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

#ifndef LINALGWRAP_EXCEPTIONS_H_
#define LINALGWRAP_EXCEPTIONS_H_

#include "linalgwrap/ExceptionSystem.hh"
#include <complex>

namespace linalgwrap {

//
// Numerics
//
/**
 * Exception denoting a division by zero.
 */
DefExceptionMsg(ExcDevideByZero, "Devision by zero encountered.");

/**
 * Exception denoting that a NaN or plus or minus infinity was encountered.
 * Use ExcDevideByZero to indicate division by zero and use ExcZero to denote
 * an unexpected zero
 *
 * The argument is of type std::complex<double>, so any real or complex scalar
 * value can be passed.
 */
DefException1(ExcNumberNotFinite, std::complex<double>,
              << "Encoutered a non-finite number, where this was not expected."
                 "(its value is"
              << arg1 << ").");

/**
* A number is zero, but it should not be here.
*/
DefExceptionMsg(ExcZero, "Encountered a zero, where this does not make sense.");

/** An overflow was detected. */
DefExceptionMsg(ExcOverflow, "Overflow detected.");

/** An underflow was detected. */
DefExceptionMsg(ExcUnderflow, "Underflow detected.");

//
// Matrix properties
//

/** A matrix that was expected to be Hermetian is not */
DefExceptionMsg(ExcMatrixNotHermitian,
                "Encountered a non-Hermetian matrix where a Hermetian matrix "
                "was expected");

/** A matrix that was expected to be Hermetian is not */
DefExceptionMsg(ExcMatrixNotSymmetric,
                "Encountered a non-Symmetric matrix where a Symmetric matrix "
                "was expected");

/** A matrix which was expected to be dense, contained sparsity */
DefExceptionMsg(
      ExcMatrixNotDense,
      "Encountered a non-dense matrix where a dense matrix was expected");

/** A matrix which was expected to be quadratic, was not */
DefExceptionMsg(
      ExcMatrixNotSquare,
      "Encountered a non-square matrix where a square matrix was expected");

//
// Range and size checking
//
/**
 * The sizes of two objects are assumed to be equal, but they were not.
 */
DefException2(ExcSizeMismatch, size_t, size_t, << "Size " << arg1
                                               << " not equal to " << arg2);

/**
 * Exception to indicate that a number is not within the expected range.
 * The constructor takes three <tt>size_t</tt> arguments
 *
 * <ol>
 * <li> the violating index
 * <li> the lower bound
 * <li> the upper bound plus one
 * </ol>
 */
template <typename T>
DefException3(ExcOutsideRange, T, T, T, << "Index " << arg1
                                        << " is not in the half-open interval ["
                                        << arg2 << "," << arg3 << ").");

/**
 * Exception to indicate that a number is larger than an upper bound.
 * The intention is that arg1 <= arg2 should have been satisfied.
 */
template <typename T>
DefException2(ExcTooLarge, T, T, << "Number " << arg1
                                 << " must be smaller or equal to " << arg2
                                 << ".");

/**
 * Exception to indicate that a number is larger or equal to an upper bound
 * The intention is that arg1 < arg2 should have been satisfied.
 */
template <typename T>
DefException2(ExcTooLargeOrEqual, T, T,
              << "Number " << arg1 << " must be smaller than " << arg2 << ".");
/**
 * Exception to indicate that two numbers, which should be exactly equal are
 * not.
 */
template <typename T>
DefException2(ExcNotEqual, T, T, << "Number " << arg1 << " must be equal to "
                                 << arg2 << ".");

//
// Program logic
//
/**
 * Exception denoting a function or functionality has not been implemented
 * either because the programmer originally thought this was too difficult
 * to do or because it was not needed.
 *
 * This should not be used to indicate that something is missing and should
 * be implemented if neccessary.
 */
DefExceptionMsg(ExcNotImplemented,
                "This functionality has not been implemented yet. "
                "Feel free to take a look and implement it.");

/**
 * This exception is used if some object is found uninitialized.
 */
DefExceptionMsg(ExcNotInitialised,
                "The object you attempt to use is not yet initialised.");

/**
 * The object is in a state not suitable for this operation.
 */
DefException1(ExcInvalidState, char *,
              << "The object you attempt to use is not in a valid state: "
              << arg1);

/**
 * The pointer is invalid (e.g. a nullptr) and should not be used.
 */
DefExceptionMsg(ExcInvalidPointer,
                "The pointer you attempt to use is not valid. Check you are "
                "not dereferencing a nullpointer or similar.");

/**
 * Internal error occurred inside a routine
 */
DefExceptionMsg(ExcInternalError,
                "An assertion inside an internal routine has failed. "
                "This is a bug and should not have happened.");

/**
 * The calling of this function was deliberately disabled for some reason (which
 * is given here).
 */
DefException1(ExcDisabled, char *,
              << "The method you attempt to call has been disabled: " << arg1);

/**
 * This is thrown if an iterator should be incremented, decremented
 * or used, but it is already at its final state.
 */
DefExceptionMsg(ExcIteratorPastEnd,
                "You are trying to use an iterator, which is pointing past "
                "the end of its range of valid elements. It is not valid to "
                "dereference or use an iterator in such a case.");

//
// IO and interaction with OS
//

/**
 * Generic IO exception, use <tt>ExcFileNotOpen</tt> to specifically indicate
 * that
 * opening a file for reading or writing has failed.
 */
DefExceptionMsg(ExcIO, "An input/output error has occurred.");
/**
 * An error occurred opening the named file.
 *
 * The constructor takes a single argument of type <tt>char*</tt> naming the
 * file.
 */
DefException1(ExcFileNotOpen, char *, << "Could not open file " << arg1);

//
// Helper macros
//
/**
 * Uses assert_dbg in order to check that a number is within a range,
 * if not raise an exception.
 *
 * The macro takes the following arguments:
 * <ol>
 * <li> The lower bound
 * <li> The number to check
 * <li> The upper bound plus one
 * </ol>
 */
#define assert_range(start, number, end)                            \
    {                                                               \
        assert_dbg((start <= number) && (number < end),             \
                   ::linalgwrap::ExcOutsideRange<decltype(number)>( \
                         number, start, end));                      \
    }

/**
 * Uses assert_dbg in order to check that a rhs is greater or equal
 * to a lhs
 *
 * Takes the following arguments
 * <ol>
 * <li> lhs number
 * <li> rhs number
 * </ol>
 */
#define assert_greater_equal(lhs, rhs)                                  \
    {                                                                   \
        assert_dbg(lhs <= rhs,                                          \
                   ::linalgwrap::ExcTooLarge<decltype(lhs)>(lhs, rhs)); \
    }

/**
 * Uses assert_dbg in order to check that a rhs is strictly greater to
 * a lhs
 *
 * Takes the following arguments
 * <ol>
 * <li> lhs number
 * <li> rhs number
 * </ol>
 */
#define assert_greater(lhs, rhs)                                               \
    {                                                                          \
        assert_dbg(lhs < rhs,                                                  \
                   ::linalgwrap::ExcTooLargeOrEqual<decltype(lhs)>(lhs, rhs)); \
    }

/**
 * Uses assert_dbg in order to check that a rhs is exactly equal to a lhs.
 *
 * Takes the following arguments
 * <ol>
 * <li> lhs number
 * <li> rhs number
 * </ol>
 */
#define assert_equal(lhs, rhs)                                          \
    {                                                                   \
        assert_dbg(lhs == rhs,                                          \
                   ::linalgwrap::ExcNotEqual<decltype(lhs)>(lhs, rhs)); \
    }

/**
 * Assert whether two sizes match
 */
#define assert_size(size1, size2) \
    { assert_dbg(size1 == size2, ::linalgwrap::ExcSizeMismatch(size1, size2)); }

/**
 * Assert whether all elements of the vector vec have the size
 */
#define assert_element_sizes(vec, vsize)                             \
    {                                                                \
        for (auto it = std::begin(vec); it != std::end(vec); ++it) { \
            assert_size(vsize, it->size());                          \
        }                                                            \
    }

/**
 * Assert that a value is finite
 */
#define assert_finite(value)                                \
    {                                                       \
        assert_dbg(std::isfinite(value),                    \
                   ::linalgwrap::ExcNumberNotFinite(value)) \
    }

}  // linalgwrap

#endif
