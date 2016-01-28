#ifndef LINALGWRAP_EXCEPTIONS_H_
#define LINALGWRAP_EXCEPTIONS_H_

#include "ExceptionSystem.hh"
#include <complex>

namespace linalgwrap {

/** Basic exceptions */
namespace exceptions {

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
                 "(its value is" << arg1 << ").");

/**
* A number is zero, but it should not be here.
*/
DefExceptionMsg(ExcZero, "Encountered a zero, where this does not make sense.");

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
 * Similar to <tt>ExcOutsideRange</tt>, the number given should actually
 * be the upper bound plus one.
 */
template <typename T>
DefException2(ExcAboveUpperBound, T, T,
              << "Number " << arg1 << " must be smaller than " << arg2 << ".");

/**
 * Exception to indicate that a number is below a lower bound.
 * Similar to <tt>ExcOutsideRange</tt>, the number given should actually
 * be the lower bound.
 */
template <typename T>
DefException2(ExcBelowLowerBound, T, T, << "Number " << arg1
                                        << " must be larger or equal " << arg2
                                        << ".");

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
DefExceptionMsg(ExcInvalidState,
                "The object you attempt to use is not in a valid state.");

/**
 * Internal error ocurred inside a routine
 */
DefExceptionMsg(ExcInternalError,
                "An assertion inside an internal routine has failed. "
                "This is a bug and should not have happened.");

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
// Stuff from deal.ii
//
//
//  /**
//   * Trying to allocate a new object failed due to lack of free memory.
//   */
//  DeclExceptionMsg (ExcOutOfMemory,
//                    "Your program tried to allocate some memory but this "
//                    "allocation failed. Typically, this either means that "
//                    "you simply do not have enough memory in your system, "
//                    "or that you are (erroneously) trying to allocate "
//                    "a chunk of memory that is simply beyond all reasonable "
//                    "size, for example because the size of the object has "
//                    "been computed incorrectly.");
//
//  /**
//   * A memory handler reached a point where all allocated objects should have
//   * been released. Since this exception is thrown, some were still allocated.
//   */
//  DeclException1 (ExcMemoryLeak, int,
//                  << "Destroying memory handler while " << arg1
//                  << " objects are still allocated");
//
//  /**
//   * The object should have been filled with something before this member
//   * function is called.
//   */
//  DeclExceptionMsg(ExcEmptyObject,
//                   "The object you are trying to access is empty but it makes
//                   "
//                   "no sense to attempt the operation you are trying on an "
//                   "empty object.");
//
//  /**
//   * This exception is thrown if the iterator you access has corrupted data.
//   * It might for instance be, that the container it refers does not have an
//   * entry at the point the iterator refers.
//   *
//   * Typically, this will be an internal error of deal.II, because the
//   * increment and decrement operators should never yield an invalid iterator.
//   */
//  DeclExceptionMsg (ExcInvalidIterator,
//                    "You are trying to use an iterator, but the iterator is "
//                    "in an invalid state. This may indicate that the iterator
//                    "
//                    "object has not been initialized, or that it has been "
//                    "moved beyond the end of the range of valid elements.");
//
//  /**
//   * This exception is thrown if the iterator you incremented or decremented
//   * was already at its final state.
//   */
//  DeclExceptionMsg (ExcIteratorPastEnd,
//                    "You are trying to use an iterator, but the iterator is "
//                    "pointing past the end of the range of valid elements. "
//                    "It is not valid to dereference the iterator in this "
//                    "case.");
//
//  /**
//   * Some of our numerical classes allow for setting all entries to zero
//   * using the assignment operator <tt>=</tt>.
//   *
//   * In many cases, this assignment operator makes sense <b>only</b> for the
//   * argument zero. In other cases, this exception is thrown.
//   */
//  DeclExceptionMsg (ExcScalarAssignmentOnlyForZeroValue,
//                    "You are trying an operation of the form 'vector=s' with "
//                    "a nonzero scalar value 's'. However, such assignments "
//                    "are only allowed if the right hand side is zero.");

}  // exceptions

// Import exceptions namespace
using namespace exceptions;

//
// Helper macros
//
/**
 * Uses assert_dbg in order to check that a number is within a range,
 * if not raise an exception.
 *
 * The macro takes the following arguments:
 * <ol>
 * <li> The number to check
 * <li> The lower bound
 * <li> The upper bound plus one
 * </ol>
 */
#define assert_range(number, start, end)                            \
    {                                                               \
        assert_dbg((start <= number) && (number < end),             \
                   ::linalgwrap::ExcOutsideRange<decltype(number)>( \
                         number, start, end));                      \
    }

/**
 * Uses assert_dbg in order to check that a lower bound is satisfied.
 *
 * Takes the following arguments
 * <ol>
 * <li> The number to check
 * <li> The lower bound
 * </ol>
 */
#define assert_lower_bound(number, start)                                      \
    {                                                                          \
        assert_dbg(start <= number,                                            \
                   ::linalgwrap::ExcBelowLowerBound<decltype(number)>(number,  \
                                                                      start)); \
    }

/**
 * Uses assert_dbg in order to check that an upper bound is satisfied.
 *
 * Takes the following arguments
 * <ol>
 * <li> The number to check
 * <li> The upper bound plus one
 * </ol>
 */
#define assert_upper_bound(number, end)                                       \
    {                                                                         \
        assert_dbg(number < end,                                              \
                   ::linalgwrap::ExcAboveUpperBound<decltype(number)>(number, \
                                                                      end));  \
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

}  // linalgwrap

#endif
