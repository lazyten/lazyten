#ifndef LINALGWRAP_EXCEPTION_SYSTEM_H_
#define LINALGWRAP_EXCEPTION_SYSTEM_H_
#include "ExceptionBase.hh"
#include <string>
#include <iostream>

namespace linalgwrap {

namespace exceptions {

/** Options how to deal with exceptions.
 *
 * THROW: Throw the exception
 * ABORT: Print it and abort the program
 */
enum class ExceptionEffect { THROW, ABORT };

/** How the assert_dbg macro deals with failing assertions.
 *  If such an assertion fails the macro will create an exception
 *  (second argument) and then process the latter as indicated with
 *  this property.
 *
 *  @note Since the assert_dbg macro is only active in DEBUG mode
 *        Changing this has no effect for RELEASE builds.
 */
static ExceptionEffect assert_dbg_effect = ExceptionEffect::ABORT;

}  // namespace exceptions

//
// Assert Macros
//
/** Assert a condition and if it fails (evaluates to false), generate an
 * exception (the 2nd argument). Deal with this exception as the current
 * value of assert_dbg_effect indicates.
 *
 * @note Active in DEBUG mode only.
 */
#ifdef DEBUG
#define assert_dbg(cond, exception)                                        \
    {                                                                      \
        using namespace linalgwrap::exceptions;                            \
        if (!(cond)) {                                                     \
            auto e = exception;                                            \
            e.add_exc_data(__FILE__, __LINE__, __PRETTY_FUNCTION__, #cond, \
                           #exception);                                    \
            if (assert_dbg_effect == ExceptionEffect::ABORT) {             \
                std::cerr << e.what() << std::endl;                        \
                std::abort();                                              \
            } else {                                                       \
                throw e;                                                   \
            }                                                              \
        }                                                                  \
    }
#else
#define assert_dbg(cond, exc) \
    {}
#endif

/** Assert a condition and if it evaluates to false, throw the exception
 *  given as the second argument.
 *
 * @note Active in DEBUG and RELEASE mode.
 */
#define assert_throw(cond, exception)                                      \
    {                                                                      \
        if (!(cond)) {                                                     \
            auto e = exception;                                            \
            e.add_exc_data(__FILE__, __LINE__, __PRETTY_FUNCTION__, #cond, \
                           #exception);                                    \
            throw e;                                                       \
        }                                                                  \
    }

/** This macro is used for actual errors that should always abort the program.
 *
 * @note Active in DEBUG and RELEASE mode.
 */
#define assert_abort(cond, exception)                                      \
    {                                                                      \
        using namespace linalgwrap::exceptions;                            \
        if (!(cond)) {                                                     \
            auto e = exception;                                            \
            e.add_exc_data(__FILE__, __LINE__, __PRETTY_FUNCTION__, #cond, \
                           #exception);                                    \
            std::cerr << e.what() << std::endl;                            \
            std::abort();                                                  \
        }                                                                  \
    }

//
// DefException Macros
//
/**
 * Define an exception with a string parameter. If no string is supplied
 * on exception construction, the default string is printed when the exception
 * is raised.
 */
#define DefExceptionMsg(Exception, defaulttext)                         \
    class Exception : public ::linalgwrap::ExceptionBase {              \
      public:                                                           \
        Exception(const std::string& msg = defaulttext) : m_msg(msg) {} \
        virtual ~Exception() noexcept {}                                \
        virtual void print_extra(std::ostream& out) const noexcept {    \
            out << m_msg << std::endl;                                  \
        }                                                               \
                                                                        \
      private:                                                          \
        const std::string m_msg;                                        \
    }

/**
 * Define an exception with one parameter.
 *
 * @ingroup Exceptions
 */
#define DefException1(Exception1, type1, outsequence)                \
    class Exception1 : public ::linalgwrap::ExceptionBase {          \
      public:                                                        \
        Exception1(const type1 a1) : arg1(a1) {}                     \
        virtual ~Exception1() noexcept {}                            \
        virtual void print_extra(std::ostream& out) const noexcept { \
            out outsequence << std::endl;                            \
        }                                                            \
                                                                     \
      private:                                                       \
        const type1 arg1;                                            \
    }

/**
 * Define an exception class derived from ExceptionBase with two additional
 * parameters.
 */
#define DefException2(Exception2, type1, type2, outsequence)               \
    class Exception2 : public ::linalgwrap::ExceptionBase {                \
      public:                                                              \
        Exception2(const type1 a1, const type2 a2) : arg1(a1), arg2(a2) {} \
        virtual ~Exception2() noexcept {}                                  \
        virtual void print_extra(std::ostream& out) const noexcept {       \
            out outsequence << std::endl;                                  \
        }                                                                  \
                                                                           \
      private:                                                             \
        const type1 arg1;                                                  \
        const type2 arg2;                                                  \
    }

/**
 * Declare an exception class derived from ExceptionBase with three additional
 * parameters.
 *
 * @ingroup Exceptions
 */
#define DefException3(Exception3, type1, type2, type3, outsequence)  \
    class Exception3 : public ::linalgwrap::ExceptionBase {          \
      public:                                                        \
        Exception3(const type1 a1, const type2 a2, const type3 a3)   \
              : arg1(a1), arg2(a2), arg3(a3) {}                      \
        virtual ~Exception3() noexcept {}                            \
        virtual void print_extra(std::ostream& out) const noexcept { \
            out outsequence << std::endl;                            \
        }                                                            \
                                                                     \
      private:                                                       \
        const type1 arg1;                                            \
        const type2 arg2;                                            \
        const type3 arg3;                                            \
    }

/**
 * Declare an exception class derived from ExceptionBase with four additional
 * parameters.
 *
 * @ingroup Exceptions
 */
#define DefException4(Exception4, type1, type2, type3, type4, outsequence) \
    class Exception4 : public ::linalgwrap::ExceptionBase {                \
      public:                                                              \
        Exception4(const type1 a1, const type2 a2, const type3 a3,         \
                   const type4 a4)                                         \
              : arg1(a1), arg2(a2), arg3(a3), arg4(a4) {}                  \
        virtual ~Exception4() noexcept {}                                  \
        virtual void print_extra(std::ostream& out) const noexcept {       \
            out outsequence << std::endl;                                  \
        }                                                                  \
                                                                           \
      private:                                                             \
        const type1 arg1;                                                  \
        const type2 arg2;                                                  \
        const type3 arg3;                                                  \
        const type4 arg4;                                                  \
    }

/**
 * Declare an exception class derived from ExceptionBase with five additional
 * parameters.
 *
 * @ingroup Exceptions
 */
#define DefException5(Exception5, type1, type2, type3, type4, type5, \
                      outsequence)                                   \
    class Exception5 : public ::linalgwrap::ExceptionBase {          \
      public:                                                        \
        Exception5(const type1 a1, const type2 a2, const type3 a3,   \
                   const type4 a4, const type5 a5)                   \
              : arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5) {}  \
        virtual ~Exception5() noexcept {}                            \
        virtual void print_extra(std::ostream& out) const noexcept { \
            out outsequence << std::endl;                            \
        }                                                            \
                                                                     \
      private:                                                       \
        const type1 arg1;                                            \
        const type2 arg2;                                            \
        const type3 arg3;                                            \
        const type4 arg4;                                            \
        const type5 arg5;                                            \
    }

}  // namespace linalgwrap

#endif  // LINALGWRAP_EXCEPTION_SYSTEM_H_
