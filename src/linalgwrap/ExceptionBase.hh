#ifndef LINALGWRAP_EXCEPTION_BASE_H_
#define LINALGWRAP_EXCEPTION_BASE_H_

#include <exception>
#include <ostream>
#include <string>

namespace linalgwrap {

/** \brief Base class for all exceptions in linalgwrap
 *
 *  Greatly inspired by the deal.ii exception system.
 */
class ExceptionBase : public std::exception {
  public:
    /** Constructor */
    ExceptionBase();

    /** Default copy constructor */
    ExceptionBase(const ExceptionBase&) = default;

    /** Default destructor */
    virtual ~ExceptionBase() noexcept = default;

    /** Add enhancing exception data */
    void add_exc_data(const char* file, int line, const char* function,
                      const char* failed_condition, const char* exception_name);

    const char* what() const noexcept;
    const char* name() const;

    /** Print exception-specific extra information to the outstream */
    virtual void print_extra(std::ostream& out) const noexcept;

    /** Print a stacktrace to the outstream */
    void print_stacktrace(std::ostream& out) const;

  protected:
    //! build the what string of this exception
    virtual std::string generate_message() const noexcept;

    virtual void print_exc_data(std::ostream& out) const noexcept;

    //! The name of the exception
    const char* m_name;

    //! The file where exception occurred
    const char* m_file;

    //! The line where exception occurred
    int m_line;

    //! The name of the function
    const char* m_function;

    //! The failed condition as a string.
    const char* m_failed_condition;

    /**
     * The number of stacktrace frames that are stored in the previous variable.
     * Zero if the system does not support stack traces.
     */
    int m_n_stacktrace_frames;

#ifdef LINALGWRAP_HAVE_GLIBC_STACKTRACE
    /**
     *   array of pointers that contains the raw stack trace
     */
    void* m_raw_stacktrace[25];
#endif

  private:
    //! The what of the exception
    mutable std::string m_what_str;
};

}  // namespace linalgwrap
#endif  // LINALGWRAP_EXCEPTION_BASE_H_
