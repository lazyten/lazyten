#ifndef LINALG_STORED_MATRIX_I_HPP_
#define LINALG_STORED_MATRIX_I_HPP_

#include "Matrix_i.hh"
#include "Subscribable.hh"
#include "Constants.hh"

namespace linalgwrap {

/** \brief Interface class for a matrix which is actually stored in memory
 * in some way
 *
 * This interface and hence all classes derived from it are subscribable using
 * the SubscriptionPointer class. This should be used very little and only when
 * other means (e.g. using shared pointers) is not possible.
 */
template <typename Scalar>
class StoredMatrix_i : public Matrix_i<Scalar>, public Subscribable {
  public:
    typedef Matrix_i<Scalar> base_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::size_type size_type;

    //
    // Assignment, construction, destruction
    //

    /** Default destructor */
    virtual ~StoredMatrix_i() = default;

    // We expect the following constructors to be implemented:
    //     - StoredMatrix_i(n_rows, n_cols, fill_zero)

    //
    // Stored matrices can have a name
    //
    /** Get the name of the matrix */
    std::string name() const { return m_name; }

    /** Set the name of the matrix */
    void name(const std::string& name) { m_name = name; }

    //
    // Modifying functions
    //
    /** Read-write access to elements */
    virtual scalar_type& operator()(size_type row, size_type col) = 0;

    /** \brief Read-write access to vectorised matrix object
     *
     * It is advisible to overload this in order to get a more performant
     * implementation.
     */
    virtual scalar_type& operator[](size_type i) {
        size_type i_row = i / this->n_rows();
        size_type i_col = i % this->n_rows();
        return (*this)(i_row, i_col);
    }

    /** Set all elements to zero
     *
     * Overload this to get a more performant implementation.
     *
     * TODO generalise with a lambda and an arbitrary scalar
     * */
    virtual void set_zero() {
        for (size_type i = 0; i < this->n_rows(); ++i) {
            for (size_type j = 0; j < this->n_cols(); ++j) {
                (*this)(i, j) = 0;
            }
        }
    }

  protected:
    //! some name identifying the matrix, or empty
    std::string m_name;
};

}  // namespace liblinalg
#endif
