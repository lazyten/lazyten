#ifndef LINALG_STORED_MATRIX_I_HPP_
#define LINALG_STORED_MATRIX_I_HPP_

#include "Matrix_i.hh"
#include "Subscribable.hh"
#include "Constants.hh"

namespace linalgwrap {

// TODO define a set of optional functions which make performance better
//      e.g. - an in-memory transpose() function
//           - transpose-add ...
//           - transpose-multiply
//           - whatever else seems sensible

// Forward-declare the interface class
template <typename Scalar>
class Matrix_i;

/** \brief Interface class for a matrix which is actually stored in memory
 * in some way
 *
 * This interface and hence all classes derived from it are subscribable using
 * the SubscriptionPointer class. This should be used very little and only when
 * other means (e.g. using shared pointers) is not possible.
 *
 * We expect any implementing class to also provide the following constructors:
 * - Construct matrix of fixed size and optionally fill with zeros or leave
 *   memory unassigned:
 *   ```
 *   StoredMatrix_i(n_rows, n_cols, fill_zero);
 *   ```
 * - Construct a matrix of the same size as a SmallMatrix and copy all entries
 *   from the SmallMatrix over, optionally providing a tolerance below which
 *   the entries are considered to be zero (The latter is useful for CRS
 *   matrices)
 *   ```
 *   StoredMatrix_i(const SmallMatrix&)
 *   StoredMatrix_i(const SmallMatrix&, scalar_type tolerance)
 *   ```
 */
template <typename Scalar>
class StoredMatrix_i : public Matrix_i<Scalar>, public Subscribable {
  public:
    typedef Matrix_i<Scalar> base_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::size_type size_type;

    //! The iterator type
    typedef DefaultMatrixIterator<StoredMatrix_i<Scalar>> iterator;

    //! The const_iterator type
    typedef typename base_type::const_iterator const_iterator;

    // Swapping:
    friend void swap(StoredMatrix_i& first, StoredMatrix_i& second) {
        using std::swap;
        swap(first.m_name, second.m_name);
    }

    //
    // Assignment, construction, destruction
    //

    /** Default destructor */
    virtual ~StoredMatrix_i() = default;

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
     * Access the element in row-major ordering (i.e. the matrix is
     * traversed row by row)
     */
    virtual scalar_type& operator[](size_type i) {
        // Check that we do not overshoot.
        assert_upper_bound(i, this->n_cols() * this->n_rows());

        const size_type i_row = i / this->n_cols();
        const size_type i_col = i % this->n_cols();
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

    //
    // Iterators
    //
    /** Return an iterator to the beginning */
    iterator begin();

    /** Return a const iterator to the beginning */
    const_iterator begin() const;

    /** Return an iterator to the end */
    iterator end();

    /** Return a const iterator to the end */
    const_iterator end() const;

    // TODO
    //   function to get stl-compatible iterator
    //   function to get number of non-zero entries

  protected:
    //! some name identifying the matrix, or empty
    std::string m_name;
};

//
// -------------------------------------------------------------
//

template <typename Scalar>
typename StoredMatrix_i<Scalar>::iterator StoredMatrix_i<Scalar>::begin() {
    return iterator(*this, {0, 0});
}

template <typename Scalar>
typename StoredMatrix_i<Scalar>::const_iterator StoredMatrix_i<Scalar>::begin()
      const {
    return base_type::cbegin();
}

template <typename Scalar>
typename StoredMatrix_i<Scalar>::iterator StoredMatrix_i<Scalar>::end() {
    return iterator(*this);
}

template <typename Scalar>
typename StoredMatrix_i<Scalar>::const_iterator StoredMatrix_i<Scalar>::end()
      const {
    return base_type::cend();
}

}  // namespace liblinalg
#endif
