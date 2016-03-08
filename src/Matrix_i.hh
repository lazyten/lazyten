#ifndef LINALG_MATRIX_I_HPP_
#define LINALG_MATRIX_I_HPP_

#include <cstddef>
#include <utility>
#include "SmallMatrix.hh"
#include <string>
#include "Exceptions.hh"
#include "Constants.hh"
#include <iostream>
#include <iomanip>
#include "MatrixIterator.hh"
#include "Subscribable.hh"

namespace linalgwrap {

// TODO make use of iterators

template <typename Scalar>
class SmallMatrix;

template <typename Scalar>
class MatrixIterator;

/** \brief Abstract matrix interface class
 *
 *  \note All classes implementing this interface should also be derived off
 *   DefaultIterators
 * */
template <typename Scalar>
class Matrix_i : public Subscribable {
  public:
    typedef size_t size_type;
    typedef Scalar scalar_type;

    typedef MatrixIterator<Scalar> iterator;
    typedef MatrixIterator<Scalar> const_iterator;

    /** \brief Destructor */
    virtual ~Matrix_i() = default;

    //
    // Matrix properties
    //
    /** \brief Number of rows of the matrix */
    virtual size_type n_rows() const = 0;

    /** \brief Number of columns of the matrix */
    virtual size_type n_cols() const = 0;

    //
    // Element access
    //
    /** \brief Extract values of this matrix partially to \p block
     * Get the values of the matrix starting at row index \p start_row
     * and column index \p start_col. As many entries as there is space
     * in \p block are extracted. So if \p block is a 2x2 matrix and
     * \p start_row is 1 and \p start_col is 3, then the entries
     * (1,3), (1,4), (2,3) and (2,4) are extracted.
     *
     * Optionally the values can be implicitly scaled by a coefficient
     * \p c_this before extracting them to \p block and the flag \p add
     * controls whether the values are added to \p block or set.
     *
     * @param start_row The row index to start from
     * @param start_col The col index to start from
     * @param block     The block to extract the values to.
     * @param add       If true add the values, else set them.
     * @param c_this    Coefficient to multiply all values of this matrix
     *                  before adding/setting them to \p block.
     */
    virtual void fill(
          size_type start_row, size_type start_col,
          SmallMatrix<scalar_type>& block, bool add = false,
          scalar_type c_this = Constants<scalar_type>::one) const = 0;

    /** \brief return an element of the matrix
     *
     * It is advisible to overload this in order to get a more performant
     * implementation.
     */
    virtual scalar_type operator()(size_type row, size_type col) const;

    /** \brief return an element of the vectorised matrix object
     *
     * Access the element in row-major ordering (i.e. the matrix is
     * traversed row by row)
     */
    virtual scalar_type operator[](size_type i) const;

    //
    // Iterators
    //
    /** Return an iterator to the beginning */
    virtual iterator begin();

    /** Return a const_iterator to the beginning */
    virtual const_iterator begin() const;

    /** Return a const_iterator to the beginning */
    virtual const_iterator cbegin() const;

    /** Return an iterator to the end */
    virtual iterator end();

    /** Return a const_iterator to the end */
    virtual const_iterator end() const;

    /** Return a const_iterator to the end */
    virtual const_iterator cend() const;

    //
    // Operations --- or maybe out of scope
    //
    /** \brief Return the inverse of the matrix
    void inverse() { assert_dbg(false, ExcNotImplemented()); }
    */

    /** \brief Return the transpose this matrix
    void transpose() { assert_dbg(false, ExcNotImplemented()); }
    */
};

/** \brief Simple output operator, that plainly shows all entries of
 *  the Matrix one by one.
 *
 *  Rows are seperated by a newline and entries by spaces.
 *  The last row is not terminated by a newline character.
 *  */
template <typename Scalar>
std::ostream& operator<<(std::ostream& o, const Matrix_i<Scalar>& m);

//
// ---------------------------------------------------------------
//

//
// Matrix_i
//
template <typename Scalar>
typename Matrix_i<Scalar>::scalar_type Matrix_i<Scalar>::operator()(
      size_type row, size_type col) const {
    // Check that we do not overshoot.
    assert_upper_bound(row, n_rows());
    assert_upper_bound(col, n_cols());

    SmallMatrix<scalar_type> block(1, 1, false);
    fill(row, col, block);
    return block(0, 0);
}

template <typename Scalar>
typename Matrix_i<Scalar>::scalar_type Matrix_i<Scalar>::operator[](
      size_type i) const {
    // Check that we do not overshoot.
    assert_upper_bound(i, n_cols() * n_rows());

    const size_type i_row = i / n_cols();
    const size_type i_col = i % n_cols();
    return (*this)(i_row, i_col);
}

template <typename Scalar>
typename Matrix_i<Scalar>::iterator Matrix_i<Scalar>::begin() {
    return iterator(*this, {0, 0});
}

template <typename Scalar>
typename Matrix_i<Scalar>::const_iterator Matrix_i<Scalar>::begin() const {
    return cbegin();
}

template <typename Scalar>
typename Matrix_i<Scalar>::const_iterator Matrix_i<Scalar>::cbegin() const {
    return const_iterator(*this, {0, 0});
}

template <typename Scalar>
typename Matrix_i<Scalar>::iterator Matrix_i<Scalar>::end() {
    return iterator(*this, iterator::invalid_pos);
}

template <typename Scalar>
typename Matrix_i<Scalar>::const_iterator Matrix_i<Scalar>::end() const {
    return cend();
}

template <typename Scalar>
typename Matrix_i<Scalar>::const_iterator Matrix_i<Scalar>::cend() const {
    return const_iterator(*this, const_iterator::invalid_pos);
}

//
// Out of scope
//
template <typename Scalar>
std::ostream& operator<<(std::ostream& o, const Matrix_i<Scalar>& m) {
    typedef typename Matrix_i<Scalar>::size_type size_type;

    for (size_type i = 0; i < m.n_rows(); ++i) {
        for (size_type j = 0; j < m.n_cols(); ++j) {
            o << std::setw(15) << m(i, j);
        }

        if (i == m.n_rows() - 1) break;
        o << std::endl;
    }

    return o;
}

}  // namespace linalg
#endif  // LINALG_MATRIX_I_HPP_
