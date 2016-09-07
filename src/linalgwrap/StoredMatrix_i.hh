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

#ifndef LINALG_STORED_MATRIX_I_HPP_
#define LINALG_STORED_MATRIX_I_HPP_

#include "linalgwrap/Constants.hh"
#include "linalgwrap/DefaultMatrixIterator.hh"
#include "linalgwrap/Matrix_i.hh"
#include "linalgwrap/type_utils.hh"

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
 *
 * All implementing classes should further provide the function
 * ```
 * stored_matrix_type extract_block(Range<size_type> row_range,
 *                                  Range<size_type> col_range) const;
 * ```
 * which should copy a block of the matrix and return it
 * (similar to the ``extract_block`` in the ``LazyMatrixExpression`` class,
 * as well as
 * ```
 * void add_block_to(stored_matrix_type& in, size_type start_row,
 *                   size_type start_col,
 *                   scalar_type c_this = Constants<scalar_type>::one) const;
 * ```
 * which --- again similar to ``add_block_to`` of ``LazyMatrixExpression`` ---
 * should add a copy of a block of the matrix to the matrix provided on in.
 *
 * Note that the operator() functions in derived classes are expected to return
 * zero even if an element is known to be zero by some sparsity pattern or
 * similar. Modification of a non-existing element should fail, however.
 */
template <typename Scalar>
class StoredMatrix_i : public Matrix_i<Scalar> {
  public:
    typedef Matrix_i<Scalar> base_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::size_type size_type;

    //! The iterator type
    typedef DefaultMatrixIterator<StoredMatrix_i<Scalar>> iterator;

    //! The const_iterator type
    typedef typename base_type::const_iterator const_iterator;

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
        assert_range(0u, i, this->n_cols() * this->n_rows());

        const size_type i_row = i / this->n_cols();
        const size_type i_col = i % this->n_cols();
        return (*this)(i_row, i_col);
    }

    /** \brief Read-only access to vectorised matrix object
     *
     * Access the element in row-major ordering (i.e. the matrix is
     * traversed row by row)
     * */
    scalar_type operator[](size_type i) const override {
        return base_type::operator[](i);
    }

    /** \name Standard operations */
    ///@{
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

    ///@}

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
    //   function to get actual number of non-zero entries
    //   function to get estimated/implicitly known number of non-zero entries
};

//@{
/** \brief struct representing a type (std::true_type, std::false_type) which
 *indicates
 *  whether T is a stored matrix
 *
 * The definition is done using SFINAE, such that even for types not having a
 *typedef
 * scalar_type this expression is valid.
 *  */
template <typename Matrix, typename = void>
struct IsStoredMatrix : public std::false_type {};

template <typename Matrix>
struct IsStoredMatrix<Matrix, void_t<typename Matrix::scalar_type>>
      : public std::is_base_of<StoredMatrix_i<typename Matrix::scalar_type>,
                               Matrix> {};
//@}

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
