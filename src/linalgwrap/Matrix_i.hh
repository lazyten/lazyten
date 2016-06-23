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

#ifndef LINALG_MATRIX_I_HPP_
#define LINALG_MATRIX_I_HPP_

#include "linalgwrap/Constants.hh"
#include "linalgwrap/DefaultMatrixIterator.hh"
#include "linalgwrap/Exceptions.hh"
#include "linalgwrap/SmallMatrix.hh"
#include "linalgwrap/Subscribable.hh"
#include "linalgwrap/io/MatrixPrinter.hh"
#include "linalgwrap/type_utils.hh"
#include <complex>
#include <cstddef>
#include <iomanip>
#include <ios>
#include <iostream>
#include <numeric>
#include <string>
#include <utility>

namespace linalgwrap {

template <typename IteratorCore>
class MatrixIterator;

template <typename Matrix, bool Constness>
class MatrixIteratorDefaultCore;

/** \brief Abstract matrix interface class
 *
 * This interface and hence all classes derived from it are subscribable using
 * the SubscriptionPointer class. This should be used very little and only when
 * other means (e.g. using shared pointers) is not possible.
 * */
template <typename Scalar>
class Matrix_i : public Subscribable {
  public:
    typedef size_t size_type;
    typedef Scalar scalar_type;

    //! The iterator type (a const iterator)
    typedef DefaultMatrixConstIterator<Matrix_i<Scalar>> iterator;

    //! The const iterator type
    typedef DefaultMatrixConstIterator<Matrix_i<Scalar>> const_iterator;

    /** \name Matrix constructor, destructor and assignment */
    ///@{
    virtual ~Matrix_i() = default;
    Matrix_i() = default;
    Matrix_i(const Matrix_i&) = default;
    Matrix_i(Matrix_i&&) = default;
    Matrix_i& operator=(const Matrix_i&) = default;
    Matrix_i& operator=(Matrix_i&&) = default;
    ///@}

    /** \name Matrix properties
     *        Access to properties common to all matrices
     */
    ///@{
    /** \brief Number of rows of the matrix */
    virtual size_type n_rows() const = 0;

    /** \brief Number of columns of the matrix */
    virtual size_type n_cols() const = 0;
    ///@}

    /** \name Data access
     *        Access to matrix data
     */
    ///@{
    /** \brief return an element of the matrix    */
    virtual scalar_type operator()(size_type row, size_type col) const = 0;

    /** \brief return an element of the vectorised matrix object
     *
     * Access the element in row-major ordering (i.e. the matrix is
     * traversed row by row)
     */
    virtual scalar_type operator[](size_type i) const;
    ///@}

    /** \name Iterators
     */
    ///@{
    /** Return an iterator to the beginning */
    iterator begin();

    /** Return a const_iterator to the beginning */
    const_iterator begin() const;

    /** Return a const_iterator to the beginning */
    const_iterator cbegin() const;

    /** Return an iterator to the end */
    iterator end();

    /** Return a const_iterator to the end */
    const_iterator end() const;

    /** Return a const_iterator to the end */
    const_iterator cend() const;
    ///@}

    /** \name Check for matrix properties
     */
    ///@{
    /** \brief Check whether the matrix is symmetric
     *
     * Loops over all elements and check wheather the difference
     * between m(i,j) and m(j,i) is less than the tolerance given
     * */
    bool is_symmetric(scalar_type tolerance =
                            Constants<scalar_type>::default_tolerance) const;

    // TODO ideas: is_real, is_hermetian, is_real_symmetric
    ///@}

    /** \name Standard operations
     */
    ///@{
    /** \brief Compute the trace of this matrix
     * Only works for quadratic matrices */
    scalar_type trace() const;

    /** Calculate the (signed) sum of all matrix entries. */
    scalar_type accumulate() const;

    /** Calculate the l1 norm (maximum of the sums over columns) */
    scalar_type norm_l1() const;

    /** Calculate the linf norm (maximum of the sums over rows) */
    scalar_type norm_linf() const;

    /** Calculate the Frobenius norm (sqrt of all matrix elements
     * squared
     *
     * \note This norm is not the matrix norm compatible to the l2 norm!
     */
    scalar_type norm_frobenius() const;

    /** Calculate the Frobenius norm squared */
    scalar_type norm_frobenius_squared() const;
    ///@}
};

//@{
/** \brief struct representing a type (std::true_type, std::false_type) which
 *  indicates whether T is a stored matrix
 *
 * The definition is done using SFINAE, such that even for types not having a
 * typedef scalar_type this expression is valid.
 *  */
template <typename T, typename = void>
struct IsMatrix : public std::false_type {};

template <typename T>
struct IsMatrix<T, void_t<typename T::scalar_type>>
      : public std::is_base_of<Matrix_i<typename T::scalar_type>, T> {};
//@}

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
typename Matrix_i<Scalar>::scalar_type Matrix_i<Scalar>::operator[](
      size_type i) const {
    // Check that we do not overshoot.
    assert_range(0, i, n_cols() * n_rows());

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
    return iterator(*this);
}

template <typename Scalar>
typename Matrix_i<Scalar>::const_iterator Matrix_i<Scalar>::end() const {
    return cend();
}

template <typename Scalar>
typename Matrix_i<Scalar>::const_iterator Matrix_i<Scalar>::cend() const {
    return const_iterator(*this);
}

template <typename Scalar>
bool Matrix_i<Scalar>::is_symmetric(scalar_type tolerance) const {
    // Check that the matrix is quadratic:
    if (n_rows() != n_cols()) return false;

    // Check if lower and upper triangle agree:
    for (auto i : range(n_rows())) {
        for (auto j : range(i + 1, n_cols())) {
            if (std::abs((*this)(i, j) - (*this)(j, i)) > tolerance)
                return false;
        }
    }
    return true;
}

template <typename Scalar>
inline typename Matrix_i<Scalar>::scalar_type Matrix_i<Scalar>::accumulate()
      const {
    return std::accumulate(begin(), end(), Constants<scalar_type>::zero);
}

template <typename Scalar>
inline typename Matrix_i<Scalar>::scalar_type Matrix_i<Scalar>::norm_l1()
      const {
    // This way is real bad for the cache and hence really slow.
    // One should do this in blocks of row indices, which fit the cache size.

    // maximum of the colsums
    //
    scalar_type res(Constants<scalar_type>::zero);
    for (auto col : range(n_cols())) {
        // sum of absolute entries of this column
        scalar_type colsum = Constants<scalar_type>::zero;
        for (auto row : range(n_rows())) {
            colsum += std::abs((*this)(row, col));
        }
        res = std::max(res, colsum);
    }
    return res;
}

template <typename Scalar>
inline typename Matrix_i<Scalar>::scalar_type Matrix_i<Scalar>::norm_linf()
      const {
    // maximum of the rowsums
    //
    scalar_type res = Constants<scalar_type>::zero;
    for (auto row : range(n_rows())) {
        // sum of absolute entries of this row
        scalar_type rowsum = Constants<scalar_type>::zero;
        for (auto col : range(n_cols())) {
            rowsum += std::abs((*this)(row, col));
        }
        res = std::max(res, rowsum);
    }
    return res;
}

template <typename Scalar>
inline typename Matrix_i<Scalar>::scalar_type Matrix_i<Scalar>::norm_frobenius()
      const {
    // sqrt of square of all matrix elements
    return std::sqrt(norm_frobenius_squared());
}

template <typename Scalar>
inline typename Matrix_i<Scalar>::scalar_type
Matrix_i<Scalar>::norm_frobenius_squared() const {
    // sum of squares of all matrix elements
    scalar_type sum = Constants<scalar_type>::zero;
    for (auto it = begin(); it != end(); ++it) {
        sum += (*it) * (*it);
    }
    return sum;
}

template <typename Scalar>
inline typename Matrix_i<Scalar>::scalar_type Matrix_i<Scalar>::trace() const {
    assert_dbg(n_rows() == n_cols(), ExcMatrixNotSquare());

    scalar_type trace{Constants<scalar_type>::zero};
    for (size_type i = 0; i < n_rows(); ++i) {
        trace += (*this)(i, i);
    }
    return trace;
}

//
// Out of scope
//
template <typename Scalar>
std::ostream& operator<<(std::ostream& o, const Matrix_i<Scalar>& m) {
    io::MatrixPrinter().print(m, o);
    return o;
}

}  // namespace linalg
#endif  // LINALG_MATRIX_I_HPP_
