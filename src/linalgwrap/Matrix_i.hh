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

#pragma once
#include "Base/Interfaces/Indexable_i.hh"
#include "Constants.hh"
#include "DefaultMatrixIterator.hh"
#include "Exceptions.hh"
#include "MultiVector.hh"
#include "PtrVector.hh"
#include "SmallMatrix.hh"
#include "io/MatrixPrinter.hh"
#include <complex>
#include <cstddef>
#include <iomanip>
#include <ios>
#include <iostream>
#include <krims/TypeUtils.hh>
#include <numeric>
#include <string>
#include <utility>

namespace linalgwrap {

// TODO subview for Matrices

template <typename IteratorCore>
class MatrixIterator;

template <typename Matrix, bool Constness>
class MatrixIteratorDefaultCore;

/** \brief Abstract matrix interface class
 *
 * This interface and hence all classes derived from it are subscribable using
 * the krims::SubscriptionPointer class. This should be used very little and
 * only when other means (e.g. using shared pointers) is not possible.
 * */
template <typename Scalar>
class Matrix_i : public Indexable_i<Scalar> {
  public:
    typedef Indexable_i<Scalar> base_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::real_type real_type;

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

    /** \brief Return the number of elements of the matrix */
    virtual size_type n_elem() const override { return n_rows() * n_cols(); }
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
    virtual scalar_type operator[](size_type i) const override;
    ///@}

    /** \name Iterators
     */
    ///@{
    /** Return an iterator to the beginning */
    iterator begin() { return iterator(*this, {0, 0}); }

    /** Return a const_iterator to the beginning */
    const_iterator begin() const { return cbegin(); }

    /** Return a const_iterator to the beginning */
    const_iterator cbegin() const { return const_iterator(*this, {0, 0}); }

    /** Return an iterator to the end */
    iterator end() { return iterator(*this); }

    /** Return a const_iterator to the end */
    const_iterator end() const { return cend(); }

    /** Return a const_iterator to the end */
    const_iterator cend() const { return const_iterator(*this); }
    ///@}

    /** \name Check for matrix properties
     */
    ///@{
    /** \brief Check whether the matrix is symmetric
     *
     * Loops over all elements and check whether the difference
     * between m(i,j) and m(j,i) is less than the tolerance given
     * */
    bool is_symmetric(
          real_type tolerance = Constants<real_type>::default_tolerance) const;

    /** \brief Check whether the matrix is Hermitian
     *
     * Loops over all elements and check whether the difference
     * between conj(m(i,j)) and m(j,i) is less than the tolerance given
     * */
    bool is_hermitian(
          real_type tolerance = Constants<real_type>::default_tolerance) const;
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
struct IsMatrix<T, krims::VoidType<typename T::scalar_type>>
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

/** Compute the trace of the matrix
 *
 * \note only sensible for square matrices
 */
template <typename Scalar>
Scalar trace(const Matrix_i<Scalar>& m);

/** Accumulate all matrix values */
template <typename Scalar>
Scalar accumulate(const Matrix_i<Scalar>& m) {
    return std::accumulate(m.begin(), m.end(), Constants<Scalar>::zero);
}

/** Compute the l1 norm (maximum of the sums over columns) */
template <typename Scalar>
Scalar norm_l1(const Matrix_i<Scalar>& m);

/** Calculate the linf norm (maximum of the sums over rows) */
template <typename Scalar>
Scalar norm_linf(const Matrix_i<Scalar>& m);

/** Calculate the Frobenius norm (sqrt of all matrix elements
 * squared
 *
 * \note This norm is not the matrix norm compatible to the l2 norm!
 */
template <typename Scalar>
Scalar norm_frobenius(const Matrix_i<Scalar>& m) {
    // sqrt of square of all matrix elements
    return std::sqrt(norm_frobenius_squared(m));
}

/** Calculate the Frobenius norm squared
 *
 * \note This norm is not the matrix norm compatible to the l2 norm!
 */
template <typename Scalar>
Scalar norm_frobenius_squared(const Matrix_i<Scalar>& m);

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
    assert_range(0u, i, n_cols() * n_rows());

    const size_type i_row = i / n_cols();
    const size_type i_col = i % n_cols();
    return (*this)(i_row, i_col);
}

template <typename Scalar>
bool Matrix_i<Scalar>::is_symmetric(real_type tolerance) const {
    // Check that the matrix is quadratic:
    if (n_rows() != n_cols()) return false;

    // Check if lower and upper triangle agree:
    for (size_type i = 0; i < n_rows(); ++i) {
        for (size_type j = 0; j < n_cols(); ++j) {
            if (std::abs((*this)(i, j) - (*this)(j, i)) > tolerance)
                return false;
        }
    }
    return true;
}

template <typename Scalar>
bool Matrix_i<Scalar>::is_hermitian(real_type tolerance) const {
    // Check that the matrix is quadratic:
    if (n_rows() != n_cols()) return false;

    // Check if lower and upper triangle agree:
    for (size_type i = 0; i < n_rows(); ++i) {
        for (size_type j = 0; j < n_cols(); ++j) {
            if (std::abs(std::conj((*this)(i, j)) - (*this)(j, i)) > tolerance)
                return false;
        }
    }
    return true;
}

//
// Out of scope
//

template <typename Scalar>
Scalar norm_l1(const Matrix_i<Scalar>& m) {
    typedef typename Matrix_i<Scalar>::size_type size_type;
    // This way is real bad for the cache and hence really slow.
    // One should do this in blocks of row indices, which fit the cache size.

    // maximum of the colsums
    //
    Scalar res(Constants<Scalar>::zero);
    for (size_type col = 0; col < m.n_cols(); ++col) {
        // sum of absolute entries of this column
        Scalar colsum = Constants<Scalar>::zero;
        for (size_type row = 0; row < m.n_rows(); ++row) {
            colsum += std::abs(m(row, col));
        }
        res = std::max(res, colsum);
    }
    return res;
}

template <typename Scalar>
Scalar norm_linf(const Matrix_i<Scalar>& m) {
    typedef typename Matrix_i<Scalar>::size_type size_type;

    Scalar res = Constants<Scalar>::zero;
    for (size_type row = 0; row < m.n_rows(); ++row) {
        // sum of absolute entries of this row
        Scalar rowsum = Constants<Scalar>::zero;
        for (size_type col = 0; col < m.n_cols(); ++col) {
            rowsum += std::abs(m(row, col));
        }
        res = std::max(res, rowsum);
    }
    return res;
}

template <typename Scalar>
Scalar norm_frobenius_squared(const Matrix_i<Scalar>& m) {
    // sum of squares of all matrix elements
    Scalar sum = Constants<Scalar>::zero;
    for (auto it = m.begin(); it != m.end(); ++it) {
        sum += (*it) * (*it);
    }
    return sum;
}

template <typename Scalar>
Scalar trace(const Matrix_i<Scalar>& m) {
    typedef typename Matrix_i<Scalar>::size_type size_type;
    assert_dbg(m.n_rows() == m.n_cols(), ExcMatrixNotSquare());

    Scalar trace{Constants<Scalar>::zero};
    for (size_type i = 0; i < m.n_rows(); ++i) {
        trace += m(i, i);
    }
    return trace;
}

template <typename Scalar>
std::ostream& operator<<(std::ostream& o, const Matrix_i<Scalar>& m) {
    io::MatrixPrinter().print(m, o);
    return o;
}

}  // namespace linalgwrap
