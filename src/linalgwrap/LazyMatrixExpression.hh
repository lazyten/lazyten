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

#ifndef LINALG_LAZY_MATRIX_EXPRESSION_HPP_
#define LINALG_LAZY_MATRIX_EXPRESSION_HPP_

#include "linalgwrap/LazyMatrixProduct.hh"
#include "linalgwrap/LazyMatrixSum.hh"
#include "linalgwrap/Matrix_i.hh"
#include "linalgwrap/Range.hh"
#include "linalgwrap/StoredMatrix_i.hh"
#include <krims/ParameterMap.hh>
#include <memory>
#include <ostream>

namespace linalgwrap {

// Forward declarations
template <typename StoredMatrix>
class LazyMatrixSum;

template <typename StoredMatrix>
class LazyMatrixProduct;

/** \brief Generic LazyMatrixExpression class
 *
 * Abstract base class for all lazy matrix expressions.
 * Implements a slightly modified paradigm of Expression templates.
 *
 * \tparam StoredMatrix The type to use for stored matrices.
 * */
template <typename StoredMatrix>
class LazyMatrixExpression
      : public Matrix_i<typename StoredMatrix::scalar_type> {
    static_assert(
          std::is_base_of<StoredMatrix_i<typename StoredMatrix::scalar_type>,
                          StoredMatrix>::value,
          "StoredMatrix is not a child of StoredMatrix_i of the same scalar "
          "type");

  public:
    typedef StoredMatrix stored_matrix_type;
    typedef Matrix_i<typename StoredMatrix::scalar_type> base_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::size_type size_type;

    /** The pointer type to use for pointers to LazyMatrixExpressions */
    typedef std::shared_ptr<LazyMatrixExpression<StoredMatrix>>
          lazy_matrix_expression_ptr_type;

    /** \name Data access
     */
    ///@{
    /** \brief Convert the expression to a stored matrix
     *
     * Achieved by calling extract_block on the whole matrix.
     */
    virtual explicit operator stored_matrix_type() const {
        return extract_block(range(this->n_rows()), range(this->n_cols()));
    }

    /** \brief Extract a block of values out of the matrix and
     *         return it as a stored matrix of the appropriate size
     *
     * The values to extract are given by the index ranges. The ranges
     * are inclusive on the lhs and exclusive on the rhs. Hence
     * choosing a row range of [2,4) and a column range of [0,2) will
     * extract the four values
     * 			(2,0), (2,1),
     * 			(3,0), (3,1)
     * into a 4x4 matrix.
     *
     * The present implementation does not preform very well given that
     * it actually uses the element-wise access operator() (row,col) on
     * each element. Deriving classes should provide better implementations
     * which apply operations block-wise. E.g. instead of adding
     * individual elements of a matrix sum, we first extract the two
     * full matrices directly and then add them up. This has the advantage
     * that the latter type of operation is done by the library which knows
     * all the internal structure of the stored matrix and (hopefully)
     * is very well optimised for such operations.
     *
     * \param row_range   The Range object representing the range of rows
     *                    to extract. Note that it is a half-open interval
     *                    i.e. the LHS is inclusive, but the RHS not.
     *                    The Range may not be empty.
     * \param col_range   The Range object representing the range of
     *                    columns to extract.
     *                    The Range may not be empty.
     */
    virtual stored_matrix_type extract_block(Range<size_type> row_range,
                                             Range<size_type> col_range) const {
        // Assertive checks:
        assert_greater(0u, row_range.length());
        assert_greater(0u, col_range.length());
        assert_greater_equal(row_range.last(), this->n_rows());
        assert_greater_equal(col_range.last(), this->n_cols());

        // Allocate storage for return:
        // TODO take care of sparsity here:
        stored_matrix_type m(row_range.size(), col_range.size(), false);

        // copy data (this is for non-sparse matrices):
        for (auto it = std::begin(m); it != std::end(m); ++it) {
            m(it.row(), it.col()) = (*this)(it.row() + row_range.first(),
                                            it.col() + col_range.first());
        }

        return m;
    }

    /** \brief Add a block of values of the matrix to the stored matrix
     *         provided by reference.
     *
     *  Extract a block of values from this matrix, where the block
     *  size is given by the size of the Stored Matrix ``in`` and the
     *  element at which we start the extraction is given by
     *  ``start_row`` and ``start_col``. In other words if ``in`` is a
     *  2x2 matrix and ``start_row == 2`` and ``start_col==1`` we extract
     *  the four elements (2,1),(2,2), (3,1) and (3,2).
     *
     *  The of a call to this function should be equivalent to
     *  ```
     *  in += c_this*extract_block({start_row, in.n_rows()},
     *                             {start_col, in.n_cols()});
     *  ```
     *  In very many cases linear algebra libraries provide quicker
     *  routes for doing this scaled addition, so this function should
     *  be used for adding blocks of data to present data instead of
     *  the code mentioned above.
     *  Note, that the LazyMatrixSum class internally makes heavy use
     *  of this function.
     *
     *  So if the derived data structure may wish to make use of such
     *  better scaled addition functions of implementing libraries
     *  this function should be overloaded.
     *
     *  \note The function assumes that the sparsity pattern in ``in``
     *  already has the correct shape or is dynamically extended such
     *  that the addition does not hit non-existing elements.
     *
     *  \param in   Matrix to add the values to. It is assumed that it
     *              already has the correct sparsity structure to take
     *              all the values. Its size defines the size of the
     *              block. May not be an empty matrix
     *  \param start_row  The row index of the first element to extract
     *  \param start_col  The column index of the first element to extract
     *  \param c_this     The coefficient to multiply this matrix with
     *                    before extracting.
     */
    virtual void add_block_to(
          stored_matrix_type& in, size_type start_row, size_type start_col,
          scalar_type c_this = Constants<scalar_type>::one) const {
        assert_greater(0u, in.n_rows());
        assert_greater(0u, in.n_cols());

        // check that we do not overshoot the indices
        assert_greater_equal(start_row + in.n_rows(), this->n_rows());
        assert_greater_equal(start_col + in.n_cols(), this->n_cols());

        // Extract the block
        const stored_matrix_type extracted =
              extract_block({start_row, start_row + in.n_rows()},
                            {start_col, start_col + in.n_cols()});

        assert_dbg(extracted.n_rows() == in.n_rows(),
                   krims::ExcInternalError());
        assert_dbg(extracted.n_cols() == in.n_cols(),
                   krims::ExcInternalError());

        // Add it to in via the equivalent function
        // in the stored matrix
        extracted.add_block_to(in, 0, 0, c_this);
    }
    ///@}

    /** \brief Update the internal data of all objects in this expression
     * given the ParameterMap
     * */
    virtual void update(const krims::ParameterMap& map) = 0;

    /** \brief Multiplication with a stored matrix */
    virtual stored_matrix_type operator*(
          const stored_matrix_type& in) const = 0;

    /** \brief Clone the expression
     *
     * Return a clone of the current object as a pointer of type
     * lazy_matrix_expression_ptr_type
     */
    virtual lazy_matrix_expression_ptr_type clone() const = 0;
};

//@{
/** \brief struct representing a type (std::true_type, std::false_type) which
 *         indicates whether Matrix is a lazy matrix
 *
 * The definition is done using SFINAE, such that even for types not having a
 * typedef stored_matrix_type this expression is valid.
 *  */
template <typename Matrix, typename = void>
struct IsLazyMatrix : public std::false_type {};

template <typename Matrix>
struct IsLazyMatrix<Matrix, void_t<typename Matrix::stored_matrix_type>>
      : public std::is_base_of<
              LazyMatrixExpression<typename Matrix::stored_matrix_type>,
              Matrix> {};
//@}

//
// CallUpadateIfAvail class
//
/** \brief Call the update function of a lazy matrix if it is
 *         available. Else give rise to an InvalidState exception
 */
template <typename Matrix, bool = std::is_const<Matrix>::value>
struct CallUpdateIfAvail {
    void operator()(Matrix&, const ParameterMap&) const {
        assert_dbg(true,
                   ExcInvalidState("Update not available for const matrix"));
    }
};

// Specialisation of the above class
template <typename Matrix>
struct CallUpdateIfAvail<Matrix, false> {
    void operator()(Matrix& matrix, const ParameterMap& map) const {
        matrix.update(map);
    }
};

//
// Multiplication
//
/** Multiply two lazy matrix expressions */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator*(
      const LazyMatrixExpression<StoredMatrix>& lhs,
      const LazyMatrixExpression<StoredMatrix>& rhs) {

    // Construct product with one factor:
    LazyMatrixProduct<StoredMatrix> prod{lhs};

    // add the other factor and return the result.
    prod.push_factor(rhs);
    return prod;
}

/** Multiply stored times lazy */
template <typename StoredMatrix>
StoredMatrix operator*(const StoredMatrix& lhs,
                       const LazyMatrixExpression<StoredMatrix>& rhs) {
    assert_dbg(
          false,
          krims::ExcDisabled(
                "The operation \"StoredMatrix * LazyMatrixExpression\" is "
                "disabled because it is usually pretty expensive. "
                "Use \"StoredMatrix * LazyMatrixExpression * StoredMatrix\" "
                "instead if possible. Otherwise you can enforce"
                "\"StoredMatrix * LazyMatrixExpression\" by explicitly "
                "converting the LazyMatrixExpression into a StoredMatrix, "
                "i.e. \"lhs * static_cast<StoredMatrix>(rhs)\""));
    return lhs * static_cast<StoredMatrix>(rhs);
}

/** Scale a lazy matrix expression */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator*(
      const LazyMatrixExpression<StoredMatrix>& m,
      typename StoredMatrix::scalar_type s) {

    // Construct product with one factor:
    LazyMatrixProduct<StoredMatrix> prod{m};

    // and scale it:
    prod *= s;
    return prod;
}

/** Scale a lazy matrix expression */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator*(
      typename StoredMatrix::scalar_type s,
      const LazyMatrixExpression<StoredMatrix>& m) {
    return m * s;
}

//
// Division by scalar
//
/** \brief Devide a lazy matrix by a scalar  */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator/(
      const LazyMatrixExpression<StoredMatrix>& m,
      typename StoredMatrix::scalar_type s) {

    typedef typename StoredMatrix::scalar_type scalar_type;
    scalar_type inverse = Constants<scalar_type>::one / s;
    return m * inverse;
}

//
// Addition
//
/** Add two lazy matrix expressions */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator+(
      const LazyMatrixExpression<StoredMatrix>& lhs,
      const LazyMatrixExpression<StoredMatrix>& rhs) {

    // Construct sum with one term:
    LazyMatrixSum<StoredMatrix> sum{lhs};

    // add the other factor and return the result.
    sum.push_term(rhs);
    return sum;
}

/** Add a lazy matrix and a stored matrix */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator+(
      const LazyMatrixExpression<StoredMatrix>& lhs, const StoredMatrix& rhs) {

    // Construct sum with one term:
    LazyMatrixSum<StoredMatrix> sum{lhs};

    // add the other factor and return the result.
    sum.push_term(rhs);
    return sum;
}

/** Add a lazy matrix and a stored matrix */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator+(
      const StoredMatrix& lhs, const LazyMatrixExpression<StoredMatrix>& rhs) {
    return rhs + lhs;
}

//
// Subtraction
//
/** Subtract two lazy matrix expressions */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator-(
      const LazyMatrixExpression<StoredMatrix>& lhs,
      const LazyMatrixExpression<StoredMatrix>& rhs) {

    // Construct sum with one term:
    LazyMatrixSum<StoredMatrix> sum{lhs};

    // typedef the scalar type
    typedef typename StoredMatrix::scalar_type scalar_type;

    // subtract the other factor and return the result.
    sum.push_term(rhs, -Constants<scalar_type>::one);
    return sum;
}

/** Subtract a stored matrix from a lazy matrix */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator-(
      const LazyMatrixExpression<StoredMatrix>& lhs, const StoredMatrix& rhs) {

    // Construct sum with one term:
    LazyMatrixSum<StoredMatrix> sum{lhs};

    // typedef the scalar type
    typedef typename StoredMatrix::scalar_type scalar_type;

    // add the other factor and return the result.
    sum.push_term(rhs, -Constants<scalar_type>::one);
    return sum;
}

/** Subtract a lazy matrix from a stored matrix */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator-(
      const StoredMatrix& lhs, const LazyMatrixExpression<StoredMatrix>& rhs) {
    // Construct sum with one term:
    LazyMatrixSum<StoredMatrix> sum{lhs};

    // typedef the scalar type
    typedef typename StoredMatrix::scalar_type scalar_type;

    // add the other factor and return the result.
    sum.push_term(rhs, -Constants<scalar_type>::one);
    return sum;
}

//
// Unary operator-
//
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator-(
      const LazyMatrixExpression<StoredMatrix>& mat) {
    typedef typename StoredMatrix::scalar_type scalar_type;
    return -Constants<scalar_type>::one * mat;
}

}  // namespace linalg
#endif  // LINALG_LAZY_MATRIX_EXPRESSION_HPP_
