#ifndef LINALG_SMALL_MATRIX_HPP_
#define LINALG_SMALL_MATRIX_HPP_

#include "Exceptions.hh"
#include "StoredMatrix_i.hh"
#include <armadillo>
#include <memory>
#include <type_traits>
#include "Constants.hh"
#include "Range.hh"

namespace linalgwrap {

// Forward-declare the interface class
template <typename Scalar>
class Matrix_i;

template <typename Scalar>
class StoredMatrix_i;

template <typename Scalar>
class SmallMatrix : public StoredMatrix_i<Scalar> {
    static_assert(std::is_same<double, Scalar>::value,
                  "SmallMatrix<Scalar> is currently only available for Scalar "
                  "== double.");

  public:
    typedef StoredMatrix_i<Scalar> base_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::size_type size_type;

    // Swapping:
    friend void swap(SmallMatrix& first, SmallMatrix& second) {
        using std::swap;
        first.m_arma.swap(second.m_arma);
        swap(static_cast<base_type&>(first), static_cast<base_type&>(second));
    }

    // TODO add constructor from tuple of initialiser list!

    //
    // Constructors, destructors and assignment
    //
    /** Construct a matrix of fixed size and optionally set the entries to zero
     */
    SmallMatrix(size_type n_rows, size_type n_cols, bool fill_zero = true)
          : m_arma(n_rows, n_cols, arma::fill::none) {
        if (fill_zero) {
            // set all elements to zero
            m_arma.zeros();
        }
    }

    /** Construct a small matrix and copy all entries from ``mat`` which are
     *  not below the tolerance threshold.
     */
    SmallMatrix(const SmallMatrix& mat, scalar_type tolerance)
          : SmallMatrix(mat.n_rows(), mat.n_cols(), false) {
        // TODO use iterator
        for (size_type i = 0; i < mat.n_rows(); ++i) {
            for (size_type j = 0; j < mat.n_cols(); ++j) {
                if (std::fabs(mat(i, j)) < tolerance) {
                    (*this)(i, j) = Constants<scalar_type>::zero();
                } else {

                    (*this)(i, j) = mat(i, j);
                }
            }
        }
    }

    explicit SmallMatrix(arma::mat inner) : m_arma(inner) {}

    //
    // Scalar operators
    //
    /** Scale matrix by a scalar value */
    SmallMatrix& operator*=(scalar_type s) {
        assert_finite(s);
        m_arma *= s;
        return *this;
    }

    /** Divide all matrix entries by a scalar value */
    SmallMatrix& operator/=(scalar_type s) {
        assert_dbg(s == 0, ExcDevideByZero());
        assert_finite(s);
        m_arma /= s;
        return *this;
    }

    /** Multiply two small matrices */
    SmallMatrix operator*(const SmallMatrix& other) const {
        assert_size(n_cols(), other.n_rows());
        arma::mat res = m_arma * other.m_arma;
        return SmallMatrix(res);
    }

    /* Add a small matrix to this one */
    SmallMatrix& operator+=(const SmallMatrix& other) {
        assert_size(n_cols(), other.n_cols());
        assert_size(n_rows(), other.n_rows());
        m_arma += other.m_arma;
        return *this;
    }

    /* Subtract a small matrix from this one */
    SmallMatrix& operator-=(const SmallMatrix& other) {
        assert_size(n_cols(), other.n_cols());
        assert_size(n_rows(), other.n_rows());
        m_arma -= other.m_arma;
        return *this;
    }

    //
    // Relational operatiors
    //
    bool operator==(const SmallMatrix& other) const {
        // does not work for some crazy arma reason
        // return (m_arma == other.m_arma);

        for (size_type i = 0; i < n_rows() * n_cols(); ++i) {
            if (m_arma[i] != other.m_arma[i]) return false;
        }
        return true;
    }

    bool operator!=(const SmallMatrix& other) const { !operator==(other); }

    //
    // matrix_i interface
    //
    /** \brief Number of rows of the matrix i
     */
    size_type n_rows() const override { return m_arma.n_rows; }

    /** \brief Number of columns of the matrix
     */
    size_type n_cols() const override { return m_arma.n_cols; }

    scalar_type operator()(size_type row, size_type col) const override {
        assert_upper_bound(row, n_rows());
        assert_upper_bound(col, n_cols());
        return m_arma.at(row, col);
    }

    // Note: operator[] is taken as the default implementation
    // since arma matrices are column-major, but our interface
    // is row-major

    //
    // StoredMatrix_i interface
    //
    /** Set all elements to zero */
    void set_zero() override { m_arma.zeros(); }

    scalar_type& operator()(size_type row, size_type col) override {
        assert_upper_bound(row, n_rows());
        assert_upper_bound(col, n_cols());
        return m_arma.at(row, col);
    }

    /** \brief Return a copy of a block of values out of the matrix and
     *         return it as a SmallMatrix of the appropriate size
     *
     * For more details of the interface see the function of the same
     * name in ``LazyMatrixExpression``.
     *
     * \param row_range   The Range object representing the range of rows
     *                    to extract. Note that it is a half-open interval
     *                    i.e. the LHS is inclusive, but the RHS not.
     * \param col_range   The Range object representing the range of
     *                    columns to extract.
     */
    virtual SmallMatrix<scalar_type> extract_block(
          Range<size_type> row_range, Range<size_type> col_range) const {
        // At least one range is empty -> no work to be done:
        if (row_range.is_empty() || col_range.is_empty()) {
            return SmallMatrix<scalar_type>{row_range.length(),
                                            col_range.length()};
        }

        // Assertive checks:
        assert_lower_bound(row_range.last(), this->n_rows() + 1);
        assert_lower_bound(col_range.last(), this->n_cols() + 1);

        // Translate ranges to armadillo spans (which are closed intervals)
        arma::span rows(row_range.first(), row_range.last() - 1);
        arma::span cols(col_range.first(), col_range.last() - 1);

        // Create a copy of the elements to extract
        arma::mat m = m_arma.submat(rows, cols);

        // Move into a now SmallMatrix:
        return SmallMatrix<scalar_type>{std::move(m)};
    }

    /** \brief Add a copy of a block of values of the matrix to
     *         the SmallMatrix provided by reference.
     *
     * For more details of the interface see the function of the same
     * name in ``LazyMatrixExpression``.
     *
     *  \param in   Matrix to add the values to. It is assumed that it
     *              already has the correct sparsity structure to take
     *              all the values. Its size defines the size of the
     *              block
     *  \param start_row  The row index of the first element to extract
     *  \param start_col  The column index of the first element to extract
     *  \param c_this     The coefficient to multiply this matrix with
     *                    before extracting.
     */
    void add_block_to(SmallMatrix<scalar_type>& in, size_type start_row,
                      size_type start_col,
                      scalar_type c_this = Constants<scalar_type>::one) const {
        // check that we do not overshoot the row index
        assert_upper_bound(start_row + in.n_rows(), this->n_rows() + 1);

        // check that we do not overshoot the column index
        assert_upper_bound(start_col + in.n_cols(), this->n_cols() + 1);

        // Do the operation:
        in.m_arma += c_this * m_arma(start_row, start_col, size(in.m_arma));
    }

    // Note: operator[] is taken as the default implementation
    // since arma matrices are column-major, but our interface
    // is row-major

  private:
    arma::mat m_arma;
};

//
// Multiply by Scalar
//
template <typename Scalar>
SmallMatrix<Scalar> operator*(Scalar s, SmallMatrix<Scalar> m) {
    m *= s;
    return m;
}

template <typename Scalar>
SmallMatrix<Scalar> operator*(SmallMatrix<Scalar> m, Scalar s) {
    return s * m;
}

template <typename Scalar>
SmallMatrix<Scalar> operator/(SmallMatrix<Scalar> m, Scalar s) {
    m /= s;
    return m;
}

//
// Add and subtract small matrices
//
template <typename Scalar>
SmallMatrix<Scalar> operator-(SmallMatrix<Scalar> lhs,
                              const SmallMatrix<Scalar>& rhs) {
    lhs -= rhs;
    return lhs;
}

template <typename Scalar>
SmallMatrix<Scalar> operator+(SmallMatrix<Scalar> lhs,
                              const SmallMatrix<Scalar>& rhs) {
    lhs += rhs;
    return lhs;
}

}  // liblinalg
#endif
