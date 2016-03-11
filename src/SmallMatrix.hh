#ifndef LINALG_SMALL_MATRIX_HPP_
#define LINALG_SMALL_MATRIX_HPP_

#include "Exceptions.hh"
#include "StoredMatrix_i.hh"
#include <armadillo>
#include <memory>
#include <type_traits>
#include "Constants.hh"

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
        arma::mat res = m_arma * other.m_arma;
        return SmallMatrix(res);
    }

    /* Add a small matrix to this one */
    SmallMatrix& operator+=(const SmallMatrix& other) {
        m_arma += other.m_arma;
        return *this;
    }

    /* Subtract a small matrix from this one */
    SmallMatrix& operator-=(const SmallMatrix& other) {
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
    /**
     * See documentation of Matrix_i function of the same name.
     */
    void fill(size_type start_row, size_type start_col,
              SmallMatrix<scalar_type>& block, bool add = false,
              scalar_type c_this = Constants<scalar_type>::one) const override {
        // check that we do not overshoot the row index
        assert_upper_bound(start_row + block.n_rows() - 1, n_rows());

        // check that we do not overshoot the column index
        assert_upper_bound(start_col + block.n_cols() - 1, n_cols());

        // TODO probably there exists a better way
        for (size_type i = 0; i < block.n_rows(); ++i) {
            for (size_type j = 0; j < block.n_cols(); ++j) {
                if (add) {
                    block(i, j) +=
                          c_this * (*this)(i + start_row, j + start_col);
                } else {
                    block(i, j) =
                          c_this * (*this)(i + start_row, j + start_col);
                }
            }
        }
    }

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
