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
    }

    //
    // Constructors, destructors and assignment
    //
    SmallMatrix(size_type n_rows, size_type n_cols, bool fill_zero = true)
          : m_arma(n_rows, n_cols, arma::fill::none) {
        if (fill_zero) {
            // set all elements to zero
            m_arma.zeros();
        }
    }

    explicit SmallMatrix(arma::mat inner) : m_arma(inner) {}

    //
    // Operators
    //
    // TODO
    SmallMatrix& operator*=(scalar_type s) {
        m_arma *= s;
        return *this;
    }

    friend SmallMatrix operator*(const SmallMatrix& lhs,
                                 const SmallMatrix& rhs) {
        arma::mat res = lhs.m_arma * rhs.m_arma;
        return SmallMatrix(res);
    }

    SmallMatrix& operator+=(const SmallMatrix& other) {
        m_arma += other.m_arma;
        return *this;
    }

    //
    // matrix_i interface
    //
    /**
     * See documentation of Matrix_i function of the same name.
     */
    virtual void fill(size_type start_row, size_type start_col,
                      SmallMatrix<scalar_type>& block, bool add = false,
                      scalar_type c_this = Constants<scalar_type>::one) const {
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
    virtual size_type n_rows() const { return m_arma.n_rows; }

    /** \brief Number of columns of the matrix
     */
    virtual size_type n_cols() const { return m_arma.n_cols; }

    virtual scalar_type operator()(size_type row, size_type col) const {
        assert_upper_bound(row, n_rows());
        assert_upper_bound(col, n_cols());
        return m_arma(row, col);
    }

    scalar_type operator[](size_type i) const {
        assert_upper_bound(i, n_rows() * n_cols());
        return m_arma[i];
    }

    //
    // StoredMatrix_i interface
    //
    /** Set all elements to zero */
    void set_zero() { m_arma.zeros(); }

    scalar_type& operator()(size_type row, size_type col) {
        assert_upper_bound(row, n_rows());
        assert_upper_bound(col, n_cols());
        return m_arma(row, col);
    }

    scalar_type& operator[](size_type i) {
        assert_upper_bound(i, n_rows() * n_cols());
        return m_arma[i];
    }

  private:
    arma::mat m_arma;
};

template <typename Scalar>
SmallMatrix<Scalar> operator*(Scalar s, SmallMatrix<Scalar> m) {
    m *= s;
    return m;
}

template <typename Scalar>
SmallMatrix<Scalar> operator*(SmallMatrix<Scalar> m, Scalar s) {
    return s * m;
}

}  // liblinalg
#endif
