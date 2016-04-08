#pragma once

#include "Constants.hh"
#include "Exceptions.hh"
#include "Range.hh"
#include "StoredMatrix_i.hh"
#include <initializer_list>
#include <memory>
#include <type_traits>

#ifdef LINALGWRAP_HAVE_ARMADILLO
#include <armadillo>
#endif

namespace linalgwrap {
#ifdef LINALGWRAP_HAVE_ARMADILLO

// Forward-declare the interface class
template <typename Scalar>
class Matrix_i;

template <typename Scalar>
class StoredMatrix_i;

/** A class for a dense stored matrix, currently implemented using armadillo.
 *
 * \note Armadillo is storing matrix data in column-major format
 *       (Fortran-style), but the default in our library is to
 *       use row-major storage (C-style) so in fact all data is
 *       stored in a transposed armadillo matrix such that memory
 *       access in a row-major fashion is contiguous.
 * */
template <typename Scalar>
class ArmadilloMatrix : public StoredMatrix_i<Scalar> {
    static_assert(
          std::is_same<double, Scalar>::value ||
                std::is_same<float, Scalar>::value ||
                std::is_same<std::complex<float>, Scalar>::value ||
                std::is_same<std::complex<double>, Scalar>::value ||
                std::is_same<short, Scalar>::value ||
                std::is_same<int, Scalar>::value ||
                std::is_same<long, Scalar>::value ||
                std::is_same<unsigned short, Scalar>::value ||
                std::is_same<unsigned int, Scalar>::value ||
                std::is_same<unsigned long, Scalar>::value,
          "ArmadilloMatrix<Scalar> is currently only available for Scalar "
          "being one of double, float, complex<double>, "
          "complex<float>,  short, int, long, unsigned short, unsigned "
          "int, unsigned long");

  public:
    typedef StoredMatrix_i<Scalar> base_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::size_type size_type;

    /** The type of the storage object used to store the data
     *  of the ArmadilloMatrix */
    typedef arma::Mat<Scalar> storage_type;

    // Swapping:
    friend void swap(ArmadilloMatrix& first, ArmadilloMatrix& second);

    /** \name Constructors
     */
    ///@{
    /** Construct a matrix of fixed size and optionally set the entries to
     * zero
     */
    ArmadilloMatrix(size_type n_rows, size_type n_cols, bool fill_zero = true);

    /** Construct a small matrix and copy all entries from ``mat`` which are
     *  not below the tolerance threshold.
     */
    ArmadilloMatrix(const ArmadilloMatrix& mat, scalar_type tolerance);

    /** \brief Construct from a nested initialiser list of scalars.
     *
     * The outermost layer gives the number of Columns, the innermost layer the
     * elements in each columns. An example would be
     * ```
     * ArmadilloMatrix mat{{1.0,2.0,0.5},{1.5,4.5,6.}})
     * ```
     * which produces a 2x3 matrix.
     */
    ArmadilloMatrix(std::initializer_list<std::initializer_list<scalar_type>>
                          list_of_lists);

    /** Construct from transposed armadillo matrix
     *
     * \note The armadillo matrix has to be the transpose of the matrix
     * which is to be constructed. This is because internally we use a
     * transpose armadillo matrix in order to store the data.
     *
     * This is due to the fact that armadillo matrices are column-major,
     * but we are row-major.
     */
    explicit ArmadilloMatrix(storage_type inner);
    ///@}

    /** \name Matrix operations */
    ///@{
    /** Scale matrix by a scalar value */
    ArmadilloMatrix& operator*=(scalar_type s) {
        assert_finite(s);
        m_arma *= s;
        return *this;
    }

    /** Divide all matrix entries by a scalar value */
    ArmadilloMatrix& operator/=(scalar_type s) {
        assert_dbg(s != 0, ExcDevideByZero());
        assert_finite(s);
        m_arma /= s;
        return *this;
    }

    /** Multiply two small matrices */
    ArmadilloMatrix operator*(const ArmadilloMatrix& other) const {
        assert_size(n_cols(), other.n_rows());

        // m_arma is the transpose of what we represent
        //   (see comment above class definition)
        // other.m_arme is the transpose of what other
        // represents
        // hence
        storage_type tres = other.m_arma * m_arma;
        // is the transpose of the result this * other
        // and this is what we need to envelope into the new ArmadilloMatrix.
        return ArmadilloMatrix(tres);
    }

    /* Add a small matrix to this one */
    ArmadilloMatrix& operator+=(const ArmadilloMatrix& other) {
        assert_size(n_cols(), other.n_cols());
        assert_size(n_rows(), other.n_rows());
        m_arma += other.m_arma;
        return *this;
    }

    /* Subtract a small matrix from this one */
    ArmadilloMatrix& operator-=(const ArmadilloMatrix& other) {
        assert_size(n_cols(), other.n_cols());
        assert_size(n_rows(), other.n_rows());
        m_arma -= other.m_arma;
        return *this;
    }
    ///@}

    //
    // Relational operatiors
    //
    bool operator==(const ArmadilloMatrix& other) const {
        // does not work for some crazy arma reason
        // return (m_arma == other.m_arma);

        for (size_type i = 0; i < n_rows() * n_cols(); ++i) {
            if (m_arma[i] != other.m_arma[i]) return false;
        }
        return true;
    }

    bool operator!=(const ArmadilloMatrix& other) const { !operator==(other); }

    //
    // matrix_i interface
    //
    /** \brief Number of rows of the matrix i
     */
    size_type n_rows() const override {
        // Note comment above the class definition why it
        // has to be this way round
        return m_arma.n_cols;
    }

    /** \brief Number of columns of the matrix
     */
    size_type n_cols() const override {
        // Note comment above the class definition why it
        // has to be this way round
        return m_arma.n_rows;
    }

    scalar_type operator()(size_type row, size_type col) const override {
        assert_greater(row, n_rows());
        assert_greater(col, n_cols());
        // Note comment above the class definition why it
        // has to be this way round
        return m_arma.at(col, row);
    }

    scalar_type operator[](size_type i) const override {
        assert_greater(i, n_cols() * n_rows());
        // Note that armadillo storage in column-major, but since
        // we store the transpose of what we represent internally
        // this give the correct interface (row-major access)
        return m_arma[i];
    }

    //
    // StoredMatrix_i interface
    //
    /** Set all elements to zero */
    void set_zero() override { m_arma.zeros(); }

    scalar_type accumulate() const override { return arma::accu(m_arma); }

    scalar_type& operator()(size_type row, size_type col) override {
        assert_greater(row, n_rows());
        assert_greater(col, n_cols());
        return m_arma.at(col, row);
    }

    /** \brief Return a copy of a block of values out of the matrix and
     *         return it as a ArmadilloMatrix of the appropriate size
     *
     * For more details of the interface see the function of the same
     * name in ``LazyMatrixExpression``.
     *
     * \param row_range   The Range object representing the range of rows
     *                    to extract. Note that it is a half-open interval
     *                    i.e. the LHS is inclusive, but the RHS not.
     *                    The Range may not be empty.
     * \param col_range   The Range object representing the range of
     *                    columns to extract.
     *                    The Range may not be empty.
     */
    virtual ArmadilloMatrix<scalar_type> extract_block(
          Range<size_type> row_range, Range<size_type> col_range) const {
        // Assertive checks:
        assert_greater(0, row_range.length());
        assert_greater(0, col_range.length());

        assert_greater_equal(row_range.last(), this->n_rows());
        assert_greater_equal(col_range.last(), this->n_cols());

        // Translate ranges to armadillo spans (which are closed intervals)
        arma::span rows(row_range.first(), row_range.last() - 1);
        arma::span cols(col_range.first(), col_range.last() - 1);

        // Create a copy of the elements to extract
        // Note comment above the class definition why it
        // has to be this way round
        storage_type m = m_arma(cols, rows);

        // Move into a now ArmadilloMatrix:
        return ArmadilloMatrix<scalar_type>{std::move(m)};
    }

    /** \brief Add a copy of a block of values of the matrix to
     *         the ArmadilloMatrix provided by reference.
     *
     * For more details of the interface see the function of the same
     * name in ``LazyMatrixExpression``.
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
    void add_block_to(ArmadilloMatrix<scalar_type>& in, size_type start_row,
                      size_type start_col,
                      scalar_type c_this = Constants<scalar_type>::one) const {
        assert_greater(0, in.n_rows());
        assert_greater(0, in.n_cols());

        // check that we do not overshoot the indices
        assert_greater_equal(start_row + in.n_rows(), this->n_rows());
        assert_greater_equal(start_col + in.n_cols(), this->n_cols());

        // Do the operation:
        // Note comment above the class definition why it
        // has to be this way round
        in.m_arma += c_this * m_arma(start_col, start_row, size(in.m_arma));
    }

    scalar_type& operator[](size_type i) override {
        assert_greater(i, n_cols() * n_rows());
        // Note that armadillo storage in column-major, but since
        // we store the transpose of what we represent internally
        // this give the correct interface (row-major access)
        return m_arma[i];
    }

    /** Read-only access to the inner storage
     *
     * \note Internally we store the data in transposed form, so
     * this will return an armadillo matrix representing the transpose
     * of the matrix this class represents.
     * */
    const storage_type& data() const { return m_arma; }

  private:
    storage_type m_arma;
};

//
// Multiply by Scalar
//
template <typename Scalar>
ArmadilloMatrix<Scalar> operator*(Scalar s, ArmadilloMatrix<Scalar> m) {
    m *= s;
    return m;
}

template <typename Scalar>
ArmadilloMatrix<Scalar> operator*(ArmadilloMatrix<Scalar> m, Scalar s) {
    return s * m;
}

template <typename Scalar>
ArmadilloMatrix<Scalar> operator/(ArmadilloMatrix<Scalar> m, Scalar s) {
    m /= s;
    return m;
}

//
// Add and subtract small matrices
//
template <typename Scalar>
ArmadilloMatrix<Scalar> operator-(ArmadilloMatrix<Scalar> lhs,
                                  const ArmadilloMatrix<Scalar>& rhs) {
    lhs -= rhs;
    return lhs;
}

template <typename Scalar>
ArmadilloMatrix<Scalar> operator+(ArmadilloMatrix<Scalar> lhs,
                                  const ArmadilloMatrix<Scalar>& rhs) {
    lhs += rhs;
    return lhs;
}

//
// ---------------------------------------------------
//

template <typename Scalar>
void swap(ArmadilloMatrix<Scalar>& first, ArmadilloMatrix<Scalar>& second) {
    using std::swap;
    typedef typename ArmadilloMatrix<Scalar>::base_type base_type;
    swap(static_cast<base_type&>(first), static_cast<base_type&>(second));
    first.m_arma.swap(second.m_arma);
}

//
// Armadillo matrix
//
template <typename Scalar>
ArmadilloMatrix<Scalar>::ArmadilloMatrix(size_type n_rows, size_type n_cols,
                                         bool fill_zero)
      : m_arma(n_cols, n_rows, arma::fill::none) {
    // Note that the armadillo matrix stores the entries in transposed
    // form
    if (fill_zero) {
        // set all elements to zero
        m_arma.zeros();
    }
}

template <typename Scalar>
ArmadilloMatrix<Scalar>::ArmadilloMatrix(const ArmadilloMatrix& mat,
                                         scalar_type tolerance)
      : ArmadilloMatrix(mat.n_rows(), mat.n_cols(), false) {
    for (const auto elem : mat) {
        if (std::fabs(*elem) < tolerance) {
            (*this)(elem.row(), elem.col()) = Constants<scalar_type>::zero();
        } else {
            (*this)(elem.row(), elem.col()) = *elem;
        }
    }
}

template <typename Scalar>
ArmadilloMatrix<Scalar>::ArmadilloMatrix(
      std::initializer_list<std::initializer_list<scalar_type>> list_of_lists)
      : ArmadilloMatrix(
              list_of_lists.size(),
              list_of_lists.size() > 0 ? list_of_lists.begin()->size() : 0,
              false) {
#ifdef DEBUG
    size_type n_rows = list_of_lists.size();
    size_type n_cols = n_rows > 0 ? list_of_lists.begin()->size() : 0;
#endif

    // Assert all columns have equal length.
    assert_element_sizes(list_of_lists, n_cols);

    size_type i = 0;
    for (auto row : list_of_lists) {
        size_type j = 0;
        for (scalar_type elem : row) {
            (*this)(i, j) = elem;
            ++j;
        }
        assert_dbg(j == n_cols, ExcInternalError());
        ++i;
    }

    assert_dbg(i == n_rows, ExcInternalError());
}

template <typename Scalar>
ArmadilloMatrix<Scalar>::ArmadilloMatrix(storage_type inner) : m_arma(inner) {}

#else
template <typename Scalar>
class ArmadilloMatrix {
    static_assert(false, "ArmadilloMatrix is not available");
};
#endif  // LINALGWRAP_HAVE_ARMADILLO

}  // namespace linalgwrap
