#ifndef LINALG_MATRIX_I_HPP_
#define LINALG_MATRIX_I_HPP_

#include <cstddef>
#include <utility>
#include "SmallMatrix.hh"
#include <string>
#include "Exceptions.hh"
#include "Constants.hh"

namespace linalgwrap {

template <typename Scalar>
class SmallMatrix;

/** \brief Abstract matrix interface class */
template <typename Scalar>
class Matrix_i {
  public:
    typedef size_t size_type;
    typedef Scalar scalar_type;

    /** \brief Destructor */
    virtual ~Matrix_i() = default;

    /** \brief Extract values of this matrix partially to \p block
     *
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
    virtual scalar_type operator()(size_type row, size_type col) const {
        SmallMatrix<scalar_type> block(1, 1, false);
        fill(row, col, block);
        return block(0, 0);
    }

    /** \brief return an element of the vectorised matrix object
     *
     * It is advisible to overload this in order to get a more performant
     * implementation.
     */
    virtual scalar_type operator[](size_type i) const {
        size_type i_row = i / n_rows();
        size_type i_col = i % n_rows();
        return (*this)(i_row, i_col);
    }

    /** \brief Number of rows of the matrix */
    virtual size_type n_rows() const = 0;

    /** \brief Number of columns of the matrix */
    virtual size_type n_cols() const = 0;

    /** \brief Return the inverse of the matrix
    void inverse() { assert_dbg(false, ExcNotImplemented()); }
    */

    /** \brief Return the transpose this matrix
    void transpose() { assert_dbg(false, ExcNotImplemented()); }
    */
};
}  // namespace linalg
#endif  // LINALG_MATRIX_I_HPP_
