#pragma once

#include "ViewBase.hh"

namespace linalgwrap {
namespace view {

template <typename Matrix>
class ScaleView : public ViewBase<Matrix> {
  public:
    /** \name Typedefs of standard types
     */
    ///@{
    typedef ViewBase<Matrix> base_type;
    typedef typename base_type::inner_matrix_type inner_matrix_type;
    typedef typename base_type::stored_matrix_type stored_matrix_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::lazy_matrix_expression_ptr_type
          lazy_matrix_expression_ptr_type;
    ///@}

    /** A swap function for ScaleView */
    friend void swap(ScaleView& first, ScaleView& second);

    /* \brief Construct from matrix and scaling factor
     *
     * \param mat      Matrix to scale
     * \param scaling  Scaling factor to scale all elements with
     */
    ScaleView(inner_matrix_type& mat, scalar_type scaling);

    //
    // Matrix_i interface
    //

    /** \brief Number of rows of the matrix */
    size_type n_rows() const override;

    /** \brief Number of columns of the matrix  */
    size_type n_cols() const override;

    /** \brief return an element of the matrix    */
    scalar_type operator()(size_type row, size_type col) const override;

    //
    // LazyMatrixExpression interface
    //
    /** \brief Extract a block of values out of the matrix and
     *         return it as a stored matrix of the appropriate size
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
    stored_matrix_type extract_block(Range<size_type> row_range,
                                     Range<size_type> col_range) const override;

    /** \brief Add a block of values of the matrix to the stored matrix
     *         provided by reference.
     *
     * For more details of the interface see the function of the same name
     * in ``LazyMatrixExpression``.
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
    void add_block_to(
          stored_matrix_type& in, size_type start_row, size_type start_col,
          scalar_type c_this = Constants<scalar_type>::one) const override;

    /** \brief Multiplication with a stored matrix */
    stored_matrix_type operator*(const stored_matrix_type& m) const override {
        assert_size(n_cols(), m.n_rows());
        return m_scaling * base_type::inner_matrix() * m;
    }

    /** \brief Clone the view */
    lazy_matrix_expression_ptr_type clone() const override;

  private:
    scalar_type m_scaling;
};

/** Convenience function to make a ScaleView
 *
 *\param m   The Matrix to scale
 *\param s   The scaling coefficient to apply to all matrix entries
 */
template <typename Matrix>
ScaleView<Matrix> scale(Matrix& m, typename Matrix::scalar_type s) {
    return ScaleView<Matrix>(m, s);
}

//
// ---------------------------------------------------
//

template <typename Matrix>
void swap(ScaleView<Matrix>& first, ScaleView<Matrix>& second) {
    using std::swap;  // enable ADL
    typedef ViewBase<Matrix> base_type;

    swap(first.m_scaling, second.m_scaling);
    swap(static_cast<base_type&>(first), static_cast<base_type&>(second));
}

template <typename Matrix>
ScaleView<Matrix>::ScaleView(inner_matrix_type& mat, scalar_type scaling)
      : base_type{mat, "ScaleView"}, m_scaling{scaling} {}

template <typename Matrix>
typename ScaleView<Matrix>::size_type ScaleView<Matrix>::n_rows() const {
    return base_type::inner_matrix().n_rows();
}

template <typename Matrix>
typename ScaleView<Matrix>::size_type ScaleView<Matrix>::n_cols() const {
    return base_type::inner_matrix().n_cols();
}

template <typename Matrix>
typename ScaleView<Matrix>::scalar_type ScaleView<Matrix>::operator()(
      size_type row, size_type col) const {
    return m_scaling * base_type::inner_matrix()(row, col);
}

template <typename Matrix>
typename ScaleView<Matrix>::stored_matrix_type ScaleView<Matrix>::extract_block(
      Range<size_type> row_range, Range<size_type> col_range) const {
    return m_scaling *
           base_type::inner_matrix().extract_block(row_range, col_range);
}

template <typename Matrix>
void ScaleView<Matrix>::add_block_to(stored_matrix_type& in,
                                     size_type start_row, size_type start_col,
                                     scalar_type c_this) const {
    base_type::inner_matrix().add_block_to(in, start_row, start_col,
                                           c_this * m_scaling);
}

template <typename Matrix>
typename ScaleView<Matrix>::lazy_matrix_expression_ptr_type
ScaleView<Matrix>::clone() const {
    return lazy_matrix_expression_ptr_type(new ScaleView(*this));
}

}  // view
}  // linalgwrap
