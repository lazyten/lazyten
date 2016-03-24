#pragma once
#include "ViewBase.hh"

namespace linalgwrap {
namespace view {

/** Implementation of fallback base functionality for all BlockViews */
template <typename Matrix>
class BlockViewBase : public ViewBase<Matrix> {
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
    friend void swap(BlockViewBase& first, BlockViewBase& second);

    /* \brief Construct from matrix and row/column Range
     *
     * Allows to select a block of entries within the provided
     * Range of rows and columns.
     *
     * \param mat      Matrix to view a block of
     * \param row_range Range of rows to keep
     * \param col_range Range of columns to keep
     */
    BlockViewBase(inner_matrix_type& mat, Range<size_type> row_range,
                  Range<size_type> col_range);

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
    stored_matrix_type operator*(const stored_matrix_type& m) const override;

  protected:
    Range<size_type> m_row_range;
    Range<size_type> m_col_range;
};

/** \brief Class implementing the type-specific BlockView functionality.
 *
 * Specialisation can make use of the second template argument.
 * This represents the matrix without const qualifier. This makes
 * it easier to provide a common specialisation for const and non-const
 * Matrices at once if wished.
 *
 * \tparam Matrix     The actual matrix type of which this is a view
 * */
template <typename Matrix, typename = typename std::remove_const<Matrix>::type>
class BlockViewSpecialise : public BlockViewBase<Matrix> {
  public:
    /** \name Typedefs of standard types
     */
    ///@{
    typedef BlockViewBase<Matrix> base_type;
    typedef typename base_type::inner_matrix_type inner_matrix_type;
    typedef typename base_type::size_type size_type;
    ///@}

    /** Constructor to pass arguments through
     *
     * \param mat      Matrix to take a block of
     * \param row_range Range of row indices in the block
     * \param col_range Range of col indices in the block
     **/
    BlockViewSpecialise(inner_matrix_type& mat, Range<size_type> row_range,
                        Range<size_type> col_range);
};

/** Top-Level BlockView class
 *
 * \tparam Matrix  The matrix type, which we view upon.
 * */
template <typename Matrix>
class BlockView : public BlockViewSpecialise<Matrix> {
  public:
    /** \name Typedefs of standard types
     */
    ///@{
    typedef BlockViewSpecialise<Matrix> base_type;
    typedef typename base_type::inner_matrix_type inner_matrix_type;
    typedef typename base_type::stored_matrix_type stored_matrix_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::lazy_matrix_expression_ptr_type
          lazy_matrix_expression_ptr_type;
    ///@}

    /* \brief Construct from matrix and row/column Range
     *
     * Allows to select a block of entries within the provided
     * Range of rows and columns.
     *
     * \param mat      Matrix to view a block of
     * \param row_range Range of rows to keep
     * \param col_range Range of columns to keep
     */
    BlockView(inner_matrix_type& mat, Range<size_type> row_range,
              Range<size_type> col_range);

    /** \brief Clone the view */
    lazy_matrix_expression_ptr_type clone() const override;
};

/** Convenience function to make a BlockView
 *
 *\param m      The Matrix to take a block of
 *\param rows   The range of row indices to consider
 *\param cols   The range of column indices to consider
 */
template <typename Matrix>
BlockView<Matrix> block(Matrix& m, Range<typename Matrix::size_type> rows,
                        Range<typename Matrix::size_type> cols);

// TODO have a special structure for a column view?
/** Convenience function to make a BlockView which selects
 *  a range of columns (i.e. a column view)
 *
 *\param m      The Matrix to take some columns of
 *\param cols   The range of columns to consider
 */
template <typename Matrix>
BlockView<Matrix> columns(Matrix& m, Range<typename Matrix::size_type> cols);

// TODO have a special structure for a row view?
/** Convenience function to make a BlockView which selects
 *  a range of rows (i.e. a row view)
 *
 *\param m      The Matrix to take some columns of
 *\param rows   The range of rows to consider
 */
template <typename Matrix>
BlockView<Matrix> rows(Matrix& m, Range<typename Matrix::size_type> rows);

//
// -------------------------------------------------------
//

template <typename Matrix>
BlockViewBase<Matrix>::BlockViewBase(inner_matrix_type& mat,
                                     Range<size_type> row_range,
                                     Range<size_type> col_range)
      : base_type(mat, "BlockView"),
        m_row_range(row_range),
        m_col_range(col_range) {

    // TODO at the moment the case of empty ranges is not allowed
    assert_dbg(!row_range.is_empty(), ExcNotImplemented());
    assert_dbg(!col_range.is_empty(), ExcNotImplemented());

    if (!row_range.is_empty()) {
        // Assert that the range makes sense:
        assert_greater_equal(row_range.last(),
                             base_type::inner_matrix().n_rows());
    }

    if (!col_range.is_empty()) {
        assert_greater_equal(col_range.last(),
                             base_type::inner_matrix().n_cols());
    }
}

template <typename Matrix>
typename BlockViewBase<Matrix>::size_type BlockViewBase<Matrix>::n_rows()
      const {
    return m_row_range.length();
}

template <typename Matrix>
typename BlockViewBase<Matrix>::size_type BlockViewBase<Matrix>::n_cols()
      const {
    return m_col_range.length();
}

template <typename Matrix>
typename BlockViewBase<Matrix>::scalar_type BlockViewBase<Matrix>::operator()(
      size_type row, size_type col) const {
    assert_range(0, row, n_rows());
    assert_range(0, col, n_cols());

    // Return the element plus the offset:
    return base_type::inner_matrix()(row + m_row_range.first(),
                                     col + m_col_range.first());
}

template <typename Matrix>
inline typename BlockViewBase<Matrix>::stored_matrix_type
BlockViewBase<Matrix>::extract_block(Range<size_type> row_range,
                                     Range<size_type> col_range) const {
    // At least one range is empty -> no work to be done:
    if (row_range.is_empty() || col_range.is_empty()) {
        return stored_matrix_type{row_range.length(), col_range.length()};
    }

    // Assertive checks:
    assert_greater_equal(row_range.last(), n_rows());
    assert_greater_equal(col_range.last(), n_cols());

    // The offsets due to the view action
    size_type row_offset = m_row_range.first();
    size_type col_offset = m_col_range.first();

    // The effective ranges
    Range<size_type> erows{row_offset + row_range.first(),
                           row_offset + row_range.last()};
    Range<size_type> ecols{col_offset + col_range.first(),
                           col_offset + col_range.last()};

    return base_type::inner_matrix().extract_block(erows, ecols);
}

template <typename Matrix>
inline void BlockViewBase<Matrix>::add_block_to(stored_matrix_type& in,
                                                size_type start_row,
                                                size_type start_col,
                                                scalar_type c_this) const {
    // check that we do not overshoot the indices
    assert_greater_equal(start_row + in.n_rows(), n_rows());
    assert_greater_equal(start_col + in.n_cols(), n_cols());

    // The offsets due to the view action
    size_type row_offset = m_row_range.first();
    size_type col_offset = m_col_range.first();

    // Perform action at inner matrix
    base_type::inner_matrix().add_block_to(in, start_row + row_offset,
                                           start_col + col_offset, c_this);
}

template <typename Matrix>
inline typename BlockViewBase<Matrix>::stored_matrix_type
      BlockViewBase<Matrix>::
      operator*(const stored_matrix_type& m) const {

    assert_size(n_cols(), m.n_rows());

    // If this structure does not actually change the number of cols, we can
    // just multiply in the inner matrix:
    if (base_type::inner_matrix().n_cols() == m.n_rows()) {
        assert_size(base_type::inner_matrix().n_cols(), n_cols());
        return base_type::inner_matrix() * m;
    }

    // For lazy matrices it is probably best to inform the user that we need
    // to do some heavy operation in the other case and disable the operation.
    assert_dbg(IsStoredMatrix<Matrix>::value,
               ExcDisabled(
                     "The operation \"BlockView<LazyMatrix> * StoredMatrix\" "
                     "is disabled because it is potentially pretty expensive. "
                     "You should either rearrange the order in which the "
                     "operations are performed or else explicitly convert a "
                     "subexpression to its stored matrix type."));

    // TODO
    // For all other cases we need to either extract the block and perform
    // the multiplication on it or else go over each element via an iterator
    // and do the multiplication manually. Both can be problematic in some
    // cases, so we bail out by throwing a not-implemented exception.
    // One can perhaps find some nice trade-off when to do what someday ...
    assert_dbg(false, ExcNotImplemented());

    // Return a dummy
    return stored_matrix_type{n_rows(), m.n_cols()};
}

//
// general BlockViewSpecialise
//
template <typename Matrix, typename MatrixBare>
BlockViewSpecialise<Matrix, MatrixBare>::BlockViewSpecialise(
      inner_matrix_type& mat, Range<size_type> row_range,
      Range<size_type> col_range)
      : base_type(mat, std::move(row_range), std::move(col_range)) {}

//
// BlockView
//
template <typename Matrix>
BlockView<Matrix>::BlockView(inner_matrix_type& mat, Range<size_type> row_range,
                             Range<size_type> col_range)
      : base_type(mat, std::move(row_range), std::move(col_range)) {}

template <typename Matrix>
typename BlockView<Matrix>::lazy_matrix_expression_ptr_type
BlockView<Matrix>::clone() const {
    return lazy_matrix_expression_ptr_type(new BlockView(*this));
}

//
// out-of-scope functions:
//
template <typename Matrix>
BlockView<Matrix> block(Matrix& m, Range<typename Matrix::size_type> rows,
                        Range<typename Matrix::size_type> cols) {
    return BlockView<Matrix>{m, rows, cols};
}

template <typename Matrix>
BlockView<Matrix> columns(Matrix& m, Range<typename Matrix::size_type> cols) {
    return BlockView<Matrix>{m, range(m.n_rows()), cols};
}

template <typename Matrix>
BlockView<Matrix> rows(Matrix& m, Range<typename Matrix::size_type> rows) {
    return BlockView<Matrix>{m, rows, range(m.n_cols())};
}

}  // view
}  // linalgwrap

// Include the specialisations for BlockView.hh
#include "BlockView.special.hh"
