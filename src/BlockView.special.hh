#pragma once
#include <SmallMatrix.hh>

/** \file BlockView.special.hh
 *  \brief Partial template specialisations for BlockView
 */

namespace linalgwrap {
namespace view {
/** \brief Class implementing the type-specific TransposeView functionality.
 *
 * Specialisation for SmallMatrices
 *
 * \tparam Matrix     The actual matrix type of which this is a view
 * */
template <typename Matrix, typename Scalar>
class BlockViewSpecialise<Matrix, SmallMatrix<Scalar>>
      : public BlockViewBase<Matrix> {
  public:
    /** \name Typedefs of standard types
     */
    ///@{
    typedef BlockViewBase<Matrix> base_type;
    typedef typename base_type::inner_matrix_type inner_matrix_type;
    typedef typename base_type::stored_matrix_type stored_matrix_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::lazy_matrix_expression_ptr_type
          lazy_matrix_expression_ptr_type;
    ///@}

    /** Constructor to pass arguments through
     *
     * \param mat      Matrix to take a block of
     * \param row_range Range of row indices in the block
     * \param col_range Range of col indices in the block
     **/
    BlockViewSpecialise(inner_matrix_type& mat, Range<size_type> row_range,
                        Range<size_type> col_range);

    /** \brief Multiplication with a stored matrix */
    stored_matrix_type operator*(const stored_matrix_type& m) const override;
};

template <typename Matrix, typename Scalar>
BlockViewSpecialise<Matrix, SmallMatrix<Scalar>>::BlockViewSpecialise(
      inner_matrix_type& mat, Range<size_type> row_range,
      Range<size_type> col_range)
      : base_type(mat, std::move(row_range), std::move(col_range)) {}

template <typename Matrix, typename Scalar>
typename BlockViewSpecialise<Matrix, SmallMatrix<Scalar>>::stored_matrix_type
      BlockViewSpecialise<Matrix, SmallMatrix<Scalar>>::
      operator*(const stored_matrix_type& m) const {
    assert_size(base_type::n_cols(), m.n_rows());

    // At least one range is empty -> no work to be done:
    if (this->m_row_range.is_empty() || this->m_col_range.is_empty()) {
        return stored_matrix_type{this->m_row_range.length(), m.n_rows()};
    }

    // Translate ranges to armadillo spans (which are closed intervals)
    arma::span rows(this->m_row_range.first(), this->m_row_range.last() - 1);
    arma::span cols(this->m_col_range.first(), this->m_col_range.last() - 1);

    // Access the inner data of both matrices:
    typedef typename stored_matrix_type::storage_type storage_type;
    const storage_type& idata = base_type::inner_matrix().data();
    const storage_type& mdata = m.data();

    // Find the product:
    storage_type res = idata(rows, cols) * mdata;

    // Return the product:
    return stored_matrix_type{std::move(res)};
}

}  // namespace view
}  // namespace linalgwrap
