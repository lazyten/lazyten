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
#include "linalgwrap/ArmadilloMatrix.hh"

/** \file BlockView.special.hh
 *  \brief Partial template specialisations for BlockView
 */

namespace linalgwrap {
namespace view {
namespace detail {
#ifdef LINALGWRAP_HAVE_ARMADILLO
/** \brief Class implementing the type-specific TransposeView functionality.
 *
 * Specialisation for SmallMatrices
 *
 * \tparam Matrix     The actual matrix type of which this is a view
 * */
template <typename Matrix, typename Scalar>
class BlockViewSpecialise<Matrix, ArmadilloMatrix<Scalar>>
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
#endif  // LINALGWRAP_HAVE_ARMADILLO

//
// ---------------------------------------------------
//

#ifdef LINALGWRAP_HAVE_ARMADILLO
template <typename Matrix, typename Scalar>
BlockViewSpecialise<Matrix, ArmadilloMatrix<Scalar>>::BlockViewSpecialise(
      inner_matrix_type& mat, Range<size_type> row_range,
      Range<size_type> col_range)
      : base_type(mat, std::move(row_range), std::move(col_range)) {}

template <typename Matrix, typename Scalar>
typename BlockViewSpecialise<Matrix,
                             ArmadilloMatrix<Scalar>>::stored_matrix_type
      BlockViewSpecialise<Matrix, ArmadilloMatrix<Scalar>>::operator*(
            const stored_matrix_type& m) const {
    assert_size(base_type::n_cols(), m.n_rows());

    assert_dbg(!this->m_row_range.empty(), ExcInternalError());
    assert_dbg(!this->m_col_range.empty(), ExcInternalError());

    // Translate ranges to armadillo spans (which are closed intervals)
    arma::span rows(this->m_row_range.first(), this->m_row_range.last() - 1);
    arma::span cols(this->m_col_range.first(), this->m_col_range.last() - 1);

    // Access the inner data of both matrices:
    typedef typename stored_matrix_type::storage_type storage_type;
    const storage_type& idata = base_type::inner_matrix().data();
    const storage_type& mdata = m.data();

    // Note that idata and mdata are in fact the transposes of the
    // matrices we want to represent (see comment above class definition
    // of ArmadilloMatrix)
    // But we also store the transpose inside the Armadillo matrix
    // we return (because there it is the same deal), hence
    // we need to determine
    // \f[ (inner_matrix * m)^T = m^T * inner_matrix^T = mdata*idata \f]
    //
    // This also explains the reverse order of cols and rows.

    // Find the product:
    storage_type res = mdata * idata(cols, rows);

    // Return the product:
    return stored_matrix_type{std::move(res)};
}
#endif  // LINALGWRAP_HAVE_ARMADILLO

}  // namespace detail
}  // namespace view
}  // namespace linalgwrap
