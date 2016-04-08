#pragma once
#include <ArmadilloMatrix.hh>

/** \file TransposeView.special.hh
 *  \brief Partial template specialisations for TransposeView.
 */

namespace linalgwrap {
namespace view {
//
// ArmadialloMatrix
//
#ifdef LINALGWRAP_HAVE_ARMADILLO
/** \brief Class implementing the type-specific TransposeView functionality.
 *
 * Specialisation for SmallMatrices
 *
 * \tparam Matrix     The actual matrix type of which this is a view
 * */
template <typename Matrix, typename Scalar>
class TransposeViewSpecialise<Matrix, ArmadilloMatrix<Scalar>>
      : public TransposeViewBase<Matrix> {
  public:
    /** \name Typedefs of standard types
     */
    ///@{
    typedef TransposeViewBase<Matrix> base_type;
    typedef typename base_type::inner_matrix_type inner_matrix_type;
    typedef typename base_type::stored_matrix_type stored_matrix_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::lazy_matrix_expression_ptr_type
          lazy_matrix_expression_ptr_type;
    ///@}

    /** Constructor to pass arguments through
     *
     * \param mat      Matrix to transpose
     **/
    explicit TransposeViewSpecialise(inner_matrix_type& mat);

    /** \brief Convert the TransposeView<SmallMatrix> to the
     *         corresponding SmallMatrix.
     *
     * Achieved by converting the inner armadillo matrix to
     * its transpose.
     */
    explicit operator stored_matrix_type() const override;

    // TODO perhaps specialise operator* between the view and an
    //      ArmadilloMatrix<Scalar> as well
};
#endif  // LINALGWRAP_HAVE_ARMADILLO

//
// ---------------------------------------------------
//

//
// ArmadilloMatrix
//
#ifdef LINALGWRAP_HAVE_ARMADILLO
template <typename Matrix, typename Scalar>
TransposeViewSpecialise<Matrix, ArmadilloMatrix<Scalar>>::
      TransposeViewSpecialise(inner_matrix_type& mat)
      : base_type(mat) {}

template <typename Matrix, typename Scalar>
TransposeViewSpecialise<Matrix, ArmadilloMatrix<Scalar>>::
operator stored_matrix_type() const {
    typedef typename inner_matrix_type::storage_type storage_type;
    const storage_type& matrix = base_type::inner_matrix().data();

    // Return a transpose copy of matrix:
    // (Armadillo-specific call)
    return SmallMatrix<Scalar>{std::move(matrix.t())};
}
#endif  // LINALGWRAP_HAVE_ARMADILLO

}  // namespace view
}  // namespace linalgwrap
