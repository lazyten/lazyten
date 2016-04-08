#pragma once
#include "linalgwrap/VectorOf.hh"
#include "linalgwrap/view/ViewBase.hh"
#include "linalgwrap/view/make_view.hh"

namespace linalgwrap {
namespace view {
namespace detail {

/** Implementation of fallback base functionality for all TransposeViews */
template <typename Matrix>
class TransposeViewBase : public ViewBase<Matrix> {
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

    /* \brief Construct from an inner matrix */
    TransposeViewBase(inner_matrix_type& mat);

    //
    // Matrix_i interface
    //
    /** \brief Number of rows of the matrix */
    size_type n_rows() const override;

    /** \brief Number of columns of the matrix */
    size_type n_cols() const override;

    /** \brief return an element of the matrix */
    scalar_type operator()(size_type row, size_type col) const override;

    //
    // LazyMatrixExpression interface
    //
    /** \brief Convert the TransposeView<Matrix> to the
     *         corresponding stored matrix.
     *
     * Achieved by converting the inner matrix to its stored
     * matrix and then taking the transpose view of the result.
     */
    explicit operator stored_matrix_type() const override;

    /** \brief Extract a block of values out of the matrix and
     *         return it as a stored matrix of the appropriate size
     *
     * For more details of the interface see the function of the same
     * name in ``LazyMatrixExpression``.
     * \param row_range   The Range object representing the range of rows
     *                    to extract. Note that it is a half-open interval
     *                    i.e. the LHS is inclusive, but the RHS not.
     *                    The Range may not be empty.
     * \param col_range   The Range object representing the range of
     *                    columns to extract.
     *                    The Range may not be empty.
     */
    stored_matrix_type extract_block(Range<size_type> row_range,
                                     Range<size_type> col_range) const override;

    /** \brief Multiplication with a stored matrix */
    stored_matrix_type operator*(const stored_matrix_type& m) const override;
};

/** \brief Class implementing the type-specific TransposeView functionality.
 *
 * Specialisation can make use of the second template argument.
 * This represents the matrix without const qualifier. This makes
 * it easier to provide a common specialisation for const and non-const
 * Matrices at once if wished.
 *
 * \tparam Matrix     The actual matrix type of which this is a view
 * */
template <typename Matrix, typename = typename std::remove_const<Matrix>::type>
class TransposeViewSpecialise : public TransposeViewBase<Matrix> {
  public:
    /** \name Typedefs of standard types
     */
    ///@{
    typedef TransposeViewBase<Matrix> base_type;
    typedef typename base_type::inner_matrix_type inner_matrix_type;
    ///@}

    /** Constructor to pass arguments through
     *
     * \param mat      Matrix to transpose
     **/
    explicit TransposeViewSpecialise(inner_matrix_type& mat);
};

/** Top-Level TransposeView class
 *
 * \tparam Matrix  The matrix type, which we view upon.
 * */
template <typename Matrix>
class TransposeView : public TransposeViewSpecialise<Matrix> {
  public:
    /** \name Typedefs of standard types
     */
    ///@{
    typedef TransposeViewSpecialise<Matrix> base_type;
    typedef typename base_type::inner_matrix_type inner_matrix_type;
    typedef typename base_type::stored_matrix_type stored_matrix_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::lazy_matrix_expression_ptr_type
          lazy_matrix_expression_ptr_type;
    ///@}

    /* \brief Construct from matrix
     *
     * \param mat      Matrix to transpose
     */
    explicit TransposeView(inner_matrix_type& mat);

    /** \brief Clone the view */
    lazy_matrix_expression_ptr_type clone() const override;
};

/** \brief Special multiplication operator for the scalar product */
template <typename MatrixType>
typename MatrixType::scalar_type operator*(
      const TransposeView<VectorOf<MatrixType>>& lhs,
      const VectorOf<MatrixType>& rhs);

//
// ---------------------------------------------------
//

//
// TransposeViewBase
//
template <typename Matrix>
TransposeViewBase<Matrix>::TransposeViewBase(inner_matrix_type& mat)
      : base_type(mat, "TransposeView") {}

template <typename Matrix>
typename TransposeViewBase<Matrix>::size_type
TransposeViewBase<Matrix>::n_rows() const {
    // Return the column number as row number
    return base_type::inner_matrix().n_cols();
}

template <typename Matrix>
typename TransposeViewBase<Matrix>::size_type
TransposeViewBase<Matrix>::n_cols() const {
    // Return the column number as row number
    return base_type::inner_matrix().n_rows();
}

template <typename Matrix>
typename TransposeViewBase<Matrix>::scalar_type TransposeViewBase<Matrix>::
operator()(size_type row, size_type col) const {
    // Return the element of the transpose part.
    return base_type::inner_matrix()(col, row);
}

template <typename Matrix>
inline TransposeViewBase<Matrix>::operator stored_matrix_type() const {
    // If we get here with a stored matrix then something
    // is wrong
    // If this assertion is removed we get a never-ending loop ...
    assert_dbg(!IsStoredMatrix<Matrix>::value,
               ExcDisabled(
                     "Conversion to stored matrix is not available for a "
                     "TransposeView of this stored matrix type. This usually "
                     "means that there has not been a proper specialisation of "
                     "TransposeViewSpecialise for this type of stored matrix"));

    // Get the converted matrix of the transposed matrix:
    stored_matrix_type mat_transposed =
          static_cast<stored_matrix_type>(base_type::inner_matrix());

    // Get a transpose view to it:
    TransposeView<stored_matrix_type> mat = view::transpose(mat_transposed);

    // Convert that view to stored matrix and return it:
    return static_cast<stored_matrix_type>(mat);
}

template <typename Matrix>
inline typename TransposeViewBase<Matrix>::stored_matrix_type
TransposeViewBase<Matrix>::extract_block(Range<size_type> row_range,
                                         Range<size_type> col_range) const {
    // Assertive checks:
    assert_greater(0, row_range.length());
    assert_greater(0, col_range.length());

    assert_greater_equal(row_range.last(), this->n_rows());
    assert_greater_equal(col_range.last(), this->n_cols());

    // Get the transpose block from the inner matrix
    stored_matrix_type transpose_extracted =
          base_type::inner_matrix().extract_block(col_range, row_range);

    // Perform the transposition of the stored matrix:
    TransposeView<stored_matrix_type> extracted =
          view::transpose(transpose_extracted);

    // Convert that view to stored matrix and return it:
    return static_cast<stored_matrix_type>(extracted);
}

template <typename Matrix>
inline typename TransposeViewBase<Matrix>::stored_matrix_type
      TransposeViewBase<Matrix>::operator*(const stored_matrix_type& m) const {
    // TODO: The default implementation done here is a really bad idea for
    // lazy matrices, since we need to calculate each lazy matrix element
    // one-by-one (Remember the philosophy of lazy matrices: This is what
    // we actually want to avoid)
    //
    // So the decision at the moment is to explicitly disable it.
    // If we really need it than we can still think of something clever.
    assert_dbg(
          base_type::is_stored_matrix_view,
          ExcDisabled(
                "The operation \"TransposeView<LazyMatrix> * "
                "StoredMatrix\" is disabled because it is pretty "
                "expensive in the current implementation. Allowed are only "
                "\"TransposeView<View_of_stored_matrix>*StoredMatrix\" or "
                "\"TransposeView<StoredMatrix>*StoredMatrix\". You have two "
                "options: 1. Rearrange your matrix expression (or simply "
                "change the evaluation order via bracketing) such that "
                "only the allowed type of operation do occur. "
                " 2. Explicitly convert a subexpression to its stored matrix "
                "type and perform the operation on the stored matrices."));

    assert_size(n_cols(), m.n_rows());

    // Allocate a matrix of zeros
    stored_matrix_type mat(n_rows(), m.n_cols(), true);

    // TODO seek an alternative for a product of two sparse matrices:

    // Explanaition of the loop below:
    //
    // We want the result:
    // res(i,col) = \sum_k inner^T(i,k) * m(k,col)
    //
    // So in fact:
    // res(i,col) = \sum_k inner(k,i) * m(k,col)
    //
    // Using iterators for inner
    //    k == iti.row()
    //    i == iti.col()
    //    inner(k,i) == *iti
    //
    // res(iti.col(),col) = \sum_{iti} *iti * m(iti.row(),col)
    //
    for (auto iti = std::begin(base_type::inner_matrix());
         iti != std::end(base_type::inner_matrix()); ++iti) {
        for (auto col : range(m.n_cols())) {
            mat(iti.col(), col) += *iti * m(iti.row(), col);
        }
    }
    return mat;
}

//
// general TransposeViewSpecialise
//

template <typename Matrix, typename MatrixBare>
TransposeViewSpecialise<Matrix, MatrixBare>::TransposeViewSpecialise(
      inner_matrix_type& mat)
      : base_type(mat) {}

//
// TransposeView
//
template <typename Matrix>
TransposeView<Matrix>::TransposeView(inner_matrix_type& mat) : base_type(mat) {}

template <typename Matrix>
typename TransposeView<Matrix>::lazy_matrix_expression_ptr_type
TransposeView<Matrix>::clone() const {
    return lazy_matrix_expression_ptr_type(new TransposeView(*this));
}

//
// out-of-scope functions:
//

template <typename MatrixType>
inline typename MatrixType::scalar_type operator*(
      const TransposeView<VectorOf<MatrixType>>& lhs,
      const VectorOf<MatrixType>& rhs) {
    return lhs.inner_matrix().dot(rhs);
}

}  // namespace detail
}  // view
}  // linalgwrap

// Include the specialisations for TransposeView.hh
#include "TransposeView.special.hh"
