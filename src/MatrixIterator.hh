#pragma once
#include "Matrix_i.hh"
#include "MatrixIteratorBase.hh"

namespace linalgwrap {
// forward declaration
template <typename Scalar>
class Matrix_i;

/** \brief Iterate over the elements of a Matrix_i object
 *         row-by-row.
 *
 * This implementation does not make any assumptions about
 * the inner matrix apart from the requirement that it should
 * satisfy the Matrix_i interface. It uses operator() of the
 * matrixto obtain values for the matrix entries and it does
 * not allow modification of the matrix entries at all.
 *
 * No attempt to skip elements based on sparsity or similar is
 * made, so a more specialised iterator should be used for
 * matrices with sparsity.
 *
 * TODO perhaps it is sensible to cache the value of
 *      operator(row,col) of the matrix and just return it on
 *      dereference of the iterator.
 */
template <typename Scalar>
class MatrixIterator
      : public FullMatrixIteratorBase<const Matrix_i<Scalar>, true> {
  public:
    typedef FullMatrixIteratorBase<const Matrix_i<Scalar>, true> base_type;
    typedef typename base_type::matrix_type matrix_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::index_type index_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::reference reference;
    typedef typename base_type::pointer pointer;

    //
    // Constructor, destructor and assignment
    //
    /** Default constructor */
    MatrixIterator();

    /** \brief Constructor from matrix reference and index of
     *  the element the iterator should point to
     *
     * In general it is assumed that the element index provided
     * is valid in the sense that it points to a proper matrix
     * element. There are exceptions however:
     * - The ``index_type`` value ``invalid_pos`` indicates an
     *   iterator-past-the-end of the matrix and should be used
     *   only in the end() function of a matrix.
     */
    MatrixIterator(matrix_type& mat, index_type start_index);

    //
    // Increment and seek
    //
    /** Prefix increment to the next non-zero element */
    MatrixIterator& operator++();

    /** Postfix increment to the next non-zero element */
    MatrixIterator operator++(int);

    /** Seek-increment to the next non-zero element after
     *  the provided index and return iterator to the new
     *  posititon
     */
    MatrixIterator& seek_to(index_type element);

    //
    // Element access
    //
    /** Return the value of the element we point to if this is a const
     *  iterator, else a reference to the element.
     */
    value_type operator*() const override;

    /** Access the members of the element we point to. */
    const pointer operator->() const override;

  private:
    //! Some dummy variable required for operator->
    mutable scalar_type dummy = Constants<scalar_type>::invalid;
};

//
// MatrixIterator
//

template <typename Scalar>
MatrixIterator<Scalar>::MatrixIterator()
      : base_type{} {}

template <typename Scalar>
MatrixIterator<Scalar>::MatrixIterator(matrix_type& mat, index_type start_index)
      : base_type{mat, start_index} {}

template <typename Scalar>
MatrixIterator<Scalar>& MatrixIterator<Scalar>::seek_to(index_type element) {
    base_type::seek_to_element(element);
    return *this;
}

template <typename Scalar>
MatrixIterator<Scalar>& MatrixIterator<Scalar>::operator++() {
    base_type::seek_next_element();
    return *this;
}

template <typename Scalar>
MatrixIterator<Scalar> MatrixIterator<Scalar>::operator++(int) {
    MatrixIterator copy{*this};
    ++(*this);
    return copy;
}

template <typename Scalar>
typename MatrixIterator<Scalar>::value_type MatrixIterator<Scalar>::operator*()
      const {
    return base_type::matrix()(base_type::row(), base_type::col());
}

template <typename Scalar>
const typename MatrixIterator<Scalar>::pointer MatrixIterator<Scalar>::
operator->() const {
    dummy = *(*this);
    return &dummy;
}

}  // linalgwrap
