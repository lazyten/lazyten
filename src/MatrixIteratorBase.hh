#pragma once
#include <iterator>
#include "Exceptions.hh"
#include "SubscriptionPointer.hh"
#include "Constants.hh"

namespace linalgwrap {
// TODO does it make sense to have a column-by-column iterator?
// Pros: makes writing matrix-matrix multiplication easier
// Cons: usually sparse matrices are CRS matrices, hence column-vise
//       traversing is expensive and we should not make the user
//       feel it is easy.

/** \brief MatrixIterator base class
 *
 * Provides the basic interface and functionality for Matrix iterators.
 *
 * \tparam  Matrix     The matrix type to iterate over
 * \tparam  constness  Is this iterator an iterator over
 *                     constant values or not.
 *
 * In addition to the given interface all MatrixIterators are expected to
 * provide the following:
 *
 * - Default constructor
 * - Copy constructor
 * - Construction from Matrix reference and index_type
 *   ```
 *   MatrixIteratorBase(Matrix&, index_type)
 *   ```
 * - Copy assignment
 * - prefix increment
 *   ```
 *   MatrixIteratorBase& operator++();
 *   ```
 * - postfix increment
 *   ```
 *   MatrixIteratorBase operator++(int);
 *   ```
 * - Seek increment
 *   ```
 *   MatrixIteratorBase& seek_to(index_type element);
 *   ```
 */
template <typename Matrix, bool Constness>
struct MatrixIteratorBase : public std::iterator<std::forward_iterator_tag,
                                                 typename Matrix::scalar_type> {
    typedef typename std::conditional<Constness, Matrix, const Matrix>::type
          matrix_type;
    typedef typename matrix_type::size_type size_type;
    typedef typename matrix_type::scalar_type scalar_type;

    typedef std::pair<size_type, size_type> index_type;

    typedef std::iterator<std::forward_iterator_tag,
                          typename Matrix::scalar_type> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::reference reference;
    typedef typename base_type::pointer pointer;

    /** Invalid iterator position */
    static constexpr std::pair<size_type, size_type> invalid_pos = {
          Constants<size_type>::invalid, Constants<size_type>::invalid};

    //
    // Information about iterator
    //
    /** Get row of the currently pointed to element */
    virtual size_type row() const = 0;

    /** Get column of currently pointed to element */
    virtual size_type col() const = 0;

    /** Return the tuple of indices of the currently pointed to element */
    index_type indices() const;

    //
    // Element access
    //
    /** Return the value of the element we point to if this is a const
     *  iterator, else a reference to the element.
     */
    virtual typename std::conditional<Constness, value_type, reference>::type
    operator*() const = 0;

    /** Access the members of the element we point to. */
    virtual typename std::conditional<Constness, const pointer, pointer>::type
    operator->() const = 0;

    //
    // Comparison
    //
    /** \brief check if two iterators are equal
     *
     * Matrix iterators are equal if they point to the same element
     */
    bool operator==(const MatrixIteratorBase& other) const {
        return other.indices() == indices();
    }

    /** \brief Check whether two iterators are unequal
     *
     * Matrix iterators are unequal if they point to different elements
     */
    bool operator!=(const MatrixIteratorBase& other) const {
        return !(*this == other);
    }
};

/** \brief Base class for matrix iterators that iterate over all
 *         entries of a matrix and do not take sparsity of any
 *         kind into account.
 *
 * This class is not yet a full iterator. It only provides some
 * helper functionality, which should make it easier to actually
 * implement the dereferencing and the element access later on.
 */
template <typename Matrix, bool Constness>
class FullMatrixIteratorBase : public MatrixIteratorBase<Matrix, Constness> {
  public:
    typedef MatrixIteratorBase<Matrix, Constness> base_type;
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
    FullMatrixIteratorBase();

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
    FullMatrixIteratorBase(matrix_type& mat, index_type start_index);

    //
    // Matrix information
    //

    /** The row currently pointed to */
    size_type row() const override;

    /** The column currently pointed to */
    size_type col() const override;

  protected:
    /** \brief Seek to the next element updating any internal
     *         state (i.e. indices) as we go along)*/
    void seek_next_element();

    /** \brief Seek to the provided  element updating any internal
     *         state (i.e. indices) as we go along)*/
    void seek_to_element(index_type element);

    /** \brief
     *  Assert that the internal state of the index and the matrix
     *  pointer makes sense. If not throw an exception in debug mode.
     *
     *  It is good practice to use this before every access to the
     *  iterator data.
     *  */
    void assert_valid_state() const;

    /** Access to the matrix */
    matrix_type& matrix() const;

  private:
    index_type m_index;
    SubscriptionPointer<matrix_type> m_matrix_ptr;
};

//
// -------------------------------------------------------------
//
// MatrixIteratorBase

template <typename Matrix, bool Constness>
typename MatrixIteratorBase<Matrix, Constness>::index_type
MatrixIteratorBase<Matrix, Constness>::indices() const {
    return make_pair(row(), col());
}

// Define invalid_pos (needs to be there, otherwise linker error)
template <typename Matrix, bool Constness>
constexpr std::pair<typename MatrixIteratorBase<Matrix, Constness>::size_type,
                    typename MatrixIteratorBase<Matrix, Constness>::size_type>
      MatrixIteratorBase<Matrix, Constness>::invalid_pos;

//
// FullMatrixIteratorBase
//
template <typename Matrix, bool Constness>
void FullMatrixIteratorBase<Matrix, Constness>::assert_valid_state() const {
    assert_dbg(m_matrix_ptr,
               ExcInvalidState("MatrixIterator does not point to any matrix"));

    assert_dbg(m_index.first != base_type::invalid_pos.first &&
                     m_index.second != base_type::invalid_pos.second,
               ExcIteratorPastEnd());
}

template <typename Matrix, bool Constness>
FullMatrixIteratorBase<Matrix, Constness>::FullMatrixIteratorBase()
      : m_index{base_type::invalid_pos},
        m_matrix_ptr{"FullMatrixIteratorBase"} {}

template <typename Matrix, bool Constness>
FullMatrixIteratorBase<Matrix, Constness>::FullMatrixIteratorBase(
      matrix_type& mat, index_type start_index)
      : m_index{start_index}, m_matrix_ptr{"FullMatrixIteratorBase", mat} {
    // if we are not anyway representing an invalid iterator,
    // assert that we are valid.
    if (m_index != base_type::invalid_pos) {
        assert_valid_state();
    }
}

template <typename Matrix, bool Constness>
typename FullMatrixIteratorBase<Matrix, Constness>::size_type
FullMatrixIteratorBase<Matrix, Constness>::row() const {
    assert_valid_state();
    return m_index.first;
}

template <typename Matrix, bool Constness>
typename FullMatrixIteratorBase<Matrix, Constness>::size_type
FullMatrixIteratorBase<Matrix, Constness>::col() const {
    assert_valid_state();
    return m_index.second;
}

template <typename Matrix, bool Constness>
void FullMatrixIteratorBase<Matrix, Constness>::seek_next_element() {
    assert_valid_state();

    // Unpack pair:
    size_type row = m_index.first;
    size_type col = m_index.second;

    // Row-wise seek to the next valid index
    if (col + 1 < m_matrix_ptr->n_cols()) {
        seek_to_element({row, col + 1});
    } else if (row + 1 < m_matrix_ptr->n_rows()) {
        seek_to_element({row + 1, 0});
    } else {
        // Make invalid:
        m_index = base_type::invalid_pos;
    }
}

template <typename Matrix, bool Constness>
void FullMatrixIteratorBase<Matrix, Constness>::seek_to_element(
      index_type element) {
    assert_valid_state();

    // Assert that the indices are not too large:
    assert_upper_bound(element.first, m_matrix_ptr->n_rows());
    assert_upper_bound(element.second, m_matrix_ptr->n_cols());

    // Assert that we make progress in the right direction:
    assert_lower_bound(element.first, m_index.first);

    if (element.first == m_index.first) {
        assert_lower_bound(element.second, m_index.second);
    }

    // Set the indices:
    m_index = element;
}

template <typename Matrix, bool Constness>
typename FullMatrixIteratorBase<Matrix, Constness>::matrix_type&
FullMatrixIteratorBase<Matrix, Constness>::matrix() const {
    return *m_matrix_ptr;
}

}  // namespace linalgwrap
