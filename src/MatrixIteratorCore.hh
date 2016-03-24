#pragma once
#include <iterator>
#include <type_traits>
#include "Exceptions.hh"
#include "SubscriptionPointer.hh"
#include "Constants.hh"
#include "EnforceReference.hh"
#include "Matrix_i.hh"
#include "StoredMatrix_i.hh"

namespace linalgwrap {

// forward declaration
template <typename Scalar>
class Matrix_i;

template <typename Scalar>
class StoredMatrix_i;

/** \brief Base class providing the required inner functionality of an iterator
 *
 * This class should be overloaded for more special (sparsity-aware) traversal
 * of the matrix and the child should be used inside the MatrixIterator class.
 */
template <typename Matrix, bool Constness>
struct MatrixIteratorCoreBase
      : public std::iterator<std::forward_iterator_tag,
                             typename Matrix::scalar_type> {
    typedef Matrix original_matrix_type;

    typedef typename std::conditional<Constness, const Matrix, Matrix>::type
          matrix_type;
    typedef typename matrix_type::size_type size_type;
    typedef typename matrix_type::scalar_type scalar_type;

    typedef std::pair<size_type, size_type> index_type;

    typedef std::iterator<std::forward_iterator_tag,
                          typename Matrix::scalar_type> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::reference reference;
    typedef typename base_type::pointer pointer;

    static_assert(std::is_base_of<Matrix_i<scalar_type>, matrix_type>::value,
                  "Matrix must be a subclass of Matrix_i");

    static_assert(Constness || std::is_base_of<StoredMatrix_i<scalar_type>,
                                               matrix_type>::value,
                  "The iterator must be a const iterator (MatrixIteratorCore "
                  "has Constness == true) or Matrix must be a subclass of "
                  "StoredMatrix_i");

    /** Invalid iterator position */
    static constexpr std::pair<size_type, size_type> invalid_pos = {
          Constants<size_type>::invalid, Constants<size_type>::invalid};

    /** It the iterator a const iterator? */
    static constexpr bool is_const_iterator = Constness;

    //@{
    /** Defaults for the big five */
    virtual ~MatrixIteratorCoreBase() = default;
    MatrixIteratorCoreBase() = default;
    MatrixIteratorCoreBase(MatrixIteratorCoreBase&&) = default;
    MatrixIteratorCoreBase(const MatrixIteratorCoreBase&) = default;
    MatrixIteratorCoreBase& operator=(MatrixIteratorCoreBase&&) = default;
    MatrixIteratorCoreBase& operator=(const MatrixIteratorCoreBase&) = default;
    //@}

    //
    // Information about the current element:
    //
    /** Get row of the currently pointed to element */
    virtual size_type row() const = 0;

    /** Get column of currently pointed to element */
    virtual size_type col() const = 0;

    /** Return the tuple of indices of the currently pointed to element */
    index_type indices() const;

  protected:
    /** Assert that the current state of the core is valid
     *
     * Throw appropriate exceptions if not */
    virtual void assert_valid_state() const = 0;

    /** Obtain value of current element. */
    virtual typename std::conditional<Constness, value_type, reference>::type
    value() const = 0;

    /** Return a pointer to the currently pointed-to value. */
    virtual typename std::conditional<Constness, const pointer, pointer>::type
    ptr_to_value() const = 0;

    //
    // Seeking
    //
    /** \brief Seek to the next element updating any internal
     *         state (i.e. indices) as we go along)*/
    virtual void seek_next_element() = 0;

    /** \brief Seek to the provided  element updating any internal
     *         state (i.e. indices) as we go along)*/
    virtual void seek_to_element(index_type element) = 0;
};

/** \brief Default MatrixIteratorCoreBase implementation
 *
 * This implementation does not make any assumptions about
 * the inner matrix apart from the requirement that it should
 * satisfy the Matrix_i interface. It uses operator() of the
 * matrix to obtain values for the matrix entries and it does
 * not allow modification of the matrix entries at all.
 *
 * No attempt to skip elements based on sparsity or similar is
 * made, so a more specialised iterator should be used for
 * matrices with sparsity.
 *
 * Either build a now implementation for MatrixIteratorCoreBase
 * for your special matrix or use (partial) template specialisation
 * of this class to provide a more suitable iterator core for
 * your matrix structure.
 */
template <typename Matrix, bool Constness>
class MatrixIteratorDefaultCore
      : public MatrixIteratorCoreBase<Matrix, Constness> {
  public:
    typedef MatrixIteratorCoreBase<Matrix, Constness> base_type;
    typedef typename base_type::original_matrix_type original_matrix_type;
    typedef typename base_type::matrix_type matrix_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::index_type index_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::reference reference;
    typedef typename base_type::pointer pointer;

    //
    // Constructor, destructor and assignment
    //
    /** Default constructor */
    MatrixIteratorDefaultCore();

    /** \brief Constructor of a Matrix iterator pointing to
     *  the past-the-end position.
     *
     *  In other words this constructor constructs an iterator
     *  in the state past-the-end. */
    explicit MatrixIteratorDefaultCore(matrix_type& mat);

    /** \brief Construct an iterator giving the initial value
     * it should point to. */
    MatrixIteratorDefaultCore(matrix_type& mat, index_type start_index);

    //
    // Matrix information
    //

    /** The row currently pointed to */
    size_type row() const override;

    /** The column currently pointed to */
    size_type col() const override;

  protected:
    /** Obtain value of current element. */
    typename std::conditional<Constness, value_type, reference>::type value()
          const override;

    /** Return a pointer to the currently pointed-to value */
    typename std::conditional<Constness, const pointer, pointer>::type
    ptr_to_value() const override;

    //
    // Seeking and increment:
    //
    /** \brief Seek to the next element updating any internal
     *         state (i.e. indices) as we go along)*/
    void seek_next_element() override;

    /** \brief Seek to the provided  element updating any internal
     *         state (i.e. indices) as we go along)*/
    void seek_to_element(index_type element) override;

    /** \brief
     *  Assert that the internal state of the index and the matrix
     *  pointer makes sense. If not throw an exception in debug mode.
     *
     *  It is good practice to use this before every access to the
     *  iterator data.
     *  */
    void assert_valid_state() const override;

  private:
    index_type m_index;
    SubscriptionPointer<matrix_type> m_matrix_ptr;

    /** The function to enforce the reference. If Constness
     *  then we can be sure, that Matrix::operator() returns a value
     *  so we expect this interface. Else we just pass the reference
     *  through.    */
    EnforceReference<value_type, Constness> make_ref;
};

//
// ----------------------------------------------------------------
//

//
// MatrixIteratorCoreBase
//
template <typename Matrix, bool Constness>
typename MatrixIteratorCoreBase<Matrix, Constness>::index_type
MatrixIteratorCoreBase<Matrix, Constness>::indices() const {
    return std::make_pair(row(), col());
}

// Define invalid_pos (needs to be there, otherwise linker error)
template <typename Matrix, bool Constness>
constexpr std::pair<
      typename MatrixIteratorCoreBase<Matrix, Constness>::size_type,
      typename MatrixIteratorCoreBase<Matrix, Constness>::size_type>
      MatrixIteratorCoreBase<Matrix, Constness>::invalid_pos;

template <typename Matrix, bool Constness>
constexpr bool MatrixIteratorCoreBase<Matrix, Constness>::is_const_iterator;

//
// MatrixIteratorCore
//
template <typename Matrix, bool Constness>
void MatrixIteratorDefaultCore<Matrix, Constness>::assert_valid_state() const {
    assert_dbg(m_matrix_ptr,
               ExcInvalidState("MatrixIterator does not point to any matrix"));

    assert_dbg(m_index.first != base_type::invalid_pos.first &&
                     m_index.second != base_type::invalid_pos.second,
               ExcIteratorPastEnd());
}

template <typename Matrix, bool Constness>
MatrixIteratorDefaultCore<Matrix, Constness>::MatrixIteratorDefaultCore()
      : m_index{base_type::invalid_pos},
        m_matrix_ptr{"MatrixIteratorDefaultCore"} {}

template <typename Matrix, bool Constness>
MatrixIteratorDefaultCore<Matrix, Constness>::MatrixIteratorDefaultCore(
      matrix_type& mat)
      : m_index{base_type::invalid_pos},
        m_matrix_ptr{"MatrixIteratorDefaultCore", mat} {}

template <typename Matrix, bool Constness>
MatrixIteratorDefaultCore<Matrix, Constness>::MatrixIteratorDefaultCore(
      matrix_type& mat, index_type start_index)
      : m_index{start_index}, m_matrix_ptr{"MatrixIteratorDefaultCore", mat} {}

template <typename Matrix, bool Constness>
typename MatrixIteratorDefaultCore<Matrix, Constness>::size_type
MatrixIteratorDefaultCore<Matrix, Constness>::row() const {
    return m_index.first;
}

template <typename Matrix, bool Constness>
typename MatrixIteratorDefaultCore<Matrix, Constness>::size_type
MatrixIteratorDefaultCore<Matrix, Constness>::col() const {
    return m_index.second;
}

template <typename Matrix, bool Constness>
typename std::conditional<
      Constness,
      typename MatrixIteratorDefaultCore<Matrix, Constness>::value_type,
      typename MatrixIteratorDefaultCore<Matrix, Constness>::reference>::type
MatrixIteratorDefaultCore<Matrix, Constness>::value() const {
    // Get the current element --- by value or by reference
    // and return it.
    return (*m_matrix_ptr)(row(), col());
}

template <typename Matrix, bool Constness>
typename std::conditional<
      Constness,
      const typename MatrixIteratorDefaultCore<Matrix, Constness>::pointer,
      typename MatrixIteratorDefaultCore<Matrix, Constness>::pointer>::type
MatrixIteratorDefaultCore<Matrix, Constness>::ptr_to_value() const {
    // Get the current element --- by value or by reference ---
    // and make a reference out of it using the EnforceReference
    // functor
    reference ref = make_ref((*m_matrix_ptr)(row(), col()));

    // Return the address this reference represents:
    return &ref;
}

template <typename Matrix, bool Constness>
void MatrixIteratorDefaultCore<Matrix, Constness>::seek_next_element() {
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
void MatrixIteratorDefaultCore<Matrix, Constness>::seek_to_element(
      index_type element) {
    // Assert that the indices are not too large:
    assert_greater(element.first, m_matrix_ptr->n_rows());
    assert_greater(element.second, m_matrix_ptr->n_cols());

    // Assert that we make progress in the right direction:
    assert_greater_equal(m_index.first, element.first);
    if (element.first == m_index.first) {
        assert_greater_equal(m_index.second, element.second);
    }

    // Set the indices:
    m_index = element;
}

}  // namespace linalgwrap
