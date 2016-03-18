#pragma once

#include "LazyMatrixExpression.hh"
#include "StoredMatrix_i.hh"
#include <type_traits>
#include "type_utils.hh"
#include "SubscriptionPointer.hh"

namespace linalgwrap {
namespace view {

/** \brief Class which contains and manages the reference
 *         to the Matrix we view upon.
 */
template <typename Matrix>
class ViewBaseMatrixContainer {
  public:
    typedef Matrix inner_matrix_type;

    /** A swap function for ViewMatrixContainer */
    friend void swap(ViewBaseMatrixContainer& first,
                     ViewBaseMatrixContainer& second);

    /** \brief Construct a container from a matrix reference and an identifier.
     *
     * The identifier is used in order to make the subscription with the matrix
     * and as an identification for the actual type of view.
     */
    ViewBaseMatrixContainer(inner_matrix_type& inner,
                            const std::string& identifier);

  protected:
    /** Access the inner matrix as a reference */
    inner_matrix_type& inner_matrix();

    /** Access the inner matrix as a const reference */
    const inner_matrix_type& inner_matrix() const;

    std::string identifier() const;

  private:
    //! Store a copy of the inner lazy matrix expression:
    SubscriptionPointer<inner_matrix_type> m_inner_ptr;
};

/** \brief View base class
 *
 * Distinguishing between the specialisation for stored or lazy matrices
 * is done by the second template argument, which is true_type for
 * stored matrices and else false_type
 */
template <typename Matrix, bool isStored = IsStoredMatrix<Matrix>::value>
class ViewBase;

/* \brief View base class for Views of lazy matrices
 * */
template <typename Matrix>
class ViewBase<Matrix, false>
      : public LazyMatrixExpression<
              typename std::remove_const<Matrix>::type::stored_matrix_type>,
        public ViewBaseMatrixContainer<Matrix> {
  public:
    //! Typedef of the container used to store the inner matrix
    typedef ViewBaseMatrixContainer<Matrix> container_type;

    //! Typedef of the type of matrix stored in this class
    typedef typename container_type::inner_matrix_type inner_matrix_type;

    //! Typedef of the base type
    typedef LazyMatrixExpression<typename std::remove_const<
          inner_matrix_type>::type::stored_matrix_type> base_type;

    //@{
    /** Default typedefs for standard types */
    typedef typename base_type::stored_matrix_type stored_matrix_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::lazy_matrix_expression_ptr_type
          lazy_matrix_expression_ptr_type;
    //@}

    static_assert(std::is_base_of<LazyMatrixExpression<stored_matrix_type>,
                                  Matrix>::value,
                  "Matrix is not a child of LazyMatrixExpression of the same "
                  "StoredMatrix type");

    /** Is this ViewBase for a view of a stored matrix */
    static constexpr bool view_of_stored_matrix = false;

    //
    // Constructors, destructors and assignment
    //
    /** Construct a ViewBase class for a LazyMatrix
     *
     * \param inner   Reference to the LazyMatrix
     * \param identifier An identifier for the type of view.
     */
    ViewBase(inner_matrix_type& inner, const std::string& identifier)
          : container_type(inner, identifier){};

    ViewBase(const ViewBase&) = default;
    ViewBase(ViewBase&&) = default;

    //
    // Partial implementation of LazyMatrixExpression interface.
    //
    /** \brief Update the internal data of all objects in this expression
     *         given the ParameterMap
     * */
    void update(const ParameterMap& map) override { do_update<Matrix>(map); }

    /** \brief Print the expression tree to this outstream
     * */
    void print_tree(std::ostream& o) const override {
        o << container_type::identifier() << " of ";
        container_type::inner_matrix().print_tree(o);
    }

  private:
    template <typename InnerMatrix>
    typename std::enable_if<std::is_const<InnerMatrix>::value>::type do_update(
          const ParameterMap&) {
        // Cannot call update for const matrix
        // TODO would be nice to have this error at compile time instead!
        assert_dbg(true,
                   ExcInvalidState("Update not available for const matrix"));
    }

    template <typename InnerMatrix>
    typename std::enable_if<!std::is_const<InnerMatrix>::value>::type do_update(
          const ParameterMap& map) {
        container_type::inner_matrix().update(map);
    }
};

/** \brief View base class for Views of stored matrices */
template <typename Matrix>
class ViewBase<Matrix, true>
      : public LazyMatrixExpression<typename std::remove_const<Matrix>::type>,
        public ViewBaseMatrixContainer<Matrix> {
  public:
    // TODO these views could in theory be made writable!

    //! Typedef of the container used to store the inner matrix
    typedef ViewBaseMatrixContainer<Matrix> container_type;

    //! Typedef of the type of matrix stored in this class
    typedef typename container_type::inner_matrix_type inner_matrix_type;

    //! Typedef of the base type
    typedef LazyMatrixExpression<typename std::remove_const<Matrix>::type>
          base_type;

    //@{
    /** Default typedefs for standard types */
    typedef typename base_type::stored_matrix_type stored_matrix_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::lazy_matrix_expression_ptr_type
          lazy_matrix_expression_ptr_type;
    //@}

    /** Is this a view of a stored matrix */
    static constexpr bool view_of_stored_matrix = true;

    //
    // Constructors, destructors and assignment
    //
    /** Construct a ViewBase class for a LazyMatrix
     *
     * \param inner   Reference to the LazyMatrix
     * \param identifier Identifier for the type of View calling this
     *                   constructor.
     */
    ViewBase(inner_matrix_type& inner, const std::string& identifier);

    //! Default copy constructor
    ViewBase(const ViewBase&) = default;

    //! Default move constructor
    ViewBase(ViewBase&&) = default;

    //
    // Partial implementation of LazyMatrixExpression interface.
    //
    /** \brief Update the internal data of all objects in this expression
     *         given the ParameterMap
     * */
    void update(const ParameterMap&) override;

    /** \brief Print the expression tree to this outstream
     * */
    void print_tree(std::ostream& o) const override;
};

//
// ------------------------------------------------------------------
//

//
// BaseViewMatrixContainer
//
template <typename Matrix>
void swap(ViewBaseMatrixContainer<Matrix>& first,
          ViewBaseMatrixContainer<Matrix>& second) {
    using std::swap;  // enable ADL
    swap(first.m_inner_ptr, second.m_inner_ptr);
}

template <typename Matrix>
ViewBaseMatrixContainer<Matrix>::ViewBaseMatrixContainer(
      inner_matrix_type& inner, const std::string& identifier)
      : m_inner_ptr{make_subscription(inner, identifier)} {}

template <typename Matrix>
typename ViewBaseMatrixContainer<Matrix>::inner_matrix_type&
ViewBaseMatrixContainer<Matrix>::inner_matrix() {
    return *m_inner_ptr;
}

template <typename Matrix>
const typename ViewBaseMatrixContainer<Matrix>::inner_matrix_type&
ViewBaseMatrixContainer<Matrix>::inner_matrix() const {
    return *m_inner_ptr;
}

template <typename Matrix>
std::string ViewBaseMatrixContainer<Matrix>::identifier() const {
    return m_inner_ptr.subscriber_id();
}

//
// ViewBase swap
//

template <typename Matrix, bool isStored>
void swap(ViewBase<Matrix, isStored>& first,
          ViewBase<Matrix, isStored>& second) {
    typedef typename ViewBase<Matrix, isStored>::container_type container_type;
    typedef typename ViewBase<Matrix, isStored>::base_type base_type;

    using std::swap;  // enable ADL

    swap(static_cast<container_type&>(first),
         static_cast<container_type&>(second));
    swap(static_cast<base_type&>(first), static_cast<base_type&>(second));
}
//
// ViewBase of lazy
//

// TODO

//
// ViewBase of stored
//

template <typename Matrix>
ViewBase<Matrix, true>::ViewBase(inner_matrix_type& inner,
                                 const std::string& identifier)
      : container_type(inner, identifier){};

template <typename Matrix>
void ViewBase<Matrix, true>::update(const ParameterMap&) {
    // Do nothing, since the underlying object is a stored matrix
}

template <typename Matrix>
void ViewBase<Matrix, true>::print_tree(std::ostream& o) const {
    o << container_type::identifier() << " of "
      << container_type::inner_matrix().name();
}

}  // view
}  // linalgwrap
