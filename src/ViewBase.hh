#pragma once

#include "LazyMatrixExpression.hh"
#include "StoredMatrix_i.hh"
#include "SubscriptionPointer.hh"
#include "type_utils.hh"
#include <type_traits>

namespace linalgwrap {
namespace view {

// Forward declaration
template <typename Matrix, typename = void>
struct ViewInnermostMatrixType;

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

    /** \brief Is the innermost matrix contained in the view structure a
     *         stored matrix?
     *
     * \note This property is evaluated in a recursive fashion, i.e. a
     *       view of a view of a stored matrix has this flag evaluating
     *       to true.
     * */
    static constexpr bool is_stored_matrix_view =
          IsStoredMatrix<typename ViewInnermostMatrixType<Matrix>::type>::value;

    /** \brief Does this view contain a matrix at innermost level,
     *         which we can write to?
     *
     * The idea is that if the innermost matrix is a stored matrix and
     * is non-const we can undo the view operations and set elements on
     * this matrix.
     * */
    static constexpr bool is_const_view =
          is_stored_matrix_view && std::is_const<Matrix>::value;

    /** Constant access to the inner matrix of the view */
    const inner_matrix_type& inner_matrix() const;

  protected:
    /** Access the inner matrix as a reference */
    inner_matrix_type& inner_matrix();

    std::string identifier() const;

  private:
    //! Store a copy of the inner lazy matrix expression:
    SubscriptionPointer<inner_matrix_type> m_inner_ptr;
};

//@{
/** \brief struct representing a type (std::true_type, std::false_type) which
 *  indicates whether T is a stored matrix
 *
 * The definition is done using SFINAE, such that even for types not having a
 * typedef inner_matrix_type this expression is valid.
 *  */
template <typename Matrix, typename = void>
struct IsView : std::false_type {};

template <typename Matrix>
struct IsView<Matrix, void_t<typename Matrix::inner_matrix_type>>
      : std::is_base_of<
              ViewBaseMatrixContainer<typename Matrix::inner_matrix_type>,
              Matrix> {};
//@}

/** \brief Struct allowing access to the type of the innermost matrix
 *         we view upon.
 *
 * Achieved by recursively checking wheather the current matrix has the
 * base type ViewBaseMatrixContainer. If this is the case we go one
 * level deeper into Matrix::inner_matrix_type, otherwise we return the
 * current type.
 */
//@{
template <typename Matrix, typename>
struct ViewInnermostMatrixType {
    // default implementation used if inner_matrix_type is not available
    // then we can be sure, that we do not have a view here ...
    typedef Matrix type;
};

template <typename Matrix>
struct ViewInnermostMatrixType<Matrix,
                               void_t<typename Matrix::inner_matrix_type>> {
    typedef typename std::conditional<IsView<Matrix>::value,
                                      typename Matrix::inner_matrix_type,
                                      Matrix>::type type;
};
//@}

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
          inner_matrix_type>::type::stored_matrix_type>
          base_type;

    static_assert(
          std::is_base_of<
                LazyMatrixExpression<typename base_type::stored_matrix_type>,
                Matrix>::value,
          "Matrix is not a child of LazyMatrixExpression of the same "
          "StoredMatrix type");

    /** Construct a ViewBase class for a LazyMatrix
     *
     * \param inner   Reference to the LazyMatrix
     * \param identifier An identifier for the type of view.
     */
    ViewBase(inner_matrix_type& inner, const std::string& identifier)
          : container_type(inner, identifier){};

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
    //      Use the is_const_view flag of ViewBaseMatrixContainer to check
    //      this before enabling the writable [] and () operators in
    //      deriving classes.

    //! Typedef of the container used to store the inner matrix
    typedef ViewBaseMatrixContainer<Matrix> container_type;

    //! Typedef of the type of matrix stored in this class
    typedef typename container_type::inner_matrix_type inner_matrix_type;

    /** Construct a ViewBase class for a LazyMatrix
     *
     * \param inner   Reference to the LazyMatrix
     * \param identifier Identifier for the type of View calling this
     *                   constructor.
     */
    ViewBase(inner_matrix_type& inner, const std::string& identifier);

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
inline void swap(ViewBaseMatrixContainer<Matrix>& first,
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
