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

#include "linalgwrap/IsView.hh"
#include "linalgwrap/LazyMatrixExpression.hh"
#include "linalgwrap/StoredMatrix_i.hh"
#include <krims/SubscriptionPointer.hh>
#include <type_traits>

namespace linalgwrap {

// Forward declaration
template <typename Matrix>
struct IsView;

namespace view {
namespace detail {

// Further forward declaration
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
    template <typename M>
    friend void swap(ViewBaseMatrixContainer<M>& first,
                     ViewBaseMatrixContainer<M>& second);

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

    /** The identifier passed upon object construction */
    std::string identifier() const { return m_inner_ptr.subscriber_id(); }

  private:
    //! Store a copy of the inner lazy matrix expression:
    krims::SubscriptionPointer<inner_matrix_type> m_inner_ptr;
};

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
struct ViewInnermostMatrixType<
      Matrix, krims::VoidType<typename Matrix::inner_matrix_type>> {
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
    void update(const krims::ParameterMap& map) override {
        do_update<Matrix>(map);
    }

  private:
    template <typename InnerMatrix>
    typename std::enable_if<std::is_const<InnerMatrix>::value>::type do_update(
          const krims::ParameterMap&) {
        // Cannot call update for const matrix
        // TODO would be nice to have this error at compile time instead!
        assert_dbg(true, krims::ExcInvalidState(
                               "Update not available for const matrix"));
    }

    template <typename InnerMatrix>
    typename std::enable_if<!std::is_const<InnerMatrix>::value>::type do_update(
          const krims::ParameterMap& map) {
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
    ViewBase(inner_matrix_type& inner, const std::string& identifier)
          : container_type(inner, identifier) {}

    /** \brief Update the internal data of all objects in this expression
     *         given the ParameterMap
     * */
    void update(const krims::ParameterMap&) override {
        // Do nothing, since the underlying object is a stored matrix
    }
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

}  // detail
}  // view
}  // linalgwrap
