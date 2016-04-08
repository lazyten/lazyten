#pragma once
#include "linalgwrap/type_utils.hh"
#include "linalgwrap/view/ViewBase.hh"

namespace linalgwrap {
namespace view {
namespace detail {

// Forward declaration
template <typename Matrix>
class ViewBaseMatrixContainer;

}  // namespace detail
}  // namespace view

namespace detail {

//! Default implementation of IsView (== false)
template <typename Matrix, typename = void>
struct IsViewImpl : std::false_type {};

//! Implementation of IsView, which is true if we have a view.
template <typename Matrix>
struct IsViewImpl<Matrix, void_t<typename Matrix::inner_matrix_type>>
      : std::is_base_of<view::detail::ViewBaseMatrixContainer<
                              typename Matrix::inner_matrix_type>,
                        Matrix> {};
}

//@{
/** \brief struct representing a type (std::true_type, std::false_type) which
 *  indicates whether T is a stored matrix
 *
 * The definition is done using SFINAE, such that even for types not having a
 * typedef inner_matrix_type this expression is valid.
 *  */
template <typename Matrix>
struct IsView : detail::IsViewImpl<Matrix> {};

//@}

}  // linalgwrap
