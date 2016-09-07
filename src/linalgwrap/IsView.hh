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
#include "linalgwrap/view/ViewBase.hh"
#include <krims/TypeUtils.hh>

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
struct IsViewImpl : public std::false_type {};

//! Implementation of IsView, which is true if we have a view.
template <typename Matrix>
struct IsViewImpl<Matrix, krims::VoidType<typename Matrix::inner_matrix_type>>
      : public std::is_base_of<view::detail::ViewBaseMatrixContainer<
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
struct IsView : public detail::IsViewImpl<Matrix> {};

//@}

}  // linalgwrap
