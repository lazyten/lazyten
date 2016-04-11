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
#include "linalgwrap/MatrixIteratorCore.hh"

namespace linalgwrap {
/** \brief The matrix iterator enwrapping a MatrixIteratorCore
 *
 * The resulting object has the required interface for an
 * STL-compatible iterator and traverses the elements of the
 * matrix in the way the core dictates.
 * */
template <typename IteratorCore>
class MatrixIterator : private IteratorCore {
  public:
    typedef IteratorCore base_type;
    typedef typename base_type::original_matrix_type original_matrix_type;
    typedef typename base_type::matrix_type matrix_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::index_type index_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::reference reference;
    typedef typename base_type::pointer pointer;

    static_assert(
          std::is_base_of<
                detail::MatrixIteratorCoreBase<original_matrix_type,
                                               IteratorCore::is_const_iterator>,
                IteratorCore>::value,
          "The IteratorCore should be a subclass of MatrixIteratorCoreBase");

    //
    // Constructor, destructor and assignment
    //
    /** Default constructor */
    MatrixIterator();

    /** \brief Constructor of a Matrix iterator pointing to
     *  the past-the-end position.
     *
     *  In other words this constructor constructs an iterator
     *  in the state past-the-end. */
    MatrixIterator(matrix_type& mat);

    /** \brief Construct an iterator giving the initial value
     * it should point to. */
    MatrixIterator(matrix_type& mat, index_type start_index);

    //
    // Iterator information
    //
    /** Get row of the currently pointed to element
     *  or invalid_pos.first if this is an iterator-past-the-end*/
    size_type row() const;

    /** Get column of currently pointed to element or invalid_pos.second
     * if this is an iterator-past-the-end.*/
    size_type col() const;

    /** Return the tuple of indices of the currently pointed to element
     * or invalid_pos if this is an iterator-past-the-end.*/
    index_type indices() const;

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
    /** Return the value of the element we point to. */
    auto operator*() const -> decltype(this->value());

    /** Access the members of the element we point to. */
    auto operator-> () const -> decltype(this->ptr_to_value());

    //
    // Comparison
    //
    /** \brief check if two iterators are equal
     *
     * Matrix iterators are equal if they point to the same element
     */
    bool operator==(const MatrixIterator& other) const;

    /** \brief Check whether two iterators are unequal
     *
     * Matrix iterators are unequal if they point to different elements
     */
    bool operator!=(const MatrixIterator& other) const;
};

//
// ----------------------------------------------------------------
//

template <typename IteratorCore>
MatrixIterator<IteratorCore>::MatrixIterator() : base_type{} {}

template <typename IteratorCore>
MatrixIterator<IteratorCore>::MatrixIterator(matrix_type& mat)
      : base_type{mat} {}

template <typename IteratorCore>
MatrixIterator<IteratorCore>::MatrixIterator(matrix_type& mat,
                                             index_type start_index)
      : base_type{mat, start_index} {
    base_type::assert_valid_state();
}

template <typename IteratorCore>
typename MatrixIterator<IteratorCore>::size_type
MatrixIterator<IteratorCore>::row() const {
    return base_type::row();
}

template <typename IteratorCore>
typename MatrixIterator<IteratorCore>::size_type
MatrixIterator<IteratorCore>::col() const {
    return base_type::col();
}

template <typename IteratorCore>
typename MatrixIterator<IteratorCore>::index_type
MatrixIterator<IteratorCore>::indices() const {
    return base_type::indices();
}

template <typename IteratorCore>
MatrixIterator<IteratorCore>& MatrixIterator<IteratorCore>::seek_to(
      index_type element) {
    base_type::assert_valid_state();
    base_type::seek_to_element(element);
    return *this;
}

template <typename IteratorCore>
MatrixIterator<IteratorCore>& MatrixIterator<IteratorCore>::operator++() {
    base_type::assert_valid_state();
    base_type::seek_next_element();
    return *this;
}

template <typename IteratorCore>
MatrixIterator<IteratorCore> MatrixIterator<IteratorCore>::operator++(int) {
    MatrixIterator copy{*this};
    ++(*this);
    return copy;
}

template <typename IteratorCore>
auto MatrixIterator<IteratorCore>::operator*() const
      -> decltype(this->value()) {
    base_type::assert_valid_state();
    return base_type::value();
}

template <typename IteratorCore>
auto MatrixIterator<IteratorCore>::operator-> () const
      -> decltype(this->ptr_to_value()) {
    base_type::assert_valid_state();
    return base_type::ptr_to_value();
}

template <typename IteratorCore>
bool MatrixIterator<IteratorCore>::operator==(
      const MatrixIterator& other) const {
    return other.indices() == indices();
}

template <typename IteratorCore>
bool MatrixIterator<IteratorCore>::operator!=(
      const MatrixIterator& other) const {
    return !((*this) == other);
}

}  // linalgwrap
