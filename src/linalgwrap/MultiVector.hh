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
#include "BaseInterfaces.hh"

namespace linalgwrap {
/** \brief Class to represent a collection of multiple vectors
 *
 * The idea is to be able to pass such a logical collection around
 * (in for example ParameterMaps) without doing this for the individual
 * objects by themselves.
 *
 * All inner objects are assumed to have the same length
 */
template <typename InnerVector>
class MultiVector : public krims::Subscribable,
                    private std::vector<InnerVector> {
  public:
    typedef InnerVector vector_type;
    typedef typename vector_type::size_type size_type;
    typedef std::vector<vector_type> base_type;

    /** Construct an empty MultiVector */
    MultiVector() {}

    /** Construct from one inner vector */
    explicit MultiVector(InnerVector i) { base_type::push_back(std::move(i)); }

    /** Construct from a vector of inner vectors */
    explicit MultiVector(std::vector<InnerVector> v)
          : base_type(std::forward<std::vector<InnerVector>>(v)) {
        assert_valid_state();
    }

    /** Construct from an iterator range of vectors */
    template <typename Iterator>
    MultiVector(Iterator begin, Iterator end) : base_type(begin, end) {
        assert_valid_state();
    }

    /** Construct by allocating n_vectors inner vectors **/
    explicit MultiVector(size_type n_vectors, const InnerVector& i)
          : base_type(n_vectors, i) {
        assert_valid_state();
    }

    /** Construct by allocating n_vectors inner vectors
     * with each n_elem elements. Optionally fill the inner vectors
     * with zeros.
     *
     * \param fill_zero  If true the inner vectors are zeroed.
     *
     * \note This constructor is only enabled for stored vectors */
    MultiVector(size_type n_vectors, size_type n_elem,
                typename std::enable_if<IsStoredVector<InnerVector>::value,
                                        bool>::type fill_zero = true)
          : MultiVector(n_vectors, InnerVector(n_elem, fill_zero)) {
        assert_valid_state();
    }

    /** Number of vectors */
    size_type n_vectors() { return base_type::size(); }

    /** Number of elements in each vector */
    size_type n_elem() { return n_vectors() > 0 ? (*this)[0].size() : 0; }

    // Forwarded things from std::vector
    // TODO assert state in these as well!
    using base_type::front;
    using base_type::back;
    using base_type::operator[];

    using base_type::begin;
    using base_type::cbegin;
    using base_type::end;
    using base_type::cend;
    using base_type::rbegin;
    using base_type::crbegin;
    using base_type::rend;
    using base_type::crend;

    using base_type::empty;
    using base_type::emplace_back;
    using base_type::push_back;
    using base_type::pop_back;
    using base_type::clear;

    using base_type::reserve;

  private:
#ifdef DEBUG
    void assert_valid_state() { assert_element_sizes(*this, n_elem()); }
#else
    void assert_valid_state() {}
#endif
};
}  // namespace linalgwrap
