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
#include "linalgwrap/BaseInterfaces/Indexable_i.hh"
#include "linalgwrap/BaseInterfaces/Vector_i.hh"
#include <cmath>
#include <complex>
#include <cstdlib>

namespace linalgwrap {

namespace detail {
/** The type returned when calling the operation on one element of the indexable
 */
template <typename Indexable, typename Operation>
using ApplyElementReturnT =
      typename std::result_of<Operation(typename Indexable::scalar_type)>::type;

/** The iterator used by ApplyElementwise */
template <typename Indexable, typename Operation>
class ApplyElementwiseIterator
      : public std::iterator<std::input_iterator_tag,
                             ApplyElementReturnT<Indexable, Operation>> {
  public:
    // TODO make this class support as much as the the Indexable iterator
    // supports

    typedef const ApplyElementReturnT<Indexable, Operation> value_type;
    typedef std::iterator<std::input_iterator_tag, value_type> base_type;
    typedef typename base_type::reference reference;
    typedef typename base_type::pointer pointer;

    typedef typename Indexable::const_iterator inner_iterator_type;

    /** Construct from inner iterator and operation */
    ApplyElementwiseIterator(inner_iterator_type inner, Operation op)
          : m_inner(std::move(inner)), m_op(std::forward<Operation>(op)) {}

    /** Default-construct */
    ApplyElementwiseIterator() : m_inner{}, m_op{} {}

    /** Compare for equality */
    //@{
    bool operator==(const inner_iterator_type& other) const {
        return m_inner == other;
    }
    bool operator!=(const inner_iterator_type& other) const {
        return !(*this == other);
    }
    bool operator==(const ApplyElementwiseIterator& other) const {
        return *this == other.m_inner;
    }
    bool operator!=(const ApplyElementwiseIterator& other) const {
        return !(*this == other);
    }
    //@}

    /** Increment iterator */
    //@{
    ApplyElementwiseIterator& operator++() {
        ++m_inner;
        return *this;
    }
    ApplyElementwiseIterator operator++(int) {
        ApplyElementwiseIterator copy(*this);
        ++m_inner;
        return copy;
    }
    //@}

    /** Return value once operation has been applied */
    value_type operator*() const { return m_op(*m_inner); }

    /* Return pointer to value once operation has been applied
     *
     * \note Very expensive and hence explicitly disabled
     */
    pointer operator->() const;

  private:
    inner_iterator_type m_inner;
    Operation m_op;
#ifndef DEBUG
    // Required for operator-> in Release mode.
    ApplyElementReturnT<Indexable, Operation> m_tmp;
#endif
};

/** Class which contains an indexable and an operation */
template <typename Indexable, typename Operation>
class ApplyElementwiseContainer {
  public:
    typedef ApplyElementwiseIterator<Indexable, Operation> iterator;
    typedef ApplyElementwiseIterator<const Indexable, Operation> const_iterator;

    ApplyElementwiseContainer(const Indexable& i, Operation op_)
          : inner_ptr{"ApplyElementwiseContainer", i},
            op(std::forward<Operation>(op_)) {}

    /** \name Iterators
     */
    ///@{
    /** Return an iterator to the beginning */
    iterator begin() { return iterator{inner_ptr->begin(), op}; }

    /** Return an const iterator to the beginning */
    const_iterator begin() const { return cbegin(); }

    /** Return a const_iterator to the beginning */
    const_iterator cbegin() const {
        return const_iterator{inner_ptr->begin(), op};
    }

    /** Return a const_iterator to the end */
    const_iterator cend() const { return const_iterator{inner_ptr->end(), op}; }

    /** Return a const_iterator to the end */
    const_iterator end() const { return cend(); }

    /** Return a iterator to the end */
    iterator end() { return iterator{inner_ptr->end(), op}; }
    ///@}

  protected:
    krims::SubscriptionPointer<const Indexable> inner_ptr;
    Operation op;
};

/** Class representing an elementwise operation to an Indexable */
template <typename Indexable, typename Operation, typename = void>
class ApplyElementwise
      : public Indexable_i<ApplyElementReturnT<Indexable, Operation>>,
        public ApplyElementwiseContainer<Indexable, Operation> {
  public:
    typedef Indexable_i<ApplyElementReturnT<Indexable, Operation>> base_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef ApplyElementwiseContainer<Indexable, Operation> container_type;

    /** Construct from Indexable and operation */
    ApplyElementwise(const Indexable& i, Operation op)
          : container_type{i, std::forward<Operation>(op)} {}

    /** \brief Return the number of elements of the indexable memory */
    size_type n_elem() const override {
        return container_type::inner_ptr->n_elem();
    }

    /** \brief Get an element */
    scalar_type operator[](size_type i) const override {
        return container_type::op((*container_type::inner_ptr)[i]);
    }
};

template <typename Vector, typename Operation>
struct ApplyElementwise<Vector, Operation,
                        typename std::enable_if<IsVector<Vector>::value>::type>
      : public Vector_i<ApplyElementReturnT<Vector, Operation>>,
        public ApplyElementwiseContainer<Vector, Operation> {
  public:
    typedef Vector_i<ApplyElementReturnT<Vector, Operation>> base_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef ApplyElementwiseContainer<Vector, Operation> container_type;

    /** Construct from Indexable and operation */
    ApplyElementwise(const Vector& v, Operation op)
          : container_type{v, std::forward<Operation>(op)} {}

    /** \brief Return the number of elements of the indexable memory */
    size_type n_elem() const override {
        return container_type::inner_ptr->n_elem();
    }

    /** \brief Return the size of the vector */
    size_type size() const override { return n_elem(); }

    /** \brief Get an element */
    scalar_type operator[](size_type i) const override {
        return container_type::op((*container_type::inner_ptr)[i]);
    }

    /** \brief Get an element */
    scalar_type operator()(size_type i) const override {
        return container_type::op((*container_type::inner_ptr)(i));
    }
};

//
// ----------------------------------------------------
//

template <typename Indexable, typename Operation>
typename ApplyElementwiseIterator<Indexable, Operation>::pointer
      ApplyElementwiseIterator<Indexable, Operation>::operator->() const {
    assert_dbg(false,
               krims::ExcDisabled("Getting a pointer to an element, which "
                                  "is not in storage is expensive and "
                                  "hence operator-> is disabled."));
#ifndef DEBUG
    // In Release make it still work even though it is expensive
    m_tmp = m_op(*m_inner);
    return &m_tmp;
#else
    // In Debug, return nullptr.
    return nullptr;
#endif
}

}  // namespace detail
}  // namespace linalgwrap
