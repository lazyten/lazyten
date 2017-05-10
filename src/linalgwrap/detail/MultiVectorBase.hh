//
// Copyright (C) 2016-17 by the linalgwrap authors
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
#include "linalgwrap/Base/Interfaces.hh"
#include <krims/IteratorUtils.hh>
#include <krims/RCPWrapper.hh>
#include <krims/TypeUtils.hh>

namespace linalgwrap {
// Forward-declare
template <typename InnerVector>
class MultiVector;

namespace detail {

/** \brief Class to represent a readonly MultiVector, i.e. a MultiVector who
 *          can only provide readonly access to the underlying elements
 *
 *  \note The class will not perform implicit copies. If it gets a reference and
 *  it cannot take ownership it will only store a subscription.
 */
template <typename InnerVector>
class MultiVectorReadonly : public krims::Subscribable {
  static_assert(IsVector<InnerVector>::value,
                "The type InnerVector needs to be a vector type");

 public:
  /** \name Forwarded types */
  //@{
  typedef InnerVector vector_type;
  typedef typename vector_type::size_type size_type;
  typedef typename vector_type::scalar_type scalar_type;
  typedef typename vector_type::real_type real_type;
  //@}

  // To allow acting on it like a std container
  typedef InnerVector value_type;

  /** Reference-counted pointer type used to point to vectors*/
  typedef krims::RCPWrapper<vector_type> vector_rcptr_type;

  /** \name Iterator types */
  typedef krims::DereferenceIterator<typename std::vector<vector_rcptr_type>::iterator>
        iterator;
  typedef krims::DereferenceIterator<
        typename std::vector<vector_rcptr_type>::const_iterator>
        const_iterator;

  /** Construct an empty MultiVectorReadonly */
  MultiVectorReadonly() : m_n_elem(0) {}

  /** \name Iterators */
  ///@{
  /** Return an iterator to the first vector */
  iterator begin() { return iterator(m_vs.begin()); }

  /** Return an iterator to the first vector */
  const_iterator begin() const { return cbegin(); }

  /** Return an iterator to the first vector */
  const_iterator cbegin() const { return const_iterator(m_vs.cbegin()); }

  /** Return an iterator to the past-the end */
  iterator end() { return iterator(m_vs.end()); }

  /** Return an iterator to the past-the end */
  const_iterator end() const { return cend(); }

  /** Return an iterator to the past-the end */
  const_iterator cend() const { return const_iterator(m_vs.cend()); }
  ///@}

  /** \name Element access */
  ///@{
  /** Obtain a const reference to the ith vector */
  const vector_type& operator[](size_type i) const { return *at_ptr(i); }

  /** Obtain a reference to the ith vector */
  vector_type& operator[](size_type i) { return *at_ptr(i); }

  /** Obtain a const reference to the ith vector */
  const vector_type& at(size_type i) const { return (*this)[i]; }

  /** Obtain a reference to the ith vector */
  vector_type& at(size_type i) { return (*this)[i]; }

  /** Obtain a pointer to the ith vector */
  krims::RCPWrapper<vector_type> at_ptr(size_type i);

  /** Obtain a const pointer to the ith vector */
  krims::RCPWrapper<const vector_type> at_ptr(size_type i) const;

  /** Obtain a reference to the front vector */
  vector_type& front() {
    assert_valid_state();
    assert_dbg(!empty(), krims::ExcInvalidState("The MultiVector needs to be non-empty"));
    return *m_vs.front();
  }

  /** Obtain a const reference to the front vector */
  const vector_type& front() const {
    assert_valid_state();
    assert_dbg(!empty(), krims::ExcInvalidState("The MultiVector needs to be non-empty"));
    return *m_vs.front();
  }

  /** Obtain a reference to the back vector */
  vector_type& back() {
    assert_valid_state();
    assert_dbg(!empty(), krims::ExcInvalidState("The MultiVector needs to be non-empty"));
    return *m_vs.back();
  }

  /** Obtain a const reference to the back vector */
  const vector_type& back() const {
    assert_valid_state();
    assert_dbg(!empty(), krims::ExcInvalidState("The MultiVector needs to be non-empty"));
    return *m_vs.back();
  }

  /** Access const memory pointers of all inner vectors
   *
   * \note The vectors do not need to occupy contiguous memory!
   * \note We use SFINAE to make sure that this method is only
   *       available for vectors which have memory pointers available
   */
  template <typename Vector = InnerVector>
  std::vector<const scalar_type*> memptrs(
        krims::enable_if_cond_same_t<IsMutableMemoryVector<Vector>::value, Vector,
                                     InnerVector>* dummy = 0) const;
  ///@}

  /** \name Modifiers */
  ///@{
  /** Internally construct a vector, taking full ownership */
  template <typename... Args>
  void emplace_back(Args&&... args);

  /** Push back a vector, not taking ownership of it */
  void push_back(vector_type& val);

  /** Push back a vector, taking ownership of it */
  void push_back(vector_type&& val);

  /** Push back a reference-counted pointer */
  void push_back(vector_rcptr_type ptr);

  /** Clear content */
  void clear() {
    m_vs.clear();
    m_n_elem = 0;
  }
  //@}

  /** \name Capacity */
  ///@{
  /** Number of vectors */
  size_type n_vectors() const { return m_vs.size(); }

  /** Number of elements in each vector */
  size_type n_elem() const {
    assert_valid_state();
    return m_n_elem;
  }

  /** Change the number of vectors.
   *
   * If we are increasing the number of vectors, then n_vectors() must be
   * greater than 0.
   * This function will add vectors and optionally set their values to zero.
   *
   * \param n   The number of vectors to resize to
   * \param fill_zero  Whether to zero added vectors
   */
  template <typename SizeT, typename = krims::enable_if_cond_convertible_t<
                                  IsStoredVector<InnerVector>::value, SizeT, size_type>>
  void resize(SizeT n, bool fill_zero = true);

  /** Request a change in capacity in the underlying storage array */
  void reserve(size_type n) { m_vs.reserve(n); }

  /** Test whether this object is empty */
  bool empty() const noexcept { return m_vs.empty(); }
  ///@}

 protected:
#ifdef DEBUG
  void assert_valid_state() const {
    for (const auto& ptr : m_vs) {
      assert_size(ptr->size(), m_n_elem);
    }
  }
#else
  void assert_valid_state() const {}
#endif

  /** The number of elements in the vectors */
  size_type m_n_elem;

  /** The object storing the vector pointers */
  std::vector<vector_rcptr_type> m_vs;
  // TODO mfh: I do not like the double indirection here:
  //      For pretty much all cases I can think of *except*
  //      stored vectors the vector_type itself only stores
  //      a pointer to some actual data and is hence cheap
  //      to copy anyway, but we wrap it with yet another
  //      pointer here. Potentially we hence pay for
  //      two cache misses here.
  //
  //      Alternatively one could introduce views for
  //      stored vectors (which need to be of the same
  //      type or of a derived type from stored vector
  //      for the templating of MultiVector to work)
  //      and just have a plain std::vector<vector_type>
  //      here.
  //
  //      The problem with this is that the virtual apply
  //      method in StoredMatrix_i and LazyMatrixExpression
  //      cannot take a template argument. Right now
  //      this is ok, since we can use an abstract
  //      vector_type as the template argument to MultiVector
  //      (we never store objects), but if we make above
  //      changes this no longer works
  //      => The whole lazy matrix system would need to be
  //      rewritten without virtual calls (which is
  //      probably a good idea anyway)
};

/** Multivector base class for readable and writable scalar types
 *
 * Provides some alternatives and extensions to MultiVectorReadonly
 **/
template <typename InnerVector>
class MultiVectorReadwrite : public MultiVectorReadonly<InnerVector> {
 public:
  typedef MultiVectorReadonly<InnerVector> base_type;

  /** Access memory pointers of all inner vectors
   *
   * \note The vectors do not need to occupy contiguous memory!
   * \note We use SFINAE to make sure that this method is only
   *       available for vectors which have memory pointers available
   */
  template <typename Vector = InnerVector>
  std::vector<typename InnerVector::scalar_type*> memptrs(
        krims::enable_if_cond_same_t<IsMutableMemoryVector<Vector>::value, Vector,
                                     InnerVector>* dummy = 0);
};

/** Multivector base class moderating between MultiVectorReadwrite and
 * MultiVectorReadonly */
template <typename InnerVector, bool ConstScalar = std::is_const<InnerVector>::value ||
                                                   !IsMutableVector<InnerVector>::value>
class MultiVectorBase : public MultiVectorReadonly<InnerVector> {};

template <typename InnerVector>
class MultiVectorBase<InnerVector, false> : public MultiVectorReadwrite<InnerVector> {};

//
// ----------------------------------------------------------------------------
//

template <typename InnerVector>
krims::RCPWrapper<InnerVector> MultiVectorReadonly<InnerVector>::at_ptr(size_type i) {
  // Assert that the range is fine and we do not dereference null
  assert_valid_state();
  assert_range(0, i, n_vectors());
  assert_internal(m_vs[i] != nullptr);
  return m_vs[i];
}

template <typename InnerVector>
krims::RCPWrapper<const InnerVector> MultiVectorReadonly<InnerVector>::at_ptr(
      size_type i) const {
  // Assert that the range is fine and we do not dereference null
  assert_valid_state();
  assert_range(0, i, n_vectors());
  assert_internal(m_vs[i] != nullptr);
  return m_vs[i];
}

template <typename InnerVector>
template <typename Vector>
std::vector<const typename InnerVector::scalar_type*>
MultiVectorReadonly<InnerVector>::memptrs(
      krims::enable_if_cond_same_t<IsMutableMemoryVector<Vector>::value, Vector,
                                   InnerVector>*) const {
  // Note: We need the dummy template argument Vector, such that one template
  // argument is only resolved once the method is called (and not on class
  // construction). Otherwise SFINAE does not work.
  // This is also why we need the dummy first argument

  std::vector<const scalar_type*> res;
  res.reserve(m_vs.size());
  for (auto& vptr : m_vs) {
    res.push_back(vptr->memptr());
  }
  return res;
}

template <typename InnerVector>
template <typename... Args>
void MultiVectorReadonly<InnerVector>::emplace_back(Args&&... args) {
  auto ptr = std::make_shared<vector_type>(std::forward<Args>(args)...);
  if (m_vs.size() == 0) {
    // Initial m_n_elem size cache
    m_n_elem = ptr->size();
  }
  assert_size(ptr->size(), m_n_elem);
  m_vs.push_back(std::move(vector_rcptr_type{ptr}));
  assert_valid_state();
}

template <typename InnerVector>
void MultiVectorReadonly<InnerVector>::push_back(vector_type& val) {
  if (m_vs.size() == 0) {
    // Initial m_n_elem size cache
    m_n_elem = val.size();
  }
  assert_size(val.size(), m_n_elem);

  auto ptr = krims::make_subscription(val, "MultiVector");
  m_vs.push_back(std::move(vector_rcptr_type{ptr}));
  assert_valid_state();
}

template <typename InnerVector>
void MultiVectorReadonly<InnerVector>::push_back(vector_type&& val) {
  if (m_vs.size() == 0) {
    // Initial m_n_elem size cache
    m_n_elem = val.size();
  }
  assert_size(val.size(), m_n_elem);

  auto ptr = std::make_shared<vector_type>(std::forward<vector_type>(val));
  m_vs.push_back(std::move(vector_rcptr_type{ptr}));
  assert_valid_state();
}

template <typename InnerVector>
void MultiVectorReadonly<InnerVector>::push_back(vector_rcptr_type ptr) {
  assert_dbg(ptr != nullptr, krims::ExcInvalidPointer());

  if (m_vs.size() == 0) {
    // Initial m_n_elem size cache
    m_n_elem = ptr->size();
  }
  assert_size(ptr->size(), m_n_elem);
  m_vs.push_back(std::move(ptr));
  assert_valid_state();
}

template <typename InnerVector>
template <typename SizeT, typename>
void MultiVectorReadonly<InnerVector>::resize(SizeT nin, bool fill_zero) {
  assert_valid_state();

  auto n = static_cast<size_type>(nin);
  if (n == 0) {
    clear();
    return;
  } else if (n <= n_vectors()) {
    m_vs.resize(n);
    return;
  } else {
    reserve(n);
    for (size_type i = n_vectors(); i < n; ++i) {
      auto ptr = std::make_shared<InnerVector>(m_n_elem, fill_zero);
      m_vs.push_back(vector_rcptr_type{ptr});
    }
  }
}

template <typename InnerVector>
template <typename Vector>
std::vector<typename InnerVector::scalar_type*>
MultiVectorReadwrite<InnerVector>::memptrs(
      krims::enable_if_cond_same_t<IsMutableMemoryVector<Vector>::value, Vector,
                                   InnerVector>*) {
  // Note: We need the dummy template argument Vector, such that one template
  // argument is only resolved once the method is called (and not on class
  // construction). Otherwise SFINAE does not work.
  // This is also why we need the dummy first argument

  std::vector<typename InnerVector::scalar_type*> res;
  res.reserve(base_type::m_vs.size());
  for (auto& vptr : base_type::m_vs) {
    res.push_back(vptr->memptr());
  }
  return res;
}

}  // namespace detail
}  // namespace linalgwrap
