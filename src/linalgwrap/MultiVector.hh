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
#include "detail/MultiVectorBase.hh"
#include "linalgwrap/Range.hh"

namespace linalgwrap {

/** \brief Class to represent a collection of multiple vectors
 *
 * The idea is to be able to pass a logical collection of vectors around
 * and to apply the same operations to them in turn.
 *
 * All inner objects are assumed to have the same length.
 *
 *  \note The class will not perform implicit copies. If it gets a reference and
 *  it cannot take ownership it will only store a subscription.
 *
 *  \note The copy-constructor and copy-assignment automatically make only
 *  *shallow* copies, i.e. they still reference the same underlying objects.
 *  If you want to make a deep copy, use ``copy_deep`` instead.
 */
template <typename InnerVector>
class MultiVector : public detail::MultiVectorBase<InnerVector> {
  public:
    typedef detail::MultiVectorBase<InnerVector> base_type;
    typedef typename base_type::vector_type vector_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::real_type real_type;

    /** \name Constructors, destructors and assignment */
    ///@{
    /** Construct an empty Multivector */
    MultiVector() : base_type() {}

    /** Construct from one inner vector, taking no ownership of the vector */
    explicit MultiVector(vector_type& v);

    /** Construct from one inner vector  taking ownership of the vector */
    explicit MultiVector(vector_type&& v);

    /** Construct a MultiVector with n_vectors vectors, which each have n_elem
     * elements */
    template <typename Boolean,
              krims::enable_if_cond_same_t<IsStoredVector<InnerVector>::value,
                                           Boolean, bool>...>
    MultiVector(size_type n_elem, size_type n_vectors,
                Boolean fill_zero = true);

    /** Default destructor */
    ~MultiVector() = default;

    /** Default move constructor */
    MultiVector(MultiVector&&) = default;

    /** Default move assignment */
    MultiVector& operator=(MultiVector&&) = default;

    /** Copy assignment making a shallow copy */
    MultiVector& operator=(const MultiVector&) = default;

    /** Copy constructor making a shallow copy */
    MultiVector(const MultiVector&) = default;

    /** conversion from different vector type
     *
     * All objects are converted via a shallow copy, so it is
     * linear in the number of vectors.
     * */
    template <typename OtherVector,
              typename = krims::enable_if_ptr_convertible_t<OtherVector,
                                                            InnerVector>>
    MultiVector(MultiVector<OtherVector>& omv);

    /** Implicit conversion from constant vector */
    // operator MultiVector<const InnerVector>();

    /** Implicit conversion from vector via move */
    template <typename OtherVector,
              typename = krims::enable_if_ptr_convertible_t<OtherVector,
                                                            InnerVector>>
    MultiVector(MultiVector<OtherVector>&& omv);
    ///@}

    /** \name MultiVector operations
     *
     * All these operations will be done on all vectors at once
     * */
    ///@{
    // TODO code and extend this section
    ///@}

    /** \name Copies and views */
    ///@{
    /** Make a deep copy */
    MultiVector copy_deep() const;

    /** Obtain a view (shallow copy) of a part of the columns */
    MultiVector subview(const Range<size_type>& colrange);

    /** Obtain a constant view of a part of the columns */
    MultiVector<const InnerVector> subview(
          const Range<size_type>& colrange) const {
        return csubview(colrange);
    }

    /** Obtain a constant view of a part of the columns */
    MultiVector<const InnerVector> csubview(
          const Range<size_type>& colrange) const;
    ///@}
};

/** \brief Simple output operator, that plainly shows all
 *   vectors.
 *
 *   Vectors are shown in a row and separated by a
 *   newline.
 *
 *   The last row is not terminated by a newline character.
 *  */
template <typename Vector>
std::ostream& operator<<(std::ostream& o, const MultiVector<Vector>& mv) {
    for (size_t i = 0; i < mv.n_vectors(); ++i) {
        o << mv[i];
        if (i + 1 < mv.n_vectors()) {
            o << std::endl;
        }
    }

    // TODO extend
    // assert_dbg(false, krims::ExcNotImplemented());
    // io::MatrixPrinter().print(m, o);
    return o;
}

//@{
/** Helper function to make a MultiVector from a single vector.
 *
 * The multivector will only contain this single vector. Both giving
 * the multivector ownership of the vector (by moving it inside) or
 * not (by just passing a reference) is possible.
 **/
template <typename Vector>
MultiVector<typename std::remove_reference<Vector>::type> as_multivector(
      Vector&& v) {
    // If we gen an rvalue
    typedef typename std::remove_reference<Vector>::type type;
    return MultiVector<type>(std::forward<type>(v));
}

template <typename Vector>
MultiVector<Vector> as_multivector(Vector& v) {
    // If we get an lvalue
    return MultiVector<Vector>(v);
}
//@}

/** Construct a multivector which contains one vector,
 * which is constructed from the arguments which are passed */
template <typename Vector, typename... Args>
MultiVector<Vector> make_as_multivector(Args&&... args) {
    return as_multivector(Vector(std::forward<Args>(args)...));
}

//
// -------------------------------------------------
//

template <typename InnerVector>
MultiVector<InnerVector>::MultiVector(vector_type& v) : MultiVector() {
    base_type::reserve(1);
    base_type::push_back(v);
}

template <typename InnerVector>
MultiVector<InnerVector>::MultiVector(vector_type&& v) : MultiVector() {
    base_type::reserve(1);
    base_type::push_back(std::move(v));
}

template <typename InnerVector>
template <typename Boolean,
          krims::enable_if_cond_same_t<IsStoredVector<InnerVector>::value,
                                       Boolean, bool>...>
MultiVector<InnerVector>::MultiVector(size_type n_elem, size_type n_vectors,
                                      Boolean fill_zero)
      : MultiVector() {
    base_type::m_n_elem = n_elem;
    base_type::resize(n_vectors, fill_zero);
    base_type::assert_valid_state();
}

template <typename InnerVector>
template <typename OtherVector, typename>
MultiVector<InnerVector>::MultiVector(MultiVector<OtherVector>& omv) {
    base_type::reserve(omv.n_vectors());
    for (size_type i = 0; i < omv.n_vectors(); ++i) {
        base_type::push_back(omv.at_ptr(i));
    }
}

template <typename InnerVector>
template <typename OtherVector, typename>
MultiVector<InnerVector>::MultiVector(MultiVector<OtherVector>&& omv) {
    base_type::reserve(omv.n_vectors());
    for (size_type i = 0; i < omv.n_vectors(); ++i) {
        base_type::push_back(std::move(omv.at_ptr(i)));
    }
}

template <typename InnerVector>
MultiVector<InnerVector> MultiVector<InnerVector>::copy_deep() const {
    base_type::assert_valid_state();

    MultiVector res;
    res.reserve(base_type::m_vs.size());
    for (const auto& vptr : base_type::m_vs) {
        res.push_back(std::move(vector_type{*vptr}));
    }
    return res;
}

template <typename InnerVector>
MultiVector<InnerVector> MultiVector<InnerVector>::subview(
      const Range<size_type>& col_range) {
    base_type::assert_valid_state();
    if (!col_range.empty()) {
        assert_greater_equal(col_range.last(), base_type::n_vectors());
    }

    MultiVector res;
    res.reserve(col_range.length());
    for (const auto& i : col_range) {
        res.push_back(base_type::m_vs[i]);
    }
    return res;
}

template <typename InnerVector>
MultiVector<const InnerVector> MultiVector<InnerVector>::csubview(
      const Range<size_type>& col_range) const {
    base_type::assert_valid_state();
    if (!col_range.empty()) {
        assert_greater_equal(col_range.last(), base_type::n_vectors());
    }

    MultiVector<const InnerVector> res;
    res.reserve(col_range.length());
    for (const auto& i : col_range) {
        res.push_back(base_type::m_vs[i]);
    }
    return res;
}

}  // namespace linalgwrap
