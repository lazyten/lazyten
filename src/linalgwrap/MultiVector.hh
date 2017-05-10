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
#include "detail/MultiVectorBase.hh"
#include <initializer_list>
#include <krims/Range.hh>

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
class MultiVector : public detail::MultiVectorBase<InnerVector>,
                    // Multivectors are cheaply copyable since this just makes a
                    // shallow copy
                    public krims::CheaplyCopyable_i {
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
  template <typename SizeType, krims::enable_if_cond_convertible_t<
                                     IsStoredVector<InnerVector>::value, SizeType,
                                     typename InnerVector::size_type>...>
  MultiVector(size_type n_elem, SizeType n_vectors, bool fill_zero = true);

  /** \brief Construct from a nested initialiser list of scalars.
   *
   * The outermost layer gives the number of vectors, the innermost layer the
   * number of elements in each vector. An example would be
   * ```
   * MultiVector<vector_type> mat{{1.0,2.0,0.5},{1.5,4.5,6.}})
   * ```
   * which produces a MultiVector with 3 Vectors a 2 elements.
   */
  template <typename Scalar, krims::enable_if_cond_convertible_t<
                                   IsStoredVector<InnerVector>::value, Scalar,
                                   typename InnerVector::scalar_type>...>
  MultiVector(std::initializer_list<std::initializer_list<Scalar>> list_of_lists);

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

  // TODO
  // TODO make some of the constructors below explicit
  // TODO

  /** conversion from different vector type
   *
   * All objects are converted via a shallow copy, so it is
   * linear in the number of vectors.
   * */
  template <typename OtherVector, typename = krims::enable_if_t<std::is_convertible<
                                        OtherVector*, InnerVector*>::value>>
  MultiVector(MultiVector<OtherVector>& omv);

  /** conversion from different vector type -- const version
   *
   * All objects are converted via a shallow copy, so it is
   * linear in the number of vectors.
   * */
  template <typename OtherVector,
            typename = krims::enable_if_t<
                  std::is_convertible<OtherVector*, InnerVector*>::value &&
                  std::is_const<InnerVector>::value>>
  MultiVector(const MultiVector<OtherVector>& omv);

  /** Implicit conversion from constant vector */
  // operator MultiVector<const InnerVector>();

  /** Implicit conversion from vector via move */
  template <typename OtherVector, typename = krims::enable_if_t<std::is_convertible<
                                        OtherVector*, InnerVector*>::value>>
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

  // TODO also have copy_shallow

  /** Obtain a view (shallow copy) of a part of the columns */
  MultiVector subview(const krims::Range<size_type>& colrange);

  /** Obtain a constant view of a part of the columns */
  MultiVector<const InnerVector> subview(const krims::Range<size_type>& colrange) const {
    return csubview(colrange);
  }

  /** Obtain a constant view of a part of the columns */
  MultiVector<const InnerVector> csubview(const krims::Range<size_type>& colrange) const;
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
template <typename Vector, typename = krims::enable_if_t<IsVector<
                                 typename std::remove_reference<Vector>::type>::value>>
MultiVector<typename std::remove_reference<Vector>::type> as_multivector(Vector&& v) {
  typedef typename std::remove_reference<Vector>::type vtype;
  return MultiVector<vtype>(std::forward<Vector>(v));
}
//@}

/** Construct a multivector which contains one vector,
 * which is constructed from the arguments which are passed */
template <typename Vector, typename... Args>
MultiVector<Vector> make_as_multivector(Args&&... args) {
  return as_multivector(Vector(std::forward<Args>(args)...));
}

/** Compute multivector norms */
///@{
/** Calculate the l1 norm of each vector of the multivector
 *  (sum of abs values of elements) */
template <typename Vector>
std::vector<typename Vector::real_type> norm_l1(const MultiVector<Vector>& mv) {
  std::vector<typename Vector::real_type> res(mv.n_vectors());
  std::transform(std::begin(mv), std::end(mv), std::begin(res),
                 [](const Vector& v) { return norm_l1(v); });
  return res;
}

/** Calculate the linf norm of the vectors of the multivector (abs. largest
 * element) */
template <typename Vector>
std::vector<typename Vector::real_type> norm_linf(const MultiVector<Vector>& mv) {
  std::vector<typename Vector::real_type> res(mv.n_vectors());
  std::transform(std::begin(mv), std::end(mv), std::begin(res),
                 [](const Vector& v) { return norm_linf(v); });
  return res;
}

/** Calculate the l2 norm squared of the vectors of the multivector. */
template <typename Vector>
std::vector<typename Vector::real_type> norm_l2_squared(const MultiVector<Vector>& mv) {
  std::vector<typename Vector::real_type> res(mv.n_vectors());
  std::transform(std::begin(mv), std::end(mv), std::begin(res),
                 [](const Vector& v) { return norm_l2_squared(v); });
  return res;
}

/** Calculate the l2 norm of the vectors of the multivector. */
template <typename Vector>
std::vector<typename Vector::real_type> norm_l2(const MultiVector<Vector>& mv) {
  std::vector<typename Vector::real_type> res(mv.n_vectors());
  std::transform(std::begin(mv), std::end(mv), std::begin(res),
                 [](const Vector& v) { return norm_l2(v); });
  return res;
}

template <typename Vector, typename Vector2,
          typename = krims::enable_if_t<
                std::is_same<typename std::remove_const<Vector>::type,
                             typename std::remove_const<Vector2>::type>::value>>
typename Vector::type_family::template matrix<typename Vector::scalar_type> dot(
      const MultiVector<Vector>& u, const linalgwrap::MultiVector<Vector2>& v) {
  typedef typename Vector::type_family::template matrix<typename Vector::scalar_type>
        matrix_type;
  if (u.n_vectors() == 0 || v.n_vectors() == 0) {
    return matrix_type(u.n_vectors(), v.n_vectors());
  }
  assert_size(u.n_elem(), v.n_elem());

  matrix_type ret(u.n_vectors(), v.n_vectors(), false);
  for (size_t ui = 0; ui < u.n_vectors(); ++ui) {
    for (size_t vi = 0; vi < v.n_vectors(); ++vi) {
      ret(ui, vi) = dot(v[vi], u[ui]);
    }
  }
  return ret;
}

// TODO later return a lazy expression here
/** Compute the sum of all outer products between the vectors,
 *  i.e. computes
 *  \f[ \sum_k  u^{(k)}_i v^{(k)}_j \f]
 * where k runs over the number of vectors in the multivector
 * and i and j runs over the respective number of elements
 * in the multivectors u and v
 */
template <typename Vector, typename Vector2,
          typename = krims::enable_if_t<
                std::is_same<typename std::remove_const<Vector>::type,
                             typename std::remove_const<Vector2>::type>::value>>
typename Vector::type_family::template matrix<typename Vector::scalar_type>
outer_prod_sum(const MultiVector<Vector>& u, const linalgwrap::MultiVector<Vector2>& v) {
  // TODO this is a temporary routine
  assert_size(u.n_vectors(), v.n_vectors());

  typedef typename Vector::type_family::template matrix<typename Vector::scalar_type>
        matrix_type;
  typedef typename Vector::size_type size_type;
  matrix_type ret(u.n_elem(), v.n_elem());

  // Later use std::inner_product and as operator* the outer_prod
  // function

  for (size_type vi = 0; vi < v.n_vectors(); ++vi) {
    const auto& uu = u[vi];
    const auto& vv = v[vi];

    for (size_type i = 0; i < u.n_elem(); ++i) {
      for (size_type j = 0; j < v.n_elem(); ++j) {
        ret(i, j) += uu(i) * vv(j);
      }  // j
    }    // i
  }      // vi

  return ret;
}

///@}

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
template <typename SizeType, krims::enable_if_cond_convertible_t<
                                   IsStoredVector<InnerVector>::value, SizeType,
                                   typename InnerVector::size_type>...>
MultiVector<InnerVector>::MultiVector(size_type n_elem, SizeType n_vectors,
                                      bool fill_zero)
      : MultiVector() {
  base_type::m_n_elem = n_elem;
  base_type::resize(n_vectors, fill_zero);
  base_type::assert_valid_state();
}

template <typename InnerVector>
template <typename Scalar,
          krims::enable_if_cond_convertible_t<IsStoredVector<InnerVector>::value, Scalar,
                                              typename InnerVector::scalar_type>...>
MultiVector<InnerVector>::MultiVector(
      std::initializer_list<std::initializer_list<Scalar>> list_of_lists)
      : MultiVector(list_of_lists.size(),
                    list_of_lists.size() > 0 ? list_of_lists.begin()->size() : 0, false) {
#ifdef DEBUG
  size_type n_elem = list_of_lists.size();
  size_type n_vectors = n_elem > 0 ? list_of_lists.begin()->size() : 0;
#endif

  // Assert all columns have equal length.
  assert_element_sizes(list_of_lists, n_vectors);

  size_type i = 0;
  for (auto row : list_of_lists) {
    size_type j = 0;
    for (scalar_type elem : row) {
      ((*this)[j])(i) = elem;
      ++j;
    }
    assert_internal(j == n_vectors);
    ++i;
  }
  assert_internal(i == n_elem);
}

template <typename InnerVector>
template <typename OtherVector, typename>
MultiVector<InnerVector>::MultiVector(const MultiVector<OtherVector>& omv) {
  base_type::reserve(omv.n_vectors());
  for (size_type i = 0; i < omv.n_vectors(); ++i) {
    base_type::push_back(omv.at_ptr(i));
  }
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
      const krims::Range<size_type>& col_range) {
  base_type::assert_valid_state();
  if (!col_range.empty()) {
    assert_greater_equal(col_range.upper_bound(), base_type::n_vectors());
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
      const krims::Range<size_type>& col_range) const {
  base_type::assert_valid_state();
  if (!col_range.empty()) {
    assert_greater_equal(col_range.upper_bound(), base_type::n_vectors());
  }

  MultiVector<const InnerVector> res;
  res.reserve(col_range.length());
  for (const auto& i : col_range) {
    res.push_back(base_type::m_vs[i]);
  }
  return res;
}

}  // namespace linalgwrap
