//
// Copyright (C) 2017 by the lazyten authors
//
// This file is part of lazyten.
//
// lazyten is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// lazyten is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with lazyten. If not, see <http://www.gnu.org/licenses/>.
//

#pragma once
#include "lazyten/config.hh"
#ifdef LAZYTEN_HAVE_BOHRIUM

#include "BohriumTypes.hh"
#include "common.hh"
#include "lazyten/Base/Interfaces.hh"
#include <bhxx/bhxx.hpp>

/** A class for exposing Bohrium multi_array objects as vectors to lazyten */
namespace lazyten {

template <typename Scalar>
class BohriumVector : public MutableMemoryVector_i<Scalar>, public Stored_i {
  static_assert(std::is_same<double, Scalar>::value ||
                      std::is_same<float, Scalar>::value ||
                      std::is_same<std::complex<float>, Scalar>::value ||
                      std::is_same<std::complex<double>, Scalar>::value,
                "BohriumVector<Scalar> is currently only available for Scalar "
                "being one of double, float, complex<double>, complex<float>");

 public:
  typedef MutableMemoryVector_i<Scalar> base_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::real_type real_type;

  /** The corresponding family of Bohrium interface types */
  typedef BohriumTypes type_family;

  typedef Scalar* iterator;
  typedef const Scalar* const_iterator;

  /** \name Constructors */
  ///@{
  /** \brief Construct vector of fixed size and optionally fill with zeros.
   *
   * \param fill_zero   If true all entries are set to zero
   */
  explicit BohriumVector(size_type size, bool /*fill_zero*/ = true)
        : BohriumVector(std::vector<Scalar>(size, 0)) {}
  // Here we fill the vector with zeros in all cases

  /** \brief Construct from initialiser list */
  explicit BohriumVector(std::initializer_list<Scalar> list)
        : BohriumVector(list.begin(), list.end()) {}

  /** \brief Construct from std::vector */
  explicit BohriumVector(std::vector<scalar_type> v)
        : BohriumVector(v.begin(), v.end()) {}

  /** \brief Construct from input iterator */
  template <class InputIterator>
  BohriumVector(InputIterator first, InputIterator last)
        : m_array(bhxx::make_base_ptr(first, last),
                  {static_cast<size_t>(std::distance(first, last))}) {
    assert_greater(0, std::distance(first, last));
  }

  /** \brief Construct from Arbitrary Indexable_i */
  template <typename Indexable,
            typename = typename std::enable_if<IsIndexable<Indexable>::value>::type>
  explicit BohriumVector(Indexable i) : BohriumVector(i.begin(), i.end()) {}

  /** \brief Construct from a one-dimensional multi-array */
  explicit BohriumVector(bhxx::BhArray<Scalar> inner) : m_array(inner.shape) {
    assert_size(inner.rank(), 1);
    bhxx::identity(m_array, inner);
  }
  ///@}

  BohriumVector(BohriumVector&& other) = default;
  BohriumVector& operator=(BohriumVector&& other) = default;
  BohriumVector(const BohriumVector& other) : m_array({other.n_elem()}) {
    bhxx::identity(m_array, other.m_array);
  }

  BohriumVector& operator=(const BohriumVector& other) {
    if (m_array.shape != other.m_array.shape) {
      // If the shape changes we need to allocate a new
      // object first and then call the identity operation
      m_array = bhxx::BhArray<Scalar>(other.m_array.shape);
    }
    bhxx::identity(m_array, other.m_array);
    return *this;
  }

  /** \name Size of the vector */
  size_type n_elem() const override { return m_array.shape[0]; }

  /** \name Size of the vector */
  size_type size() const override { return m_array.shape[0]; }

  /** \name Data access
   */
  ///@{
  /** \brief return an element of the vector    */
  scalar_type operator()(size_type i) const override {
    assert_range(0, i, n_elem());
    return memptr()[i];
  }

  /** \brief return an element of the vector */
  scalar_type operator[](size_type i) const override {
    assert_range(0, i, n_elem());
    return memptr()[i];
  }

  /** \brief return an element of the vector    */
  scalar_type& operator()(size_type i) override {
    assert_range(0, i, n_elem());
    return memptr()[i];
  }

  /** \brief return an element of the vector */
  scalar_type& operator[](size_type i) override {
    assert_range(0, i, n_elem());
    return memptr()[i];
  }
  ///@}

  /** Set all elements of the vector to zero */
  void set_zero() override { bhxx::identity(m_array, 0); }

  /** \name Iterators
   */
  ///@{
  /** Return an iterator to the beginning */
  iterator begin() { return memptr(); }

  /** Return a const_iterator to the beginning */
  const_iterator begin() const { return memptr(); }

  /** Return a const_iterator to the beginning */
  const_iterator cbegin() const { return memptr(); }

  /** Return an iterator to the end */
  iterator end() { return memptr() + size(); }

  /** Return a const_iterator to the end */
  const_iterator end() const { return memptr() + size(); }

  /** Return a const_iterator to the end */
  const_iterator cend() const { return memptr() + size(); }
  ///@}

  /** \name Vector operations */
  ///@{
  /** Scale vector by a scalar value */
  BohriumVector& operator*=(scalar_type s) {
    assert_finite(s);
    bhxx::multiply(m_array, m_array, s);
    return *this;
  }

  /** Divide all vector entries by a scalar value */
  BohriumVector& operator/=(scalar_type s) {
    assert_nonzero(s);
    assert_finite(s);
    bhxx::divide(m_array, m_array, s);
    return *this;
  }

  /* Add a vector to this one */
  BohriumVector& operator+=(const BohriumVector& other) {
    assert_size(n_elem(), other.n_elem());
    bhxx::add(m_array, m_array, other.m_array);
    return *this;
  }

  /* Add a vector to this one */
  BohriumVector& operator-=(const BohriumVector& other) {
    assert_size(n_elem(), other.n_elem());
    bhxx::subtract(m_array, m_array, other.m_array);
    return *this;
  }

  bool operator==(const BohriumVector& other) const {
    using bhxx::BhArray;
    if (m_array.n_elem() != other.m_array.n_elem()) return false;

    // Final result
    BhArray<bool> is_equal({1});
    {
      // Compute elementwise equality and reduce using &&
      BhArray<bool> eq({m_array.n_elem()});
      bhxx::equal(eq, m_array, other.m_array);
      logical_and_reduce(is_equal, eq, 0);
    }

    return bhxx::as_scalar(is_equal);
  }

  bool operator!=(const BohriumVector& other) const {
    using bhxx::BhArray;
    if (m_array.n_elem() != other.m_array.n_elem()) return true;

    // Final result
    BhArray<bool> is_unequal({1});
    {
      // Compute elementwise unequality and reduce using ||
      BhArray<bool> ueq({m_array.n_elem()});
      bhxx::not_equal(ueq, m_array, other.m_array);
      logical_or_reduce(is_unequal, ueq, 0);
    }

    return bhxx::as_scalar(is_unequal);
  }
  ///@}

  /** Read-only access to the raw memory */
  const scalar_type* memptr() const override {
    prepare_array_memory_for_use();
    return m_array.data();
  }

  /** Access to the raw memory */
  scalar_type* memptr() override {
    prepare_array_memory_for_use();
    return m_array.data();
  }

  //@{
  /** Access to the inner BhArray */
  const bhxx::BhArray<Scalar>& bh_array() const { return m_array; }
  bhxx::BhArray<Scalar>& bh_array() { return m_array; }
  //@}

 private:
  bhxx::BhArray<Scalar> m_array;

  /** Initialise the m_array object and/or make it contiguous */
  void prepare_array_memory_for_use() const {
    // Cast the const away ...
    auto& array = const_cast<bhxx::BhArray<Scalar>&>(m_array);

    if (!array.is_contiguous()) {
      // Get contiguous representation of the array:
      bhxx::BhArray<Scalar> contiguous{array.shape};
      identity(contiguous, array);
      array = std::move(contiguous);
    }
    sync(array);
    bhxx::Runtime::instance().flush();

    assert_internal(m_array.is_contiguous());
    assert_internal(m_array.data() != nullptr);
  }
};

//
// Standard operations
//
template <typename Scalar>
BohriumVector<Scalar> operator*(Scalar s, BohriumVector<Scalar> m) {
  m *= s;
  return m;
}

template <typename Scalar>
BohriumVector<Scalar> operator*(BohriumVector<Scalar> m, Scalar s) {
  return s * m;
}

template <typename Scalar>
BohriumVector<Scalar> operator/(BohriumVector<Scalar> m, Scalar s) {
  m /= s;
  return m;
}

template <typename Scalar>
BohriumVector<Scalar> operator-(BohriumVector<Scalar> v) {
  return Scalar(-1) * v;
}

//
// Add and subtract small matrices
//
template <typename Scalar>
BohriumVector<Scalar> operator-(BohriumVector<Scalar> lhs,
                                const BohriumVector<Scalar>& rhs) {
  lhs -= rhs;
  return lhs;
}

template <typename Scalar>
BohriumVector<Scalar> operator+(BohriumVector<Scalar> lhs,
                                const BohriumVector<Scalar>& rhs) {
  lhs += rhs;
  return lhs;
}

//
// Specialisations of fallback operations
//

// -- dot
template <typename Scalar1, typename Scalar2>
typename std::common_type<Scalar1, Scalar2>::type dot(const BohriumVector<Scalar1>& A,
                                                      const BohriumVector<Scalar2>& B) {
  assert_size(A.size(), B.size());
  typedef typename std::common_type<Scalar1, Scalar2>::type ctype;
  bhxx::BhArray<ctype> cA({A.size()});
  bhxx::identity(cA, A.bh_array());
  bhxx::BhArray<ctype> cB({B.size()});
  bhxx::identity(cB, B.bh_array());

  bhxx::BhArray<ctype> res({A.size()});
  bhxx::multiply(res, cA, cB);

  bhxx::BhArray<ctype> dot({1});
  bhxx::add_reduce(dot, res, 0);

  return bhxx::as_scalar(dot);
}

// TODO Not there yet in Bohrium.
// template <typename Scalar1, typename Scalar2>
// typename std::common_type<Scalar1, Scalar2>::type cdot(
//      const BohriumVector<Scalar1>& A, const BohriumVector<Scalar2>& B) {
//  // arma does the complex conjugate in the first argument as well.
//  return arma::cdot(A.data(), B.data());
//}

// -- accumulate
template <typename Scalar>
Scalar accumulate(const BohriumVector<Scalar>& A) {
  bhxx::BhArray<Scalar> out({1});
  bhxx::add_reduce(out, A.bh_array(), 0);
  return bhxx::as_scalar(out);
}

// -- minmax
template <typename Scalar>
Scalar min(const BohriumVector<Scalar>& A) {
  bhxx::BhArray<Scalar> out({1});
  bhxx::minimum_reduce(out, A.bh_array(), 0);
  return bhxx::as_scalar(out);
}

template <typename Scalar>
Scalar max(const BohriumVector<Scalar>& A) {
  bhxx::BhArray<Scalar> out({1});
  bhxx::maximum_reduce(out, A.bh_array(), 0);
  return bhxx::as_scalar(out);
}

// -- norms
// Use default implementation for all norms ( linf, l1 and l2  )

template <typename Scalar,
          krims::enable_if_t<krims::IsComplexNumber<Scalar>::value, int> = 0>
typename BohriumVector<Scalar>::real_type norm_l2_squared(
      const BohriumVector<Scalar>& v) {
  typedef typename krims::RealTypeOf<Scalar>::type real_type;
  bhxx::BhArray<real_type> re({v.size()});
  bhxx::BhArray<real_type> im({v.size()});
  bhxx::real(re, v.bh_array());
  bhxx::imag(im, v.bh_array());

  bhxx::BhArray<real_type> resq({v.size()});
  bhxx::BhArray<real_type> imsq({v.size()});
  bhxx::power(resq, re, real_type(2));
  bhxx::power(imsq, im, real_type(2));

  bhxx::BhArray<real_type> sum({v.size()});
  bhxx::add(sum, resq, imsq);

  bhxx::BhArray<real_type> res({1});
  bhxx::add_reduce(res, sum, 0);

  return bhxx::as_scalar(res);
}

template <typename Scalar,
          krims::enable_if_t<!krims::IsComplexNumber<Scalar>::value, int> = 0>
Scalar norm_l2_squared(const BohriumVector<Scalar>& v) {
  return dot(v, v);
}

// -- elementwise
template <typename Scalar>
BohriumVector<typename krims::RealTypeOf<Scalar>::type> abs(
      const BohriumVector<Scalar>& v) {
  typedef typename krims::RealTypeOf<Scalar>::type real_type;
  bhxx::BhArray<real_type> out({v.size()});
  bhxx::absolute(out, v.bh_array());
  return BohriumVector<real_type>(std::move(out));
}

/*  TODO Does not exist in Bohrium yet.
template <typename Scalar>
BohriumVector<Scalar> conj(const BohriumVector<Scalar>& v) {
  return BohriumVector<Scalar>(arma::conj(v.data()));
}
*/

template <typename Scalar>
BohriumVector<Scalar> sqrt(const BohriumVector<Scalar>& v) {
  bhxx::BhArray<Scalar> out({v.size()});
  bhxx::sqrt(out, v.bh_array());
  return BohriumVector<Scalar>(std::move(out));
}

template <typename Scalar>
BohriumVector<Scalar> square(const BohriumVector<Scalar>& v) {
  bhxx::BhArray<Scalar> out({v.size()});
  bhxx::power(out, v.bh_array(), static_cast<Scalar>(2));
  return BohriumVector<Scalar>(std::move(out));
}

}  // end namespace lazyten
#endif  // LAZYTEN_HAVE_BOHRIUM
