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
#include "Base/Interfaces.hh"
#include "StoredMatrix_i.hh"
#include <chrono>
#include <krims/TypeUtils.hh>
#include <random>

namespace lazyten {

/** Functor to generate a random real scalar value
 *
 * \note The values are only pseudorandom and *not* cryptographically
 *       safe in any way.
 * */
template <typename T, bool Complex = krims::IsComplexNumber<T>::value>
struct RandomScalar {
  static_assert(std::is_floating_point<T>::value,
                "T needs to be a floating point value here.");

  /** Construct a RandomScalar functor and make sure the generated
   *  values stay in the range [min,max)
   */
  RandomScalar(T min, T max) : distribution(min, max) {}

  /** Construct a RandomScalar functor, which produces
   * random values in the range [-100,100) */
  RandomScalar() : RandomScalar(-100, 100) {}

  /** Return a random scalar */
  T operator()() { return distribution(engine); }

  static std::default_random_engine engine;
  std::uniform_real_distribution<T> distribution;
};

/** Functor to generate a random complex scalar value
 *
 * Essentially uses the real version above to generate both the
 * value of the real and the imaginary type */
template <typename T, bool Complex>
std::default_random_engine RandomScalar<T, Complex>::engine{
      static_cast<typename std::default_random_engine::result_type>(
            std::chrono::system_clock::now().time_since_epoch().count())};

template <typename T>
struct RandomScalar<T, true> : public RandomScalar<typename krims::RealTypeOf<T>::type> {
  using real_type = typename krims::RealTypeOf<T>::type;
  using base_type = RandomScalar<real_type>;

  /** Construct a RandomScalar functor, which produces
   * random values with real and imaginary part in the range [-100,100) */
  RandomScalar() : base_type() {}

  /** Construct a RandomScalar functor and make sure the generated
   *  real and imaginary parts stay in the range [min,max)
   */
  RandomScalar(real_type min, real_type max) : base_type(min, max) {}

  /** Return a random scalar */
  T operator()() {
    const real_type real = base_type::operator()();
    const real_type imag = base_type::operator()();
    return T(real, imag);
  }
};

//@{
/** Return a random scalar using a RandomScalar functor in the default range */
template <typename T,
          typename Enabled = krims::enable_if_t<std::is_arithmetic<T>::value ||
                                                krims::IsComplexNumber<T>::value>>
T random() {
  return RandomScalar<T>{}();
}
//@}

/** Return a random indexable constructed from the passed argument
 *  and filled with random values generated using a RandomScalar functor
 *  in the default range */
template <typename Stored, typename... Args,
          typename Enabled = krims::enable_if_t<IsStoredVector<Stored>::value ||
                                                IsStoredMatrix<Stored>::value>>
Stored random(Args&&... args) {
  Stored stored(std::forward<Args>(args)...);
  std::generate(stored.begin(), stored.end(),
                RandomScalar<typename Stored::scalar_type>{});
  return stored;
}
}  // namespace lazyten
