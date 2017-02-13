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
#include <cmath>
#include <complex>
#include <cstdlib>
#include <krims/TypeUtils.hh>

namespace linalgwrap {
namespace detail {
struct PlusFctr {
  template <typename T1, typename T2>
  auto operator()(const T1& t1, const T2& t2) -> decltype(t1 + t2) {
    return t1 + t2;
  }
};

struct MinusFctr {
  template <typename T1, typename T2>
  auto operator()(const T1& t1, const T2& t2) -> decltype(t1 - t2) {
    return t1 - t2;
  }
};

struct MultipliesFctr {
  template <typename T1, typename T2>
  auto operator()(const T1& t1, const T2& t2) -> decltype(t1 * t2) {
    return t1 * t2;
  }
};

struct DividesFctr {
  template <typename T1, typename T2>
  auto operator()(const T1& t1, const T2& t2) -> decltype(t1 / t2) {
    return t1 / t2;
  }
};

struct NegateFctr {
  template <typename T>
  auto operator()(const T& t) -> decltype(-t) {
    return t;
  }
};

template <typename S>
struct ScaleByFctr {
  S sfac;
  ScaleByFctr(S sfac_) : sfac(sfac_) {}

  template <typename T>
  auto operator()(const T& t) -> decltype(sfac * t) {
    return sfac * t;
  }
};

template <typename S>
struct DivideByFctr {
  S sfac;
  DivideByFctr(S sfac_) : sfac(sfac_) {}

  template <typename T>
  auto operator()(const T& t) -> decltype(t / sfac) {
    return t / sfac;
  }
};

struct AbsFctr {
  template <typename Scalar>
  auto operator()(const Scalar& s) const -> decltype(std::abs(s)) {
    return std::abs(s);
  }
};

struct SqrtFctr {
  template <typename Scalar>
  auto operator()(const Scalar& s) const -> decltype(std::sqrt(s)) {
    return std::sqrt(s);
  }
};

struct SquareFctr {
  template <typename Scalar>
  auto operator()(const Scalar& s) const -> decltype(s * s) {
    return s * s;
  }
};

struct ConjFctr {
  template <typename Scalar>
  auto operator()(const std::complex<Scalar>& s) const -> decltype(std::conj(s)) {
    return std::conj(s);
  }

  template <typename Scalar>
  typename std::enable_if<!krims::IsComplexNumber<Scalar>::value, Scalar>::type
  operator()(const Scalar& s) const {
    return s;
  }
};

struct RealFctr {
  template <typename Scalar>
  auto operator()(const Scalar& s) const -> decltype(std::real(s)) {
    return std::real(s);
  }
};

struct ImagFctr {
  template <typename Scalar>
  auto operator()(const Scalar& s) const -> decltype(std::imag(s)) {
    return std::imag(s);
  }
};

}  // namespace detail
}  // namespace linalgwrap
