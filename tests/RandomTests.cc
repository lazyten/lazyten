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

#include <catch.hpp>
#include <lazyten/SmallMatrix.hh>
#include <lazyten/SmallVector.hh>
#include <lazyten/random.hh>

namespace lazyten {
namespace tests {
struct RandomConstants {
  // Minimal and maximal value of the default range
  double min = -100;
  double max = 100;

  size_t runs = 1000;

  virtual bool payload_assertion() const = 0;

  bool operator()() const {
    for (size_t i = 0; i < runs; ++i) {
      if (!payload_assertion()) return false;
    }
    return true;
  }

  virtual ~RandomConstants() = default;
  RandomConstants() = default;
  RandomConstants(RandomConstants&&) = default;
  RandomConstants(const RandomConstants&) = default;
  RandomConstants& operator=(RandomConstants&&) = default;
  RandomConstants& operator=(const RandomConstants&) = default;
};

template <typename T>
struct TestableReal : public RandomConstants {
  bool payload_assertion() const override {
    T res = random<T>();
    INFO(res);
    return min <= res && res < max;
  }
};

template <typename T>
struct TestableComplex : public RandomConstants {
  bool payload_assertion() const override {
    T res = random<T>();
    INFO(res);
    const bool real_ok = min <= res.real() && res.real() < max;
    const bool imag_ok = min <= res.imag() && res.imag() < max;
    return real_ok && imag_ok;
  }
};

template <typename T>
struct TestableIndexableBase : public RandomConstants {
  virtual T generate_random() const = 0;

  bool payload_assertion() const override {
    T res = generate_random();
    INFO(res);

    return std::all_of(
          std::begin(res), std::end(res),
          [this](const typename T::scalar_type& t) { return min <= t && t < max; });
  }
};

template <typename T>
struct TestableMatrix : public TestableIndexableBase<T> {
  T generate_random() const { return random<T>(100, 100); }
};

template <typename T>
struct TestableVector : public TestableIndexableBase<T> {
  T generate_random() const { return random<T>(100); }
};

TEST_CASE("Random number/vector/matrix generator", "[random]") {
  SECTION("Test double numbers in range [12,15)") {
    RandomScalar<double> random(12, 15);

    bool ok = true;
    for (size_t i = 0; i < 1000; ++i) {
      const double val = random();
      if (!(12 <= val && val < 15)) {
        ok = false;
        break;
      }
    }
    CHECK(ok);
  }

  SECTION("Test float numbers") { CHECK(TestableReal<float>{}()); }
  SECTION("Test double numbers") { CHECK(TestableReal<double>{}()); }
  SECTION("Test complex<double> numbers") {
    CHECK(TestableComplex<std::complex<double>>{}());
  }

  typedef SmallMatrix<double> matrix_type;
  typedef SmallVector<double> vector_type;
  SECTION("Test vector generation") {
    TestableVector<vector_type> testable{};
    testable.runs = 50;
    CHECK(testable());
  }
  SECTION("Test matrix generation") {
    TestableMatrix<matrix_type> testable{};
    testable.runs = 10;
    CHECK(testable());
  }

}  // Random object generator

}  // namespace tests
}  // namespace lazyten
