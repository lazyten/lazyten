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

#include "mutable_vector_tests.hh"
#include <linalgwrap/IteratorVector.hh>

namespace linalgwrap {

// TODO
// The fact that we have to badly overload this internal class shows
// that the testing library still has a fundamental flaw:
// - Objects, which are needed inside the tests may need to be generated
//   on the fly, but this is not possible if they are not stored objects,
//   which may well happen.
namespace gen {
namespace detail {

static std::vector<double> rdata = {1.};
template <>
struct NumericTensor<VectorMemoryWrapper<typename std::vector<double>::iterator>> {
  typedef std::vector<double> container_type;
  typedef VectorMemoryWrapper<typename container_type::iterator> vector_type;

  static rc::Gen<vector_type> numeric_tensor(typename vector_type::size_type count) {
    if (rdata.size() < 100) {
      rdata = *gen::scale(10., gen::numeric_container<container_type>(100));
    }

    return rc::gen::exec([&] {
      RC_PRE(count < rdata.size());
      size_t start = *rc::gen::inRange(static_cast<size_t>(0u), rdata.size());
      size_t end = *rc::gen::inRange(start, rdata.size());
      return make_vector_mem_wrap(rdata.begin() + start, rdata.begin() + end);
    });
  }

  static rc::Gen<vector_type> numeric_tensor() {
    return rc::gen::exec([] { return *numeric_tensor(*gen::numeric_size<1>()); });
  }
};

static std::vector<std::complex<double>> cdata = {1.};
template <>
struct NumericTensor<
      VectorMemoryWrapper<typename std::vector<std::complex<double>>::iterator>> {
  typedef std::vector<std::complex<double>> container_type;
  typedef VectorMemoryWrapper<typename container_type::iterator> vector_type;

  static rc::Gen<vector_type> numeric_tensor(typename vector_type::size_type count) {
    if (cdata.size() < 100) {
      cdata = *gen::scale(10., gen::numeric_container<container_type>(100));
    }

    return rc::gen::exec([&] {
      RC_PRE(count < cdata.size());
      size_t start = *rc::gen::inRange(static_cast<size_t>(0u), cdata.size());
      size_t end = *rc::gen::inRange(start, cdata.size());
      return make_vector_mem_wrap(cdata.begin() + start, cdata.begin() + end);
    });
  }

  static rc::Gen<vector_type> numeric_tensor() {
    return rc::gen::exec([] { return *numeric_tensor(*gen::numeric_size<1>()); });
  }
};

}  // namespace detail
}  // namespace gen

namespace tests {
namespace vector_memory_wrapper_tests {
using namespace rc;
using namespace vector_tests;

template <typename S>
struct vectorgen {
  // Assume the args_type is std::vector<S>
  typedef typename std::vector<S>::iterator iterator_type;
  typedef VectorMemoryWrapper<iterator_type> return_type;
  typedef GeneratorLibrary<return_type> genlib;
  typedef typename genlib::args_type args_type;
  static_assert(std::is_same<std::vector<S>, args_type>::value,
                "Assumption about args_type was wrong.");

  return_type operator()(args_type args) {
    m_data = args;
    return make_vector_mem_wrap(m_data.begin(), m_data.end());
  }

 private:
  args_type m_data;
};

template <typename S>
using genlib = typename vectorgen<S>::genlib;

TEST_CASE("VectorMemoryWrapper class", "[VectorMemoryWrapper]") {
  SECTION("Default mutable vector tests") {
    // Decrease tolerance to require a more accurate numerical agreement
    // for passing.
    auto lowertol = NumCompConstants::change_temporary(
          0.1 * krims::NumCompConstants::default_tolerance_factor);

    typedef double realtt;
    mutable_vector_tests::run_with_generator(
          genlib<realtt>::testgenerator(vectorgen<realtt>{}),
          "VectorMemoryWrapper<double>: ");

    typedef std::complex<double> cplxtt;
    mutable_vector_tests::run_with_generator(
          genlib<cplxtt>::testgenerator(vectorgen<cplxtt>{}),
          "VectorMemoryWrapper<complex double>: ");
  }
}

}  // builtin_vector_tests
}  // tests
}  // linalgwrap
